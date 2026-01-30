#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

import os, re
import h5py
import pandas as pd
import numpy as np
import endf

from data import PERIODIC_TABLE, SECTIONS_ABBREVS, SUBSHELL_LABELS
from function import parse_mf26_mt525, linear_interpolation, build_pdf, small_angle_scattering_cosine

def extract_sections(mat, mf, mt):
    """
    Extracts section (mf,mt) from an endf.Material and returns a pandas DataFrame.
    
    Parameters:
        mat    : endf.Material         # loaded ENDF file object
        mf,mt  : int, int             # material function and reaction type
    Returns:
        DataFrame                    # one of two formats depending on mf
    """

    key = (mf, mt)
    if key not in SECTIONS_ABBREVS:
        raise KeyError(f"No abbreviation defined for MF={mf}, MT={mt}")
    
    abbrev = SECTIONS_ABBREVS[key]
    sec = mat.section_data[key]

    # --- MF=23: Cross Sections ---
    if mf == 23:
        data = sec['sigma']   # Tabulated1D: .x: energies, .y: xs values
        df_y = pd.DataFrame({
            'energy_eV':      data.x,
            'cross_section':  data.y
        })

        df_y['section'] = abbrev

        df_y.attrs = {
            'breakpoints': data.breakpoints,
            'interpolation': data.interpolation
        }

        return None, df_y
    
    # --- MF=26: Angular & Energy Distributions ---
    elif mf == 26:
        prod = sec['products'][0]

        # Metadata
        ZAP   = prod.get('ZAP')
        AWI   = prod.get('AWI')
        LAW   = prod.get('LAW')
        y_tab = prod.get('y')

        # Primary distribution dict
        dist = prod.get('distribution') or {}

        ET = None
        if 'ET' in dist and dist['ET'] is not None:
            ET = dist['ET']

        LANG = None
        if 'LANG' in dist and dist['LANG'] is not None:
            LANG = dist['LANG']

        LEP = None
        if 'LEP' in dist and dist['LEP'] is not None:
            LEP = dist['LEP']

        NR = None
        if 'NR' in dist and dist['NR'] is not None:
            NR = dist['NR']

        NE = None
        if 'NE' in dist and dist['NE'] is not None:
            NE = dist['NE']

        E_inc = None
        if 'E' in dist and dist['E'] is not None:
            E_inc = dist['E']

        # 525 and 528
        if mt == 525:
            text = mat.section_text[key]
            data_26525 = parse_mf26_mt525(text)

            rows = []
            for reg in data_26525:
                E = reg['E_in']
                for mu, prob in reg['pairs']:
                    rows.append({
                    'inc_energy': E,
                    'mu':        mu,
                    'prob':      prob,
                    'section':   abbrev
                })

            df = pd.DataFrame(rows)
            df['section'] = abbrev

            df_y = pd.DataFrame({
                'y_inc_energy': y_tab.x,
                'y_yield': y_tab.y
            })

            df_y['section'] = abbrev

            df_y.attrs = {
                'breakpoints': y_tab.breakpoints,
                'interpolation': y_tab.interpolation
            }
            return df, df_y

        elif mt == 528:
            df_y = pd.DataFrame({
                'energy_eV':   ET.x,
                'avg_loss_eV': ET.y,

            })
            df_y['section'] = abbrev
            df_y.attrs = {
                'y_inc_energy': y_tab.x,
                'y_yield': y_tab.y,
                'breakpoints': y_tab.breakpoints,
                'interpolation': y_tab.interpolation
            }
            return None, df_y
        
        elif mt == 527:
            ph = next((p for p in sec['products'] if p.get('ZAP') == 0), None)
            el = next((p for p in sec['products'] if p.get('ZAP') == 11), None)

            df_y = None
            if el:
                el_dist = el.get('distribution') or {}
                ET = el_dist.get('ET')
                if ET is not None:
                    df_y = pd.DataFrame({
                        'energy_eV':   ET.x,
                        'avg_loss_eV': ET.y,
                        'section':     abbrev
                    })

            if not ph:
                return None, df_y
            
            dist = ph.get('distribution') or {}
            E_inc = dist.get('E', [])
            sub_list = dist.get('distribution', [])

            records = []
            for idx, sub in enumerate(sub_list):
                E_out = sub.get("E'", [])
                b_arr = sub.get('b')
                for eo, bb in zip(E_out, b_arr):
                    records.append({
                        'inc_energy_eV':  E_inc[idx],
                        'out_energy_eV':  eo,
                        'b':              bb,
                        'section':       abbrev
                    })
            df = pd.DataFrame(records)
            return df, df_y

        else: 
            df_y = pd.DataFrame({
                'y_inc_energy': y_tab.x,
                'y_yield': y_tab.y,
                'section': abbrev
            })

            df_y.attrs = {
                'breakpoints': y_tab.breakpoints,
                'interpolation': y_tab.interpolation
            }

        # Sub‐distribution list: one dict per incident energy

        sub_list = dist.get('distribution', [])

        # Flatten into records
        records = []
        for idx, sub in enumerate(sub_list):
            ND     = sub.get('ND')
            NA     = sub.get('NA')
            NW     = sub.get('NW')
            NEP    = sub.get('NEP')
            E_out  = sub.get("E'", [])
            b_arr  = sub.get('b')
            b_flat = b_arr.flatten() if b_arr is not None else []
            
            for Eo, b in zip(E_out, b_flat):
                records.append({
                    'section':        abbrev,
                    'ZAP':            ZAP,
                    'AWI':            AWI,
                    'LAW':            LAW,
                    'ET':             ET,
                    'LANG':           LANG,
                    'LEP':            LEP,
                    'NR':             NR,
                    'NE':             NE,
                    'inc_energy_eV':  E_inc[idx],
                    'out_energy_eV':  Eo,
                    'b' :             b,
                    'ND':             ND,
                    'NA':             NA,
                    'NW':             NW,
                    'NEP':            NEP,
                })

        df = pd.DataFrame.from_records(records)

        return df, df_y

    # Other MFs
    else:
        raise NotImplementedError(f"MF={mf} is not currently supported")


def save_element_h5(mat_path, out_dir):
    """
    Load one ENDF file, extract all sections, and write to <symbol>.h5.
    """

    fname = os.path.basename(mat_path)
    m = re.search(r'ZA(\d{3})000', fname)
    mat = endf.Material(mat_path)
    Z = int(m.group(1))
    entry = PERIODIC_TABLE.get(Z, {})
    sym = entry.get("symbol", f"Z{Z:03d}")

    sec0 = None
    for s in mat.section_data.values():
        if isinstance(s, dict) and all(k in s for k in ("ZA","AWR","NK")):
            sec0 = s
            break

    h5_path = os.path.join(out_dir, f"{sym}.h5")
    with h5py.File(h5_path, "w") as h5f:
        # metadata group
        meta = h5f.create_group("metadata")
        meta.create_dataset("Z", data=np.int64(Z))
        meta.create_dataset("Sym", data=sym)
        if sec0:
            meta.create_dataset("ZA", data=sec0["ZA"])
            meta.create_dataset("AWR", data=sec0["AWR"])
            meta.create_dataset("NK", data=sec0["NK"])

        grp_total   = h5f.create_group("total_xs")
        grp_scatter = h5f.create_group("elastic_scatter")
        grp_ion     = h5f.create_group("ionization")
        grp_brems   = h5f.create_group("bremsstrahlung")
        grp_exc     = h5f.create_group("excitation")

        subshell_mts = set(SUBSHELL_LABELS.keys())

        for (mf, mt), abbrev in SECTIONS_ABBREVS.items():
            if mf != 23 or (mf, mt) not in mat.section_data:
                continue
            df, df_y = extract_sections(mat, mf, mt)
            if df_y is None or df_y.empty:
                continue

            if mt == 501:
                g = grp_total.create_group("cross_section")
            elif mt == 522:
                g_s = grp_ion.create_group("cross_section")
                g = g_s.create_group("total")
            elif mt == 525:
                g_s = grp_scatter.create_group("cross_section")
                g = g_s.create_group("large_angle")
            elif mt == 526:
                g_s = grp_scatter.require_group("cross_section")
                g = g_s.create_group("total")
            elif mt == 527:
                g = grp_brems.create_group("cross_section")
            elif mt == 528:
                g = grp_exc.create_group("cross_section")
            elif mt in subshell_mts:
                xs_root = grp_ion.require_group("cross_section")
                label = SUBSHELL_LABELS[mt]
                g = xs_root.require_group(label)
            else:
                continue

            bpts  = getattr(df_y, "attrs", {}).get("breakpoints", None)
            interp = getattr(df_y, "attrs", {}).get("interpolation", None)
            if bpts is not None:
                g.create_dataset("breakpoints", data=np.asarray(bpts, dtype="f8"), dtype="f8")
            if interp is not None:
                g.create_dataset("interpolation", data=np.asarray(interp, dtype="f8"), dtype="f8")

            g.create_dataset("energy",        
                             data=df_y["energy_eV"].to_numpy(dtype="f8"),     
                             dtype="f8")
            g.create_dataset("cross_section", 
                             data=df_y["cross_section"].to_numpy(dtype="f8"), 
                             dtype="f8")

        for (mf, mt), _ in SECTIONS_ABBREVS.items():
            if mf != 26 or (mf, mt) not in mat.section_data:
                continue
            df, df_y = extract_sections(mat, mf, mt)

            if mt == 525:
                g_s = grp_scatter.require_group("distributions")
                g = g_s.require_group("large_angle")
                if df_y is not None:
                    g.create_dataset("y_inc_energy",
                                     data=df_y["y_inc_energy"].to_numpy(dtype="f8"),
                                     dtype="f8")
                    g.create_dataset("y_yield",
                                     data=df_y["y_yield"].to_numpy(dtype="f8"),
                                     dtype="f8")
                if df is not None and not df.empty:
                    g.create_dataset("inc_energy",
                                     data=df["inc_energy"].to_numpy(dtype="f8"),
                                     dtype="f8")
                    g.create_dataset("mu",
                                     data=df["mu"].to_numpy(dtype="f8"),         
                                     dtype="f8")
                    g.create_dataset("probability",
                                     data=df["prob"].to_numpy(dtype="f8"),
                                     dtype="f8")
                continue

            if mt == 528:
                g = grp_exc.require_group("distributions")
                if df_y is not None:
                    g.create_dataset("loss_inc_energy",
                                    data=df_y["energy_eV"].to_numpy(dtype="f8"),   dtype="f8")
                    g.create_dataset("avg_loss",
                                    data=df_y["avg_loss_eV"].to_numpy(dtype="f8"),
                                    dtype="f8")
                    if hasattr(df_y, "attrs"):
                        bpts  = df_y.attrs.get("breakpoints", None)
                        interp = df_y.attrs.get("interpolation", None)
                        if bpts is not None:
                            g.create_dataset("y_breakpoints",
                                            data=np.asarray(bpts, dtype="f8"),
                                            dtype="f8")
                        if interp is not None:
                            g.create_dataset("y_interpolation", 
                                            data=np.asarray(interp, dtype="f8"), 
                                            dtype="f8")
                continue

            if mt == 527:
                g = grp_brems.require_group("distributions")
                if df is not None and not df.empty:
                    g.create_dataset("inc_energy", 
                                    data=df["inc_energy_eV"].to_numpy(dtype="f8"),
                                    dtype="f8")
                    g.create_dataset("out_energy", 
                                    data=df["out_energy_eV"].to_numpy(dtype="f8"), 
                                    dtype="f8")
                    g.create_dataset("b",
                                    data=df["b"].to_numpy(dtype="f8"),
                                    dtype="f8")
                if df_y is not None:
                    g.create_dataset("loss_inc_energy",
                                    data=df_y["energy_eV"].to_numpy(dtype="f8"),
                                    dtype="f8")
                    g.create_dataset("avg_loss",
                                    data=df_y["avg_loss_eV"].to_numpy(dtype="f8"), 
                                    dtype="f8")
                    if hasattr(df_y, "attrs"):
                        bpts  = df_y.attrs.get("breakpoints", None)
                        interp = df_y.attrs.get("interpolation", None)
                        if bpts is not None:
                            g.create_dataset("y_breakpoints", 
                                            data=np.asarray(bpts, dtype="f8"), 
                                            dtype="f8")
                        if interp is not None:
                            g.create_dataset("y_interpolation", 
                                            data=np.asarray(interp, dtype="f8"), 
                                            dtype="f8")
                continue

            if mt in SUBSHELL_LABELS:
                shell_root = grp_ion.require_group("distributions")
                g = shell_root.require_group(SUBSHELL_LABELS[mt])
                if df is not None and not df.empty:
                    g.create_dataset("inc_energy",
                                    data=df["inc_energy_eV"].to_numpy(dtype="f8"),
                                    dtype="f8")
                    g.create_dataset("out_energy",
                                    data=df["out_energy_eV"].to_numpy(dtype="f8"),
                                    dtype="f8")
                    g.create_dataset("b",
                                    data=df["b"].to_numpy(dtype="f8"),
                                    dtype="f8")
                if df_y is not None:
                    if "y_inc_energy" in df_y:
                        g.create_dataset("y_inc_energy",
                                        data=df_y["y_inc_energy"].to_numpy(dtype="f8"),
                                        dtype="f8")
                        g.create_dataset("binding_energy",
                                        data=df_y["y_inc_energy"][0],
                                        dtype="f8")
                    if "y_yield" in df_y:
                        g.create_dataset("y_yield",
                                        data=df_y["y_yield"].to_numpy(dtype="f8"),
                                        dtype="f8")

                continue

    print(f"Saved! {h5_path}")

def create_mcdc_file(in_path: str, out_dir: str) -> str:
    with h5py.File(in_path, "r") as src:
        Z   = int(src["/metadata/Z"][()])
        AWR = float(src["/metadata/AWR"][()])
        Sym = src["/metadata/Sym"].asstr()[()]

        total_energy = src["/total_xs/cross_section/energy"][:]
        total_xs     = src["/total_xs/cross_section/cross_section"][:]

        total_scat_xs_energy = src["/elastic_scatter/cross_section/total/energy"][:]
        total_scat_xs        = src["/elastic_scatter/cross_section/total/cross_section"][:]
        la_scat_xs_energy    = src["/elastic_scatter/cross_section/large_angle/energy"][:]
        la_scat_xs           = src["/elastic_scatter/cross_section/large_angle/cross_section"][:]

        la_scat_dist_mu     = src["/elastic_scatter/distributions/large_angle/mu"][:]
        la_scat_dist_energy = src["/elastic_scatter/distributions/large_angle/inc_energy"][:]
        la_scat_dist_prob   = src["/elastic_scatter/distributions/large_angle/probability"][:]

        brem_xs_energy = src["/bremsstrahlung/cross_section/energy"][:]
        brem_xs        = src["/bremsstrahlung/cross_section/cross_section"][:]
        brem_avg_loss_energy = src["/bremsstrahlung/distributions/loss_inc_energy"][:]
        brem_avg_loss        = src["/bremsstrahlung/distributions/avg_loss"][:]

        exc_xs_energy = src["/excitation/cross_section/energy"][:]
        exc_xs        = src["/excitation/cross_section/cross_section"][:]
        exc_avg_loss_energy = src["/excitation/distributions/loss_inc_energy"][:]
        exc_avg_loss        = src["/excitation/distributions/avg_loss"][:]

        ion_total_xs_energy = src["/ionization/cross_section/total/energy"][:]
        ion_total_xs        = src["/ionization/cross_section/total/cross_section"][:]

        subshell_names = []
        xs_root = src["/ionization/cross_section"]
        if isinstance(xs_root, h5py.Group):
            for name in xs_root.keys():
                if name != "total":
                    subshell_names.append(name)

        subshell_xs_data = {}
        for name in subshell_names:
            xs = src[f"/ionization/cross_section/{name}/cross_section"][:]
            xe = src[f"/ionization/cross_section/{name}/energy"][:]
            subshell_xs_data[name] = {f"xs_{name}": xs, f"xs_energy_{name}": xe}

        subshell_dist_data = {}
        for name in subshell_names:
            inc  = src[f"/ionization/distributions/{name}/inc_energy"][:]
            out  = src[f"/ionization/distributions/{name}/out_energy"][:]
            prob = src[f"/ionization/distributions/{name}/b"][:]
            be   = src[f"/ionization/distributions/{name}/binding_energy"][()]
            subshell_dist_data[name] = {
                f"dist_inc_energy_{name}": inc,
                f"dist_out_energy_{name}": out,
                f"dist_prob_{name}": prob,
                f"binding_energy_{name}": be,
            }

    xs_energy_grid = total_energy
    xs_sc_total    = linear_interpolation(xs_energy_grid, total_scat_xs_energy, total_scat_xs)
    xs_sc_la       = linear_interpolation(xs_energy_grid, la_scat_xs_energy,    la_scat_xs)
    xs_brem        = linear_interpolation(xs_energy_grid, brem_xs_energy,       brem_xs)
    xs_exc         = linear_interpolation(xs_energy_grid, exc_xs_energy,        exc_xs)
    xs_ion_total   = linear_interpolation(xs_energy_grid, ion_total_xs_energy,  ion_total_xs)
    xs_sc_sa       = xs_sc_total - xs_sc_la

    xs_total     = np.asarray(total_xs,    "f8")
    xs_sc_total  = np.asarray(xs_sc_total, "f8")
    xs_sc_la     = np.asarray(xs_sc_la,    "f8")
    xs_sc_sa     = np.asarray(xs_sc_sa,    "f8")
    xs_brem      = np.asarray(xs_brem,     "f8")
    xs_exc       = np.asarray(xs_exc,      "f8")
    xs_ion_total = np.asarray(xs_ion_total,"f8")

    mask_sa = xs_sc_sa > 0.0
    eg_sa, off_sa, val_sa, pdf_sa = small_angle_scattering_cosine(int(Z), np.array(xs_energy_grid)[mask_sa], n_mu=200)

    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{Sym}.h5")

    with h5py.File(out_path, "w") as h5f:
        h5f.create_dataset("atomic_weight_ratio", data=AWR)
        h5f.create_dataset("atomic_number", data=Z)
        h5f.create_dataset("element_name", data=Sym)
        h5f.create_group("electron_reactions")
        h5f.create_dataset("electron_reactions/xs_energy_grid", data=xs_energy_grid)

        h5f.create_group("electron_reactions/total")
        h5f.create_dataset("electron_reactions/total/xs", data=xs_total)

        h5f.create_group("electron_reactions/elastic_scattering")
        h5f.create_dataset("electron_reactions/elastic_scattering/xs", data=xs_sc_total)

        h5f.create_group("electron_reactions/elastic_scattering/large_angle")
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/xs", data=xs_sc_la)
        h5f.create_group("electron_reactions/elastic_scattering/large_angle/scattering_cosine")
        eg, off, val, PDF = build_pdf(la_scat_dist_energy, la_scat_dist_mu, la_scat_dist_prob)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/energy_grid",  data=eg)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/energy_offset", data=off)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/value",        data=val)
        h5f.create_dataset("electron_reactions/elastic_scattering/large_angle/scattering_cosine/PDF",          data=PDF)

        h5f.create_group("electron_reactions/elastic_scattering/small_angle")
        h5f.create_dataset("electron_reactions/elastic_scattering/small_angle/xs", data=xs_sc_sa)
        h5f.create_group("electron_reactions/elastic_scattering/small_angle/scattering_cosine")
        h5f.create_dataset("electron_reactions/elastic_scattering/small_angle/scattering_cosine/energy_grid",  data=eg_sa)
        h5f.create_dataset("electron_reactions/elastic_scattering/small_angle/scattering_cosine/energy_offset", data=off_sa)
        h5f.create_dataset("electron_reactions/elastic_scattering/small_angle/scattering_cosine/value",        data=val_sa)
        h5f.create_dataset("electron_reactions/elastic_scattering/small_angle/scattering_cosine/PDF",          data=pdf_sa)

        h5f.create_group("electron_reactions/bremsstrahlung")
        h5f.create_dataset("electron_reactions/bremsstrahlung/xs", data=xs_brem)
        h5f.create_group("electron_reactions/bremsstrahlung/energy_loss")
        h5f.create_dataset("electron_reactions/bremsstrahlung/energy_loss/energy", data=brem_avg_loss_energy)
        h5f.create_dataset("electron_reactions/bremsstrahlung/energy_loss/value",  data=brem_avg_loss)

        h5f.create_group("electron_reactions/excitation")
        h5f.create_dataset("electron_reactions/excitation/xs", data=xs_exc)
        h5f.create_group("electron_reactions/excitation/energy_loss")
        h5f.create_dataset("electron_reactions/excitation/energy_loss/energy", data=exc_avg_loss_energy)
        h5f.create_dataset("electron_reactions/excitation/energy_loss/value",  data=exc_avg_loss)

        h5f.create_group("electron_reactions/ionization")
        h5f.create_dataset("electron_reactions/ionization/xs", data=xs_ion_total)
        subs = h5f.create_group("electron_reactions/ionization/subshells")
        for shell in subshell_names:
            g = subs.create_group(shell)
            xs_e = subshell_xs_data[shell][f"xs_energy_{shell}"]
            xs_v = subshell_xs_data[shell][f"xs_{shell}"]
            xs_on_grid = linear_interpolation(xs_energy_grid, xs_e, xs_v)
            g.create_dataset("xs", data=xs_on_grid)
            inc = subshell_dist_data[shell][f"dist_inc_energy_{shell}"]
            out = subshell_dist_data[shell][f"dist_out_energy_{shell}"]
            prob= subshell_dist_data[shell][f"dist_prob_{shell}"]
            egp, offp, valp, PDFp = build_pdf(inc, out, prob)
            dg = g.create_group("product")
            dg.create_dataset("energy_grid",  data=egp)
            dg.create_dataset("energy_offset",data=offp)
            dg.create_dataset("value",        data=valp)
            dg.create_dataset("PDF",          data=PDFp)
            g.create_dataset("binding_energy", data=subshell_dist_data[shell][f"binding_energy_{shell}"])

        print(f"Saved! {out_path}")
