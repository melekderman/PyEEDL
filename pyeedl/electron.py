#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

import os
import re
import h5py
import pandas as pd
import numpy as np
import endf

from data import *
from function import parse_mf26_mt525

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

        else: 
            df_y = pd.DataFrame({
                'y_inc_energy': y_tab.x,
                'y_yield': y_tab.y
            })
            df_y['section'] = abbrev

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

        # Create DataFrame from list of dicts

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
    m = re.search(r'Z(\d{3})000', fname)
    mat = endf.Material(mat_path)
    Z = int(m.group(1))
    entry = PERIODIC_TABLE.get(Z, {})
    sym = entry.get("symbol", f"Z{Z:03d}")

    sec0 = None
    for sec in mat.section_data.values():
        if all(k in sec for k in ("ZA","AWR","NK")):
            sec0 = sec
            break

    h5_path = os.path.join(out_dir, f"{sym}.h5")
    with h5py.File(h5_path, "w") as h5f:
        # metadata group
        meta = h5f.create_group("metadata")
        meta.attrs["Z"]   = Z
        meta.attrs["Sym"] = sym
        meta.attrs["ZA"]  = sec0["ZA"]
        meta.attrs["AWR"] = sec0["AWR"]
        meta.attrs["NK"]  = sec0["NK"]

        # cross_sections group
        xs_grp = h5f.create_group("cross_sections")
        for (mf, mt), abbrev in SECTIONS_ABBREVS.items():
            if mf != 23 or (mf, mt) not in mat.section_data:
                continue
            df, df_y = extract_sections(mat, mf, mt)
            g  = xs_grp.create_group(abbrev)
            g.create_dataset(
                "energy",
                data=df_y["energy_eV"].values,
                dtype="f8"
            )

            g.create_dataset(
                "cross_section",
                data=df_y["cross_section"].values,
                dtype="f8"
            )

            if hasattr(df_y, "attrs"):
                bpts = df_y.attrs.get("breakpoints", None)
                interp = df_y.attrs.get("interpolation", None)

                if bpts is not None:
                    g.create_dataset("breakpoints", data=bpts, dtype="f8")

                if interp is not None:
                    g.create_dataset("interpolation", data=interp, dtype="f8")

        # distributions (MF=26)
        dist_grp = h5f.create_group("distributions")
        for (mf, mt), abbrev in SECTIONS_ABBREVS.items():
            if mf != 26 or (mf, mt) not in mat.section_data:
                continue

            df, df_y = extract_sections(mat, mf, mt)
            g  = dist_grp.create_group(abbrev)

            if hasattr(df_y, "attrs"):
                bpts = df_y.attrs.get("breakpoints", None)
                interp = df_y.attrs.get("interpolation", None)

                if bpts is not None:
                    g.create_dataset("y_breakpoints", data=bpts, dtype="f8")

                if interp is not None:
                    g.create_dataset("y_interpolation", data=interp, dtype="f8")

            if mt == 525:
                g.create_dataset(
                    "y_inc_energy",
                    data=df_y["y_inc_energy"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "y_yield",
                    data=df_y["y_yield"].values,
                    dtype="f8"
                )

                g.create_dataset(
                    "inc_energy",
                    data=df["inc_energy"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "mu",
                    data=df["mu"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "probability",
                    data=df["prob"].values,
                    dtype="f8"
                )
                continue
        
            elif mt == 528:
                g.create_dataset(
                    "energy",
                    data=df_y["energy_eV"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "avg_loss",
                    data=df_y["avg_loss_eV"].values,
                    dtype="f8"
                )
                continue


            # full energy–energy distribution
            if not df.empty and 'inc_energy_eV' in df.columns:
                g.create_dataset(
                    "inc_energy",
                    data=df["inc_energy_eV"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "out_energy",
                    data=df["out_energy_eV"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "b",
                    data=df["b"].values,
                    dtype="f8"
                )
                g.create_dataset(
                    "y_yield",
                    data=np.array(df_y["y_yield"].tolist()),
                    dtype="f8"
                )
                g.create_dataset(
                    "y_inc_energy",
                    data=np.array(df_y["y_inc_energy"].tolist()),
                    dtype="f8"
                )

                for key in ("ZAP","AWI","LAW","ET","LANG","LEP","NR","NE"):
                    if key in df.columns:
                        val = df[key].iloc[0]
                        if pd.notnull(val):
                            g.attrs[key] = val
            else:
                dist_grp.pop(abbrev)

    print(f"Saved! {h5_path}")

