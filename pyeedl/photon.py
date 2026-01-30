#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Photon data processing module for EPDL (Evaluated Photon Data Library).

This module provides functions to extract and process photon interaction
cross sections from EPDL data files in ENDF format.
"""

import os
import re
import h5py
import pandas as pd
import numpy as np
import endf

from data import (
    PERIODIC_TABLE,
    PHOTON_SECTIONS_ABBREVS,
    PHOTON_MF_MT,
    SUBSHELL_LABELS,
)
from function import linear_interpolation, build_pdf


def extract_photon_sections(mat, mf, mt):
    """
    Extracts photon section (mf, mt) from an endf.Material and returns a pandas DataFrame.

    Parameters:
        mat    : endf.Material         # loaded ENDF file object
        mf, mt : int, int              # material function and reaction type

    Returns:
        DataFrame or tuple(DataFrame, DataFrame)
    """
    key = (mf, mt)
    if key not in PHOTON_SECTIONS_ABBREVS:
        raise KeyError(f"No abbreviation defined for photon MF={mf}, MT={mt}")

    abbrev = PHOTON_SECTIONS_ABBREVS[key]
    sec = mat.section_data[key]

    # --- MF=23: Photon Cross Sections ---
    if mf == 23:
        data = sec['sigma']
        df_y = pd.DataFrame({
            'energy_eV': data.x,
            'cross_section': data.y
        })
        df_y['section'] = abbrev
        df_y.attrs = {
            'breakpoints': data.breakpoints,
            'interpolation': data.interpolation
        }
        return None, df_y

    # --- MF=27: Form Factors and Scattering Functions ---
    elif mf == 27:
        data = sec['sigma']
        df_y = pd.DataFrame({
            'x': data.x,  # momentum transfer or energy depending on MT
            'y': data.y   # form factor or scattering function
        })
        df_y['section'] = abbrev
        df_y.attrs = {
            'breakpoints': data.breakpoints,
            'interpolation': data.interpolation
        }
        return None, df_y

    else:
        raise NotImplementedError(f"MF={mf} is not currently supported for photon data")


def save_photon_element_h5(mat_path, out_dir):
    """
    Load one EPDL ENDF file, extract all photon sections, and write to <symbol>.h5.

    Parameters:
        mat_path: Path to the EPDL ENDF file
        out_dir: Output directory for HDF5 files
    """
    fname = os.path.basename(mat_path)
    m = re.search(r'ZA(\d{3})000', fname)
    mat = endf.Material(mat_path)
    Z = int(m.group(1))
    entry = PERIODIC_TABLE.get(Z, {})
    sym = entry.get("symbol", f"Z{Z:03d}")

    # Get metadata from first section
    sec0 = None
    for s in mat.section_data.values():
        if isinstance(s, dict) and all(k in s for k in ("ZA", "AWR")):
            sec0 = s
            break

    os.makedirs(out_dir, exist_ok=True)
    h5_path = os.path.join(out_dir, f"{sym}.h5")

    with h5py.File(h5_path, "w") as h5f:
        # Metadata group
        meta = h5f.create_group("metadata")
        meta.create_dataset("Z", data=np.int64(Z))
        meta.create_dataset("Sym", data=sym)
        if sec0:
            meta.create_dataset("ZA", data=sec0["ZA"])
            meta.create_dataset("AWR", data=sec0["AWR"])

        # Create main groups
        grp_total = h5f.create_group("total_xs")
        grp_coherent = h5f.create_group("coherent_scattering")
        grp_incoherent = h5f.create_group("incoherent_scattering")
        grp_photoelectric = h5f.create_group("photoelectric")
        grp_pair = h5f.create_group("pair_production")
        grp_form_factors = h5f.create_group("form_factors")

        # Process MF=23 cross sections
        for (mf, mt), abbrev in PHOTON_SECTIONS_ABBREVS.items():
            if mf != 23 or (mf, mt) not in mat.section_data:
                continue

            df, df_y = extract_photon_sections(mat, mf, mt)
            if df_y is None or df_y.empty:
                continue

            # Determine target group based on MT
            if mt == 501:
                g = grp_total.create_group("cross_section")
            elif mt == 502:
                g = grp_coherent.create_group("cross_section")
            elif mt == 504:
                g = grp_incoherent.create_group("cross_section")
            elif mt == 516:
                g_s = grp_pair.create_group("cross_section")
                g = g_s.create_group("total")
            elif mt == 517:
                g_s = grp_pair.require_group("cross_section")
                g = g_s.create_group("nuclear")
            elif mt == 518:
                g_s = grp_pair.require_group("cross_section")
                g = g_s.create_group("electron")
            elif mt == 522:
                g_s = grp_photoelectric.create_group("cross_section")
                g = g_s.create_group("total")
            elif mt in SUBSHELL_LABELS:
                # Subshell photoelectric
                xs_root = grp_photoelectric.require_group("cross_section")
                label = SUBSHELL_LABELS[mt]
                g = xs_root.require_group(label)
            else:
                continue

            # Save cross section data
            bpts = getattr(df_y, "attrs", {}).get("breakpoints", None)
            interp = getattr(df_y, "attrs", {}).get("interpolation", None)
            if bpts is not None:
                g.create_dataset("breakpoints", data=np.asarray(bpts, dtype="f8"))
            if interp is not None:
                g.create_dataset("interpolation", data=np.asarray(interp, dtype="f8"))

            g.create_dataset("energy", data=df_y["energy_eV"].to_numpy(dtype="f8"))
            g.create_dataset("cross_section", data=df_y["cross_section"].to_numpy(dtype="f8"))

        # Process MF=27 form factors
        for (mf, mt), abbrev in PHOTON_SECTIONS_ABBREVS.items():
            if mf != 27 or (mf, mt) not in mat.section_data:
                continue

            df, df_y = extract_photon_sections(mat, mf, mt)
            if df_y is None or df_y.empty:
                continue

            if mt == 502:
                g = grp_form_factors.create_group("coherent")
                g.create_dataset("momentum_transfer", data=df_y["x"].to_numpy(dtype="f8"))
                g.create_dataset("form_factor", data=df_y["y"].to_numpy(dtype="f8"))
            elif mt == 504:
                g = grp_form_factors.create_group("incoherent")
                g.create_dataset("momentum_transfer", data=df_y["x"].to_numpy(dtype="f8"))
                g.create_dataset("scattering_function", data=df_y["y"].to_numpy(dtype="f8"))
            elif mt == 505:
                g = grp_form_factors.require_group("anomalous")
                g.create_dataset("energy", data=df_y["x"].to_numpy(dtype="f8"))
                g.create_dataset("imaginary", data=df_y["y"].to_numpy(dtype="f8"))
            elif mt == 506:
                g = grp_form_factors.require_group("anomalous")
                if "energy" not in g:
                    g.create_dataset("energy", data=df_y["x"].to_numpy(dtype="f8"))
                g.create_dataset("real", data=df_y["y"].to_numpy(dtype="f8"))

    print(f"Saved! {h5_path}")


def create_photon_mcdc_file(in_path: str, out_dir: str) -> str:
    """
    Create MCDC-format HDF5 file from raw photon data.

    Parameters:
        in_path: Path to raw HDF5 file
        out_dir: Output directory for MCDC file

    Returns:
        Path to created MCDC file
    """
    with h5py.File(in_path, "r") as src:
        Z = int(src["/metadata/Z"][()])
        AWR = float(src["/metadata/AWR"][()])
        Sym = src["/metadata/Sym"].asstr()[()]

        # Total cross section
        total_energy = src["/total_xs/cross_section/energy"][:]
        total_xs = src["/total_xs/cross_section/cross_section"][:]

        # Coherent scattering
        coherent_xs_energy = src["/coherent_scattering/cross_section/energy"][:]
        coherent_xs = src["/coherent_scattering/cross_section/cross_section"][:]

        # Incoherent scattering
        incoherent_xs_energy = src["/incoherent_scattering/cross_section/energy"][:]
        incoherent_xs = src["/incoherent_scattering/cross_section/cross_section"][:]

        # Photoelectric
        pe_total_xs_energy = src["/photoelectric/cross_section/total/energy"][:]
        pe_total_xs = src["/photoelectric/cross_section/total/cross_section"][:]

        # Pair production (if available)
        pair_total_xs = None
        pair_nuclear_xs = None
        pair_electron_xs = None
        if "/pair_production/cross_section/total" in src:
            pair_total_energy = src["/pair_production/cross_section/total/energy"][:]
            pair_total_xs = src["/pair_production/cross_section/total/cross_section"][:]
        if "/pair_production/cross_section/nuclear" in src:
            pair_nuclear_energy = src["/pair_production/cross_section/nuclear/energy"][:]
            pair_nuclear_xs = src["/pair_production/cross_section/nuclear/cross_section"][:]
        if "/pair_production/cross_section/electron" in src:
            pair_electron_energy = src["/pair_production/cross_section/electron/energy"][:]
            pair_electron_xs = src["/pair_production/cross_section/electron/cross_section"][:]

        # Form factors
        ff_coherent_x = None
        ff_coherent_y = None
        sf_incoherent_x = None
        sf_incoherent_y = None
        if "/form_factors/coherent" in src:
            ff_coherent_x = src["/form_factors/coherent/momentum_transfer"][:]
            ff_coherent_y = src["/form_factors/coherent/form_factor"][:]
        if "/form_factors/incoherent" in src:
            sf_incoherent_x = src["/form_factors/incoherent/momentum_transfer"][:]
            sf_incoherent_y = src["/form_factors/incoherent/scattering_function"][:]

        # Subshell photoelectric data
        subshell_names = []
        xs_root = src.get("/photoelectric/cross_section")
        if xs_root is not None:
            for name in xs_root.keys():
                if name != "total":
                    subshell_names.append(name)

        subshell_xs_data = {}
        for name in subshell_names:
            xs = src[f"/photoelectric/cross_section/{name}/cross_section"][:]
            xe = src[f"/photoelectric/cross_section/{name}/energy"][:]
            subshell_xs_data[name] = {"xs": xs, "energy": xe}

    # Use total energy as common grid
    xs_energy_grid = total_energy

    # Interpolate all cross sections to common grid
    xs_coherent = linear_interpolation(xs_energy_grid, coherent_xs_energy, coherent_xs)
    xs_incoherent = linear_interpolation(xs_energy_grid, incoherent_xs_energy, incoherent_xs)
    xs_pe_total = linear_interpolation(xs_energy_grid, pe_total_xs_energy, pe_total_xs)

    if pair_total_xs is not None:
        xs_pair_total = linear_interpolation(xs_energy_grid, pair_total_energy, pair_total_xs)
    else:
        xs_pair_total = np.zeros_like(xs_energy_grid)

    if pair_nuclear_xs is not None:
        xs_pair_nuclear = linear_interpolation(xs_energy_grid, pair_nuclear_energy, pair_nuclear_xs)
    else:
        xs_pair_nuclear = np.zeros_like(xs_energy_grid)

    if pair_electron_xs is not None:
        xs_pair_electron = linear_interpolation(xs_energy_grid, pair_electron_energy, pair_electron_xs)
    else:
        xs_pair_electron = np.zeros_like(xs_energy_grid)

    # Create output file
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{Sym}.h5")

    with h5py.File(out_path, "w") as h5f:
        h5f.create_dataset("atomic_weight_ratio", data=AWR)
        h5f.create_dataset("atomic_number", data=Z)
        h5f.create_dataset("element_name", data=Sym)

        # Photon reactions group
        h5f.create_group("photon_reactions")
        h5f.create_dataset("photon_reactions/xs_energy_grid", data=xs_energy_grid)

        # Total
        h5f.create_group("photon_reactions/total")
        h5f.create_dataset("photon_reactions/total/xs", data=total_xs)

        # Coherent scattering
        h5f.create_group("photon_reactions/coherent_scattering")
        h5f.create_dataset("photon_reactions/coherent_scattering/xs", data=xs_coherent)
        if ff_coherent_x is not None:
            h5f.create_group("photon_reactions/coherent_scattering/form_factor")
            h5f.create_dataset("photon_reactions/coherent_scattering/form_factor/momentum_transfer", data=ff_coherent_x)
            h5f.create_dataset("photon_reactions/coherent_scattering/form_factor/value", data=ff_coherent_y)

        # Incoherent scattering
        h5f.create_group("photon_reactions/incoherent_scattering")
        h5f.create_dataset("photon_reactions/incoherent_scattering/xs", data=xs_incoherent)
        if sf_incoherent_x is not None:
            h5f.create_group("photon_reactions/incoherent_scattering/scattering_function")
            h5f.create_dataset("photon_reactions/incoherent_scattering/scattering_function/momentum_transfer",
                               data=sf_incoherent_x)
            h5f.create_dataset("photon_reactions/incoherent_scattering/scattering_function/value",
                               data=sf_incoherent_y)

        # Photoelectric
        h5f.create_group("photon_reactions/photoelectric")
        h5f.create_dataset("photon_reactions/photoelectric/xs", data=xs_pe_total)
        subs = h5f.create_group("photon_reactions/photoelectric/subshells")
        for shell in subshell_names:
            g = subs.create_group(shell)
            xs_e = subshell_xs_data[shell]["energy"]
            xs_v = subshell_xs_data[shell]["xs"]
            xs_on_grid = linear_interpolation(xs_energy_grid, xs_e, xs_v)
            g.create_dataset("xs", data=xs_on_grid)

        # Pair production
        h5f.create_group("photon_reactions/pair_production")
        h5f.create_dataset("photon_reactions/pair_production/xs", data=xs_pair_total)
        h5f.create_group("photon_reactions/pair_production/nuclear")
        h5f.create_dataset("photon_reactions/pair_production/nuclear/xs", data=xs_pair_nuclear)
        h5f.create_group("photon_reactions/pair_production/electron")
        h5f.create_dataset("photon_reactions/pair_production/electron/xs", data=xs_pair_electron)

    print(f"Saved! {out_path}")
    return out_path
