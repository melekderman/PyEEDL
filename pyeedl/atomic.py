#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Atomic relaxation data processing module for EADL (Evaluated Atomic Data Library).

This module provides functions to extract and process atomic relaxation data
including binding energies, fluorescence yields, and transition probabilities
from EADL data files in ENDF format.
"""

import os
import re
import h5py
import pandas as pd
import numpy as np
import endf

from data import (
    PERIODIC_TABLE,
    SUBSHELL_DESIGNATORS,
)


def extract_atomic_relaxation(mat):
    """
    Extracts atomic relaxation data (MF=28, MT=533) from an endf.Material.

    Parameters:
        mat: endf.Material - loaded ENDF file object

    Returns:
        dict: Dictionary containing:
            - 'subshells': dict of subshell data with binding energies,
                          electron numbers, and transitions
    """
    key = (28, 533)
    if key not in mat.section_data:
        return None

    sec = mat.section_data[key]

    # EADL atomic relaxation data structure
    result = {
        'NSS': sec.get('NSS', 0),  # Number of subshells
        'subshells': {}
    }

    # Process each subshell
    subshells = sec.get('subshells', [])
    for subshell in subshells:
        subi = int(subshell.get('SUBI', 0))
        shell_name = SUBSHELL_DESIGNATORS.get(subi, f"S{subi}")

        shell_data = {
            'designator': subi,
            'name': shell_name,
            'binding_energy_eV': subshell.get('EBI', 0.0),
            'n_electrons': subshell.get('ELN', 0.0),
            'n_transitions': int(subshell.get('NTR', 0)),
            'transitions': []
        }

        # Process transitions
        transitions = subshell.get('transitions', [])
        for trans in transitions:
            subj = int(trans.get('SUBJ', 0))
            subk = int(trans.get('SUBK', 0))

            trans_data = {
                'origin_subshell': SUBSHELL_DESIGNATORS.get(subj, f"S{subj}"),
                'origin_designator': subj,
                'secondary_subshell': SUBSHELL_DESIGNATORS.get(subk, f"S{subk}") if subk > 0 else "radiative",
                'secondary_designator': subk,
                'energy_eV': trans.get('ETR', 0.0),
                'probability': trans.get('FTR', 0.0),
                'is_radiative': subk == 0,  # SUBK=0 means radiative (X-ray)
                'is_auger': subk > 0,  # SUBK>0 means Auger/Coster-Kronig
            }
            shell_data['transitions'].append(trans_data)

        result['subshells'][shell_name] = shell_data

    return result


def save_atomic_element_h5(mat_path, out_dir):
    """
    Load one EADL ENDF file, extract atomic relaxation data, and write to <symbol>.h5.

    Parameters:
        mat_path: Path to the EADL ENDF file
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

    # Extract atomic relaxation data
    relax_data = extract_atomic_relaxation(mat)

    with h5py.File(h5_path, "w") as h5f:
        # Metadata group
        meta = h5f.create_group("metadata")
        meta.create_dataset("Z", data=np.int64(Z))
        meta.create_dataset("Sym", data=sym)
        if sec0:
            meta.create_dataset("ZA", data=sec0["ZA"])
            meta.create_dataset("AWR", data=sec0["AWR"])

        if relax_data is None:
            print(f"No atomic relaxation data found for {sym}")
            return

        # Atomic relaxation group
        relax_grp = h5f.create_group("atomic_relaxation")
        relax_grp.create_dataset("n_subshells", data=relax_data['NSS'])

        # Create subshells group
        subshells_grp = relax_grp.create_group("subshells")

        # Store binding energies summary
        shell_names = []
        binding_energies = []
        n_electrons = []

        for shell_name, shell_data in relax_data['subshells'].items():
            shell_names.append(shell_name)
            binding_energies.append(shell_data['binding_energy_eV'])
            n_electrons.append(shell_data['n_electrons'])

            # Create group for each subshell
            shell_grp = subshells_grp.create_group(shell_name)
            shell_grp.create_dataset("designator", data=shell_data['designator'])
            shell_grp.create_dataset("binding_energy_eV", data=shell_data['binding_energy_eV'])
            shell_grp.create_dataset("n_electrons", data=shell_data['n_electrons'])
            shell_grp.create_dataset("n_transitions", data=shell_data['n_transitions'])

            # Transitions
            if shell_data['transitions']:
                trans_grp = shell_grp.create_group("transitions")

                # Store transition data as arrays
                origin_shells = []
                secondary_shells = []
                energies = []
                probabilities = []
                is_radiative = []

                for trans in shell_data['transitions']:
                    origin_shells.append(trans['origin_designator'])
                    secondary_shells.append(trans['secondary_designator'])
                    energies.append(trans['energy_eV'])
                    probabilities.append(trans['probability'])
                    is_radiative.append(trans['is_radiative'])

                trans_grp.create_dataset("origin_designator", data=np.array(origin_shells, dtype='i4'))
                trans_grp.create_dataset("secondary_designator", data=np.array(secondary_shells, dtype='i4'))
                trans_grp.create_dataset("energy_eV", data=np.array(energies, dtype='f8'))
                trans_grp.create_dataset("probability", data=np.array(probabilities, dtype='f8'))
                trans_grp.create_dataset("is_radiative", data=np.array(is_radiative, dtype='bool'))

        # Summary datasets
        relax_grp.create_dataset("shell_names", data=np.array(shell_names, dtype='S8'))
        relax_grp.create_dataset("binding_energies_eV", data=np.array(binding_energies, dtype='f8'))
        relax_grp.create_dataset("n_electrons", data=np.array(n_electrons, dtype='f8'))

    print(f"Saved! {h5_path}")


def create_atomic_mcdc_file(in_path: str, out_dir: str) -> str:
    """
    Create MCDC-format HDF5 file from raw atomic relaxation data.

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

        # Check if atomic relaxation data exists
        if "/atomic_relaxation" not in src:
            print(f"No atomic relaxation data in {in_path}")
            return None

        n_subshells = int(src["/atomic_relaxation/n_subshells"][()])
        shell_names = [s.decode() for s in src["/atomic_relaxation/shell_names"][:]]
        binding_energies = src["/atomic_relaxation/binding_energies_eV"][:]
        n_electrons_arr = src["/atomic_relaxation/n_electrons"][:]

        # Read transition data for each subshell
        subshell_data = {}
        for shell_name in shell_names:
            shell_path = f"/atomic_relaxation/subshells/{shell_name}"
            if shell_path not in src:
                continue

            shell_data = {
                'binding_energy_eV': float(src[f"{shell_path}/binding_energy_eV"][()]),
                'n_electrons': float(src[f"{shell_path}/n_electrons"][()]),
            }

            trans_path = f"{shell_path}/transitions"
            if trans_path in src:
                shell_data['transitions'] = {
                    'origin_designator': src[f"{trans_path}/origin_designator"][:],
                    'secondary_designator': src[f"{trans_path}/secondary_designator"][:],
                    'energy_eV': src[f"{trans_path}/energy_eV"][:],
                    'probability': src[f"{trans_path}/probability"][:],
                    'is_radiative': src[f"{trans_path}/is_radiative"][:],
                }

            subshell_data[shell_name] = shell_data

    # Create output file
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{Sym}.h5")

    with h5py.File(out_path, "w") as h5f:
        h5f.create_dataset("atomic_weight_ratio", data=AWR)
        h5f.create_dataset("atomic_number", data=Z)
        h5f.create_dataset("element_name", data=Sym)

        # Atomic relaxation group
        relax = h5f.create_group("atomic_relaxation")
        relax.create_dataset("n_subshells", data=n_subshells)
        relax.create_dataset("binding_energies_eV", data=binding_energies)
        relax.create_dataset("n_electrons", data=n_electrons_arr)

        # Subshells group with detailed data
        subs = relax.create_group("subshells")
        for i, shell_name in enumerate(shell_names):
            if shell_name not in subshell_data:
                continue

            g = subs.create_group(shell_name)
            shell = subshell_data[shell_name]
            g.create_dataset("binding_energy_eV", data=shell['binding_energy_eV'])
            g.create_dataset("n_electrons", data=shell['n_electrons'])

            if 'transitions' in shell:
                trans = shell['transitions']

                # Separate radiative and non-radiative transitions
                is_rad = trans['is_radiative']

                # Radiative transitions (X-ray emission)
                rad_mask = is_rad
                if np.any(rad_mask):
                    rad_grp = g.create_group("radiative")
                    rad_grp.create_dataset("origin_designator", data=trans['origin_designator'][rad_mask])
                    rad_grp.create_dataset("energy_eV", data=trans['energy_eV'][rad_mask])
                    rad_grp.create_dataset("probability", data=trans['probability'][rad_mask])

                    # Calculate fluorescence yield
                    fluorescence_yield = np.sum(trans['probability'][rad_mask])
                    rad_grp.create_dataset("fluorescence_yield", data=fluorescence_yield)

                # Non-radiative transitions (Auger/Coster-Kronig)
                auger_mask = ~is_rad
                if np.any(auger_mask):
                    auger_grp = g.create_group("non_radiative")
                    auger_grp.create_dataset("origin_designator", data=trans['origin_designator'][auger_mask])
                    auger_grp.create_dataset("secondary_designator", data=trans['secondary_designator'][auger_mask])
                    auger_grp.create_dataset("energy_eV", data=trans['energy_eV'][auger_mask])
                    auger_grp.create_dataset("probability", data=trans['probability'][auger_mask])

                    # Calculate Auger yield
                    auger_yield = np.sum(trans['probability'][auger_mask])
                    auger_grp.create_dataset("auger_yield", data=auger_yield)

    print(f"Saved! {out_path}")
    return out_path


def get_binding_energies(h5_path: str) -> pd.DataFrame:
    """
    Extract binding energies for all subshells from an EADL HDF5 file.

    Parameters:
        h5_path: Path to EADL HDF5 file

    Returns:
        DataFrame with columns: shell, binding_energy_eV, n_electrons
    """
    with h5py.File(h5_path, "r") as f:
        if "/atomic_relaxation" not in f:
            return pd.DataFrame()

        shell_names = [s.decode() for s in f["/atomic_relaxation/shell_names"][:]]
        binding_energies = f["/atomic_relaxation/binding_energies_eV"][:]
        n_electrons = f["/atomic_relaxation/n_electrons"][:]

    return pd.DataFrame({
        'shell': shell_names,
        'binding_energy_eV': binding_energies,
        'n_electrons': n_electrons
    })


def get_fluorescence_yield(h5_path: str, shell: str) -> float:
    """
    Get fluorescence yield for a specific subshell.

    Parameters:
        h5_path: Path to EADL MCDC HDF5 file
        shell: Subshell name (e.g., 'K', 'L1', 'M1')

    Returns:
        Fluorescence yield (probability of radiative transition)
    """
    with h5py.File(h5_path, "r") as f:
        path = f"/atomic_relaxation/subshells/{shell}/radiative/fluorescence_yield"
        if path in f:
            return float(f[path][()])
    return 0.0
