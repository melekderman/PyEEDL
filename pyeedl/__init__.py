#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
PyEEDL - Python library for reading EPICS (Electron Photon Interaction Cross Sections) data.

This package provides tools to download, parse, and process:
- EEDL: Evaluated Electron Data Library
- EPDL: Evaluated Photon Data Library
- EADL: Evaluated Atomic Data Library

All data is sourced from IAEA Nuclear Data Services in ENDF format.
"""

__version__ = "0.2.0"
__author__ = "CEMeNT"

# Electron data (EEDL)
from .electron import (
    extract_sections,
    save_element_h5,
    create_mcdc_file,
)

# Photon data (EPDL)
from .photon import (
    extract_photon_sections,
    save_photon_element_h5,
    create_photon_mcdc_file,
)

# Atomic relaxation data (EADL)
from .atomic import (
    extract_atomic_relaxation,
    save_atomic_element_h5,
    create_atomic_mcdc_file,
    get_binding_energies,
    get_fluorescence_yield,
)

# Data constants and mappings
from .data import (
    # Electron
    MF_MT,
    SECTIONS_ABBREVS,
    MF23,
    MF26,
    SUBSHELL_LABELS,
    # Photon
    PHOTON_MF_MT,
    PHOTON_SECTIONS_ABBREVS,
    MF27,
    # Atomic
    ATOMIC_MF_MT,
    ATOMIC_SECTIONS_ABBREVS,
    MF28,
    SUBSHELL_DESIGNATORS,
    SUBSHELL_DESIGNATORS_INV,
    # Common
    PERIODIC_TABLE,
    # Constants
    FINE_STRUCTURE,
    ELECTRON_MASS,
    BARN_TO_CM2,
    PLANCK_CONSTANT,
    SPEED_OF_LIGHT,
    ELECTRON_CHARGE,
)

# Utility functions
from .function import (
    float_endf,
    int_endf,
    linear_interpolation,
    build_pdf,
    small_angle_eta,
    small_angle_scattering_cosine,
)

__all__ = [
    # Version
    "__version__",
    # Electron functions
    "extract_sections",
    "save_element_h5",
    "create_mcdc_file",
    # Photon functions
    "extract_photon_sections",
    "save_photon_element_h5",
    "create_photon_mcdc_file",
    # Atomic functions
    "extract_atomic_relaxation",
    "save_atomic_element_h5",
    "create_atomic_mcdc_file",
    "get_binding_energies",
    "get_fluorescence_yield",
    # Data mappings - Electron
    "MF_MT",
    "SECTIONS_ABBREVS",
    "MF23",
    "MF26",
    "SUBSHELL_LABELS",
    # Data mappings - Photon
    "PHOTON_MF_MT",
    "PHOTON_SECTIONS_ABBREVS",
    "MF27",
    # Data mappings - Atomic
    "ATOMIC_MF_MT",
    "ATOMIC_SECTIONS_ABBREVS",
    "MF28",
    "SUBSHELL_DESIGNATORS",
    "SUBSHELL_DESIGNATORS_INV",
    # Common
    "PERIODIC_TABLE",
    # Constants
    "FINE_STRUCTURE",
    "ELECTRON_MASS",
    "BARN_TO_CM2",
    "PLANCK_CONSTANT",
    "SPEED_OF_LIGHT",
    "ELECTRON_CHARGE",
    # Utility functions
    "float_endf",
    "int_endf",
    "linear_interpolation",
    "build_pdf",
    "small_angle_eta",
    "small_angle_scattering_cosine",
]
