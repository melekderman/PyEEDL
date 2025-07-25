#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 Melek Derman (@melekderman)
#
# SPDX-License-Identifier: MIT
#
# Released under the MIT License:
#   https://opensource.org/licenses/MIT
#
# Acknowledgments:
#   This work was supported by the Center for Exascale Monte Carlo Neutron
#   Transport (CEMeNT), a PSAAP‑III project funded by the U.S. Department of
#   Energy under grant DE‑NA003967.

# Note:
#   All external code snippets and referenced sources used in this code are
#   properly cited at their respective locations, and license information for
#   any incorporated third‑party code has been included. If you notice any
#   missing references or citations, please contact me at
#   derman1@oregonstate.edu.
# -----------------------------------------------------------------------------

import glob
import pandas as pd
import endf

import re
from typing import Tuple

__all__ = [
    "MF_MT",
    "SECTIONS_ABBREVS",
    "MF26",
    "PERIODIC_TABLE",
]

MF_MT = {
    # MF=23  : Electron Cross Sections
    (23, 501): "Total Electron Cross Sections",
    (23, 522): "Ionization (sum of subshells)",
    (23, 525): "Large Angle Elastic Scattering Cross Section",
    (23, 526): "Elastic Scatter (Total) Cross Sections",
    (23, 527): "Bremsstrahlung Cross Sections",
    (23, 528): "Excitation Cross Sections",
    (23, 534): "K (1S1/2) Electroionization Subshell Cross Sections",
    (23, 535): "L1 (2s1/2) Electroionization Subshell Cross Sections",
    (23, 536): "L2 (2p1/2) Electroionization Subshell Cross Sections",
    (23, 537): "L3 (2p3/2) Electroionization Subshell Cross Sections",
    (23, 538): "M1 (3s1/2) Electroionization Subshell Cross Sections",
    (23, 539): "M2 (3p1/2) Electroionization Subshell Cross Sections",
    (23, 540): "M3 (3p3/2) Electroionization Subshell Cross Sections",
    (23, 541): "M4 (3d3/2) Electroionization Subshell Cross Sections",
    (23, 542): "M5 (3d5/2) Electroionization Subshell Cross Sections",
    (23, 543): "N1 (4s1/2) Electroionization Subshell Cross Sections",
    (23, 544): "N2 (4p1/2) Electroionization Subshell Cross Sections",
    (23, 545): "N3 (4p3/2) Electroionization Subshell Cross Sections",
    (23, 546): "N4 (4d3/2) Electroionization Subshell Cross Sections",
    (23, 547): "N5 (4d5/2) Electroionization Subshell Cross Sections",
    (23, 548): "N6 (4f5/2) Electroionization Subshell Cross Sections",
    (23, 549): "N7 (4f7/2) Electroionization Subshell Cross Sections",
    (23, 550): "O1 (5s1/2) Electroionization Subshell Cross Sections",
    (23, 551): "O2 (5p1/2) Electroionization Subshell Cross Sections",
    (23, 552): "O3 (5p3/2) Electroionization Subshell Cross Sections",
    (23, 553): "O4 (5d3/2) Electroionization Subshell Cross Sections",
    (23, 554): "O5 (5d5/2) Electroionization Subshell Cross Sections",
    (23, 555): "O6 (5f5/2) Electroionization Subshell Cross Sections",
    (23, 556): "O7 (5f7/2) Electroionization Subshell Cross Sections",
    (23, 559): "P1 (6s1/2) Electroionization Subshell Cross Sections",
    (23, 560): "P2 (6p1/2) Electroionization Subshell Cross Sections",
    (23, 561): "P3 (6p3/2) Electroionization Subshell Cross Sections",
    (23, 570): "Q1 (7s1/2) Electroionization Subshell Cross Sections",

    # MF=26  : Angular and Energy Distributions
    (26, 525): "Large Angle Elastic Angular Distributions",
    (26, 527): "Bremsstrahlung Photon Energy Spectra and Electron Average Energy Loss",
    (26, 528): "Excitation Electron Average Energy Loss",
    (26, 534): "K (1S1/2) Electroionization Subshell Energy Spectra",
    (26, 535): "L1 (2s1/2) Electroionization Subshell Energy Spectra",
    (26, 536): "L2 (2p1/2) Electroionization Subshell Energy Spectra",
    (26, 537): "L3 (2p3/2) Electroionization Subshell Energy Spectra",
    (26, 538): "M1 (3s1/2) Electroionization Subshell Energy Spectra",
    (26, 539): "M2 (3p1/2) Electroionization Subshell Energy Spectra",
    (26, 540): "M3 (3p3/2) Electroionization Subshell Energy Spectra",
    (26, 541): "M4 (3d3/2) Electroionization Subshell Energy Spectra",
    (26, 542): "M5 (3d5/2) Electroionization Subshell Energy Spectra",
    (26, 543): "N1 (4s1/2) Electroionization Subshell Energy Spectra",
    (26, 544): "N2 (4p1/2) Electroionization Subshell Energy Spectra",
    (26, 545): "N3 (4p3/2) Electroionization Subshell Energy Spectra",
    (26, 546): "N4 (4d3/2) Electroionization Subshell Energy Spectra",
    (26, 547): "N5 (4d5/2) Electroionization Subshell Energy Spectra",
    (26, 548): "N6 (4f5/2) Electroionization Subshell Energy Spectra",
    (26, 549): "N7 (4f7/2) Electroionization Subshell Energy Spectra",
    (26, 550): "O1 (5s1/2) Electroionization Subshell Energy Spectra",
    (26, 551): "O2 (5p1/2) Electroionization Subshell Energy Spectra",
    (26, 552): "O3 (5p3/2) Electroionization Subshell Energy Spectra",
    (26, 553): "O4 (5d3/2) Electroionization Subshell Energy Spectra",
    (26, 554): "O5 (5d5/2) Electroionization Subshell Energy Spectra",
    (26, 555): "O6 (5f5/2) Electroionization Subshell Energy Spectra",
    (26, 556): "O7 (5f7/2) Electroionization Subshell Energy Spectra",
    (26, 559): "P1 (6s1/2) Electroionization Subshell Energy Spectra",
    (26, 560): "P2 (6p1/2) Electroionization Subshell Energy Spectra",
    (26, 561): "P3 (6p3/2) Electroionization Subshell Energy Spectra",
    (26, 570): "Q1 (7s1/2) Electroionization Subshell Energy Spectra",
}

SECTIONS_ABBREVS = {
    # MF=23 : Electron Cross Sections
    (23, 501): "xs_tot",
    (23, 522): "xs_ion",
    (23, 525): "xs_lge",       # large-angle elastic
    (23, 526): "xs_el",        # elastic total
    (23, 527): "xs_brem",
    (23, 528): "xs_exc",
    (23, 534): "xs_K",
    (23, 535): "xs_L1",
    (23, 536): "xs_L2",
    (23, 537): "xs_L3",
    (23, 538): "xs_M1",
    (23, 539): "xs_M2",
    (23, 540): "xs_M3",
    (23, 541): "xs_M4",
    (23, 542): "xs_M5",
    (23, 543): "xs_N1",
    (23, 544): "xs_N2",
    (23, 545): "xs_N3",
    (23, 546): "xs_N4",
    (23, 547): "xs_N5",
    (23, 548): "xs_N6",
    (23, 549): "xs_N7",
    (23, 550): "xs_O1",
    (23, 551): "xs_O2",
    (23, 552): "xs_O3",
    (23, 553): "xs_O4",
    (23, 554): "xs_O5",
    (23, 555): "xs_O6",
    (23, 556): "xs_O7",
    (23, 559): "xs_P1",
    (23, 560): "xs_P2",
    (23, 561): "xs_P3",
    (23, 570): "xs_Q1",

    # MF=26 : Angular and Energy Distributions
    (26, 525): "ang_lge",
    (26, 527): "loss_brem_spec",
    (26, 528): "loss_exc",
    (26, 534): "spec_K",
    (26, 535): "spec_L1",
    (26, 536): "spec_L2",
    (26, 537): "spec_L3",
    (26, 538): "spec_M1",
    (26, 539): "spec_M2",
    (26, 540): "spec_M3",
    (26, 541): "spec_M4",
    (26, 542): "spec_M5",
    (26, 543): "spec_N1",
    (26, 544): "spec_N2",
    (26, 545): "spec_N3",
    (26, 546): "spec_N4",
    (26, 547): "spec_N5",
    (26, 548): "spec_N6",
    (26, 549): "spec_N7",
    (26, 550): "spec_O1",
    (26, 551): "spec_O2",
    (26, 552): "spec_O3",
    (26, 553): "spec_O4",
    (26, 554): "spec_O5",
    (26, 555): "spec_O6",
    (26, 556): "spec_O7",
    (26, 559): "spec_P1",
    (26, 560): "spec_P2",
    (26, 561): "spec_P3",
    (26, 570): "spec_Q1",
}

# Mapping of MF=23 cross-section fields to their ENDF‑6 descriptions
# Refs:
# 1) ENDF‑6 Formats Manual (ENDF‑102, 2023), Section 23.2 (“File 23: Electron Interaction Data”)
# 2) ENDF‑6 Formats Manual (ENDF‑102, 2023), Section 2.2.1 (“TAB1: Tabulated Functions”)
MF23 = {
    "sigma.x":              "Array of incident-energy grid points (eV)",
    "sigma.y":              "Array of cross-section values (e.g. barns)",
    "sigma.breakpoints":    "NBT: number of points in each interpolation region",
    "sigma.interpolation":  "INT: interpolation law code for each region",
    "ZA":                   "ZA identifier of the target (atomic number ×1000 + mass number)",
    "AWR":                  "Atomic weight ratio of the target",
    # optional:
    "LRF":                  "Resonance/interpolation flag (if present)",
}

# Mapping of MF=26 distribution fields to their ENDF‑6 descriptions
# Refs:
# 1) ENDF‑6 Formats Manual (ENDF‑102, 2023), Section 26.2 (“File 26: Secondary Energy Distributions”)
# 2) ENDF‑6 Formats Manual (ENDF‑102, 2023), Table 6.2 (“TAB2: Energy–Energy Distributions”)
# 3) https://t2.lanl.gov/nis/endf/law1for6.html

MF26 = {
    "ZAP":            "Product identifier (e.g. 11 for electrons, 0 for photons)",  
    "AWI":            "Atomic weight ratio of the incident particle",  
    "LAW":            "Representation law (1=continuum energy, 2=angular, 8=energy transfer)",  
    "y":              "Yield y(E) vs incident energy (Tabulated1D, typically unity)",  

    # Nested distribution parameters
    "distribution.LANG":   "Indicator that selects the angular representation to be used: (e.g. 1=Legendre, 2=Kalbach)",  
    "distribution.LEP":    "Interpolation scheme for secondary energy: (1=histogram, 2=linear-linear)",  
    "distribution.NR":     "Number of interpolation regions in incident-energy axis",  
    "distribution.NE":     "Number of incident-energy points",  
    "distribution.E_int":  "Tabulated2D full energy–energy distribution object",  
    "distribution.E":      "Array of incident-energy grid points (eV)",  

    # Per-incident-energy entries (one dict per E)
    "distribution.distribution.ND":   "Number of discrete outgoing-energy points",  
    "distribution.distribution.NA":   "Number of angular parameters (0=Isotropic, 1=Kalbach (with LANG=2))",  
    "distribution.distribution.NW":   "Total words in the LIST record (NEP * (NA+2))",  
    "distribution.distribution.NEP":  "Number of secondary energy points in the distribution.",  
    "distribution.distribution.E'":   "Array of outgoing energies E′ (eV)",  
    "distribution.distribution.b":    "PDF or coefficient array for each outgoing-energy point",  
}

PERIODIC_TABLE = {
    1:   {"name": "Hydrogen",      "symbol": "H"},
    2:   {"name": "Helium",        "symbol": "He"},
    3:   {"name": "Lithium",       "symbol": "Li"},
    4:   {"name": "Beryllium",     "symbol": "Be"},
    5:   {"name": "Boron",         "symbol": "B"},
    6:   {"name": "Carbon",        "symbol": "C"},
    7:   {"name": "Nitrogen",      "symbol": "N"},
    8:   {"name": "Oxygen",        "symbol": "O"},
    9:   {"name": "Fluorine",      "symbol": "F"},
    10:  {"name": "Neon",          "symbol": "Ne"},
    11:  {"name": "Sodium",        "symbol": "Na"},
    12:  {"name": "Magnesium",     "symbol": "Mg"},
    13:  {"name": "Aluminium",     "symbol": "Al"},
    14:  {"name": "Silicon",       "symbol": "Si"},
    15:  {"name": "Phosphorus",    "symbol": "P"},
    16:  {"name": "Sulfur",        "symbol": "S"},
    17:  {"name": "Chlorine",      "symbol": "Cl"},
    18:  {"name": "Argon",         "symbol": "Ar"},
    19:  {"name": "Potassium",     "symbol": "K"},
    20:  {"name": "Calcium",       "symbol": "Ca"},
    21:  {"name": "Scandium",      "symbol": "Sc"},
    22:  {"name": "Titanium",      "symbol": "Ti"},
    23:  {"name": "Vanadium",      "symbol": "V"},
    24:  {"name": "Chromium",      "symbol": "Cr"},
    25:  {"name": "Manganese",     "symbol": "Mn"},
    26:  {"name": "Iron",          "symbol": "Fe"},
    27:  {"name": "Cobalt",        "symbol": "Co"},
    28:  {"name": "Nickel",        "symbol": "Ni"},
    29:  {"name": "Copper",        "symbol": "Cu"},
    30:  {"name": "Zinc",          "symbol": "Zn"},
    31:  {"name": "Gallium",       "symbol": "Ga"},
    32:  {"name": "Germanium",     "symbol": "Ge"},
    33:  {"name": "Arsenic",       "symbol": "As"},
    34:  {"name": "Selenium",      "symbol": "Se"},
    35:  {"name": "Bromine",       "symbol": "Br"},
    36:  {"name": "Krypton",       "symbol": "Kr"},
    37:  {"name": "Rubidium",      "symbol": "Rb"},
    38:  {"name": "Strontium",     "symbol": "Sr"},
    39:  {"name": "Yttrium",       "symbol": "Y"},
    40:  {"name": "Zirconium",     "symbol": "Zr"},
    41:  {"name": "Niobium",       "symbol": "Nb"},
    42:  {"name": "Molybdenum",    "symbol": "Mo"},
    43:  {"name": "Technetium",    "symbol": "Tc"},
    44:  {"name": "Ruthenium",     "symbol": "Ru"},
    45:  {"name": "Rhodium",       "symbol": "Rh"},
    46:  {"name": "Palladium",     "symbol": "Pd"},
    47:  {"name": "Silver",        "symbol": "Ag"},
    48:  {"name": "Cadmium",       "symbol": "Cd"},
    49:  {"name": "Indium",        "symbol": "In"},
    50:  {"name": "Tin",           "symbol": "Sn"},
    51:  {"name": "Antimony",      "symbol": "Sb"},
    52:  {"name": "Tellurium",     "symbol": "Te"},
    53:  {"name": "Iodine",        "symbol": "I"},
    54:  {"name": "Xenon",         "symbol": "Xe"},
    55:  {"name": "Caesium",       "symbol": "Cs"},
    56:  {"name": "Barium",        "symbol": "Ba"},
    57:  {"name": "Lanthanum",     "symbol": "La"},
    58:  {"name": "Cerium",        "symbol": "Ce"},
    59:  {"name": "Praseodymium",  "symbol": "Pr"},
    60:  {"name": "Neodymium",     "symbol": "Nd"},
    61:  {"name": "Promethium",    "symbol": "Pm"},
    62:  {"name": "Samarium",      "symbol": "Sm"},
    63:  {"name": "Europium",      "symbol": "Eu"},
    64:  {"name": "Gadolinium",    "symbol": "Gd"},
    65:  {"name": "Terbium",       "symbol": "Tb"},
    66:  {"name": "Dysprosium",    "symbol": "Dy"},
    67:  {"name": "Holmium",       "symbol": "Ho"},
    68:  {"name": "Erbium",        "symbol": "Er"},
    69:  {"name": "Thulium",       "symbol": "Tm"},
    70:  {"name": "Ytterbium",     "symbol": "Yb"},
    71:  {"name": "Lutetium",      "symbol": "Lu"},
    72:  {"name": "Hafnium",       "symbol": "Hf"},
    73:  {"name": "Tantalum",      "symbol": "Ta"},
    74:  {"name": "Tungsten",      "symbol": "W"},
    75:  {"name": "Rhenium",       "symbol": "Re"},
    76:  {"name": "Osmium",        "symbol": "Os"},
    77:  {"name": "Iridium",       "symbol": "Ir"},
    78:  {"name": "Platinum",      "symbol": "Pt"},
    79:  {"name": "Gold",          "symbol": "Au"},
    80:  {"name": "Mercury",       "symbol": "Hg"},
    81:  {"name": "Thallium",      "symbol": "Tl"},
    82:  {"name": "Lead",          "symbol": "Pb"},
    83:  {"name": "Bismuth",       "symbol": "Bi"},
    84:  {"name": "Polonium",      "symbol": "Po"},
    85:  {"name": "Astatine",      "symbol": "At"},
    86:  {"name": "Radon",         "symbol": "Rn"},
    87:  {"name": "Francium",      "symbol": "Fr"},
    88:  {"name": "Radium",        "symbol": "Ra"},
    89:  {"name": "Actinium",      "symbol": "Ac"},
    90:  {"name": "Thorium",       "symbol": "Th"},
    91:  {"name": "Protactinium",  "symbol": "Pa"},
    92:  {"name": "Uranium",       "symbol": "U"},
    93:  {"name": "Neptunium",     "symbol": "Np"},
    94:  {"name": "Plutonium",     "symbol": "Pu"},
    95:  {"name": "Americium",     "symbol": "Am"},
    96:  {"name": "Curium",        "symbol": "Cm"},
    97:  {"name": "Berkelium",     "symbol": "Bk"},
    98:  {"name": "Californium",   "symbol": "Cf"},
    99:  {"name": "Einsteinium",   "symbol": "Es"},
    100: {"name": "Fermium",       "symbol": "Fm"},
    101: {"name": "Mendelevium",   "symbol": "Md"},
    102: {"name": "Nobelium",      "symbol": "No"},
    103: {"name": "Lawrencium",    "symbol": "Lr"},
    104: {"name": "Rutherfordium", "symbol": "Rf"},
    105: {"name": "Dubnium",       "symbol": "Db"},
    106: {"name": "Seaborgium",    "symbol": "Sg"},
    107: {"name": "Bohrium",       "symbol": "Bh"},
    108: {"name": "Hassium",       "symbol": "Hs"},
    109: {"name": "Meitnerium",    "symbol": "Mt"},
    110: {"name": "Darmstadtium",  "symbol": "Ds"},
    111: {"name": "Roentgenium",   "symbol": "Rg"},
    112: {"name": "Copernicium",   "symbol": "Cn"},
    113: {"name": "Nihonium",      "symbol": "Nh"},
    114: {"name": "Flerovium",     "symbol": "Fl"},
    115: {"name": "Moscovium",     "symbol": "Mc"},
    116: {"name": "Livermorium",   "symbol": "Lv"},
    117: {"name": "Tennessine",    "symbol": "Ts"},
    118: {"name": "Oganesson",     "symbol": "Og"},
}
