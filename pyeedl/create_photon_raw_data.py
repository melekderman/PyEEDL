#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Script to create raw HDF5 files from EPDL (photon) ENDF data.
"""

import os
from photon import save_photon_element_h5


def main():
    base_dir = os.path.dirname(__file__)
    in_dir = os.path.join(base_dir, "../epdl")
    out_dir = os.path.join(base_dir, "../raw_data_photon")
    os.makedirs(out_dir, exist_ok=True)

    for Z in range(1, 101):
        fn = os.path.join(in_dir, f"EPDL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_photon_element_h5(fn, out_dir)


if __name__ == "__main__":
    main()
