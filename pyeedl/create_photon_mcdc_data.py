#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Script to create MCDC-format HDF5 files from raw photon data.
"""

import os
import glob
from photon import create_photon_mcdc_file


def main():
    base_dir = os.path.dirname(__file__)
    raw_dir = os.path.join(base_dir, "../raw_data_photon")
    out_dir = os.path.join(base_dir, "../mcdc_data_photon")
    os.makedirs(out_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))
    for p in files:
        create_photon_mcdc_file(p, out_dir)


if __name__ == "__main__":
    main()
