#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Script to create raw HDF5 files from EADL (atomic relaxation) ENDF data.
"""

import os
from atomic import save_atomic_element_h5


def main():
    base_dir = os.path.dirname(__file__)
    in_dir = os.path.join(base_dir, "../eadl")
    out_dir = os.path.join(base_dir, "../raw_data_atomic")
    os.makedirs(out_dir, exist_ok=True)

    for Z in range(1, 101):
        fn = os.path.join(in_dir, f"EADL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_atomic_element_h5(fn, out_dir)


if __name__ == "__main__":
    main()
