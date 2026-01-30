#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Master script to create all EPICS data files (electron, photon, atomic).

This script runs the complete pipeline for all three EPICS libraries:
1. EEDL (Electron Data Library)
2. EPDL (Photon Data Library)
3. EADL (Atomic Data Library)

Usage:
    python create_all_data.py              # Process all libraries
    python create_all_data.py electron     # Process only electron data
    python create_all_data.py photon       # Process only photon data
    python create_all_data.py atomic       # Process only atomic data
"""

import os
import glob
import argparse


def process_electron_data():
    """Process EEDL electron data."""
    from electron import save_element_h5, create_mcdc_file

    base_dir = os.path.dirname(__file__)
    in_dir = os.path.join(base_dir, "../eedl")
    raw_dir = os.path.join(base_dir, "../raw_data")
    mcdc_dir = os.path.join(base_dir, "../mcdc_data")

    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(mcdc_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print("Processing EEDL (Electron Data Library)")
    print("=" * 60)

    # Create raw data
    print("\nStep 1: Creating raw HDF5 files...")
    for Z in range(1, 101):
        fn = os.path.join(in_dir, f"EEDL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_element_h5(fn, raw_dir)

    # Create MCDC data
    print("\nStep 2: Creating MCDC HDF5 files...")
    files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))
    for p in files:
        create_mcdc_file(p, mcdc_dir)


def process_photon_data():
    """Process EPDL photon data."""
    from photon import save_photon_element_h5, create_photon_mcdc_file

    base_dir = os.path.dirname(__file__)
    in_dir = os.path.join(base_dir, "../epdl")
    raw_dir = os.path.join(base_dir, "../raw_data_photon")
    mcdc_dir = os.path.join(base_dir, "../mcdc_data_photon")

    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(mcdc_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print("Processing EPDL (Photon Data Library)")
    print("=" * 60)

    # Create raw data
    print("\nStep 1: Creating raw HDF5 files...")
    for Z in range(1, 101):
        fn = os.path.join(in_dir, f"EPDL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_photon_element_h5(fn, raw_dir)

    # Create MCDC data
    print("\nStep 2: Creating MCDC HDF5 files...")
    files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))
    for p in files:
        create_photon_mcdc_file(p, mcdc_dir)


def process_atomic_data():
    """Process EADL atomic relaxation data."""
    from atomic import save_atomic_element_h5, create_atomic_mcdc_file

    base_dir = os.path.dirname(__file__)
    in_dir = os.path.join(base_dir, "../eadl")
    raw_dir = os.path.join(base_dir, "../raw_data_atomic")
    mcdc_dir = os.path.join(base_dir, "../mcdc_data_atomic")

    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(mcdc_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print("Processing EADL (Atomic Data Library)")
    print("=" * 60)

    # Create raw data
    print("\nStep 1: Creating raw HDF5 files...")
    for Z in range(1, 101):
        fn = os.path.join(in_dir, f"EADL.ZA{Z:03d}000.endf")
        if os.path.exists(fn):
            save_atomic_element_h5(fn, raw_dir)

    # Create MCDC data
    print("\nStep 2: Creating MCDC HDF5 files...")
    files = sorted(glob.glob(os.path.join(raw_dir, "*.h5")))
    for p in files:
        create_atomic_mcdc_file(p, mcdc_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Process EPICS data libraries",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python create_all_data.py              # Process all libraries
    python create_all_data.py electron     # Process only EEDL
    python create_all_data.py photon       # Process only EPDL
    python create_all_data.py atomic       # Process only EADL
    python create_all_data.py electron photon  # Process EEDL and EPDL
        """
    )
    parser.add_argument(
        "libraries",
        nargs="*",
        choices=["electron", "photon", "atomic"],
        default=[],
        help="Libraries to process (default: all)"
    )

    args = parser.parse_args()

    processors = {
        "electron": process_electron_data,
        "photon": process_photon_data,
        "atomic": process_atomic_data,
    }

    if not args.libraries:
        # Process all
        for name, processor in processors.items():
            processor()
    else:
        for lib in args.libraries:
            processors[lib]()

    print("\n" + "=" * 60)
    print("Processing complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
