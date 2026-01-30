#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
# Maintainer: @melekderman
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

"""
Download EPICS (Electron Photon Interaction Cross Sections) data from IAEA.

This script downloads three libraries:
- EEDL: Evaluated Electron Data Library
- EPDL: Evaluated Photon Data Library
- EADL: Evaluated Atomic Data Library
"""

import os
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup

# EPICS data libraries configuration
LIBRARIES = {
    "eedl": {
        "url": "https://www-nds.iaea.org/epics/ENDF2023/EEDL.ELEMENTS/getza.htm",
        "prefix": "EEDL",
        "description": "Evaluated Electron Data Library"
    },
    "epdl": {
        "url": "https://www-nds.iaea.org/epics/ENDF2023/EPDL.ELEMENTS/getza.htm",
        "prefix": "EPDL",
        "description": "Evaluated Photon Data Library"
    },
    "eadl": {
        "url": "https://www-nds.iaea.org/epics/ENDF2023/EADL.ELEMENTS/getza.htm",
        "prefix": "EADL",
        "description": "Evaluated Atomic Data Library"
    }
}


def download_library(library_name: str, out_dir: str = None):
    """
    Download a specific EPICS library from IAEA.

    Parameters:
        library_name: One of 'eedl', 'epdl', or 'eadl'
        out_dir: Output directory (defaults to library name)
    """
    if library_name not in LIBRARIES:
        raise ValueError(f"Unknown library: {library_name}. Choose from: {list(LIBRARIES.keys())}")

    config = LIBRARIES[library_name]
    base_url = config["url"]
    prefix = config["prefix"]

    if out_dir is None:
        out_dir = library_name

    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Downloading {config['description']} ({prefix})")
    print(f"{'='*60}")

    # Fetch the main page containing element links
    resp = requests.get(base_url)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")

    # Find individual elements and download them
    for a in soup.find_all("a", href=True):
        href = a["href"]
        url = urljoin(base_url, href)
        dst = os.path.join(out_dir, f"{prefix}.{os.path.basename(href)}.endf")

        print(f"Downloading {url} → {dst}")
        with requests.get(url, stream=True) as r, open(dst, "wb") as f:
            r.raise_for_status()
            for chunk in r.iter_content(8192):
                f.write(chunk)

    print(f"Completed downloading {prefix} to {out_dir}/")


def download_all():
    """Download all EPICS libraries (EEDL, EPDL, EADL)."""
    for library_name in LIBRARIES:
        download_library(library_name)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Download EPICS data libraries from IAEA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python download_data.py              # Download all libraries
    python download_data.py eedl         # Download only EEDL (electron data)
    python download_data.py epdl         # Download only EPDL (photon data)
    python download_data.py eadl         # Download only EADL (atomic data)
    python download_data.py eedl epdl    # Download EEDL and EPDL
        """
    )
    parser.add_argument(
        "libraries",
        nargs="*",
        choices=["eedl", "epdl", "eadl"],
        default=[],
        help="Libraries to download (default: all)"
    )

    args = parser.parse_args()

    if not args.libraries:
        download_all()
    else:
        for lib in args.libraries:
            download_library(lib)
