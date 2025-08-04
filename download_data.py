#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
# Maintainer: @melekderman
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

import os
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup

BASE_URL = "https://www-nds.iaea.org/epics/ENDF2023/EEDL.ELEMENTS/getza.htm"
OUT_DIR  = "eedl"

os.makedirs(OUT_DIR, exist_ok=True)

# fetch the main page containing EEDL links
resp = requests.get(BASE_URL)
resp.raise_for_status()
soup = BeautifulSoup(resp.text, "html.parser")

# find individual elements linking to EEDL files and download them
for a in soup.find_all("a", href=True):
    href = a["href"]
    url  = urljoin(BASE_URL, href)
    dst  = os.path.join(OUT_DIR, f"EEDL.{os.path.basename(href)}.endf")

    print(f"Downloading {url} → {dst}")
    with requests.get(url, stream=True) as r, open(dst, "wb") as f:
        r.raise_for_status()
        for chunk in r.iter_content(8192):
            f.write(chunk)
