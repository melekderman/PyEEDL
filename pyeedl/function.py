#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------


def float_endf(s: str) -> float:
    """Convert an ENDF to Python float."""
    t = s.replace('D','E').strip()

    # insert missing exponent marker if needed (e.g. "1.2345-03")
    if t and ('+' in t[1:] or '-' in t[1:]) and 'E' not in t.upper():
        # find final + or - (beyond first char)
        idx = max(t.rfind('+',1), t.rfind('-',1))
        if idx > 0:
            t = t[:idx] + 'E' + t[idx:]
    return float(t) if t else 0.0


def int_endf(s: str) -> int:
    """Convert ENDF int to Python int."""
    t = s.strip()
    return int(t) if t.isdigit() else 0


def parse_mf26_mt525(raw: str):
    """
    Parse MF=26, MT=525 data from line 9 of the ENDF block.
    Returns list of dicts with E_in, E_out, NW, NL, and ang-PDF pairs.
    """
    lines = raw.splitlines()[8:]  # drop header lines 1–8
    groups = []
    i = 0
    while i < len(lines):
        # --- CONT line ---
        cont = lines[i]
        E_in  = float_endf(cont[ 0:11])
        E_out = float_endf(cont[11:22])
        NW = int_endf(cont[44:55])
        NL = int_endf(cont[55:66])
        i += 1

        # --- read NW mu–p pairs ---
        pairs = []
        while len(pairs) < NL and i < len(lines):
            ln = lines[i][:66]  # only cols 0–65 contain data
            vals = [float_endf(ln[j:j+11]) for j in range(0, 66, 11)]
            for j in range(0, len(vals), 2):
                if len(pairs) < NL:
                    pairs.append((vals[j], vals[j+1]))
            i += 1

        groups.append({
            'E_loss':  E_in,
            'E_in': E_out,
            'NW':    NW,
            'NL':    NL,
            'pairs': pairs
        })

    return groups

