#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------
import numpy as np
from data import FINE_STRUCTURE, ELECTRON_MASS


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
    lines = raw.splitlines()[6:]  # drop header lines 1–6
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

        # --- read NW mu–p pairs (NW stores 2*NL numeric entries) ---
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


def linear_interpolation(xs_energy_grid, energy_ref, xs_ref):
    xs_new = np.interp(np.array(xs_energy_grid), np.array(energy_ref), np.array(xs_ref))
    if np.any (xs_energy_grid > energy_ref[-1]):
        print("Warning: xs_energy_grid has values larger than energy grid max.")
    return xs_new


def build_pdf(inc_energy, value, probability):
    inc = np.array(inc_energy)
    val = np.array(value)
    pdf = np.array(probability)

    energy_vals = []
    offsets = [0]

    if inc.size:
        cur = inc[0]
        cnt = 1
        for e in inc[1:]:
            if e == cur:
                cnt += 1
            else:
                energy_vals.append(cur)
                offsets.append(offsets[-1] + cnt)
                cur = e
                cnt = 1
        energy_vals.append(cur)
        offsets.append(offsets[-1] + cnt)
    else:
        offsets = [0]

    energy_grid = np.asarray(energy_vals, dtype="f8")
    energy_offset = np.asarray(offsets[:-1], dtype="i8")

    return energy_grid, energy_offset, val, pdf


def small_angle_eta(Z, energy_eV):
    alpha = FINE_STRUCTURE
    mec2 = ELECTRON_MASS                              # MeV
    T  = np.array(energy_eV, dtype="f8") / 1e6        # MeV
    pc  = np.sqrt(T * (T + 2.0*mec2))                 # MeV
    E  = T + mec2                                     # MeV
    beta = pc / E
    tau  = T / mec2
    term = (alpha * mec2 / (0.885 * pc))**2
    corr = 1.13 + 3.76 * (alpha * Z / beta)**2
    return 0.25 * term * (Z**(2.0/3.0)) * corr * np.sqrt(tau/(tau+1.0))


def small_angle_scattering_cosine2(Z, energy_eV, n_mu):
    energy_grid = np.array(energy_eV, dtype="f8").ravel()

    mu = np.linspace(0.999999, 1.0, n_mu, dtype="f8")
    eta = small_angle_eta(Z, energy_grid)
    M = mu.size
    value = np.tile(mu, energy_grid.size)
    PDF = np.empty(energy_grid.size * M, dtype="f8")
    for i, et in enumerate(eta):
        f = 1.0 / (et + (1.0 - mu))**2
        s = f.sum()
        PDF[i*M:(i+1)*M] = f / s if s > 0 else 1.0 / M
    energy_offset = np.arange(0, (energy_grid.size + 1) * M, M, dtype="i8")
    return energy_grid, energy_offset[:-1], value, PDF

def small_angle_scattering_cosine(Z, energy_eV, n_mu):
    energy_grid = np.array(energy_eV, dtype="f8").ravel()
    if energy_grid.size == 0:
        return np.array([], dtype="f8"), np.array([0], dtype="i8"), np.array([], dtype="f8"), np.array([], dtype="f8")

    mu = np.linspace(0.999999, 1.0, n_mu, endpoint=False, dtype="f8")
    eta = small_angle_eta(Z, energy_grid)

    N, M = energy_grid.size, mu.size
    value = np.empty(N*M, dtype="f8")
    PDF   = np.empty(N*M, dtype="f8")

    dmu = np.diff(mu)
    widths = np.empty(M, dtype="f8")
    if dmu.size:
        widths[:-1] = dmu
        widths[-1] = dmu[-1]
    else:
        widths[:] = 1.0

    for i, et in enumerate(eta):
        s = slice(i*M, (i+1)*M)
        value[s] = mu
        f = 1.0 / (et + (1.0 - mu))**2
        denom = np.dot(f, widths)
        PDF[s] = f / denom if denom > 0 else 0.0

    energy_offset = np.arange(0, (N+1)*M, M, dtype="i8")
    return energy_grid, energy_offset, value, PDF

def build_coupled_scattering_cosine(
    inc_energy,
    mu_arr,
    prob_arr,
    Z,
    mu_cut=0.999999,
    n_mu_tail=200,
    tail_x_floor=1.0e-12,
):
    """
    Build a coupled elastic angular distribution by stitching the tabulated
    large-angle body to an analytical Screened Rutherford tail and then
    renormalizing the combined PDF.

    The tabulated part is used for mu <= mu_cut. The tail is scaled to match
    the tabulated PDF at mu_cut, following the continuity-based stitching used
    in earlier coupled-elastic experiments.
    """
    energy_grid = np.unique(np.asarray(inc_energy, dtype="f8"))

    out_inc = []
    out_mu = []
    out_prob = []

    x_cut = max(1.0 - mu_cut, tail_x_floor)

    for E in energy_grid:
        mask = inc_energy == E
        mu_body = np.asarray(mu_arr[mask], dtype="f8")
        pdf_body = np.asarray(prob_arr[mask], dtype="f8")

        order = np.argsort(mu_body)
        mu_body = mu_body[order]
        pdf_body = pdf_body[order]

        body_mask = mu_body <= mu_cut
        mu_body = mu_body[body_mask]
        pdf_body = pdf_body[body_mask]

        if mu_body.size == 0:
            mu_body = np.array([-1.0, mu_cut], dtype="f8")
            pdf_body = np.array([0.5, 0.5], dtype="f8")

        if mu_body[-1] < mu_cut:
            pdf_cut = np.interp(mu_cut, mu_body, pdf_body)
            mu_body = np.append(mu_body, mu_cut)
            pdf_body = np.append(pdf_body, pdf_cut)

        body_area = np.trapezoid(pdf_body, mu_body)
        if body_area > 0.0:
            pdf_body = pdf_body / body_area

        x_floor = min(tail_x_floor, 0.1 * x_cut)
        x_floor = max(x_floor, 1.0e-15)
        x_tail = np.geomspace(x_cut, x_floor, max(2, int(n_mu_tail)), dtype="f8")
        mu_tail = 1.0 - x_tail

        eta = float(small_angle_eta(Z, [E])[0])
        tail_shape = 1.0 / (eta + (1.0 - mu_tail)) ** 2

        cutoff_pdf = float(np.interp(mu_cut, mu_body, pdf_body))
        cutoff_shape = 1.0 / (eta + (1.0 - mu_cut)) ** 2
        scale = cutoff_pdf / cutoff_shape if cutoff_shape > 0.0 else 1.0
        pdf_tail = tail_shape * scale

        mu_combined = np.concatenate((mu_body, mu_tail[1:]))
        pdf_combined = np.concatenate((pdf_body, pdf_tail[1:]))

        area = np.trapezoid(pdf_combined, mu_combined)
        if area > 0.0:
            pdf_combined = pdf_combined / area

        out_inc.extend([E] * len(mu_combined))
        out_mu.extend(mu_combined)
        out_prob.extend(pdf_combined)

    return (
        np.array(out_inc, dtype="f8"),
        np.array(out_mu, dtype="f8"),
        np.array(out_prob, dtype="f8"),
    )

# Helpers to densify angular grid
def densify_angular_grid(inc_energy, mu_arr, prob_arr, max_gap_ratio=1.1):
    """
    Fill large gaps in the angular-distribution energy grid by log-interpolating
    the PDF tables.  FRENSIE V&V (Kersting et al., NSE 2019) showed log-log
    grid policies best match experimental results for EEDL data.
    Any two adjacent energies with E_hi / E_lo > *max_gap_ratio*
    get extra points inserted at geometric midpoints until the ratio criterion is met.

    Returns new (inc_energy, mu, probability) arrays with the extra tables spliced in.
    """
    # Group data by incident energy
    unique_E = np.unique(inc_energy)
    tables = {}
    for E in unique_E:
        mask = inc_energy == E
        mu = mu_arr[mask]
        prob = prob_arr[mask]
        order = np.argsort(mu)
        tables[E] = (mu[order], prob[order])

    # Find gaps and insert interpolated tables
    sorted_E = np.sort(list(tables.keys()))
    new_tables = dict(tables)  # keep originals

    for i in range(len(sorted_E) - 1):
        E_lo, E_hi = sorted_E[i], sorted_E[i + 1]
        ratio = E_hi / E_lo
        if ratio <= max_gap_ratio:
            continue

        # How many midpoints needed
        n_insert = int(np.ceil(np.log(ratio) / np.log(max_gap_ratio))) - 1
        if n_insert < 1:
            n_insert = 1

        log_lo, log_hi = np.log(E_lo), np.log(E_hi)
        insert_energies = np.exp(np.linspace(log_lo, log_hi, n_insert + 2)[1:-1])

        mu_lo, pdf_lo = tables[E_lo]
        mu_hi, pdf_hi = tables[E_hi]

        # Common mu grid (union of both tables' mu grids)
        mu_common = np.union1d(mu_lo, mu_hi)
        pdf_lo_interp = np.interp(mu_common, mu_lo, pdf_lo)
        pdf_hi_interp = np.interp(mu_common, mu_hi, pdf_hi)

        # Use log-interpolation on PDF (avoids negative values)
        log_pdf_lo = np.log(np.maximum(pdf_lo_interp, 1e-30))
        log_pdf_hi = np.log(np.maximum(pdf_hi_interp, 1e-30))

        for E_new in insert_energies:
            f = (np.log(E_new) - log_lo) / (log_hi - log_lo)
            log_pdf_new = log_pdf_lo + f * (log_pdf_hi - log_pdf_lo)
            pdf_new = np.exp(log_pdf_new)
            # Normalize
            area = np.trapezoid(pdf_new, mu_common)
            if area > 0:
                pdf_new /= area
            new_tables[E_new] = (mu_common.copy(), pdf_new)

    # Rebuild flat arrays sorted by energy
    all_E = np.sort(list(new_tables.keys()))
    out_inc, out_mu, out_prob = [], [], []
    for E in all_E:
        mu, prob = new_tables[E]
        out_inc.extend([E] * len(mu))
        out_mu.extend(mu)
        out_prob.extend(prob)

    return np.array(out_inc), np.array(out_mu), np.array(out_prob)


def _build_cdf(mu_orig, pdf_orig):
    dmu = np.diff(mu_orig)
    cdf = np.zeros(len(mu_orig), dtype="f8")
    for i in range(len(dmu)):
        cdf[i + 1] = cdf[i] + 0.5 * (pdf_orig[i] + pdf_orig[i + 1]) * dmu[i]

    if cdf[-1] > 0:
        cdf /= cdf[-1]
    return cdf


def _recover_pdf(mu_common, cdf_common):
    n_mu = len(mu_common)
    pdf = np.zeros(n_mu, dtype="f8")
    dmu_common = np.diff(mu_common)
    for i in range(n_mu - 1):
        pdf[i] = (cdf_common[i + 1] - cdf_common[i]) / dmu_common[i]
    if n_mu > 1:
        pdf[-1] = pdf[-2]

    area = np.trapezoid(pdf, mu_common)
    if area > 0:
        pdf /= area
    return pdf


def _unify_mu_grid_on_support(inc_energy, mu_arr, prob_arr, mu_common):
    unique_E = np.unique(inc_energy)
    mu_common = np.asarray(mu_common, dtype="f8")

    out_inc = []
    out_mu = []
    out_prob = []

    for E in unique_E:
        mask = inc_energy == E
        mu_orig = mu_arr[mask]
        pdf_orig = prob_arr[mask]
        order = np.argsort(mu_orig)
        mu_orig = mu_orig[order]
        pdf_orig = pdf_orig[order]

        cdf_orig = _build_cdf(mu_orig, pdf_orig)
        cdf_new = np.interp(mu_common, mu_orig, cdf_orig)
        pdf_new = _recover_pdf(mu_common, cdf_new)

        out_inc.extend([E] * len(mu_common))
        out_mu.extend(mu_common)
        out_prob.extend(pdf_new)

    return np.array(out_inc), np.array(out_mu), np.array(out_prob)


def unify_mu_grid(inc_energy, mu_arr, prob_arr, n_mu=200):
    """
    Re-interpolate all angular distribution tables onto a common mu grid
    using CDF-based interpolation (preserves distribution shape better than
    direct PDF interpolation).

    Steps per energy table:
      1. Build CDF from original (mu, PDF) via trapezoidal integration
      2. Interpolate CDF onto common mu grid
      3. Recover PDF as finite-difference derivative of interpolated CDF
      4. Normalize

    Parameters
    ----------
    inc_energy, mu_arr, prob_arr : flat arrays (output of densify_angular_grid)
    n_mu : int
        Number of equally-spaced mu points in [-1, 0.999999].

    Returns new (inc_energy, mu, probability) flat arrays with uniform table sizes.
    """
    mu_common = np.linspace(-1.0, 0.999999, n_mu, dtype="f8")
    return _unify_mu_grid_on_support(inc_energy, mu_arr, prob_arr, mu_common)


def unify_mu_grid_log1m(inc_energy, mu_arr, prob_arr, n_mu=200, delta=1.0e-10):
    """
    Re-interpolate all angular tables onto a common support that is uniform in
    log((1 - mu) + delta), following the elastic-scattering change of variables
    used for forward-peaked EEDL interpolation.

    This keeps far more resolution near mu -> 1 than a uniform-mu grid.
    """
    mu_min = float(np.min(mu_arr))
    mu_max = float(np.max(mu_arr))

    x_max = (1.0 - mu_min) + delta
    x_min = (1.0 - mu_max) + delta
    x_common = np.geomspace(x_max, x_min, n_mu, dtype="f8")

    mu_common = 1.0 + delta - x_common
    mu_common[0] = mu_min
    mu_common[-1] = mu_max

    return _unify_mu_grid_on_support(inc_energy, mu_arr, prob_arr, mu_common)

# GFP2 Functions


def _legendre_p1(mu):
    return mu


def _legendre_p2(mu):
    return 0.5 * (3.0 * mu * mu - 1.0)


def compute_legendre_moments_tabulated(mu_tab, pdf_tab):
    mu = np.asarray(mu_tab, dtype="f8")
    pdf = np.asarray(pdf_tab, dtype="f8")

    if mu.size < 2:
        return np.array([1.0, 0.0, 0.0], dtype="f8")

    order = np.argsort(mu)
    mu = mu[order]
    pdf = pdf[order]

    norm = np.trapezoid(pdf, mu)
    if norm <= 0.0:
        return np.array([1.0, 0.0, 0.0], dtype="f8")

    pdf = pdf / norm

    g = np.empty(3, dtype="f8")
    g[0] = 1.0
    g[1] = np.trapezoid(_legendre_p1(mu) * pdf, mu)
    g[2] = np.trapezoid(_legendre_p2(mu) * pdf, mu)
    return g


def compute_legendre_moments_screened_rutherford(eta, n_quad=200):
    eta = max(float(eta), 1.0e-30)

    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    f_vals = 1.0 / (eta + 1.0 - nodes) ** 2
    norm = 2.0 / (eta * (eta + 2.0))
    pdf = f_vals / norm

    g = np.empty(3, dtype="f8")
    g[0] = 1.0
    g[1] = np.dot(weights, _legendre_p1(nodes) * pdf)
    g[2] = np.dot(weights, _legendre_p2(nodes) * pdf)
    return g


def _precompute_legendre_anchors(inc_energy, mu_arr, prob_arr):
    energy_grid = np.unique(np.asarray(inc_energy, dtype="f8"))
    anchors = np.empty((energy_grid.size, 3), dtype="f8")

    for i, E in enumerate(energy_grid):
        mask = inc_energy == E
        anchors[i, :] = compute_legendre_moments_tabulated(
            np.asarray(mu_arr[mask], dtype="f8"),
            np.asarray(prob_arr[mask], dtype="f8"),
        )

    return energy_grid, anchors


def _interpolate_legendre_anchors(energy_eV, anchor_energy, anchor_moments):
    M = anchor_energy.size

    if M == 0:
        return np.array([1.0, 0.0, 0.0], dtype="f8")

    if M == 1 or energy_eV <= anchor_energy[0]:
        return anchor_moments[0, :].copy()

    if energy_eV >= anchor_energy[-1]:
        return anchor_moments[-1, :].copy()

    j = np.searchsorted(anchor_energy, energy_eV, side="right") - 1
    j = max(0, min(j, M - 2))

    e0 = anchor_energy[j]
    e1 = anchor_energy[j + 1]
    f = (energy_eV - e0) / (e1 - e0)
    return (1.0 - f) * anchor_moments[j, :] + f * anchor_moments[j + 1, :]


def compute_gfp2_params(
    Sigma_s0,
    Sigma_s1,
    Sigma_s2,
    mu_star=0.9,
    beta_max=100.0,
):
    Sigma_a1 = float(Sigma_s0 - Sigma_s1)
    Sigma_a2 = float(Sigma_s0 - Sigma_s2)
    Sigma_tr = Sigma_a1

    if Sigma_s0 <= 0.0:
        return {
            "alpha": 0.0,
            "beta": 0.0,
            "beta_raw": 0.0,
            "mu_star": float(mu_star),
            "Sigma_delta0": 0.0,
            "transition_rate": 0.0,
            "Sigma_tr": 0.0,
            "Sigma_a1": 0.0,
            "Sigma_a2": 0.0,
            "regime": "boltzmann",
            "success": False,
            "warning": "Sigma_s0 <= 0",
        }

    denom = Sigma_a2 - Sigma_a1
    if abs(denom) < 1.0e-30 * max(float(Sigma_s0), 1.0):
        beta_raw = np.inf
    else:
        beta_raw = (3.0 * Sigma_a1 - Sigma_a2) / (6.0 * denom)

    if beta_raw < 0.0 or Sigma_a1 <= 0.0:
        if Sigma_a1 > 0.0:
            ratio = Sigma_a2 / Sigma_a1
            warning = (
                "GFP2 inadmissible: require Sigma_a1 < Sigma_a2 < 3*Sigma_a1, "
                f"got Sigma_a2/Sigma_a1 = {ratio:.6e}"
            )
        else:
            warning = "GFP2 inadmissible: Sigma_a1 <= 0"

        return {
            "alpha": 0.0,
            "beta": float(beta_raw),
            "beta_raw": float(beta_raw),
            "mu_star": float(mu_star),
            "Sigma_delta0": 0.0,
            "transition_rate": 0.0,
            "Sigma_tr": Sigma_tr,
            "Sigma_a1": Sigma_a1,
            "Sigma_a2": Sigma_a2,
            "regime": "boltzmann",
            "success": False,
            "warning": warning,
        }

    if np.isfinite(beta_raw) and beta_raw <= beta_max:
        regime = "gfp2"
        beta_val = float(beta_raw)
        warning = ""
    else:
        regime = "gfp2_capped"
        beta_val = float(beta_max)
        warning = f"beta_raw = {beta_raw:.6e} > beta_max = {beta_max:.6e}; capped"

    one_minus_mu_star = 1.0 - float(mu_star)
    if one_minus_mu_star <= 0.0:
        warning = (warning + "; " if warning else "") + "mu_star >= 1"
        one_minus_mu_star = 1.0e-12

    alpha_val = 0.5 * Sigma_tr * (1.0 + 2.0 * beta_val)
    Sigma_delta0 = Sigma_tr * (1.0 + 2.0 * beta_val) / one_minus_mu_star
    if beta_val <= 0.0:
        transition_rate = np.inf
    else:
        transition_rate = Sigma_tr * (1.0 + 2.0 * beta_val) / (2.0 * beta_val)

    return {
        "alpha": alpha_val,
        "beta": beta_val,
        "beta_raw": float(beta_raw),
        "mu_star": float(mu_star),
        "Sigma_delta0": Sigma_delta0,
        "transition_rate": transition_rate,
        "Sigma_tr": Sigma_tr,
        "Sigma_a1": Sigma_a1,
        "Sigma_a2": Sigma_a2,
        "regime": regime,
        "success": True,
        "warning": warning,
    }


def build_gfp2_elastic_data(
    xs_energy_grid,
    xs_total,
    xs_large,
    inc_energy,
    mu_arr,
    prob_arr,
    Z,
    mu_star=0.9,
    beta_max=100.0,
    dcs_source="eedl",
):
    if dcs_source not in {"eedl", "sr"}:
        raise ValueError("dcs_source must be 'eedl' or 'sr'")

    xs_energy_grid = np.asarray(xs_energy_grid, dtype="f8")
    xs_total = np.asarray(xs_total, dtype="f8")
    xs_large = np.asarray(xs_large, dtype="f8")
    inc_energy = np.asarray(inc_energy, dtype="f8")
    mu_arr = np.asarray(mu_arr, dtype="f8")
    prob_arr = np.asarray(prob_arr, dtype="f8")

    N = xs_energy_grid.size
    Sigma_s0 = np.empty(N, dtype="f8")
    Sigma_s1 = np.empty(N, dtype="f8")
    Sigma_s2 = np.empty(N, dtype="f8")
    alpha = np.empty(N, dtype="f8")
    beta = np.empty(N, dtype="f8")
    beta_raw = np.empty(N, dtype="f8")
    Sigma_delta0 = np.empty(N, dtype="f8")
    transition_rate = np.empty(N, dtype="f8")
    Sigma_tr = np.empty(N, dtype="f8")
    success = np.empty(N, dtype=bool)
    regime = np.empty(N, dtype="U12")

    if dcs_source == "eedl" and inc_energy.size > 0:
        anchor_energy, anchor_moments = _precompute_legendre_anchors(
            inc_energy, mu_arr, prob_arr
        )
    else:
        anchor_energy = np.array([], dtype="f8")
        anchor_moments = np.empty((0, 3), dtype="f8")

    for i, E in enumerate(xs_energy_grid):
        sigma_el = float(xs_total[i])
        sigma_la = float(xs_large[i])
        sigma_sa = max(sigma_el - sigma_la, 0.0)

        if sigma_el <= 0.0:
            Sigma_s0[i] = 0.0
            Sigma_s1[i] = 0.0
            Sigma_s2[i] = 0.0
            alpha[i] = 0.0
            beta[i] = 0.0
            beta_raw[i] = 0.0
            Sigma_delta0[i] = 0.0
            transition_rate[i] = 0.0
            Sigma_tr[i] = 0.0
            success[i] = False
            regime[i] = "boltzmann"
            continue

        if dcs_source == "sr":
            eta = float(small_angle_eta(Z, [E])[0])
            g_full = compute_legendre_moments_screened_rutherford(eta)
        else:
            if sigma_la > 0.0 and anchor_energy.size > 0:
                g_la = _interpolate_legendre_anchors(E, anchor_energy, anchor_moments)
            else:
                g_la = np.array([1.0, 0.0, 0.0], dtype="f8")

            if sigma_sa > 0.0:
                eta = float(small_angle_eta(Z, [E])[0])
                g_sa = compute_legendre_moments_screened_rutherford(eta)
            else:
                g_sa = np.array([1.0, 0.0, 0.0], dtype="f8")

            g_full = (sigma_la * g_la + sigma_sa * g_sa) / sigma_el

        g_full[0] = 1.0
        Sigma_sl = sigma_el * g_full

        res = compute_gfp2_params(
            Sigma_sl[0],
            Sigma_sl[1],
            Sigma_sl[2],
            mu_star=mu_star,
            beta_max=beta_max,
        )

        Sigma_s0[i] = Sigma_sl[0]
        Sigma_s1[i] = Sigma_sl[1]
        Sigma_s2[i] = Sigma_sl[2]
        alpha[i] = res["alpha"]
        beta[i] = res["beta"]
        beta_raw[i] = res["beta_raw"]
        Sigma_delta0[i] = res["Sigma_delta0"]
        transition_rate[i] = res["transition_rate"]
        Sigma_tr[i] = res["Sigma_tr"]
        success[i] = res["success"]
        regime[i] = res["regime"]

    return {
        "energy_grid": xs_energy_grid,
        "regime": regime,
        "Sigma_s0": Sigma_s0,
        "Sigma_s1": Sigma_s1,
        "Sigma_s2": Sigma_s2,
        "alpha": alpha,
        "beta": beta,
        "beta_raw": beta_raw,
        "Sigma_delta0": Sigma_delta0,
        "transition_rate": transition_rate,
        "Sigma_tr": Sigma_tr,
        "success": success,
        "mu_star": float(mu_star),
        "beta_max": float(beta_max),
        "dcs_source": dcs_source,
        "n_angular_anchors": int(anchor_energy.size),
    }

def densify_angular_grid(inc_energy, mu_arr, prob_arr, max_gap_ratio=1.1):
    """
    Fill large gaps in the angular-distribution energy grid by log-interpolating
    the PDF tables.  FRENSIE V&V (Kersting et al., NSE 2019) showed log-log
    grid policies best match experimental results for EEDL data.
    Any two adjacent energies with E_hi / E_lo > *max_gap_ratio*
    get extra points inserted at geometric midpoints until the ratio criterion is met.

    Returns new (inc_energy, mu, probability) arrays with the extra tables spliced in.
    """
    # Group data by incident energy
    unique_E = np.unique(inc_energy)
    tables = {}
    for E in unique_E:
        mask = inc_energy == E
        mu = mu_arr[mask]
        prob = prob_arr[mask]
        order = np.argsort(mu)
        tables[E] = (mu[order], prob[order])

    # Find gaps and insert interpolated tables
    sorted_E = np.sort(list(tables.keys()))
    new_tables = dict(tables)  # keep originals

    for i in range(len(sorted_E) - 1):
        E_lo, E_hi = sorted_E[i], sorted_E[i + 1]
        ratio = E_hi / E_lo
        if ratio <= max_gap_ratio:
            continue

        # How many midpoints needed
        n_insert = int(np.ceil(np.log(ratio) / np.log(max_gap_ratio))) - 1
        if n_insert < 1:
            n_insert = 1

        log_lo, log_hi = np.log(E_lo), np.log(E_hi)
        insert_energies = np.exp(np.linspace(log_lo, log_hi, n_insert + 2)[1:-1])

        mu_lo, pdf_lo = tables[E_lo]
        mu_hi, pdf_hi = tables[E_hi]

        # Common mu grid (union of both tables' mu grids)
        mu_common = np.union1d(mu_lo, mu_hi)
        pdf_lo_interp = np.interp(mu_common, mu_lo, pdf_lo)
        pdf_hi_interp = np.interp(mu_common, mu_hi, pdf_hi)

        # Use log-interpolation on PDF (avoids negative values)
        log_pdf_lo = np.log(np.maximum(pdf_lo_interp, 1e-30))
        log_pdf_hi = np.log(np.maximum(pdf_hi_interp, 1e-30))

        for E_new in insert_energies:
            f = (np.log(E_new) - log_lo) / (log_hi - log_lo)
            log_pdf_new = log_pdf_lo + f * (log_pdf_hi - log_pdf_lo)
            pdf_new = np.exp(log_pdf_new)
            # Normalize
            area = np.trapezoid(pdf_new, mu_common)
            if area > 0:
                pdf_new /= area
            new_tables[E_new] = (mu_common.copy(), pdf_new)

    # Rebuild flat arrays sorted by energy
    all_E = np.sort(list(new_tables.keys()))
    out_inc, out_mu, out_prob = [], [], []
    for E in all_E:
        mu, prob = new_tables[E]
        out_inc.extend([E] * len(mu))
        out_mu.extend(mu)
        out_prob.extend(prob)

    return np.array(out_inc), np.array(out_mu), np.array(out_prob)


def _build_cdf(mu_orig, pdf_orig):
    dmu = np.diff(mu_orig)
    cdf = np.zeros(len(mu_orig), dtype="f8")
    for i in range(len(dmu)):
        cdf[i + 1] = cdf[i] + 0.5 * (pdf_orig[i] + pdf_orig[i + 1]) * dmu[i]

    if cdf[-1] > 0:
        cdf /= cdf[-1]
    return cdf


def _recover_pdf(mu_common, cdf_common):
    n_mu = len(mu_common)
    pdf = np.zeros(n_mu, dtype="f8")
    dmu_common = np.diff(mu_common)
    for i in range(n_mu - 1):
        pdf[i] = (cdf_common[i + 1] - cdf_common[i]) / dmu_common[i]
    if n_mu > 1:
        pdf[-1] = pdf[-2]

    area = np.trapezoid(pdf, mu_common)
    if area > 0:
        pdf /= area
    return pdf


def _unify_mu_grid_on_support(inc_energy, mu_arr, prob_arr, mu_common):
    unique_E = np.unique(inc_energy)
    mu_common = np.asarray(mu_common, dtype="f8")

    out_inc = []
    out_mu = []
    out_prob = []

    for E in unique_E:
        mask = inc_energy == E
        mu_orig = mu_arr[mask]
        pdf_orig = prob_arr[mask]
        order = np.argsort(mu_orig)
        mu_orig = mu_orig[order]
        pdf_orig = pdf_orig[order]

        cdf_orig = _build_cdf(mu_orig, pdf_orig)
        cdf_new = np.interp(mu_common, mu_orig, cdf_orig)
        pdf_new = _recover_pdf(mu_common, cdf_new)

        out_inc.extend([E] * len(mu_common))
        out_mu.extend(mu_common)
        out_prob.extend(pdf_new)

    return np.array(out_inc), np.array(out_mu), np.array(out_prob)


def unify_mu_grid(inc_energy, mu_arr, prob_arr, n_mu=200):
    """
    Re-interpolate all angular distribution tables onto a common mu grid
    using CDF-based interpolation (preserves distribution shape better than
    direct PDF interpolation).

    Steps per energy table:
      1. Build CDF from original (mu, PDF) via trapezoidal integration
      2. Interpolate CDF onto common mu grid
      3. Recover PDF as finite-difference derivative of interpolated CDF
      4. Normalize

    Parameters
    ----------
    inc_energy, mu_arr, prob_arr : flat arrays (output of densify_angular_grid)
    n_mu : int
        Number of equally-spaced mu points in [-1, 0.999999].

    Returns new (inc_energy, mu, probability) flat arrays with uniform table sizes.
    """
    mu_common = np.linspace(-1.0, 0.999999, n_mu, dtype="f8")
    return _unify_mu_grid_on_support(inc_energy, mu_arr, prob_arr, mu_common)


def unify_mu_grid_log1m(inc_energy, mu_arr, prob_arr, n_mu=200, delta=1.0e-10):
    """
    Re-interpolate all angular tables onto a common support that is uniform in
    log((1 - mu) + delta), following the elastic-scattering change of variables
    used for forward-peaked EEDL interpolation.

    This keeps far more resolution near mu -> 1 than a uniform-mu grid.
    """
    mu_min = float(np.min(mu_arr))
    mu_max = float(np.max(mu_arr))

    x_max = (1.0 - mu_min) + delta
    x_min = (1.0 - mu_max) + delta
    x_common = np.geomspace(x_max, x_min, n_mu, dtype="f8")

    mu_common = 1.0 + delta - x_common
    mu_common[0] = mu_min
    mu_common[-1] = mu_max

    return _unify_mu_grid_on_support(inc_energy, mu_arr, prob_arr, mu_common)
