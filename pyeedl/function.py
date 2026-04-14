#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------
import numpy as np
from .data import FINE_STRUCTURE, ELECTRON_MASS


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
    if not t:
        return 0
    try:
        return int(t)
    except ValueError:
        return 0


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


def _legendre_p3(mu):
    return 0.5 * (5.0 * mu * mu * mu - 3.0 * mu)


def compute_legendre_moments_tabulated(mu_tab, pdf_tab, n_moments=4):
    """Compute Legendre moments g_0 .. g_{n_moments-1} from tabulated PDF."""
    mu = np.asarray(mu_tab, dtype="f8")
    pdf = np.asarray(pdf_tab, dtype="f8")

    g = np.zeros(n_moments, dtype="f8")
    g[0] = 1.0

    if mu.size < 2:
        return g

    order = np.argsort(mu)
    mu = mu[order]
    pdf = pdf[order]

    norm = np.trapezoid(pdf, mu)
    if norm <= 0.0:
        return g

    pdf = pdf / norm

    _P = [None, _legendre_p1, _legendre_p2, _legendre_p3]
    for ell in range(1, min(n_moments, len(_P))):
        g[ell] = np.trapezoid(_P[ell](mu) * pdf, mu)
    return g


def compute_legendre_moments_screened_rutherford(eta, n_moments=4, n_quad=200):
    """Compute Legendre moments g_0 .. g_{n_moments-1} for screened Rutherford."""
    eta = max(float(eta), 1.0e-30)

    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    f_vals = 1.0 / (eta + 1.0 - nodes) ** 2
    norm = 2.0 / (eta * (eta + 2.0))
    pdf = f_vals / norm

    g = np.zeros(n_moments, dtype="f8")
    g[0] = 1.0

    _P = [None, _legendre_p1, _legendre_p2, _legendre_p3]
    for ell in range(1, min(n_moments, len(_P))):
        g[ell] = np.dot(weights, _P[ell](nodes) * pdf)
    return g


def _precompute_legendre_anchors(inc_energy, mu_arr, prob_arr, n_moments=4):
    energy_grid = np.unique(np.asarray(inc_energy, dtype="f8"))
    anchors = np.empty((energy_grid.size, n_moments), dtype="f8")

    for i, E in enumerate(energy_grid):
        mask = inc_energy == E
        anchors[i, :] = compute_legendre_moments_tabulated(
            np.asarray(mu_arr[mask], dtype="f8"),
            np.asarray(prob_arr[mask], dtype="f8"),
            n_moments=n_moments,
        )

    return energy_grid, anchors


def _interpolate_legendre_anchors(energy_eV, anchor_energy, anchor_moments):
    n_moments = anchor_moments.shape[1] if anchor_moments.ndim == 2 else 3
    M = anchor_energy.size

    if M == 0:
        g = np.zeros(n_moments, dtype="f8")
        g[0] = 1.0
        return g

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
    mu_star_check = _validate_mu_star(mu_star, context="GFP2 mu_star")

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

    alpha_val = 0.5 * Sigma_tr * (1.0 + 2.0 * beta_val)
    if beta_val <= 0.0:
        transition_rate = np.inf
    else:
        transition_rate = Sigma_tr * (1.0 + 2.0 * beta_val) / (2.0 * beta_val)

    if mu_star_check["success"]:
        one_minus_mu_star = 1.0 - float(mu_star_check["mu_star"])
        Sigma_delta0 = Sigma_tr * (1.0 + 2.0 * beta_val) / one_minus_mu_star
        success = True
    else:
        Sigma_delta0 = np.nan
        success = False
        warning = (
            warning + "; " if warning else ""
        ) + mu_star_check["warning"]

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
        "success": success,
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
    Sigma_s3 = np.empty(N, dtype="f8")
    alpha = np.empty(N, dtype="f8")
    beta = np.empty(N, dtype="f8")
    beta_raw = np.empty(N, dtype="f8")
    Sigma_delta0 = np.empty(N, dtype="f8")
    transition_rate = np.empty(N, dtype="f8")
    Sigma_tr = np.empty(N, dtype="f8")
    success = np.empty(N, dtype=bool)
    regime = np.empty(N, dtype="U12")
    warning = np.empty(N, dtype="U256")

    if dcs_source == "eedl" and inc_energy.size > 0:
        anchor_energy, anchor_moments = _precompute_legendre_anchors(
            inc_energy, mu_arr, prob_arr, n_moments=4
        )
    else:
        anchor_energy = np.array([], dtype="f8")
        anchor_moments = np.empty((0, 4), dtype="f8")

    for i, E in enumerate(xs_energy_grid):
        sigma_el = float(xs_total[i])
        sigma_la = float(xs_large[i])
        sigma_sa = max(sigma_el - sigma_la, 0.0)

        if sigma_el <= 0.0:
            Sigma_s0[i] = 0.0
            Sigma_s1[i] = 0.0
            Sigma_s2[i] = 0.0
            Sigma_s3[i] = 0.0
            alpha[i] = 0.0
            beta[i] = 0.0
            beta_raw[i] = 0.0
            Sigma_delta0[i] = 0.0
            transition_rate[i] = 0.0
            Sigma_tr[i] = 0.0
            success[i] = False
            regime[i] = "boltzmann"
            warning[i] = "sigma_el <= 0"
            continue

        if dcs_source == "sr":
            eta = float(small_angle_eta(Z, [E])[0])
            g_full = compute_legendre_moments_screened_rutherford(eta)
        else:
            if sigma_la > 0.0 and anchor_energy.size > 0:
                g_la = _interpolate_legendre_anchors(E, anchor_energy, anchor_moments)
            else:
                g_la = np.array([1.0, 0.0, 0.0, 0.0], dtype="f8")

            if sigma_sa > 0.0:
                eta = float(small_angle_eta(Z, [E])[0])
                g_sa = compute_legendre_moments_screened_rutherford(eta, n_moments=4)
            else:
                g_sa = np.array([1.0, 0.0, 0.0, 0.0], dtype="f8")

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
        Sigma_s3[i] = Sigma_sl[3] if len(Sigma_sl) > 3 else 0.0
        alpha[i] = res["alpha"]
        beta[i] = res["beta"]
        beta_raw[i] = res["beta_raw"]
        Sigma_delta0[i] = res["Sigma_delta0"]
        transition_rate[i] = res["transition_rate"]
        Sigma_tr[i] = res["Sigma_tr"]
        success[i] = res["success"]
        regime[i] = res["regime"]
        warning[i] = res["warning"]

    return {
        "energy_grid": xs_energy_grid,
        "regime": regime,
        "Sigma_s0": Sigma_s0,
        "Sigma_s1": Sigma_s1,
        "Sigma_s2": Sigma_s2,
        "Sigma_s3": Sigma_s3,
        "alpha": alpha,
        "beta": beta,
        "beta_raw": beta_raw,
        "Sigma_delta0": Sigma_delta0,
        "transition_rate": transition_rate,
        "Sigma_tr": Sigma_tr,
        "success": success,
        "warning": warning,
        "mu_star": float(mu_star),
        "beta_max": float(beta_max),
        "dcs_source": dcs_source,
        "n_angular_anchors": int(anchor_energy.size),
    }


# =============================================================================
# Optimal mu_star Selection  (Prinja, CARRE slides §12)
# =============================================================================


def _mu_star_result(mu_star, success, warning="", mu_star_raw=None):
    raw = mu_star if mu_star_raw is None else mu_star_raw
    raw = float(raw) if np.isfinite(raw) else np.nan
    value = float(mu_star) if success else np.nan
    return {
        "mu_star": value,
        "mu_star_raw": raw,
        "success": bool(success),
        "warning": warning,
    }


def _validate_mu_star(mu_star, context="mu_star"):
    if mu_star is None:
        return _mu_star_result(np.nan, False, f"{context} is required")

    try:
        mu_star_raw = float(mu_star)
    except (TypeError, ValueError):
        return _mu_star_result(np.nan, False, f"{context} is not a real number")

    if not np.isfinite(mu_star_raw):
        return _mu_star_result(
            mu_star_raw, False, f"{context} must be finite", mu_star_raw=mu_star_raw
        )
    if mu_star_raw < -1.0:
        return _mu_star_result(
            mu_star_raw,
            False,
            f"{context} = {mu_star_raw:.6e} is below the Eq.(17) domain [-1, 1]",
            mu_star_raw=mu_star_raw,
        )
    if mu_star_raw >= 1.0:
        return _mu_star_result(
            mu_star_raw,
            False,
            f"{context} = {mu_star_raw:.6e} makes 1-mu_star non-positive in Eqs.(18),(38)",
            mu_star_raw=mu_star_raw,
        )

    return _mu_star_result(mu_star_raw, True, mu_star_raw=mu_star_raw)


def compute_optimal_mu_star_fp(Sigma_a1, Sigma_a2):
    """
    Optimal mu_star for the delta-kernel FP model by matching the ell=2 eigenvalue.

    From Prinja Eq.(22), the delta-kernel eigenvalue at ell=2 is:
        lambda_{delta,2} = -Sigma_tr/(1-mu*) * [1 - P_2(mu*)]

    Using [1 - P_2(mu)]/(1-mu) = 3(1+mu)/2 and setting lambda_{delta,2} = -Sigma_a2:
        mu* = 2*Sigma_a2 / (3*Sigma_a1) - 1

    Admissibility: Sigma_a2 <= 3*Sigma_a1  (same as GFP2).
    """
    if Sigma_a1 <= 0.0:
        return _mu_star_result(
            np.nan, False, "Sigma_a1 <= 0, so Eq.(22) cannot define mu_star"
        )

    mu_star = 2.0 * Sigma_a2 / (3.0 * Sigma_a1) - 1.0
    result = _validate_mu_star(mu_star, context="FP optimal mu_star")
    if result["success"]:
        return result

    result["warning"] = (
        f"FP ell=2 constraint gives mu_star = {mu_star:.6e}, outside the admissible "
        "delta-kernel range"
    )
    return result


def compute_optimal_mu_star_gfp2(alpha, beta, Sigma_a3):
    """
    Optimal mu_star for the delta-kernel GFP2 model by matching the ell=3 eigenvalue.

    From the delta-kernel GFP2 two-state system (Prinja Eqs. 37-40), eliminating
    the state-1 scattering term gives the modal eigenvalue

        lambda_{delta,GFP2,ell} = -2*alpha*x_ell / (1 + 2*beta*x_ell),

    where x_ell = [1 - P_ell(mu*)] / (1 - mu*).

    For ell=3, using  [1-P_3(mu)]/(1-mu) = (5*mu^2 + 5*mu + 2)/2, we set
    the eigenvalue equal to -Sigma_a3.  This yields a closed-form quadratic:
        5*mu*^2 + 5*mu* + (2 - C) = 0,   C = Sigma_a3 / (alpha - beta*Sigma_a3)
        mu* = (-5 + sqrt(20*C - 15)) / 10
    """
    if alpha <= 0.0:
        return _mu_star_result(
            np.nan, False, "alpha <= 0, so the GFP2 ell=3 constraint is undefined"
        )
    if (alpha - beta * Sigma_a3) <= 0.0:
        return _mu_star_result(
            np.nan,
            False,
            "alpha - beta*Sigma_a3 <= 0, so Eq.(ell=3) has no admissible delta-kernel root",
        )
    x3_target = Sigma_a3 / (2.0 * (alpha - beta * Sigma_a3))
    C = 2.0 * x3_target
    disc = 20.0 * C - 15.0
    if disc < 0.0:
        return _mu_star_result(
            np.nan,
            False,
            f"GFP2 ell=3 constraint has no real mu_star root: discriminant={disc:.6e}",
        )
    mu_star = (-5.0 + np.sqrt(disc)) / 10.0
    result = _validate_mu_star(mu_star, context="GFP2 optimal mu_star")
    if result["success"]:
        return result

    result["warning"] = (
        f"GFP2 ell=3 constraint gives mu_star = {mu_star:.6e}, outside the admissible "
        "delta-kernel range"
    )
    return result


# =============================================================================
# Delta-Kernel Fokker-Planck  (Prinja Eqs. 17-18, 25-27)
# =============================================================================


def compute_delta_kernel_fp_params(Sigma_s0, Sigma_s1, mu_star=None):
    """
    Compute delta-kernel Fokker-Planck parameters.

    Prinja Eq.(17):  Sigma_delta(mu0) = Sigma_tr / [2*pi*(1-mu*)] * delta(mu0-mu*)
    Prinja Eq.(18):  Sigma_delta0     = Sigma_tr / (1-mu*)

    A valid mu_star must be supplied explicitly.
    """
    Sigma_a1 = float(Sigma_s0 - Sigma_s1)
    if Sigma_a1 <= 0.0:
        return {
            "Sigma_tr": Sigma_a1,
            "Sigma_delta0": np.nan,
            "mu_star": np.nan,
            "mu_star_raw": np.nan,
            "success": False,
            "warning": "Sigma_a1 <= 0, so the delta-kernel FP model is undefined",
        }

    mu_star_check = _validate_mu_star(mu_star, context="FP mu_star")
    if not mu_star_check["success"]:
        return {
            "Sigma_tr": Sigma_a1,
            "Sigma_delta0": np.nan,
            "mu_star": np.nan,
            "mu_star_raw": mu_star_check["mu_star_raw"],
            "success": False,
            "warning": mu_star_check["warning"],
        }

    oms = 1.0 - float(mu_star_check["mu_star"])
    return {
        "Sigma_tr": Sigma_a1,
        "Sigma_delta0": Sigma_a1 / oms,
        "mu_star": float(mu_star_check["mu_star"]),
        "mu_star_raw": float(mu_star_check["mu_star_raw"]),
        "success": True,
        "warning": "",
    }


def build_delta_kernel_fp_data(
    xs_energy_grid, xs_total, xs_large, inc_energy, mu_arr, prob_arr, Z,
    mu_star_policy="optimal", dcs_source="eedl", mu_star_fixed=None,
):
    """
    Build per-energy delta-kernel FP data with optional optimal mu_star.

    Parameters
    ----------
    mu_star_policy : str
        "optimal" - compute mu* from ell=2 eigenvalue matching per energy point
        "fixed"   - use the explicitly provided mu_star_fixed everywhere
    """
    if mu_star_policy not in {"optimal", "fixed"}:
        raise ValueError("mu_star_policy must be 'optimal' or 'fixed'")

    xs_energy_grid = np.asarray(xs_energy_grid, dtype="f8")
    xs_total = np.asarray(xs_total, dtype="f8")
    xs_large = np.asarray(xs_large, dtype="f8")
    inc_energy = np.asarray(inc_energy, dtype="f8")
    mu_arr = np.asarray(mu_arr, dtype="f8")
    prob_arr = np.asarray(prob_arr, dtype="f8")

    N = xs_energy_grid.size
    Sigma_tr_arr = np.empty(N, dtype="f8")
    Sigma_delta0_arr = np.empty(N, dtype="f8")
    mu_star_arr = np.empty(N, dtype="f8")
    mu_star_raw_arr = np.empty(N, dtype="f8")
    success_arr = np.empty(N, dtype=bool)
    warning_arr = np.empty(N, dtype="U256")

    if dcs_source == "eedl" and inc_energy.size > 0:
        anchor_energy, anchor_moments = _precompute_legendre_anchors(
            inc_energy, mu_arr, prob_arr, n_moments=4
        )
    else:
        anchor_energy = np.array([], dtype="f8")
        anchor_moments = np.empty((0, 4), dtype="f8")

    for i, E in enumerate(xs_energy_grid):
        sigma_el = float(xs_total[i])
        sigma_la = float(xs_large[i])
        sigma_sa = max(sigma_el - sigma_la, 0.0)

        if sigma_el <= 0.0:
            Sigma_tr_arr[i] = 0.0
            Sigma_delta0_arr[i] = np.nan
            mu_star_arr[i] = np.nan
            mu_star_raw_arr[i] = np.nan
            success_arr[i] = False
            warning_arr[i] = "sigma_el <= 0"
            continue

        if dcs_source == "sr":
            eta = float(small_angle_eta(Z, [E])[0])
            g = compute_legendre_moments_screened_rutherford(eta, n_moments=4)
        else:
            if sigma_la > 0.0 and anchor_energy.size > 0:
                g_la = _interpolate_legendre_anchors(E, anchor_energy, anchor_moments)
            else:
                g_la = np.array([1.0, 0.0, 0.0, 0.0], dtype="f8")
            if sigma_sa > 0.0:
                eta = float(small_angle_eta(Z, [E])[0])
                g_sa = compute_legendre_moments_screened_rutherford(eta, n_moments=4)
            else:
                g_sa = np.array([1.0, 0.0, 0.0, 0.0], dtype="f8")
            g = (sigma_la * g_la + sigma_sa * g_sa) / sigma_el

        g[0] = 1.0
        Sigma_sl = sigma_el * g
        Sigma_a1 = Sigma_sl[0] - Sigma_sl[1]
        Sigma_a2 = Sigma_sl[0] - Sigma_sl[2]

        if mu_star_policy == "optimal":
            mu_star_sel = compute_optimal_mu_star_fp(Sigma_a1, Sigma_a2)
        else:
            mu_star_sel = _validate_mu_star(mu_star_fixed, context="fixed FP mu_star")

        if not mu_star_sel["success"]:
            Sigma_tr_arr[i] = Sigma_a1
            Sigma_delta0_arr[i] = np.nan
            mu_star_arr[i] = np.nan
            mu_star_raw_arr[i] = mu_star_sel["mu_star_raw"]
            success_arr[i] = False
            warning_arr[i] = mu_star_sel["warning"]
            continue

        res = compute_delta_kernel_fp_params(
            Sigma_sl[0], Sigma_sl[1], mu_star=mu_star_sel["mu_star"]
        )
        Sigma_tr_arr[i] = res["Sigma_tr"]
        Sigma_delta0_arr[i] = res["Sigma_delta0"]
        mu_star_arr[i] = res["mu_star"]
        mu_star_raw_arr[i] = res["mu_star_raw"]
        success_arr[i] = res["success"]
        warning_arr[i] = res["warning"]

    return {
        "energy_grid": xs_energy_grid,
        "Sigma_tr": Sigma_tr_arr,
        "Sigma_delta0": Sigma_delta0_arr,
        "mu_star": mu_star_arr,
        "mu_star_raw": mu_star_raw_arr,
        "success": success_arr,
        "warning": warning_arr,
    }


# =============================================================================
# GFP3 Parameters  (Prinja §8, extended to 3rd order)
# =============================================================================


def compute_gfp3_params(Sigma_s0, Sigma_s1, Sigma_s2, Sigma_s3, mu_star=0.9):
    """
    Compute GFP3 parameters (alpha, beta1, beta2) by matching ell=1,2,3
    Boltzmann eigenvalues.

    GFP3 operator:  L_{GFP3} = alpha * L * (1 - beta1*L)^{-1} * (1 - beta2*L)^{-1}

    Eigenvalues:
        lambda_{GFP3,ell} = -alpha * ell(ell+1)
                             / [(1 + beta1*ell(ell+1)) * (1 + beta2*ell(ell+1))]

    Using p = beta1+beta2 and q = beta1*beta2, the matching conditions at
    ell=1,2,3 reduce to a 2x2 linear system for (p, q), then alpha is obtained
    from the ell=1 condition.
    """
    Sigma_a1 = float(Sigma_s0 - Sigma_s1)
    Sigma_a2 = float(Sigma_s0 - Sigma_s2)
    Sigma_a3 = float(Sigma_s0 - Sigma_s3)

    fail = {
        "alpha": 0.0, "beta1": 0.0, "beta2": 0.0,
        "Sigma_delta0": np.nan,
        "transition_rate_02": np.nan,
        "transition_rate_21": np.nan,
        "Sigma_a1": Sigma_a1, "Sigma_a2": Sigma_a2, "Sigma_a3": Sigma_a3,
        "mu_star": float(mu_star), "success": False, "warning": "",
    }
    mu_star_check = _validate_mu_star(mu_star, context="GFP3 mu_star")

    if Sigma_a1 <= 0.0:
        fail["warning"] = "Sigma_a1 <= 0"
        return fail

    R2 = Sigma_a2 / Sigma_a1
    R3 = Sigma_a3 / Sigma_a1

    # 2x2 linear system:  a*p + b*q = c
    a11 = 6.0 * (R2 - 1.0)
    a12 = 12.0 * (3.0 * R2 - 1.0)
    b1 = 3.0 - R2
    a21 = 12.0 * (R3 - 1.0)
    a22 = 24.0 * (6.0 * R3 - 1.0)
    b2 = 6.0 - R3

    det = a11 * a22 - a12 * a21
    if abs(det) < 1e-30:
        fail["warning"] = "singular system"
        return fail

    p = (b1 * a22 - b2 * a12) / det
    q = (a11 * b2 - a21 * b1) / det

    disc = p * p - 4.0 * q
    if disc < 0.0:
        fail["warning"] = f"complex roots: disc={disc:.6e}"
        return fail

    sqrt_disc = np.sqrt(disc)
    beta1 = (p + sqrt_disc) / 2.0
    beta2 = (p - sqrt_disc) / 2.0

    if beta1 < 0.0 or beta2 < 0.0:
        fail["warning"] = f"negative beta: beta1={beta1:.6e}, beta2={beta2:.6e}"
        return fail

    alpha = Sigma_a1 * (1.0 + 2.0 * p + 4.0 * q) / 2.0
    if alpha <= 0.0:
        fail["warning"] = f"alpha <= 0: {alpha:.6e}"
        return fail

    if not mu_star_check["success"]:
        fail["alpha"] = alpha
        fail["beta1"] = beta1
        fail["beta2"] = beta2
        fail["warning"] = mu_star_check["warning"]
        fail["Sigma_delta0"] = np.nan
        fail["transition_rate_02"] = np.nan
        fail["transition_rate_21"] = np.nan
        return fail

    # Delta-kernel transport quantities for a 3-state system
    oms = 1.0 - float(mu_star_check["mu_star"])
    Sigma_delta0 = alpha * 2.0 / oms

    # Transition rates between states (GFP3 has state0 <-> state2 <-> state1)
    if beta2 > 0.0:
        transition_rate_02 = 1.0 / beta2   # state-0 <-> state-2
    else:
        transition_rate_02 = np.inf
    if beta1 > 0.0:
        transition_rate_21 = 1.0 / beta1   # state-2 <-> state-1
    else:
        transition_rate_21 = np.inf

    return {
        "alpha": alpha,
        "beta1": beta1,
        "beta2": beta2,
        "p": p,
        "q": q,
        "Sigma_delta0": Sigma_delta0,
        "transition_rate_02": transition_rate_02,
        "transition_rate_21": transition_rate_21,
        "Sigma_a1": Sigma_a1,
        "Sigma_a2": Sigma_a2,
        "Sigma_a3": Sigma_a3,
        "mu_star": float(mu_star_check["mu_star"]),
        "success": True,
        "warning": "",
    }


def build_gfp3_elastic_data(
    xs_energy_grid, xs_total, xs_large,
    inc_energy, mu_arr, prob_arr, Z,
    mu_star=0.9, dcs_source="eedl",
):
    """Build per-energy GFP3 data arrays."""
    xs_energy_grid = np.asarray(xs_energy_grid, dtype="f8")
    xs_total = np.asarray(xs_total, dtype="f8")
    xs_large = np.asarray(xs_large, dtype="f8")
    inc_energy = np.asarray(inc_energy, dtype="f8")
    mu_arr = np.asarray(mu_arr, dtype="f8")
    prob_arr = np.asarray(prob_arr, dtype="f8")

    N = xs_energy_grid.size
    fields = ["alpha", "beta1", "beta2", "Sigma_delta0",
              "transition_rate_02", "transition_rate_21",
              "Sigma_s0", "Sigma_s1", "Sigma_s2", "Sigma_s3"]
    arrays = {f: np.empty(N, dtype="f8") for f in fields}
    success = np.empty(N, dtype=bool)
    warning = np.empty(N, dtype="U256")

    if dcs_source == "eedl" and inc_energy.size > 0:
        anchor_energy, anchor_moments = _precompute_legendre_anchors(
            inc_energy, mu_arr, prob_arr, n_moments=4
        )
    else:
        anchor_energy = np.array([], dtype="f8")
        anchor_moments = np.empty((0, 4), dtype="f8")

    for i, E in enumerate(xs_energy_grid):
        sigma_el = float(xs_total[i])
        sigma_la = float(xs_large[i])
        sigma_sa = max(sigma_el - sigma_la, 0.0)

        if sigma_el <= 0.0:
            for f in fields:
                arrays[f][i] = 0.0
            success[i] = False
            warning[i] = "sigma_el <= 0"
            continue

        if dcs_source == "sr":
            eta = float(small_angle_eta(Z, [E])[0])
            g = compute_legendre_moments_screened_rutherford(eta, n_moments=4)
        else:
            if sigma_la > 0.0 and anchor_energy.size > 0:
                g_la = _interpolate_legendre_anchors(E, anchor_energy, anchor_moments)
            else:
                g_la = np.array([1.0, 0.0, 0.0, 0.0], dtype="f8")
            if sigma_sa > 0.0:
                eta = float(small_angle_eta(Z, [E])[0])
                g_sa = compute_legendre_moments_screened_rutherford(eta, n_moments=4)
            else:
                g_sa = np.array([1.0, 0.0, 0.0, 0.0], dtype="f8")
            g = (sigma_la * g_la + sigma_sa * g_sa) / sigma_el

        g[0] = 1.0
        Sigma_sl = sigma_el * g

        res = compute_gfp3_params(
            Sigma_sl[0], Sigma_sl[1], Sigma_sl[2], Sigma_sl[3], mu_star=mu_star,
        )

        arrays["Sigma_s0"][i] = Sigma_sl[0]
        arrays["Sigma_s1"][i] = Sigma_sl[1]
        arrays["Sigma_s2"][i] = Sigma_sl[2]
        arrays["Sigma_s3"][i] = Sigma_sl[3]
        arrays["alpha"][i] = res["alpha"]
        arrays["beta1"][i] = res["beta1"]
        arrays["beta2"][i] = res["beta2"]
        arrays["Sigma_delta0"][i] = res["Sigma_delta0"]
        arrays["transition_rate_02"][i] = res["transition_rate_02"]
        arrays["transition_rate_21"][i] = res["transition_rate_21"]
        success[i] = res["success"]
        warning[i] = res["warning"]

    arrays["energy_grid"] = xs_energy_grid
    arrays["success"] = success
    arrays["warning"] = warning
    arrays["mu_star"] = float(mu_star)
    arrays["dcs_source"] = dcs_source
    return arrays


# =============================================================================
# mu_star -> 1 Convergence Validation  (Prinja §7, §12)
# =============================================================================


def validate_mu_star_convergence(
    Sigma_s0, Sigma_s1, Sigma_s2,
    mu_star_values=None, model="gfp2", beta_max=100.0,
):
    """
    Verify that delta-kernel eigenvalues converge to continuous FP/GFP2
    eigenvalues as mu_star -> 1.

    Returns a dict of {mu_star: {ell: (delta_eigenvalue, continuous_eigenvalue)}}
    for ell = 1..3.

    The continuous eigenvalues are:
      FP:   lambda_{FP,ell}   = -(Sigma_tr/2) * ell*(ell+1)
      GFP2: lambda_{GFP2,ell} = -alpha * ell*(ell+1) / (1 + beta*ell*(ell+1))

    The delta-kernel GFP2 eigenvalues follow from the two-state system
    (Prinja Eqs. 37-40):
      lambda_{delta,GFP2,ell} = -2*alpha*x_ell / (1 + 2*beta*x_ell),
      x_ell = [1 - P_ell(mu*)] / (1 - mu*).
    """
    if mu_star_values is None:
        mu_star_values = [0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999]

    Sigma_a1 = Sigma_s0 - Sigma_s1
    Sigma_a2 = Sigma_s0 - Sigma_s2
    Sigma_tr = Sigma_a1

    # Continuous eigenvalues
    if model == "fp":
        cont = {}
        for ell in range(1, 4):
            cont[ell] = -Sigma_tr / 2.0 * ell * (ell + 1)
    else:
        res = compute_gfp2_params(Sigma_s0, Sigma_s1, Sigma_s2, beta_max=beta_max)
        alpha_c = res["alpha"]
        beta_c = res["beta"]
        cont = {}
        for ell in range(1, 4):
            ll1 = ell * (ell + 1)
            cont[ell] = -alpha_c * ll1 / (1.0 + beta_c * ll1)

    # Legendre polynomials at mu*
    def P_ell(ell, mu):
        if ell == 0:
            return 1.0
        if ell == 1:
            return mu
        if ell == 2:
            return 0.5 * (3.0 * mu**2 - 1.0)
        if ell == 3:
            return 0.5 * (5.0 * mu**3 - 3.0 * mu)
        return 0.0

    results = {}
    for ms in mu_star_values:
        oms = 1.0 - ms
        if oms <= 0.0:
            oms = 1e-15

        row = {}
        for ell in range(1, 4):
            x_ell = (1.0 - P_ell(ell, ms)) / oms

            if model == "fp":
                lam_delta = -Sigma_tr * x_ell
            else:
                lam_delta = -2.0 * alpha_c * x_ell / (1.0 + 2.0 * beta_c * x_ell)

            row[ell] = {
                "delta": lam_delta,
                "continuous": cont[ell],
                "rel_error": abs(lam_delta - cont[ell]) / max(abs(cont[ell]), 1e-30),
            }
        results[ms] = row

    return results
