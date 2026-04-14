#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025 CEMeNT
#
# SPDX-License-Identifier: MIT
# -----------------------------------------------------------------------------

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

from ..function import (
    build_pdf,
    compute_delta_kernel_fp_params,
    compute_gfp2_params,
    compute_gfp3_params,
    compute_legendre_moments_screened_rutherford,
    compute_optimal_mu_star_fp,
    compute_optimal_mu_star_gfp2,
)


LL_TWO_REGION_CASES = {
    "LL_R1": {
        "sigma_t": 100.0,
        "eta": 1.0e-4,
        "atomic_number": 1,
        "atomic_weight_ratio": 1.0,
    },
    "LL_R2": {
        "sigma_t": 50.0,
        "eta": 1.0e-2,
        "atomic_number": 2,
        "atomic_weight_ratio": 2.0,
    },
}

DEFAULT_ENERGY_GRID_EV = np.array([1.0e5, 2.0e6], dtype="f8")


def _validate_scalar(name, value, *, positive=False):
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if positive and value <= 0.0:
        raise ValueError(f"{name} must be > 0")
    return value


def _validate_energy_grid(energy_grid_eV):
    energy_grid = np.asarray(energy_grid_eV, dtype="f8")
    if energy_grid.ndim != 1:
        raise ValueError("energy_grid_eV must be one-dimensional")
    if energy_grid.size < 2:
        raise ValueError("energy_grid_eV must contain at least two points")
    if not np.all(np.isfinite(energy_grid)):
        raise ValueError("energy_grid_eV must be finite")
    if np.any(np.diff(energy_grid) <= 0.0):
        raise ValueError("energy_grid_eV must be strictly increasing")
    return energy_grid


def _validate_mu_star(mu_star):
    mu_star = float(mu_star)
    if not np.isfinite(mu_star):
        raise ValueError("mu_star must be finite")
    if mu_star < -1.0 or mu_star >= 1.0:
        raise ValueError("mu_star must satisfy -1 <= mu_star < 1")
    return mu_star


def build_log1m_mu_grid(n_mu=1200, mu_min=-1.0, mu_max=1.0, delta=1.0e-12):
    """
    Build a monotone mu grid with dense resolution near mu -> 1.

    The grid is uniform in log((1-mu)+delta), which is much better suited to the
    screened-Rutherford forward peak than a uniform-mu discretization.
    """
    n_mu = int(n_mu)
    if n_mu < 2:
        raise ValueError("n_mu must be at least 2")

    mu_min = float(mu_min)
    mu_max = float(mu_max)
    delta = _validate_scalar("delta", delta, positive=True)

    if mu_min < -1.0 or mu_min >= mu_max or mu_max > 1.0:
        raise ValueError("require -1 <= mu_min < mu_max <= 1")

    x_max = (1.0 - mu_min) + delta
    x_min = (1.0 - mu_max) + delta
    x_grid = np.geomspace(x_max, x_min, n_mu, dtype="f8")

    mu_grid = 1.0 + delta - x_grid
    mu_grid[0] = mu_min
    mu_grid[-1] = mu_max
    return mu_grid


def screened_rutherford_pdf(mu, eta):
    """
    Return the normalized screened-Rutherford PDF on the supplied mu grid.

    The exact shape is proportional to
        1 / (eta + 1 - mu)^2.
    We normalize numerically over the tabulated support so the resulting HDF5
    table is self-consistent even when the support is discretized.
    """
    mu = np.asarray(mu, dtype="f8")
    eta = _validate_scalar("eta", eta, positive=True)

    if mu.ndim != 1 or mu.size < 2:
        raise ValueError("mu must be a one-dimensional array with at least two points")
    if np.any(np.diff(mu) <= 0.0):
        raise ValueError("mu must be strictly increasing")

    pdf = 1.0 / (eta + 1.0 - mu) ** 2
    area = np.trapezoid(pdf, mu)
    if area <= 0.0 or not np.isfinite(area):
        raise ValueError("screened-Rutherford PDF could not be normalized")
    return pdf / area


def _constant_multitable(energy_grid_eV, value_grid, pdf_grid):
    inc_energy = np.repeat(np.asarray(energy_grid_eV, dtype="f8"), value_grid.size)
    value = np.tile(np.asarray(value_grid, dtype="f8"), len(energy_grid_eV))
    pdf = np.tile(np.asarray(pdf_grid, dtype="f8"), len(energy_grid_eV))
    return build_pdf(inc_energy, value, pdf)


def _string_dataset(group, name, text):
    group.create_dataset(name, data=np.bytes_(text))


def _write_optional_zero_reaction_group(parent, name):
    parent.create_group(name)


def _write_gfp2_group(group, energy_grid, sigma_t, g_moments, mu_star, beta_max):
    Sigma_sl = sigma_t * g_moments
    gfp2_res = compute_gfp2_params(
        Sigma_sl[0],
        Sigma_sl[1],
        Sigma_sl[2],
        mu_star=mu_star,
        beta_max=beta_max,
    )

    energy_grid = np.asarray(energy_grid, dtype="f8")
    n_energy = energy_grid.size

    g2 = group.create_group("gfp2")
    g2.create_dataset("mu_star", data=float(mu_star))
    g2.create_dataset("beta_max", data=float(beta_max))
    g2.attrs["dcs_source"] = "synthetic_sr"
    g2.attrs["n_angular_anchors"] = 0
    g2.attrs["mu_grid_policy"] = "log1m_exact"
    g2.create_dataset("energy_grid", data=energy_grid)
    g2.create_dataset(
        "regime",
        data=np.full(n_energy, gfp2_res["regime"], dtype="U12").astype("S12"),
    )
    g2.create_dataset("Sigma_s0", data=np.full(n_energy, Sigma_sl[0], dtype="f8"))
    g2.create_dataset("Sigma_s1", data=np.full(n_energy, Sigma_sl[1], dtype="f8"))
    g2.create_dataset("Sigma_s2", data=np.full(n_energy, Sigma_sl[2], dtype="f8"))
    g2.create_dataset("Sigma_s3", data=np.full(n_energy, Sigma_sl[3], dtype="f8"))
    g2.create_dataset("alpha", data=np.full(n_energy, gfp2_res["alpha"], dtype="f8"))
    g2.create_dataset("beta", data=np.full(n_energy, gfp2_res["beta"], dtype="f8"))
    g2.create_dataset(
        "beta_raw", data=np.full(n_energy, gfp2_res["beta_raw"], dtype="f8")
    )
    g2.create_dataset(
        "Sigma_delta0",
        data=np.full(n_energy, gfp2_res["Sigma_delta0"], dtype="f8"),
    )
    g2.create_dataset(
        "transition_rate",
        data=np.full(n_energy, gfp2_res["transition_rate"], dtype="f8"),
    )
    g2.create_dataset(
        "Sigma_tr", data=np.full(n_energy, gfp2_res["Sigma_tr"], dtype="f8")
    )
    g2.create_dataset(
        "success", data=np.full(n_energy, gfp2_res["success"], dtype=bool)
    )
    g2.create_dataset(
        "warning",
        data=np.full(n_energy, gfp2_res["warning"], dtype="U256").astype("S256"),
    )

    Sigma_a3 = Sigma_sl[0] - Sigma_sl[3]
    mu_opt = compute_optimal_mu_star_gfp2(
        gfp2_res["alpha"], gfp2_res["beta"], Sigma_a3
    )
    g2.create_dataset(
        "mu_star_optimal",
        data=np.full(n_energy, mu_opt["mu_star"], dtype="f8"),
    )
    g2.create_dataset(
        "mu_star_optimal_raw",
        data=np.full(n_energy, mu_opt["mu_star_raw"], dtype="f8"),
    )
    g2.create_dataset(
        "mu_star_optimal_success",
        data=np.full(n_energy, mu_opt["success"], dtype=bool),
    )
    g2.create_dataset(
        "mu_star_optimal_warning",
        data=np.full(n_energy, mu_opt["warning"], dtype="U256").astype("S256"),
    )


def _write_delta_kernel_fp_group(group, energy_grid, sigma_t, g_moments):
    Sigma_sl = sigma_t * g_moments
    Sigma_a1 = Sigma_sl[0] - Sigma_sl[1]
    Sigma_a2 = Sigma_sl[0] - Sigma_sl[2]
    mu_opt = compute_optimal_mu_star_fp(Sigma_a1, Sigma_a2)

    if mu_opt["success"]:
        fp_res = compute_delta_kernel_fp_params(
            Sigma_sl[0], Sigma_sl[1], mu_star=mu_opt["mu_star"]
        )
    else:
        fp_res = {
            "Sigma_tr": Sigma_a1,
            "Sigma_delta0": np.nan,
            "mu_star": np.nan,
            "mu_star_raw": mu_opt["mu_star_raw"],
            "success": False,
            "warning": mu_opt["warning"],
        }

    energy_grid = np.asarray(energy_grid, dtype="f8")
    n_energy = energy_grid.size

    fp = group.create_group("delta_kernel_fp")
    fp.create_dataset("energy_grid", data=energy_grid)
    fp.create_dataset("Sigma_tr", data=np.full(n_energy, fp_res["Sigma_tr"], dtype="f8"))
    fp.create_dataset(
        "Sigma_delta0", data=np.full(n_energy, fp_res["Sigma_delta0"], dtype="f8")
    )
    fp.create_dataset("mu_star", data=np.full(n_energy, fp_res["mu_star"], dtype="f8"))
    fp.create_dataset(
        "mu_star_raw", data=np.full(n_energy, fp_res["mu_star_raw"], dtype="f8")
    )
    fp.create_dataset("success", data=np.full(n_energy, fp_res["success"], dtype=bool))
    fp.create_dataset(
        "warning",
        data=np.full(n_energy, fp_res["warning"], dtype="U256").astype("S256"),
    )


def _write_gfp3_group(group, energy_grid, sigma_t, g_moments, mu_star):
    Sigma_sl = sigma_t * g_moments
    gfp3_res = compute_gfp3_params(
        Sigma_sl[0],
        Sigma_sl[1],
        Sigma_sl[2],
        Sigma_sl[3],
        mu_star=mu_star,
    )

    energy_grid = np.asarray(energy_grid, dtype="f8")
    n_energy = energy_grid.size

    g3 = group.create_group("gfp3")
    g3.create_dataset("energy_grid", data=energy_grid)
    g3.create_dataset("alpha", data=np.full(n_energy, gfp3_res["alpha"], dtype="f8"))
    g3.create_dataset("beta1", data=np.full(n_energy, gfp3_res["beta1"], dtype="f8"))
    g3.create_dataset("beta2", data=np.full(n_energy, gfp3_res["beta2"], dtype="f8"))
    g3.create_dataset(
        "Sigma_delta0", data=np.full(n_energy, gfp3_res["Sigma_delta0"], dtype="f8")
    )
    g3.create_dataset(
        "transition_rate_02",
        data=np.full(n_energy, gfp3_res["transition_rate_02"], dtype="f8"),
    )
    g3.create_dataset(
        "transition_rate_21",
        data=np.full(n_energy, gfp3_res["transition_rate_21"], dtype="f8"),
    )
    g3.create_dataset("Sigma_s0", data=np.full(n_energy, Sigma_sl[0], dtype="f8"))
    g3.create_dataset("Sigma_s1", data=np.full(n_energy, Sigma_sl[1], dtype="f8"))
    g3.create_dataset("Sigma_s2", data=np.full(n_energy, Sigma_sl[2], dtype="f8"))
    g3.create_dataset("Sigma_s3", data=np.full(n_energy, Sigma_sl[3], dtype="f8"))
    g3.create_dataset("success", data=np.full(n_energy, gfp3_res["success"], dtype=bool))
    g3.create_dataset(
        "warning",
        data=np.full(n_energy, gfp3_res["warning"], dtype="U256").astype("S256"),
    )
    g3.create_dataset("mu_star", data=float(mu_star))


def write_screened_rutherford_element(
    out_dir,
    element_name,
    sigma_t,
    eta,
    *,
    mu_star=0.9,
    beta_max=100.0,
    energy_grid_eV=None,
    n_mu=1200,
    mu_min=-1.0,
    mu_max=1.0,
    atomic_number=1,
    atomic_weight_ratio=1.0,
):
    """
    Write a synthetic elastic-only electron library in MCDC HDF5 format.

    The output is a continuous-energy file whose cross sections and angular
    distributions are constant in energy. When used with a monoenergetic source
    and no nonelastic reactions, MCDC transport is effectively one-speed.
    """
    sigma_t = _validate_scalar("sigma_t", sigma_t, positive=True)
    eta = _validate_scalar("eta", eta, positive=True)
    mu_star = _validate_mu_star(mu_star)
    beta_max = _validate_scalar("beta_max", beta_max, positive=True)
    atomic_number = int(atomic_number)
    atomic_weight_ratio = _validate_scalar(
        "atomic_weight_ratio", atomic_weight_ratio, positive=True
    )

    if atomic_number <= 0:
        raise ValueError("atomic_number must be a positive integer")

    if energy_grid_eV is None:
        energy_grid = DEFAULT_ENERGY_GRID_EV.copy()
    else:
        energy_grid = _validate_energy_grid(energy_grid_eV)

    mu_grid = build_log1m_mu_grid(n_mu=n_mu, mu_min=mu_min, mu_max=mu_max)
    sr_pdf = screened_rutherford_pdf(mu_grid, eta)
    ang_energy, ang_offset, ang_value, ang_pdf = _constant_multitable(
        energy_grid, mu_grid, sr_pdf
    )

    g_moments = compute_legendre_moments_screened_rutherford(eta, n_moments=4)
    xs_total = np.full(energy_grid.size, sigma_t, dtype="f8")

    out_path = Path(out_dir).expanduser().resolve() / f"{element_name}.h5"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_path, "w") as h5f:
        h5f.create_dataset("atomic_weight_ratio", data=atomic_weight_ratio)
        h5f.create_dataset("atomic_number", data=atomic_number)
        _string_dataset(h5f, "element_name", str(element_name))

        er = h5f.create_group("electron_reactions")
        er.create_dataset("xs_energy_grid", data=energy_grid)

        total_group = er.create_group("total")
        total_group.create_dataset("xs", data=xs_total)

        elastic_root = er.create_group("elastic_scattering")
        elastic = elastic_root.create_group("MT525")
        elastic.attrs["MT"] = 525
        elastic.attrs["synthetic_case"] = "screened_rutherford_one_speed"
        elastic.attrs["synthetic_eta"] = eta
        elastic.attrs["synthetic_sigma_t"] = sigma_t

        xs_ds = elastic.create_dataset("xs", data=xs_total)
        xs_ds.attrs["offset"] = 0
        _string_dataset(elastic, "reference_frame", "LAB")
        elastic.create_dataset("xs_energy", data=energy_grid)
        # For the synthetic exact benchmark we tabulate the full elastic PDF,
        # so xs_large equals the total elastic cross section and no analytical
        # small-angle branch is needed at runtime.
        elastic.create_dataset("xs_large", data=xs_total)

        sc = elastic.create_group("scattering_cosine")
        sc.attrs["mu_grid_policy"] = "log1m_exact"
        sc.create_dataset("energy_grid", data=ang_energy)
        sc.create_dataset("energy_offset", data=ang_offset)
        sc.create_dataset("value", data=ang_value)
        sc.create_dataset("PDF", data=ang_pdf)

        sc_c = elastic.create_group("scattering_cosine_coupled")
        sc_c.attrs["mu_grid_policy"] = "log1m_exact"
        sc_c.create_dataset("energy_grid", data=ang_energy)
        sc_c.create_dataset("energy_offset", data=ang_offset)
        sc_c.create_dataset("value", data=ang_value)
        sc_c.create_dataset("PDF", data=ang_pdf)

        _write_gfp2_group(elastic, energy_grid, sigma_t, g_moments, mu_star, beta_max)
        _write_delta_kernel_fp_group(elastic, energy_grid, sigma_t, g_moments)
        _write_gfp3_group(elastic, energy_grid, sigma_t, g_moments, mu_star)

        # Keep the remaining sections explicitly empty so the file layout still
        # mirrors the MCDC-readable element structure without introducing any
        # energy-loss physics.
        _write_optional_zero_reaction_group(er, "excitation")
        _write_optional_zero_reaction_group(er, "bremsstrahlung")
        _write_optional_zero_reaction_group(er, "ionization")

    return str(out_path)


def write_ll_two_region_library(
    out_dir,
    *,
    mu_star=0.9,
    beta_max=100.0,
    energy_grid_eV=None,
    n_mu=1200,
):
    """
    Write the two-region Leakeas-Larsen one-speed screened-Rutherford library.
    """
    paths = {}
    for name, spec in LL_TWO_REGION_CASES.items():
        paths[name] = write_screened_rutherford_element(
            out_dir,
            name,
            spec["sigma_t"],
            spec["eta"],
            mu_star=mu_star,
            beta_max=beta_max,
            energy_grid_eV=energy_grid_eV,
            n_mu=n_mu,
            atomic_number=spec["atomic_number"],
            atomic_weight_ratio=spec["atomic_weight_ratio"],
        )
    return paths


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Write synthetic one-speed screened-Rutherford electron libraries "
            "in the MCDC HDF5 format."
        )
    )
    parser.add_argument(
        "--out-dir",
        default=str(Path(__file__).resolve().parents[2] / "synthetic_mcdc_data"),
        help="Output directory for the generated HDF5 files.",
    )
    parser.add_argument(
        "--mu-star",
        type=float,
        default=0.9,
        help="Fixed mu_star written into the GFP2/GFP3 groups.",
    )
    parser.add_argument(
        "--beta-max",
        type=float,
        default=100.0,
        help="Maximum beta used by the continuous GFP2 fit.",
    )
    parser.add_argument(
        "--n-mu",
        type=int,
        default=1200,
        help="Number of mu points in each tabulated screened-Rutherford PDF.",
    )
    parser.add_argument(
        "--energy-min",
        type=float,
        default=float(DEFAULT_ENERGY_GRID_EV[0]),
        help="Minimum energy grid point in eV.",
    )
    parser.add_argument(
        "--energy-max",
        type=float,
        default=float(DEFAULT_ENERGY_GRID_EV[1]),
        help="Maximum energy grid point in eV.",
    )
    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    energy_grid = np.array([args.energy_min, args.energy_max], dtype="f8")
    paths = write_ll_two_region_library(
        args.out_dir,
        mu_star=args.mu_star,
        beta_max=args.beta_max,
        energy_grid_eV=energy_grid,
        n_mu=args.n_mu,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
