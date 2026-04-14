import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyeedl.function import (  # noqa: E402
    compute_delta_kernel_fp_params,
    compute_gfp2_params,
    compute_gfp3_params,
    compute_optimal_mu_star_fp,
    compute_optimal_mu_star_gfp2,
    validate_mu_star_convergence,
)


def _p3(mu):
    return 0.5 * (5.0 * mu**3 - 3.0 * mu)


class TestDeltaKernelGFP2(unittest.TestCase):
    def test_optimal_mu_star_matches_ell3_constraint(self):
        res = compute_gfp2_params(10.0, 8.0, 5.0)
        Sigma_a3 = 7.0

        mu_res = compute_optimal_mu_star_gfp2(res["alpha"], res["beta"], Sigma_a3)
        self.assertTrue(mu_res["success"])
        mu_star = mu_res["mu_star"]
        self.assertGreater(mu_star, -1.0)
        self.assertLess(mu_star, 1.0)

        x3 = (1.0 - _p3(mu_star)) / (1.0 - mu_star)
        lam_delta = -2.0 * res["alpha"] * x3 / (1.0 + 2.0 * res["beta"] * x3)
        self.assertAlmostEqual(lam_delta, -Sigma_a3, places=12)

    def test_convergence_matches_continuous_gfp2(self):
        results = validate_mu_star_convergence(
            10.0, 8.0, 5.0, mu_star_values=[0.9, 0.99, 0.9999, 0.99999]
        )

        self.assertAlmostEqual(results[0.99999][1]["rel_error"], 0.0, places=12)
        self.assertLess(results[0.99999][2]["rel_error"], 1.0e-5)
        self.assertLess(results[0.99999][3]["rel_error"], 1.0e-5)
        self.assertLess(
            results[0.99999][3]["rel_error"],
            results[0.99][3]["rel_error"],
        )

    def test_fp_optimal_mu_star_boundary_is_reported_not_clipped(self):
        mu_res = compute_optimal_mu_star_fp(1.0, 3.0)
        self.assertFalse(mu_res["success"])
        self.assertTrue(np.isnan(mu_res["mu_star"]))
        self.assertAlmostEqual(mu_res["mu_star_raw"], 1.0, places=12)

    def test_gfp2_optimal_mu_star_boundary_is_reported_not_clipped(self):
        mu_res = compute_optimal_mu_star_gfp2(1.0, 0.0, 12.0)
        self.assertFalse(mu_res["success"])
        self.assertTrue(np.isnan(mu_res["mu_star"]))
        self.assertAlmostEqual(mu_res["mu_star_raw"], 1.0, places=12)

    def test_delta_kernel_fp_requires_explicit_valid_mu_star(self):
        res_missing = compute_delta_kernel_fp_params(10.0, 8.0, mu_star=None)
        self.assertFalse(res_missing["success"])
        self.assertTrue(np.isnan(res_missing["mu_star"]))
        self.assertTrue(np.isnan(res_missing["Sigma_delta0"]))

        res_invalid = compute_delta_kernel_fp_params(10.0, 8.0, mu_star=1.0)
        self.assertFalse(res_invalid["success"])
        self.assertTrue(np.isnan(res_invalid["mu_star"]))
        self.assertAlmostEqual(res_invalid["mu_star_raw"], 1.0, places=12)

    def test_gfp2_invalid_fixed_mu_star_is_reported(self):
        res = compute_gfp2_params(10.0, 8.0, 5.0, mu_star=1.0)
        self.assertFalse(res["success"])
        self.assertTrue(np.isnan(res["Sigma_delta0"]))
        self.assertIn("1-mu_star", res["warning"])

    def test_gfp3_invalid_fixed_mu_star_is_reported(self):
        res = compute_gfp3_params(10.0, 9.5, 9.0, 8.7, mu_star=1.0)
        self.assertFalse(res["success"])
        self.assertTrue(np.isnan(res["Sigma_delta0"]))
        self.assertIn("1-mu_star", res["warning"])


if __name__ == "__main__":
    unittest.main()
