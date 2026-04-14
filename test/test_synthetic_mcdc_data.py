import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyeedl.synthetic import (  # noqa: E402
    LL_TWO_REGION_CASES,
    write_ll_two_region_library,
    write_screened_rutherford_element,
)


class TestSyntheticScreenedRutherfordLibrary(unittest.TestCase):
    def test_invalid_mu_star_raises_instead_of_being_forced(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                write_screened_rutherford_element(
                    tmpdir,
                    "BAD",
                    sigma_t=100.0,
                    eta=1.0e-4,
                    mu_star=1.0,
                )

    def test_single_element_writer_creates_expected_mcdc_groups(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_screened_rutherford_element(
                tmpdir,
                "SR_TEST",
                sigma_t=100.0,
                eta=1.0e-4,
                mu_star=0.9,
                n_mu=256,
            )

            with h5py.File(path, "r") as h5f:
                self.assertEqual(h5f["atomic_number"][()], 1)
                self.assertAlmostEqual(h5f["atomic_weight_ratio"][()], 1.0)
                np.testing.assert_allclose(
                    h5f["electron_reactions/xs_energy_grid"][:],
                    np.array([1.0e5, 2.0e6]),
                )
                np.testing.assert_allclose(
                    h5f["electron_reactions/elastic_scattering/MT525/xs"][:],
                    np.array([100.0, 100.0]),
                )
                np.testing.assert_allclose(
                    h5f["electron_reactions/elastic_scattering/MT525/xs_large"][:],
                    np.array([100.0, 100.0]),
                )

                gfp2 = h5f["electron_reactions/elastic_scattering/MT525/gfp2"]
                self.assertAlmostEqual(gfp2["mu_star"][()], 0.9)
                self.assertTrue(np.all(gfp2["success"][:]))
                self.assertTrue(np.all(~gfp2["mu_star_optimal_success"][:]))

                sc = h5f[
                    "electron_reactions/elastic_scattering/MT525/scattering_cosine_coupled"
                ]
                offset = sc["energy_offset"][:]
                value = sc["value"][:]
                pdf = sc["PDF"][:]
                start = int(offset[0])
                end = int(offset[1])
                area = np.trapezoid(pdf[start:end], value[start:end])
                self.assertAlmostEqual(area, 1.0, places=10)

    def test_ll_two_region_writer_creates_both_region_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = write_ll_two_region_library(tmpdir, mu_star=0.9, n_mu=256)
            self.assertEqual(set(paths.keys()), set(LL_TWO_REGION_CASES.keys()))

            for name, spec in LL_TWO_REGION_CASES.items():
                with h5py.File(paths[name], "r") as h5f:
                    elastic = h5f["electron_reactions/elastic_scattering/MT525"]
                    np.testing.assert_allclose(
                        elastic["xs"][:],
                        np.full(2, spec["sigma_t"]),
                    )
                    self.assertAlmostEqual(elastic.attrs["synthetic_eta"], spec["eta"])
                    self.assertAlmostEqual(
                        elastic.attrs["synthetic_sigma_t"], spec["sigma_t"]
                    )


if __name__ == "__main__":
    unittest.main()
