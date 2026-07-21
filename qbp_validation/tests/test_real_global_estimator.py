from __future__ import annotations

import unittest

import numpy as np

from qbp_validation.cases import observables, regular_magnitudes
from qbp_validation.circuits import probabilities, real_global_measurement_circuit
from qbp_validation.decoders import (
    decode_balanced_magnitude_gradient,
    global_moments_direct,
    global_moments_fwht,
)
from qbp_validation.reference import real_gradient, real_tree_data


class RealGlobalEstimatorTests(unittest.TestCase):
    def test_exact_distribution_returns_complete_gradient(self) -> None:
        for n in range(1, 5):
            theta = regular_magnitudes(n)
            data = real_tree_data(theta)
            for obs_index, observable in enumerate(observables(n)):
                with self.subTest(n=n, observable=obs_index):
                    probs = probabilities(real_global_measurement_circuit(theta, observable))
                    self.assertAlmostEqual(float(probs.sum()), 1.0, places=13)
                    np.testing.assert_allclose(
                        global_moments_fwht(probs),
                        global_moments_direct(probs),
                        atol=2e-12,
                        rtol=0.0,
                    )
                    decoded = decode_balanced_magnitude_gradient(probs, data.metric, n)
                    exact = real_gradient(theta, observable)
                    np.testing.assert_allclose(decoded, exact, atol=3e-12, rtol=0.0)

    def test_per_record_coordinate_range_is_at_most_two(self) -> None:
        for n in range(1, 6):
            metric = real_tree_data(regular_magnitudes(n)).metric
            self.assertLessEqual(float(np.max(2.0 * np.sqrt(metric))), 2.0 + 1e-13)


if __name__ == "__main__":
    unittest.main()
