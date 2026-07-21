from __future__ import annotations

import unittest

import numpy as np

from qbp_validation.cases import complex_magnitudes, observables, phases
from qbp_validation.circuits import complex_phase_measurement_circuit, probabilities
from qbp_validation.decoders import decode_phase_gradient, phase_record
from qbp_validation.reference import complex_phase_gradient


class ComplexPhaseTests(unittest.TestCase):
    def test_signed_one_hot_phase_estimator(self) -> None:
        for n in range(1, 5):
            magnitudes = complex_magnitudes(n)
            leaf_phases = phases(n)
            for obs_index, observable in enumerate(observables(n)):
                with self.subTest(n=n, observable=obs_index):
                    probs = probabilities(complex_phase_measurement_circuit(magnitudes, leaf_phases, observable))
                    decoded = decode_phase_gradient(probs)
                    exact = complex_phase_gradient(magnitudes, leaf_phases, observable)
                    np.testing.assert_allclose(decoded, exact, atol=3e-12, rtol=0.0)
                    self.assertAlmostEqual(float(decoded.sum()), 0.0, places=11)

    def test_every_phase_record_has_fixed_norm_two(self) -> None:
        for N in (2, 4, 8, 16):
            for ancilla in (0, 1):
                for leaf in range(N):
                    self.assertAlmostEqual(float(np.linalg.norm(phase_record(ancilla, leaf, N))), 2.0)

    def test_zero_amplitude_leaf_phase_derivatives_vanish(self) -> None:
        n = 3
        magnitudes = complex_magnitudes(n)
        magnitudes[0] = 0.0  # The entire left half has nonzero mass; right half vanishes.
        leaf_phases = phases(n)
        observable = observables(n)[-1]
        decoded = decode_phase_gradient(
            probabilities(complex_phase_measurement_circuit(magnitudes, leaf_phases, observable))
        )
        exact = complex_phase_gradient(magnitudes, leaf_phases, observable)
        np.testing.assert_allclose(decoded, exact, atol=3e-12, rtol=0.0)
        np.testing.assert_allclose(decoded[4:], 0.0, atol=3e-12, rtol=0.0)


if __name__ == "__main__":
    unittest.main()
