from __future__ import annotations

import unittest

import numpy as np
from qibo import models

from qbp_validation.cases import complex_magnitudes, observables, phases
from qbp_validation.circuits import (
    add_phase_layer,
    add_real_frame,
    complex_magnitude_measurement_circuit,
    probabilities,
)
from qbp_validation.conventions import marker_label
from qbp_validation.decoders import decode_balanced_magnitude_gradient
from qbp_validation.reference import (
    complex_magnitude_derivatives,
    complex_magnitude_gradient,
    complex_state,
    real_tree_data,
)


class ComplexMagnitudeTests(unittest.TestCase):
    def test_phase_dressed_frame_columns(self) -> None:
        for n in range(1, 5):
            magnitudes = complex_magnitudes(n)
            leaf_phases = phases(n)
            circuit = models.Circuit(n)
            system = tuple(range(n))
            add_real_frame(circuit, magnitudes, system)
            add_phase_layer(circuit, leaf_phases, system)
            unitary = np.asarray(circuit.unitary())
            np.testing.assert_allclose(unitary[:, 0], complex_state(magnitudes, leaf_phases), atol=2e-12, rtol=0.0)
            derivatives = complex_magnitude_derivatives(magnitudes, leaf_phases)
            incoming = real_tree_data(magnitudes).incoming
            for node, derivative in enumerate(derivatives, start=1):
                expected_column = derivative / incoming[node]
                np.testing.assert_allclose(
                    unitary[:, marker_label(node, n)], expected_column, atol=2e-12, rtol=0.0
                )

    def test_exact_magnitude_estimator(self) -> None:
        for n in range(1, 5):
            magnitudes = complex_magnitudes(n)
            leaf_phases = phases(n)
            metric = real_tree_data(magnitudes).metric
            for obs_index, observable in enumerate(observables(n)):
                with self.subTest(n=n, observable=obs_index):
                    probs = probabilities(
                        complex_magnitude_measurement_circuit(magnitudes, leaf_phases, observable)
                    )
                    decoded = decode_balanced_magnitude_gradient(probs, metric, n)
                    exact = complex_magnitude_gradient(magnitudes, leaf_phases, observable)
                    np.testing.assert_allclose(decoded, exact, atol=3e-12, rtol=0.0)


if __name__ == "__main__":
    unittest.main()
