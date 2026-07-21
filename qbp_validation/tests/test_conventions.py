from __future__ import annotations

import unittest

import numpy as np
from qibo import gates, models

from qbp_validation.cases import observables
from qbp_validation.circuits import add_y_basis_rotation, probabilities
from qbp_validation.conventions import bit_at, frame_gate_specs, marker_label, parity
from qbp_validation.reference import is_hermitian_unitary


class ConventionTests(unittest.TestCase):
    def test_marker_map_is_bijection_onto_nonzero_labels(self) -> None:
        for n in range(1, 6):
            labels = [marker_label(node, n) for node in range(1, 1 << n)]
            self.assertEqual(sorted(labels), list(range(1, 1 << n)))

    def test_big_endian_bit_convention(self) -> None:
        label = 0b1011
        self.assertEqual([bit_at(label, q, 4) for q in range(4)], [1, 0, 1, 1])
        self.assertEqual(parity(0b0110, 0b1011), 1)

    def test_qibo_qubit_zero_is_most_significant(self) -> None:
        circuit = models.Circuit(3)
        circuit.add(gates.X(0))
        state = np.asarray(circuit().state())
        self.assertEqual(int(np.argmax(np.abs(state))), 0b100)

    def test_y_basis_rotation_has_manuscript_sign(self) -> None:
        # Prepare |+i> = S H |0>.  S^dagger followed by H must map it to |0>.
        circuit = models.Circuit(1)
        circuit.add(gates.H(0))
        circuit.add(gates.Unitary(np.diag([1.0, 1.0j]), 0))
        add_y_basis_rotation(circuit, 0)
        probs = probabilities(circuit)
        np.testing.assert_allclose(probs, [1.0, 0.0], atol=1e-13, rtol=0.0)

    def test_frame_specs_match_active_pairs(self) -> None:
        for n in range(1, 5):
            for spec in frame_gate_specs(n):
                self.assertEqual(spec.anchor ^ spec.marker, 1 << (n - 1 - spec.target))
                self.assertEqual(bit_at(spec.anchor, spec.target, n), 0)
                self.assertEqual(bit_at(spec.marker, spec.target, n), 1)

    def test_validation_observables_satisfy_access_model(self) -> None:
        for n in range(1, 5):
            for observable in observables(n):
                self.assertTrue(is_hermitian_unitary(observable))


if __name__ == "__main__":
    unittest.main()
