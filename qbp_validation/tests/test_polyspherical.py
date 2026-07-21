from __future__ import annotations

import unittest

import numpy as np

from qbp_validation.cases import (
    balanced_poly_tree,
    observables,
    poly_angles,
    unbalanced_poly_tree,
)
from qbp_validation.circuits import (
    polyspherical_frame_circuit,
    polyspherical_global_measurement_circuit,
    probabilities,
)
from qbp_validation.conventions import poly_relative_markers
from qbp_validation.decoders import decode_polyspherical_gradient
from qbp_validation.reference import (
    polyspherical_gradient,
    polyspherical_shifted_frame_matrix,
    polyspherical_tree_data,
)


class PolysphericalTests(unittest.TestCase):
    def _cases(self):
        return (
            (balanced_poly_tree(list(range(8))), poly_angles(balanced_poly_tree(list(range(8))), 6100)),
            (unbalanced_poly_tree(), poly_angles(unbalanced_poly_tree(), 6200)),
        )

    def test_frame_and_local_complements(self) -> None:
        n = 3
        for case_index, (tree, angles) in enumerate(self._cases()):
            with self.subTest(case=case_index):
                qibo_frame = np.asarray(polyspherical_frame_circuit(tree, angles, n).unitary())
                reference = polyspherical_shifted_frame_matrix(tree, angles, n)
                np.testing.assert_allclose(qibo_frame, reference, atol=3e-12, rtol=0.0)
                np.testing.assert_allclose(
                    qibo_frame.conj().T @ qibo_frame, np.eye(1 << n), atol=3e-12, rtol=0.0
                )

    def test_parity_estimator_for_representative_trees(self) -> None:
        n = 3
        for case_index, (tree, angles) in enumerate(self._cases()):
            data = polyspherical_tree_data(tree, angles, n)
            for obs_index, observable in enumerate(observables(n)):
                with self.subTest(case=case_index, observable=obs_index):
                    probs = probabilities(
                        polyspherical_global_measurement_circuit(tree, angles, n, observable)
                    )
                    decoded = decode_polyspherical_gradient(
                        probs, data.incoming, poly_relative_markers(tree)
                    )
                    exact = polyspherical_gradient(tree, angles, n, observable)
                    for key in exact:
                        self.assertAlmostEqual(decoded[key], exact[key], places=11)

    def test_derivative_factorization_and_diagonal_metric(self) -> None:
        tree = unbalanced_poly_tree()
        angles = poly_angles(tree, 6300)
        data = polyspherical_tree_data(tree, angles, 3)
        keys = list(data.derivatives)
        gram = np.asarray(
            [
                [np.real(np.vdot(data.derivatives[a], data.derivatives[b])) for b in keys]
                for a in keys
            ]
        )
        expected = np.diag([data.incoming[key] ** 2 for key in keys])
        np.testing.assert_allclose(gram, expected, atol=3e-12, rtol=0.0)


if __name__ == "__main__":
    unittest.main()
