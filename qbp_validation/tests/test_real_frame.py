from __future__ import annotations

import unittest

import numpy as np

from qbp_validation.cases import regular_magnitudes, singular_magnitudes
from qbp_validation.circuits import depth_preparation_circuit, frame_circuit, statevector
from qbp_validation.conventions import marker_label
from qbp_validation.reference import real_frame_matrix, real_tree_data


class RealFrameTests(unittest.TestCase):
    def _check_case(self, theta: np.ndarray) -> None:
        data = real_tree_data(theta)
        N = 1 << data.n
        qibo_unitary = np.asarray(frame_circuit(theta).unitary())
        reference = real_frame_matrix(theta)
        np.testing.assert_allclose(qibo_unitary, reference, atol=2e-12, rtol=0.0)
        np.testing.assert_allclose(qibo_unitary.conj().T @ qibo_unitary, np.eye(N), atol=2e-12, rtol=0.0)
        np.testing.assert_allclose(qibo_unitary[:, 0], data.state, atol=2e-12, rtol=0.0)

        for node in range(1, N):
            marker = marker_label(node, data.n)
            complement = data.complements[node]
            self.assertIsNotNone(complement)
            np.testing.assert_allclose(qibo_unitary[:, marker], complement, atol=2e-12, rtol=0.0)
            np.testing.assert_allclose(
                data.derivatives[node - 1],
                data.incoming[node] * complement,
                atol=2e-12,
                rtol=0.0,
            )

        depth_state = statevector(depth_preparation_circuit(theta))
        np.testing.assert_allclose(depth_state, data.state, atol=2e-12, rtol=0.0)

    def test_regular_frames_n1_to_n4(self) -> None:
        for n in range(1, 5):
            with self.subTest(n=n):
                self._check_case(regular_magnitudes(n))

    def test_singular_frame_completion_remains_unitary(self) -> None:
        for n in range(1, 5):
            with self.subTest(n=n):
                self._check_case(singular_magnitudes(n))


if __name__ == "__main__":
    unittest.main()
