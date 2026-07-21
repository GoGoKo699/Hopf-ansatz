from __future__ import annotations

import unittest

import numpy as np

from qbp_validation.decoders import (
    decode_checkpoint_gradient,
    decode_phase_gradient,
    fwht,
    global_moments_direct,
    global_moments_fwht,
)


class DecoderTests(unittest.TestCase):
    def test_fwht_matches_dense_walsh_matrix(self) -> None:
        for n in range(1, 7):
            N = 1 << n
            rng = np.random.default_rng(7000 + n)
            values = rng.normal(size=N)
            dense = np.asarray(
                [
                    [(-1.0) ** ((row & column).bit_count() & 1) for column in range(N)]
                    for row in range(N)
                ]
            )
            np.testing.assert_allclose(fwht(values), dense @ values, atol=2e-12, rtol=0.0)

    def test_global_direct_and_fwht_decoders_agree(self) -> None:
        for n in range(1, 6):
            rng = np.random.default_rng(7100 + n)
            probs = rng.uniform(size=2 * (1 << n))
            probs /= probs.sum()
            np.testing.assert_allclose(
                global_moments_direct(probs), global_moments_fwht(probs), atol=2e-12, rtol=0.0
            )

    def test_phase_decoder_is_signed_leaf_histogram(self) -> None:
        probs = np.asarray([[0.10, 0.20, 0.05, 0.15], [0.05, 0.10, 0.25, 0.10]])
        np.testing.assert_allclose(decode_phase_gradient(probs.reshape(-1)), 2.0 * (probs[0] - probs[1]))

    def test_checkpoint_decoder_discards_suffix_bits(self) -> None:
        n = 3
        depth = 1
        probs = np.zeros((2, 1 << n), dtype=float)
        # Same ancilla/prefix/target data but two different suffix values.
        probs[0, 0b010] = 0.2
        probs[0, 0b011] = 0.3
        probs[1, 0b110] = 0.1
        probs[1, 0b111] = 0.4
        decoded = decode_checkpoint_gradient(probs.reshape(-1), n, depth)
        np.testing.assert_allclose(decoded, [1.0, -1.0], atol=1e-13, rtol=0.0)


if __name__ == "__main__":
    unittest.main()
