from __future__ import annotations

import unittest

from qbp_validation.conventions import (
    checkpoint_cnot_charge_without_observable,
    controlled_ry_cnot_charge,
    depth_layer_cnot_charge,
    depth_preparation_cnot_charge,
    frame_cnot_charge,
)


class ResourceLedgerTests(unittest.TestCase):
    def test_assigned_controlled_ry_charges(self) -> None:
        expected = {0: 0, 1: 2, 2: 6, 3: 14, 4: 30, 5: 56, 6: 72, 7: 88}
        self.assertEqual({m: controlled_ry_cnot_charge(m) for m in expected}, expected)

    def test_four_qubit_frame_and_depth_totals(self) -> None:
        self.assertEqual(frame_cnot_charge(4), 210)
        self.assertEqual([depth_layer_cnot_charge(d) for d in range(4)], [0, 4, 24, 112])
        self.assertEqual(depth_preparation_cnot_charge(4), 140)

    def test_four_qubit_checkpoint_totals(self) -> None:
        self.assertEqual(
            [checkpoint_cnot_charge_without_observable(4, d) for d in range(4)],
            [280, 276, 252, 140],
        )


if __name__ == "__main__":
    unittest.main()
