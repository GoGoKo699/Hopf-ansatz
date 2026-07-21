"""Qibo circuit builders for deterministic exact-logical validation.

The builders emit actual ``qibo.models.Circuit`` objects.  Analytic state and
coordinate-derivative formulas live exclusively in :mod:`reference`.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from qibo import gates, models, set_backend

from .conventions import (
    PolyBranch,
    PolyTree,
    bit_at,
    depth_gate_specs,
    frame_gate_specs,
    infer_n_from_magnitudes,
    poly_anchor,
    poly_marker,
    poly_preorder,
    validate_poly_tree,
)

set_backend("numpy")

SDG = np.asarray([[1.0, 0.0], [0.0, -1.0j]], dtype=complex)


def _controlled_ry_with_pattern(
    circuit: models.Circuit,
    target: int,
    controls: Sequence[int],
    control_values: Sequence[int],
    angle: float,
) -> None:
    if len(controls) != len(control_values):
        raise ValueError("Control qubits and values have different lengths.")
    open_controls = [q for q, value in zip(controls, control_values) if value == 0]
    for q in open_controls:
        circuit.add(gates.X(q))
    gate = gates.RY(target, 2.0 * float(angle))
    if controls:
        gate = gate.controlled_by(*controls)
    circuit.add(gate)
    for q in reversed(open_controls):
        circuit.add(gates.X(q))


def add_real_frame(
    circuit: models.Circuit,
    theta: np.ndarray,
    system_qubits: Sequence[int],
    *,
    inverse: bool = False,
) -> None:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    n = infer_n_from_magnitudes(theta)
    if len(system_qubits) != n:
        raise ValueError("System-qubit list length does not match theta.")
    specs = list(frame_gate_specs(n))
    if inverse:
        specs.reverse()
    for spec in specs:
        target = system_qubits[spec.target]
        controls = [q for index, q in enumerate(system_qubits) if index != spec.target]
        values = [bit_at(spec.anchor, index, n) for index in range(n) if index != spec.target]
        angle = -theta[spec.node - 1] if inverse else theta[spec.node - 1]
        _controlled_ry_with_pattern(circuit, target, controls, values, angle)


def add_depth_layer(
    circuit: models.Circuit,
    theta: np.ndarray,
    system_qubits: Sequence[int],
    depth: int,
    *,
    inverse: bool = False,
) -> None:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    n = infer_n_from_magnitudes(theta)
    if len(system_qubits) != n:
        raise ValueError("System-qubit list length does not match theta.")
    for spec in depth_gate_specs(n, depth):
        controls = list(system_qubits[:depth])
        values = [((spec.position >> (depth - 1 - pos)) & 1) for pos in range(depth)]
        angle = -theta[spec.node - 1] if inverse else theta[spec.node - 1]
        _controlled_ry_with_pattern(circuit, system_qubits[depth], controls, values, angle)


def add_depth_preparation(
    circuit: models.Circuit,
    theta: np.ndarray,
    system_qubits: Sequence[int],
) -> None:
    n = infer_n_from_magnitudes(theta)
    for depth in range(n):
        add_depth_layer(circuit, theta, system_qubits, depth)


def add_inverse_depth_suffix(
    circuit: models.Circuit,
    theta: np.ndarray,
    system_qubits: Sequence[int],
    selected_depth: int,
) -> None:
    n = infer_n_from_magnitudes(theta)
    for depth in range(n - 1, selected_depth, -1):
        add_depth_layer(circuit, theta, system_qubits, depth, inverse=True)


def add_phase_layer(
    circuit: models.Circuit,
    phases: np.ndarray,
    system_qubits: Sequence[int],
    *,
    inverse: bool = False,
) -> None:
    phases = np.asarray(phases, dtype=float).reshape(-1)
    if phases.size != 1 << len(system_qubits):
        raise ValueError("Phase vector length does not match the system register.")
    sign = -1.0 if inverse else 1.0
    matrix = np.diag(np.exp(1j * sign * phases))
    circuit.add(gates.Unitary(matrix, *system_qubits))


def add_controlled_observable(
    circuit: models.Circuit,
    observable: np.ndarray,
    ancilla: int,
    system_qubits: Sequence[int],
) -> None:
    observable = np.asarray(observable, dtype=complex)
    N = 1 << len(system_qubits)
    if observable.shape != (N, N):
        raise ValueError("Observable dimension does not match system register.")
    circuit.add(gates.Unitary(observable, *system_qubits).controlled_by(ancilla))


def add_x_basis_rotation(circuit: models.Circuit, qubits: Sequence[int]) -> None:
    for q in qubits:
        circuit.add(gates.H(q))


def add_y_basis_rotation(circuit: models.Circuit, qubit: int) -> None:
    # State-update order is S^dagger followed by H, i.e. net H S^dagger.
    circuit.add(gates.Unitary(SDG, qubit))
    circuit.add(gates.H(qubit))


def frame_circuit(theta: np.ndarray, *, inverse: bool = False) -> models.Circuit:
    n = infer_n_from_magnitudes(theta)
    circuit = models.Circuit(n)
    add_real_frame(circuit, theta, tuple(range(n)), inverse=inverse)
    return circuit


def depth_preparation_circuit(theta: np.ndarray) -> models.Circuit:
    n = infer_n_from_magnitudes(theta)
    circuit = models.Circuit(n)
    add_depth_preparation(circuit, theta, tuple(range(n)))
    return circuit


def real_global_measurement_circuit(theta: np.ndarray, observable: np.ndarray) -> models.Circuit:
    n = infer_n_from_magnitudes(theta)
    ancilla = 0
    system = tuple(range(1, n + 1))
    circuit = models.Circuit(n + 1)
    circuit.add(gates.H(ancilla))
    add_real_frame(circuit, theta, system)
    add_controlled_observable(circuit, observable, ancilla, system)
    add_real_frame(circuit, theta, system, inverse=True)
    add_x_basis_rotation(circuit, (ancilla, *system))
    return circuit


def complex_magnitude_measurement_circuit(
    magnitudes: np.ndarray,
    phases: np.ndarray,
    observable: np.ndarray,
) -> models.Circuit:
    n = infer_n_from_magnitudes(magnitudes)
    ancilla = 0
    system = tuple(range(1, n + 1))
    circuit = models.Circuit(n + 1)
    circuit.add(gates.H(ancilla))
    add_real_frame(circuit, magnitudes, system)
    add_phase_layer(circuit, phases, system)
    add_controlled_observable(circuit, observable, ancilla, system)
    add_phase_layer(circuit, phases, system, inverse=True)
    add_real_frame(circuit, magnitudes, system, inverse=True)
    add_x_basis_rotation(circuit, (ancilla, *system))
    return circuit


def complex_phase_measurement_circuit(
    magnitudes: np.ndarray,
    phases: np.ndarray,
    observable: np.ndarray,
) -> models.Circuit:
    n = infer_n_from_magnitudes(magnitudes)
    ancilla = 0
    system = tuple(range(1, n + 1))
    circuit = models.Circuit(n + 1)
    circuit.add(gates.H(ancilla))
    add_real_frame(circuit, magnitudes, system)
    add_phase_layer(circuit, phases, system)
    add_controlled_observable(circuit, observable, ancilla, system)
    add_y_basis_rotation(circuit, ancilla)
    return circuit


def real_checkpoint_measurement_circuit(
    theta: np.ndarray,
    observable: np.ndarray,
    depth: int,
) -> models.Circuit:
    n = infer_n_from_magnitudes(theta)
    if not 0 <= depth < n:
        raise ValueError("Depth must lie in 0, ..., n-1.")
    ancilla = 0
    system = tuple(range(1, n + 1))
    circuit = models.Circuit(n + 1)
    circuit.add(gates.H(ancilla))
    add_depth_preparation(circuit, theta, system)
    add_controlled_observable(circuit, observable, ancilla, system)
    add_inverse_depth_suffix(circuit, theta, system, depth)
    add_y_basis_rotation(circuit, ancilla)
    add_y_basis_rotation(circuit, system[depth])
    return circuit


def complex_checkpoint_measurement_circuit(
    magnitudes: np.ndarray,
    phases: np.ndarray,
    observable: np.ndarray,
    depth: int,
) -> models.Circuit:
    n = infer_n_from_magnitudes(magnitudes)
    if not 0 <= depth < n:
        raise ValueError("Depth must lie in 0, ..., n-1.")
    ancilla = 0
    system = tuple(range(1, n + 1))
    circuit = models.Circuit(n + 1)
    circuit.add(gates.H(ancilla))
    add_depth_preparation(circuit, magnitudes, system)
    add_phase_layer(circuit, phases, system)
    add_controlled_observable(circuit, observable, ancilla, system)
    add_phase_layer(circuit, phases, system, inverse=True)
    add_inverse_depth_suffix(circuit, magnitudes, system, depth)
    add_y_basis_rotation(circuit, ancilla)
    add_y_basis_rotation(circuit, system[depth])
    return circuit


def _two_level_rotation(dimension: int, first: int, second: int, angle: float) -> np.ndarray:
    matrix = np.eye(dimension, dtype=complex)
    c = float(np.cos(angle))
    s = float(np.sin(angle))
    matrix[first, first] = c
    matrix[second, first] = s
    matrix[first, second] = -s
    matrix[second, second] = c
    return matrix


def _add_poly_translation(circuit: models.Circuit, system_qubits: Sequence[int], label: int) -> None:
    n = len(system_qubits)
    for index, qubit in enumerate(system_qubits):
        if bit_at(label, index, n):
            circuit.add(gates.X(qubit))


def add_polyspherical_frame(
    circuit: models.Circuit,
    tree: PolyTree,
    angles: Mapping[str, float],
    system_qubits: Sequence[int],
    *,
    shifted: bool = True,
    inverse: bool = False,
) -> None:
    n = len(system_qubits)
    validate_poly_tree(tree, n, angles)
    nodes = list(poly_preorder(tree))
    dimension = 1 << n
    if not inverse:
        if shifted:
            _add_poly_translation(circuit, system_qubits, poly_anchor(tree))
        for node in nodes:
            matrix = _two_level_rotation(
                dimension,
                poly_anchor(node.left),
                poly_anchor(node.right),
                float(angles[node.key]),
            )
            circuit.add(gates.Unitary(matrix, *system_qubits))
    else:
        for node in reversed(nodes):
            matrix = _two_level_rotation(
                dimension,
                poly_anchor(node.left),
                poly_anchor(node.right),
                -float(angles[node.key]),
            )
            circuit.add(gates.Unitary(matrix, *system_qubits))
        if shifted:
            _add_poly_translation(circuit, system_qubits, poly_anchor(tree))


def polyspherical_frame_circuit(
    tree: PolyTree,
    angles: Mapping[str, float],
    n: int,
    *,
    shifted: bool = True,
) -> models.Circuit:
    circuit = models.Circuit(n)
    add_polyspherical_frame(circuit, tree, angles, tuple(range(n)), shifted=shifted)
    return circuit


def polyspherical_global_measurement_circuit(
    tree: PolyTree,
    angles: Mapping[str, float],
    n: int,
    observable: np.ndarray,
) -> models.Circuit:
    ancilla = 0
    system = tuple(range(1, n + 1))
    circuit = models.Circuit(n + 1)
    circuit.add(gates.H(ancilla))
    add_polyspherical_frame(circuit, tree, angles, system, shifted=True)
    add_controlled_observable(circuit, observable, ancilla, system)
    add_polyspherical_frame(circuit, tree, angles, system, shifted=True, inverse=True)
    add_x_basis_rotation(circuit, (ancilla, *system))
    return circuit


def statevector(circuit: models.Circuit) -> np.ndarray:
    result = circuit()
    return np.asarray(result.state(), dtype=complex).reshape(-1)


def probabilities(circuit: models.Circuit) -> np.ndarray:
    state = statevector(circuit)
    probs = np.abs(state) ** 2
    probs /= probs.sum()
    return probs
