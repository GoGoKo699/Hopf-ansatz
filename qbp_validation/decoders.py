"""Exact-distribution and record-level decoders for Hopf-QBP validation."""
from __future__ import annotations

import numpy as np

from .conventions import marker_label, parity


def fwht(values: np.ndarray) -> np.ndarray:
    """Unnormalized in-place-style Walsh-Hadamard transform on a copy."""

    out = np.asarray(values, dtype=float).reshape(-1).copy()
    n = out.size
    if n == 0 or n & (n - 1):
        raise ValueError("FWHT input length must be a positive power of two.")
    width = 1
    while width < n:
        for start in range(0, n, 2 * width):
            a = out[start : start + width].copy()
            b = out[start + width : start + 2 * width].copy()
            out[start : start + width] = a + b
            out[start + width : start + 2 * width] = a - b
        width *= 2
    return out


def reshape_ancilla_system(probabilities: np.ndarray) -> np.ndarray:
    probs = np.asarray(probabilities, dtype=float).reshape(-1)
    if probs.size < 4 or probs.size & (probs.size - 1):
        raise ValueError("Probability vector must describe one ancilla and at least one system qubit.")
    N = probs.size // 2
    return probs.reshape(2, N)


def signed_system_histogram(probabilities: np.ndarray) -> np.ndarray:
    probs = reshape_ancilla_system(probabilities)
    return probs[0] - probs[1]


def global_moments_direct(probabilities: np.ndarray) -> np.ndarray:
    signed = signed_system_histogram(probabilities)
    N = signed.size
    moments = np.zeros(N, dtype=float)
    for label in range(N):
        moments[label] = sum(
            value * (-1.0 if parity(label, outcome) else 1.0)
            for outcome, value in enumerate(signed)
        )
    return moments


def global_moments_fwht(probabilities: np.ndarray) -> np.ndarray:
    return fwht(signed_system_histogram(probabilities))


def decode_balanced_magnitude_gradient(
    probabilities: np.ndarray,
    metric: np.ndarray,
    n: int,
    *,
    use_fwht: bool = True,
) -> np.ndarray:
    metric = np.asarray(metric, dtype=float).reshape(-1)
    if metric.size != (1 << n) - 1:
        raise ValueError("Metric length does not match n.")
    moments = global_moments_fwht(probabilities) if use_fwht else global_moments_direct(probabilities)
    return np.asarray(
        [2.0 * np.sqrt(max(metric[node - 1], 0.0)) * moments[marker_label(node, n)] for node in range(1, 1 << n)],
        dtype=float,
    )


def decode_phase_gradient(probabilities: np.ndarray) -> np.ndarray:
    probs = reshape_ancilla_system(probabilities)
    return 2.0 * (probs[0] - probs[1])


def decode_checkpoint_gradient(probabilities: np.ndarray, n: int, depth: int) -> np.ndarray:
    if not 0 <= depth < n:
        raise ValueError("Depth must lie in 0, ..., n-1.")
    probs = reshape_ancilla_system(probabilities)
    width = 1 << depth
    result = np.zeros(width, dtype=float)
    target_shift = n - depth - 1
    for ancilla in (0, 1):
        for system_label in range(1 << n):
            prefix = system_label >> (n - depth) if depth else 0
            target = (system_label >> target_shift) & 1
            sign = -2.0 * (-1.0 if (ancilla + target) & 1 else 1.0)
            result[prefix] += sign * probs[ancilla, system_label]
    return result


def decode_polyspherical_gradient(
    probabilities: np.ndarray,
    incoming: dict[str, float],
    relative_markers: dict[str, int],
) -> dict[str, float]:
    moments = global_moments_fwht(probabilities)
    return {
        key: 2.0 * incoming[key] * moments[relative_markers[key]]
        for key in relative_markers
    }


def phase_record(ancilla_bit: int, leaf: int, N: int) -> np.ndarray:
    record = np.zeros(N, dtype=float)
    record[leaf] = 2.0 * (-1.0 if ancilla_bit else 1.0)
    return record


def checkpoint_record(ancilla_bit: int, target_bit: int, prefix: int, width: int) -> np.ndarray:
    record = np.zeros(width, dtype=float)
    record[prefix] = -2.0 * (-1.0 if (ancilla_bit + target_bit) & 1 else 1.0)
    return record
