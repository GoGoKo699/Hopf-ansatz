"""Deterministic validation cases and Hermitian-unitary observables."""
from __future__ import annotations

import numpy as np

from .conventions import PolyBranch, PolyLeaf, PolyTree


def regular_magnitudes(n: int, seed: int = 1000) -> np.ndarray:
    rng = np.random.default_rng(seed + n)
    N = 1 << n
    theta = rng.uniform(0.17, 1.31, size=N - 1)
    # Exercise the real final-layer sign convention without approaching zeros.
    final_start = (1 << (n - 1)) - 1
    theta[final_start:] = rng.uniform(0.2, 2.0 * np.pi - 0.2, size=N // 2)
    return theta


def complex_magnitudes(n: int, seed: int = 2000) -> np.ndarray:
    rng = np.random.default_rng(seed + n)
    return rng.uniform(0.14, 1.40, size=(1 << n) - 1)


def phases(n: int, seed: int = 3000) -> np.ndarray:
    rng = np.random.default_rng(seed + n)
    return rng.uniform(0.0, 2.0 * np.pi, size=1 << n)


def singular_magnitudes(n: int) -> np.ndarray:
    theta = complex_magnitudes(n, seed=4000)
    if n == 1:
        theta[0] = 0.0
    else:
        theta[0] = 0.0
        theta[1] = np.pi / 2.0
    return theta


def pauli_matrix(letter: str) -> np.ndarray:
    matrices = {
        "I": np.eye(2, dtype=complex),
        "X": np.asarray([[0, 1], [1, 0]], dtype=complex),
        "Y": np.asarray([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.asarray([[1, 0], [0, -1]], dtype=complex),
    }
    return matrices[letter]


def pauli_string(word: str) -> np.ndarray:
    matrix = np.asarray([[1.0]], dtype=complex)
    for letter in word:
        matrix = np.kron(matrix, pauli_matrix(letter))
    return matrix


def diagonal_reflection(n: int) -> np.ndarray:
    N = 1 << n
    signs = np.asarray([1.0 if ((3 * index + index.bit_count()) % 5) < 3 else -1.0 for index in range(N)])
    return np.diag(signs.astype(complex))


def householder_reflection(n: int, seed: int = 5000) -> np.ndarray:
    rng = np.random.default_rng(seed + n)
    N = 1 << n
    vector = rng.normal(size=N) + 1j * rng.normal(size=N)
    vector /= np.linalg.norm(vector)
    return np.eye(N, dtype=complex) - 2.0 * np.outer(vector, np.conjugate(vector))


def observables(n: int) -> tuple[np.ndarray, ...]:
    if n == 1:
        word = "Y"
    else:
        letters = "XYZI"
        word = "".join(letters[index % len(letters)] for index in range(n))
    return (pauli_string(word), diagonal_reflection(n), householder_reflection(n))


def balanced_poly_tree(labels: list[int], prefix: str = "v") -> PolyTree:
    if len(labels) == 1:
        return PolyLeaf(labels[0])
    midpoint = len(labels) // 2
    return PolyBranch(
        f"{prefix}_{len(labels)}_{labels[0]}",
        balanced_poly_tree(labels[:midpoint], prefix + "L"),
        balanced_poly_tree(labels[midpoint:], prefix + "R"),
    )


def unbalanced_poly_tree() -> PolyTree:
    # Eight leaves, deliberately permuted so the root anchor is nonzero and the
    # Pauli-X translation in the theorem is exercised.
    leaves = [PolyLeaf(label) for label in (3, 0, 6, 1, 7, 2, 5, 4)]
    left = PolyBranch(
        "u_left",
        leaves[0],
        PolyBranch("u_left_tail", leaves[1], PolyBranch("u_left_deep", leaves[2], leaves[3])),
    )
    right = PolyBranch(
        "u_right",
        PolyBranch("u_right_a", leaves[4], leaves[5]),
        PolyBranch("u_right_b", leaves[6], leaves[7]),
    )
    return PolyBranch("u_root", left, right)


def poly_angles(tree: PolyTree, seed: int = 6000) -> dict[str, float]:
    from .conventions import poly_preorder

    rng = np.random.default_rng(seed)
    return {node.key: float(rng.uniform(-1.2, 1.2)) for node in poly_preorder(tree)}
