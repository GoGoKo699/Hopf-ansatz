"""VQE_qibo.py
================================================================================
Layerwise-circuit VQE safeguard for two local n=4 Hopf toy models.

Scope:
    This is a small functional circuit-realizability demo. It explicitly builds
    n=4 label-controlled derivative/branch circuits so the layerwise energy
    estimator can be checked against the exact Hopf gradient. It is not the
    asymptotically optimized indexed-gradient circuit construction or a scaling
    benchmark for the gate-count claims in the paper.

This script replaces the per-parameter Hadamard-test gradient with the
layerwise signed-energy estimator used in the draft:

    * one baseline energy E_psi,
    * one indexed tangent-energy layer per Hopf magnitude depth,
    * one indexed signed-branch layer per Hopf magnitude depth,
    * and, for the complex Hopf ansatz, one additional indexed phase layer.

For each parameter i, the sampled layer energies are combined as

    dE/dtheta_i = 2 sqrt(g_ii) s
        (E_phi_i^(s) - 0.5 * (E_psi + E_partial_i)),

averaged over the sampled s=+ and s=- branch outcomes when both are present.

The two examples are:

    * real Hopf:    local real X/Z + XX/YY/ZZ nearest-neighbor chain;
    * complex Hopf: local chiral X/Y/Z + XX/YY/ZZ + XY-YX chain.

The default initialization is fully random.  Adam uses a fixed learning rate;
there is no line search.

Default run:

    MPLBACKEND=Agg python VQE_qibo.py

Output:

    VQE_qibo.png

The default sampler is "auto": it uses explicit Qibo indexed-preparation
circuits when Qibo is importable, and otherwise falls back to an equivalent
statevector sampler for the ideal layer states.  The fallback is included so
that the non-Qibo path can still be tested on systems without Qibo.
================================================================================
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

try:
    from qibo import models, gates, set_backend

    set_backend("numpy")
    _QIBO_AVAILABLE = True
    _QIBO_IMPORT_ERROR = None
except Exception as exc:  # pragma: no cover - used only without qibo
    models = None
    gates = None
    set_backend = None
    _QIBO_AVAILABLE = False
    _QIBO_IMPORT_ERROR = exc

from hopf_utils import (
    gates_order,
    theta_from_vector,
    vector_from_theta,
    theta_hopf_tangent_state,
    jacobian,
    metric_diagonal,
    clip_theta_hopf_real,
)


# =============================================================================
# Pauli utilities
# =============================================================================

I2 = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
H1 = (1.0 / math.sqrt(2.0)) * np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex)
SDG = np.array([[1.0, 0.0], [0.0, -1.0j]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}

# A Pauli string is represented as (("X", 0), ("Y", 1), ...), where qubit 0 is
# the leftmost Kronecker factor.
PauliString = Tuple[Tuple[str, int], ...]


@dataclass(frozen=True)
class PauliTerm:
    coeff: float
    pauli: PauliString

    @property
    def label(self) -> str:
        if not self.pauli:
            return "I"
        return "".join(f"{p}{q}" for p, q in self.pauli)


def pauli_on(n: int, local_ops: Dict[int, str]) -> np.ndarray:
    """Return an n-qubit Pauli matrix."""
    out = np.array([[1.0]], dtype=complex)
    for q in range(n):
        out = np.kron(out, PAULI[local_ops.get(q, "I")])
    return out


def pauli_string_matrix(n: int, pauli: PauliString) -> np.ndarray:
    return pauli_on(n, {q: p for p, q in pauli})


def hamiltonian_from_terms(n: int, terms: Sequence[PauliTerm], *, real_output: bool) -> np.ndarray:
    H = np.zeros((1 << n, 1 << n), dtype=complex)
    for term in terms:
        H += float(term.coeff) * pauli_string_matrix(n, term.pauli)
    H = 0.5 * (H + H.conj().T)
    if real_output:
        return np.real_if_close(H, tol=1000).real
    return H


# =============================================================================
# Local n=4 toy Hamiltonians
# =============================================================================


def real_local_terms(n: int = 4) -> List[PauliTerm]:
    """Local real-Hopf toy: X/Z fields plus nearest-neighbor XX/YY/ZZ."""
    if n != 4:
        raise ValueError("The fixed real toy coefficients are for n=4 only.")

    hx = [0.132751, 0.470071, -0.415367, -0.480177]
    hz = [0.285754, 0.159707, 0.125592, 0.177667]
    Jx = [0.996826, -0.960835, 0.679046]
    Jy = [-0.022817, 0.930282, -0.403418]
    Jz = [0.866037, -0.663238, -0.213498]

    terms: List[PauliTerm] = []
    for q in range(n):
        terms.append(PauliTerm(hx[q], (("X", q),)))
        terms.append(PauliTerm(hz[q], (("Z", q),)))

    for q in range(n - 1):
        r = q + 1
        terms.append(PauliTerm(Jx[q], (("X", q), ("X", r))))
        terms.append(PauliTerm(Jy[q], (("Y", q), ("Y", r))))
        terms.append(PauliTerm(Jz[q], (("Z", q), ("Z", r))))

    return terms


def complex_local_terms(n: int = 4) -> List[PauliTerm]:
    """Local complex-Hopf toy: chiral nearest-neighbor spin chain."""
    if n != 4:
        raise ValueError("The fixed complex toy coefficients are for n=4 only.")

    hx = [0.31, -0.27, 0.43, -0.35]
    hy = [0.22, -0.18, 0.29, -0.24]
    hz = [-0.47, 0.39, -0.33, 0.28]
    Jx = [-0.83, 0.71, -0.62]
    Jy = [0.54, -0.49, 0.58]
    Jz = [0.42, -0.37, 0.31]
    D = [0.36, -0.32, 0.27]

    terms: List[PauliTerm] = []
    for q in range(n):
        terms.append(PauliTerm(hx[q], (("X", q),)))
        terms.append(PauliTerm(hy[q], (("Y", q),)))
        terms.append(PauliTerm(hz[q], (("Z", q),)))

    for q in range(n - 1):
        r = q + 1
        terms.append(PauliTerm(Jx[q], (("X", q), ("X", r))))
        terms.append(PauliTerm(Jy[q], (("Y", q), ("Y", r))))
        terms.append(PauliTerm(Jz[q], (("Z", q), ("Z", r))))
        terms.append(PauliTerm(D[q], (("X", q), ("Y", r))))
        terms.append(PauliTerm(-D[q], (("Y", q), ("X", r))))

    return terms


# =============================================================================
# Hopf clipping, initialization, and exact quantities
# =============================================================================


def clip_theta_hopf_complex(theta: np.ndarray, n: int, eps: float = 1e-6) -> np.ndarray:
    """Clip complex-Hopf magnitude angles and wrap phases."""
    theta = np.asarray(theta, dtype=float).copy()
    d = 1 << n
    num_magnitudes = d - 1
    theta[:num_magnitudes] = np.clip(theta[:num_magnitudes], eps, (math.pi / 2.0) - eps)
    theta[num_magnitudes:] = np.mod(theta[num_magnitudes:], 2.0 * math.pi)
    return theta


def clip_for_case(theta: np.ndarray, case: str, n: int) -> np.ndarray:
    if case == "real":
        return clip_theta_hopf_real(theta, n=n)
    if case == "complex":
        return clip_theta_hopf_complex(theta, n=n)
    raise ValueError("case must be 'real' or 'complex'.")


def random_initial_theta(case: str, n: int, seed: int) -> np.ndarray:
    """Fully random normalized state, converted into the corresponding Hopf chart."""
    rng = np.random.default_rng(seed)
    d = 1 << n
    if case == "real":
        psi = rng.normal(size=d)
        psi /= np.linalg.norm(psi)
        # Fix only the irrelevant global sign before converting to coordinates.
        if psi[np.argmax(np.abs(psi))] < 0.0:
            psi *= -1.0
    elif case == "complex":
        psi = rng.normal(size=d) + 1.0j * rng.normal(size=d)
        psi /= np.linalg.norm(psi)
    else:
        raise ValueError("case must be 'real' or 'complex'.")
    return clip_for_case(theta_from_vector(psi, case), case, n)


def energy(theta: np.ndarray, H: np.ndarray, case: str) -> float:
    psi = vector_from_theta(theta, case)
    return float(np.real(np.vdot(psi, H @ psi)))


def exact_hopf_gradient(theta: np.ndarray, H: np.ndarray, case: str) -> np.ndarray:
    psi = vector_from_theta(theta, case)
    Hpsi = H @ psi
    J = jacobian(theta, case)
    return 2.0 * np.real(np.conjugate(J).T @ Hpsi)


def exact_spectrum(H: np.ndarray) -> Tuple[float, float]:
    eigvals = np.linalg.eigvalsh(H)
    return float(np.real(eigvals[0])), float(np.real(eigvals[1] - eigvals[0]))


# =============================================================================
# Layer specifications
# =============================================================================


@dataclass(frozen=True)
class LayerSpec:
    name: str
    kind: str  # "magnitude" or "phase"
    depth: int | None
    index_qubits: int
    param_indices: Tuple[int, ...]  # one-based Hopf parameter indices

    @property
    def size(self) -> int:
        return len(self.param_indices)


def layer_specs(case: str, n: int) -> List[LayerSpec]:
    specs: List[LayerSpec] = []
    for depth in range(n):
        indices = tuple((1 << depth) + r for r in range(1 << depth))
        specs.append(
            LayerSpec(
                name=f"mag-depth-{depth}",
                kind="magnitude",
                depth=depth,
                index_qubits=depth,
                param_indices=indices,
            )
        )

    if case == "complex":
        first_phase = 1 << n  # one-based index of theta_{2^n+0}
        indices = tuple(first_phase + ell for ell in range(1 << n))
        specs.append(
            LayerSpec(
                name="phase-layer",
                kind="phase",
                depth=None,
                index_qubits=n,
                param_indices=indices,
            )
        )

    return specs


# =============================================================================
# Statevector layer states and Pauli readout sampling
# =============================================================================


def readout_unitary_for_pauli(n: int, pauli: PauliString) -> np.ndarray:
    """Unitary applied before computational-basis measurement of a Pauli string."""
    local: Dict[int, np.ndarray] = {}
    for letter, q in pauli:
        if letter == "X":
            local[q] = H1
        elif letter == "Y":
            # Apply S^dagger then H; total state update is H S^dagger.
            local[q] = H1 @ SDG
        elif letter == "Z":
            local[q] = I2
        else:
            raise ValueError(f"Unsupported Pauli letter: {letter!r}")

    U = np.array([[1.0]], dtype=complex)
    for q in range(n):
        U = np.kron(U, local.get(q, I2))
    return U


def pauli_eigenvalues_for_outcomes(n: int, pauli: PauliString) -> np.ndarray:
    """Eigenvalue of a Pauli string after readout rotations for each bitstring."""
    eigs = np.ones(1 << n, dtype=float)
    for outcome in range(1 << n):
        val = 1.0
        for _letter, q in pauli:
            bit = (outcome >> (n - 1 - q)) & 1
            if bit:
                val *= -1.0
        eigs[outcome] = val
    return eigs


def multinomial_counts(rng: np.random.Generator, probs: np.ndarray, shots: int) -> np.ndarray:
    probs = np.asarray(probs, dtype=float).reshape(-1)
    probs = np.maximum(probs, 0.0)
    total = float(probs.sum())
    if total <= 0.0:
        raise ValueError("Probability vector has zero mass.")
    probs = probs / total
    return rng.multinomial(int(shots), probs).reshape(probs.shape)


def sample_baseline_term_statevector(
    theta: np.ndarray,
    pauli: PauliString,
    *,
    case: str,
    n: int,
    shots: int,
    rng: np.random.Generator,
) -> float:
    psi = vector_from_theta(theta, case)
    U = readout_unitary_for_pauli(n, pauli)
    probs = np.abs(U @ psi) ** 2
    counts = multinomial_counts(rng, probs, shots)
    eigs = pauli_eigenvalues_for_outcomes(n, pauli)
    return float(np.dot(counts, eigs) / max(int(shots), 1))


def sample_derivative_layer_term_statevector(
    theta: np.ndarray,
    spec: LayerSpec,
    pauli: PauliString,
    *,
    case: str,
    n: int,
    shots_per_label: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-label Pauli means and per-label shot counts for |Phi_layer>."""
    m = spec.size
    total_shots = max(1, int(shots_per_label) * m)
    U = readout_unitary_for_pauli(n, pauli)
    probs = np.zeros((m, 1 << n), dtype=float)

    for label, param_index in enumerate(spec.param_indices):
        tangent_theta = theta_hopf_tangent_state(theta, param_index, case)
        tangent = vector_from_theta(tangent_theta, case)
        probs[label, :] = (np.abs(U @ tangent) ** 2) / float(m)

    flat_counts = rng.multinomial(total_shots, np.maximum(probs.reshape(-1), 0.0) / probs.sum())
    counts = flat_counts.reshape(m, 1 << n)
    eigs = pauli_eigenvalues_for_outcomes(n, pauli)

    means = np.zeros(m, dtype=float)
    label_counts = counts.sum(axis=1).astype(int)
    for label in range(m):
        if label_counts[label] > 0:
            means[label] = float(np.dot(counts[label], eigs) / label_counts[label])
        else:
            means[label] = 0.0
    return means, label_counts


def sample_branch_layer_term_statevector(
    theta: np.ndarray,
    spec: LayerSpec,
    pauli: PauliString,
    *,
    case: str,
    n: int,
    shots_per_label: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-label/per-sign Pauli means and counts for |Omega_layer>.

    Sign axis 0 is s=+1, sign axis 1 is s=-1.
    """
    m = spec.size
    total_shots = max(1, 2 * int(shots_per_label) * m)
    U = readout_unitary_for_pauli(n, pauli)
    psi = vector_from_theta(theta, case)
    probs = np.zeros((m, 2, 1 << n), dtype=float)

    for label, param_index in enumerate(spec.param_indices):
        tangent_theta = theta_hopf_tangent_state(theta, param_index, case)
        tangent = vector_from_theta(tangent_theta, case)
        plus = psi + tangent
        minus = psi - tangent
        probs[label, 0, :] = np.abs(U @ plus) ** 2 / (4.0 * float(m))
        probs[label, 1, :] = np.abs(U @ minus) ** 2 / (4.0 * float(m))

    flat = probs.reshape(-1)
    flat_counts = rng.multinomial(total_shots, np.maximum(flat, 0.0) / flat.sum())
    counts = flat_counts.reshape(m, 2, 1 << n)
    eigs = pauli_eigenvalues_for_outcomes(n, pauli)

    means = np.zeros((m, 2), dtype=float)
    sign_counts = counts.sum(axis=2).astype(int)
    for label in range(m):
        for sign_index in range(2):
            if sign_counts[label, sign_index] > 0:
                means[label, sign_index] = float(
                    np.dot(counts[label, sign_index], eigs) / sign_counts[label, sign_index]
                )
            else:
                means[label, sign_index] = 0.0
    return means, sign_counts




def get_readout_data(n: int, pauli: PauliString) -> Tuple[np.ndarray, np.ndarray]:
    """Cached readout rotation and eigenvalue table for a Pauli string."""
    key = (n, tuple(pauli))
    cache = getattr(get_readout_data, "_cache", {})
    if key not in cache:
        cache[key] = (readout_unitary_for_pauli(n, pauli), pauli_eigenvalues_for_outcomes(n, pauli))
        setattr(get_readout_data, "_cache", cache)
    return cache[key]


def sample_pauli_mean_from_probs(
    probs: np.ndarray,
    eigs: np.ndarray,
    shots: int,
    rng: np.random.Generator,
) -> float:
    counts = rng.multinomial(max(1, int(shots)), np.maximum(probs, 0.0) / np.sum(probs))
    return float(np.dot(counts, eigs) / max(1, int(shots)))


def sampled_layerwise_gradient_statevector_fast(
    theta: np.ndarray,
    terms: Sequence[PauliTerm],
    *,
    case: str,
    n: int,
    shots_per_label: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Fast statevector sampler for the layerwise signed-energy estimator."""
    theta = np.asarray(theta, dtype=float)
    L = theta.size
    specs = layer_specs(case, n)
    metric = metric_diagonal(theta, case)

    psi = vector_from_theta(theta, case)
    tangent_vectors: List[np.ndarray] = [None] * L  # type: ignore[list-item]
    for spec in specs:
        for param_index in spec.param_indices:
            j = param_index - 1
            if tangent_vectors[j] is None:
                tangent_theta = theta_hopf_tangent_state(theta, param_index, case)
                tangent_vectors[j] = vector_from_theta(tangent_theta, case)

    Epsi = 0.0
    Epartial = np.zeros(L, dtype=float)
    Ebranch_plus = np.zeros(L, dtype=float)
    Ebranch_minus = np.zeros(L, dtype=float)
    counts_partial = np.zeros(L, dtype=int)
    counts_plus = np.zeros(L, dtype=int)
    counts_minus = np.zeros(L, dtype=int)

    for term in terms:
        U, eigs = get_readout_data(n, term.pauli)
        c = float(term.coeff)

        probs_psi = np.abs(U @ psi) ** 2
        Epsi += c * sample_pauli_mean_from_probs(probs_psi, eigs, shots_per_label, rng)

        for spec in specs:
            m = spec.size
            # Tangent-energy layer |Phi>.
            d_probs = np.zeros((m, 1 << n), dtype=float)
            for label, param_index in enumerate(spec.param_indices):
                u = tangent_vectors[param_index - 1]
                d_probs[label, :] = np.abs(U @ u) ** 2 / float(m)
            d_total = max(1, int(shots_per_label) * m)
            d_counts = rng.multinomial(d_total, d_probs.reshape(-1) / d_probs.sum()).reshape(m, 1 << n)
            d_label_counts = d_counts.sum(axis=1).astype(int)

            # Signed branch layer |Omega> followed by X-basis ancilla readout.
            b_probs = np.zeros((m, 2, 1 << n), dtype=float)
            for label, param_index in enumerate(spec.param_indices):
                u = tangent_vectors[param_index - 1]
                b_probs[label, 0, :] = np.abs(U @ (psi + u)) ** 2 / (4.0 * float(m))
                b_probs[label, 1, :] = np.abs(U @ (psi - u)) ** 2 / (4.0 * float(m))
            b_total = max(1, 2 * int(shots_per_label) * m)
            b_counts = rng.multinomial(b_total, b_probs.reshape(-1) / b_probs.sum()).reshape(m, 2, 1 << n)
            b_sign_counts = b_counts.sum(axis=2).astype(int)

            for label, param_index in enumerate(spec.param_indices):
                j = param_index - 1
                if d_label_counts[label] > 0:
                    Epartial[j] += c * float(np.dot(d_counts[label], eigs) / d_label_counts[label])
                counts_partial[j] += int(d_label_counts[label])

                if b_sign_counts[label, 0] > 0:
                    Ebranch_plus[j] += c * float(np.dot(b_counts[label, 0], eigs) / b_sign_counts[label, 0])
                if b_sign_counts[label, 1] > 0:
                    Ebranch_minus[j] += c * float(np.dot(b_counts[label, 1], eigs) / b_sign_counts[label, 1])
                counts_plus[j] += int(b_sign_counts[label, 0])
                counts_minus[j] += int(b_sign_counts[label, 1])

    grad = np.zeros(L, dtype=float)
    for j in range(L):
        sqrt_g = math.sqrt(max(float(metric[j]), 0.0))
        if sqrt_g < 1e-14:
            continue
        center = 0.5 * (Epsi + Epartial[j])
        estimates: List[float] = []
        if counts_plus[j] > 0:
            estimates.append(2.0 * sqrt_g * (Ebranch_plus[j] - center))
        if counts_minus[j] > 0:
            estimates.append(-2.0 * sqrt_g * (Ebranch_minus[j] - center))
        grad[j] = float(np.mean(estimates)) if estimates else 0.0
    return grad


# =============================================================================
# Explicit Qibo indexed-preparation circuits
# =============================================================================


def _require_qibo() -> None:
    if not _QIBO_AVAILABLE:
        raise RuntimeError(f"Qibo is not available: {_QIBO_IMPORT_ERROR!r}")


def _mask_bitpos(mask: int) -> int:
    if mask <= 0 or (mask & (mask - 1)) != 0:
        raise ValueError(f"Expected power-of-two mask, got {mask}.")
    return mask.bit_length() - 1


def _bitpos_to_qubit(bitpos: int, n: int) -> int:
    """Map integer bit position 0=LSB to Hopf-system qubit index 0=MSB."""
    return n - 1 - bitpos


def _mask_to_system_qubits(mask: int, n: int, system_offset: int) -> List[int]:
    qubits: List[int] = []
    for bitpos in range(n):
        if (mask >> bitpos) & 1:
            qubits.append(system_offset + _bitpos_to_qubit(bitpos, n))
    return qubits


def _controlled_gate(circ, gate, controls: Sequence[int]) -> None:
    if controls:
        circ.add(gate.controlled_by(*controls))
    else:
        circ.add(gate)


def add_hopf_block_qibo(
    circ,
    theta: np.ndarray,
    *,
    case: str,
    n: int,
    system_offset: int,
    extra_controls: Sequence[int] = (),
) -> None:
    """Add a Hopf block on system qubits, optionally controlled by extra qubits.

    This follows the same Ctrl/Anti/Targ/Index schedule used by hopf_utils.
    """
    _require_qibo()
    ctrl_list, anti_list, targ_list, idx_list = gates_order(n, case)

    for ctrl_mask, anti_mask, targ_mask, idx in zip(ctrl_list, anti_list, targ_list, idx_list):
        targ_mask = int(targ_mask)
        t_bitpos = _mask_bitpos(targ_mask)
        target = system_offset + _bitpos_to_qubit(t_bitpos, n)

        ctrl_mask = int(ctrl_mask) & (~targ_mask)
        anti_mask = int(anti_mask) & (~targ_mask)
        controls_mask = ctrl_mask | anti_mask

        system_controls = _mask_to_system_qubits(controls_mask, n, system_offset)
        system_anticontrols = _mask_to_system_qubits(anti_mask, n, system_offset)
        controls = list(extra_controls) + system_controls

        for q in system_anticontrols:
            circ.add(gates.X(q))

        if case == "real":
            theta_gate = float(theta[int(idx) - 1])
            _controlled_gate(circ, gates.RY(target, 2.0 * theta_gate), controls)
        elif case == "complex":
            if isinstance(idx, (list, tuple, np.ndarray)) and len(idx) == 3:
                a = float(theta[int(idx[0]) - 1])
                b = float(theta[int(idx[1]) - 1])
                c = float(theta[int(idx[2]) - 1])
                U = np.array(
                    [
                        [np.exp(1.0j * b) * np.cos(a), -np.exp(-1.0j * c) * np.sin(a)],
                        [np.exp(1.0j * c) * np.sin(a), np.exp(-1.0j * b) * np.cos(a)],
                    ],
                    dtype=complex,
                )
                _controlled_gate(circ, gates.Unitary(U, target), controls)
            else:
                theta_gate = float(theta[int(idx) - 1])
                _controlled_gate(circ, gates.RY(target, 2.0 * theta_gate), controls)
        else:
            raise ValueError("case must be 'real' or 'complex'.")

        for q in system_anticontrols:
            circ.add(gates.X(q))


def _index_bits(label: int, width: int) -> List[int]:
    return [(label >> (width - 1 - pos)) & 1 for pos in range(width)]


def _apply_label_anticontrol_x(circ, index_qubits: Sequence[int], label: int) -> List[int]:
    flipped: List[int] = []
    bits = _index_bits(label, len(index_qubits))
    for q, bit in zip(index_qubits, bits):
        if bit == 0:
            circ.add(gates.X(q))
            flipped.append(q)
    return flipped


def _undo_flips(circ, flipped: Sequence[int]) -> None:
    for q in reversed(list(flipped)):
        circ.add(gates.X(q))


def add_pauli_readout_rotations_qibo(circ, pauli: PauliString, *, system_offset: int) -> None:
    _require_qibo()
    for letter, q in pauli:
        target = system_offset + q
        if letter == "X":
            circ.add(gates.H(target))
        elif letter == "Y":
            circ.add(gates.Unitary(SDG, target))
            circ.add(gates.H(target))
        elif letter == "Z":
            pass
        else:
            raise ValueError(f"Unsupported Pauli letter: {letter!r}")


def _samples_array(result) -> np.ndarray:
    samples = np.asarray(result.samples())
    if samples.ndim != 2:
        samples = samples.reshape((samples.shape[0], -1))
    return samples.astype(int, copy=False)


def _labels_from_bits(bits: np.ndarray) -> np.ndarray:
    if bits.size == 0:
        return np.zeros(bits.shape[0], dtype=int)
    labels = np.zeros(bits.shape[0], dtype=int)
    width = bits.shape[1]
    for pos in range(width):
        labels = (labels << 1) | bits[:, pos].astype(int)
    return labels


def _pauli_values_from_system_bits(system_bits: np.ndarray, pauli: PauliString, n: int) -> np.ndarray:
    values = np.ones(system_bits.shape[0], dtype=float)
    for _letter, q in pauli:
        values *= np.where(system_bits[:, q] == 0, 1.0, -1.0)
    return values


def sample_baseline_term_qibo(
    theta: np.ndarray,
    pauli: PauliString,
    *,
    case: str,
    n: int,
    shots: int,
) -> float:
    _require_qibo()
    circ = models.Circuit(n)
    add_hopf_block_qibo(circ, theta, case=case, n=n, system_offset=0)
    add_pauli_readout_rotations_qibo(circ, pauli, system_offset=0)
    circ.add(gates.M(*range(n)))
    result = circ(nshots=int(shots))
    samples = _samples_array(result)
    values = _pauli_values_from_system_bits(samples[:, :n], pauli, n)
    return float(values.mean())


def sample_derivative_layer_term_qibo(
    theta: np.ndarray,
    spec: LayerSpec,
    pauli: PauliString,
    *,
    case: str,
    n: int,
    shots_per_label: int,
) -> Tuple[np.ndarray, np.ndarray]:
    _require_qibo()
    num_index = spec.index_qubits
    m = spec.size
    total_shots = max(1, int(shots_per_label) * m)
    index_qubits = list(range(num_index))
    system_offset = num_index
    system_qubits = list(range(system_offset, system_offset + n))

    circ = models.Circuit(num_index + n)
    for q in index_qubits:
        circ.add(gates.H(q))

    for label, param_index in enumerate(spec.param_indices):
        tangent_theta = theta_hopf_tangent_state(theta, param_index, case)
        flipped = _apply_label_anticontrol_x(circ, index_qubits, label)
        add_hopf_block_qibo(
            circ,
            tangent_theta,
            case=case,
            n=n,
            system_offset=system_offset,
            extra_controls=index_qubits,
        )
        _undo_flips(circ, flipped)

    add_pauli_readout_rotations_qibo(circ, pauli, system_offset=system_offset)
    circ.add(gates.M(*(index_qubits + system_qubits)))
    result = circ(nshots=total_shots)
    samples = _samples_array(result)

    if num_index > 0:
        labels = _labels_from_bits(samples[:, :num_index])
    else:
        labels = np.zeros(samples.shape[0], dtype=int)
    system_bits = samples[:, num_index : num_index + n]
    values = _pauli_values_from_system_bits(system_bits, pauli, n)

    means = np.zeros(m, dtype=float)
    counts = np.zeros(m, dtype=int)
    for label in range(m):
        mask = labels == label
        counts[label] = int(np.count_nonzero(mask))
        if counts[label] > 0:
            means[label] = float(values[mask].mean())
    return means, counts


def sample_branch_layer_term_qibo(
    theta: np.ndarray,
    spec: LayerSpec,
    pauli: PauliString,
    *,
    case: str,
    n: int,
    shots_per_label: int,
) -> Tuple[np.ndarray, np.ndarray]:
    _require_qibo()
    num_index = spec.index_qubits
    m = spec.size
    total_shots = max(1, 2 * int(shots_per_label) * m)
    index_qubits = list(range(num_index))
    anc = num_index
    system_offset = num_index + 1
    system_qubits = list(range(system_offset, system_offset + n))

    circ = models.Circuit(num_index + 1 + n)
    for q in index_qubits:
        circ.add(gates.H(q))
    circ.add(gates.H(anc))

    # Baseline branch on ancilla |0>.
    circ.add(gates.X(anc))
    add_hopf_block_qibo(circ, theta, case=case, n=n, system_offset=system_offset, extra_controls=[anc])
    circ.add(gates.X(anc))

    # Tangent branch on ancilla |1>, indexed by the layer register.
    for label, param_index in enumerate(spec.param_indices):
        tangent_theta = theta_hopf_tangent_state(theta, param_index, case)
        flipped = _apply_label_anticontrol_x(circ, index_qubits, label)
        add_hopf_block_qibo(
            circ,
            tangent_theta,
            case=case,
            n=n,
            system_offset=system_offset,
            extra_controls=[anc] + index_qubits,
        )
        _undo_flips(circ, flipped)

    # X-basis ancilla measurement gives branch sign s=(-1)^b.
    circ.add(gates.H(anc))
    add_pauli_readout_rotations_qibo(circ, pauli, system_offset=system_offset)
    circ.add(gates.M(*(index_qubits + [anc] + system_qubits)))
    result = circ(nshots=total_shots)
    samples = _samples_array(result)

    if num_index > 0:
        labels = _labels_from_bits(samples[:, :num_index])
    else:
        labels = np.zeros(samples.shape[0], dtype=int)
    anc_bits = samples[:, num_index]
    system_bits = samples[:, num_index + 1 : num_index + 1 + n]
    values = _pauli_values_from_system_bits(system_bits, pauli, n)

    means = np.zeros((m, 2), dtype=float)
    counts = np.zeros((m, 2), dtype=int)
    for label in range(m):
        for sign_index, anc_value in enumerate((0, 1)):
            mask = (labels == label) & (anc_bits == anc_value)
            counts[label, sign_index] = int(np.count_nonzero(mask))
            if counts[label, sign_index] > 0:
                means[label, sign_index] = float(values[mask].mean())
    return means, counts


# =============================================================================
# Layerwise sampled gradient
# =============================================================================


def resolve_sampler(sampler: str) -> str:
    if sampler == "auto":
        return "qibo-explicit" if _QIBO_AVAILABLE else "statevector"
    if sampler in {"qibo-explicit", "statevector"}:
        if sampler == "qibo-explicit":
            _require_qibo()
        return sampler
    raise ValueError("sampler must be 'auto', 'qibo-explicit', or 'statevector'.")


def sampled_layerwise_gradient(
    theta: np.ndarray,
    terms: Sequence[PauliTerm],
    *,
    case: str,
    n: int,
    shots_per_label: int,
    rng: np.random.Generator,
    sampler: str,
) -> np.ndarray:
    """Estimate the Euclidean Hopf-coordinate gradient using layerwise circuits."""
    sampler = resolve_sampler(sampler)
    if sampler == "statevector":
        return sampled_layerwise_gradient_statevector_fast(
            theta,
            terms,
            case=case,
            n=n,
            shots_per_label=shots_per_label,
            rng=rng,
        )

    theta = np.asarray(theta, dtype=float)
    L = theta.size
    metric = metric_diagonal(theta, case)

    Epsi = 0.0
    Epartial = np.zeros(L, dtype=float)
    Ebranch_plus = np.zeros(L, dtype=float)
    Ebranch_minus = np.zeros(L, dtype=float)
    counts_partial = np.zeros(L, dtype=int)
    counts_plus = np.zeros(L, dtype=int)
    counts_minus = np.zeros(L, dtype=int)

    # Baseline energy: one readout setting per Pauli term.
    for term in terms:
        if sampler == "qibo-explicit":
            mean = sample_baseline_term_qibo(theta, term.pauli, case=case, n=n, shots=shots_per_label)
        else:
            mean = sample_baseline_term_statevector(
                theta,
                term.pauli,
                case=case,
                n=n,
                shots=shots_per_label,
                rng=rng,
            )
        Epsi += float(term.coeff) * mean

    # Layerwise tangent and branch energies.
    for spec in layer_specs(case, n):
        for term in terms:
            if sampler == "qibo-explicit":
                d_means, d_counts = sample_derivative_layer_term_qibo(
                    theta,
                    spec,
                    term.pauli,
                    case=case,
                    n=n,
                    shots_per_label=shots_per_label,
                )
                b_means, b_counts = sample_branch_layer_term_qibo(
                    theta,
                    spec,
                    term.pauli,
                    case=case,
                    n=n,
                    shots_per_label=shots_per_label,
                )
            else:
                d_means, d_counts = sample_derivative_layer_term_statevector(
                    theta,
                    spec,
                    term.pauli,
                    case=case,
                    n=n,
                    shots_per_label=shots_per_label,
                    rng=rng,
                )
                b_means, b_counts = sample_branch_layer_term_statevector(
                    theta,
                    spec,
                    term.pauli,
                    case=case,
                    n=n,
                    shots_per_label=shots_per_label,
                    rng=rng,
                )

            for label, param_index in enumerate(spec.param_indices):
                j = param_index - 1
                c = float(term.coeff)
                Epartial[j] += c * d_means[label]
                Ebranch_plus[j] += c * b_means[label, 0]
                Ebranch_minus[j] += c * b_means[label, 1]
                counts_partial[j] += int(d_counts[label])
                counts_plus[j] += int(b_counts[label, 0])
                counts_minus[j] += int(b_counts[label, 1])

    grad = np.zeros(L, dtype=float)
    for j in range(L):
        sqrt_g = math.sqrt(max(float(metric[j]), 0.0))
        if sqrt_g < 1e-14:
            grad[j] = 0.0
            continue

        center = 0.5 * (Epsi + Epartial[j])
        estimates: List[float] = []
        if counts_plus[j] > 0:
            estimates.append(2.0 * sqrt_g * (Ebranch_plus[j] - center))
        if counts_minus[j] > 0:
            estimates.append(-2.0 * sqrt_g * (Ebranch_minus[j] - center))
        grad[j] = float(np.mean(estimates)) if estimates else 0.0

    return grad


def exact_layerwise_gradient(theta: np.ndarray, H: np.ndarray, case: str, n: int) -> np.ndarray:
    """Deterministic check of the signed-energy layerwise identity."""
    psi = vector_from_theta(theta, case)
    Epsi = float(np.real(np.vdot(psi, H @ psi)))
    metric = metric_diagonal(theta, case)
    grad = np.zeros(theta.size, dtype=float)

    for spec in layer_specs(case, n):
        for param_index in spec.param_indices:
            j = param_index - 1
            sqrt_g = math.sqrt(max(float(metric[j]), 0.0))
            if sqrt_g < 1e-14:
                continue
            tangent_theta = theta_hopf_tangent_state(theta, param_index, case)
            tangent = vector_from_theta(tangent_theta, case)
            Epartial = float(np.real(np.vdot(tangent, H @ tangent)))
            plus = (psi + tangent) / math.sqrt(2.0)
            minus = (psi - tangent) / math.sqrt(2.0)
            Eplus = float(np.real(np.vdot(plus, H @ plus)))
            Eminus = float(np.real(np.vdot(minus, H @ minus)))
            g_plus = 2.0 * sqrt_g * (Eplus - 0.5 * (Epsi + Epartial))
            g_minus = -2.0 * sqrt_g * (Eminus - 0.5 * (Epsi + Epartial))
            grad[j] = 0.5 * (g_plus + g_minus)
    return grad


# =============================================================================
# Optimizer and VQE loop
# =============================================================================


class Adam:
    def __init__(self, dim: int, lr: float = 0.05, beta1: float = 0.9, beta2: float = 0.999, eps: float = 1e-8):
        self.lr = float(lr)
        self.beta1 = float(beta1)
        self.beta2 = float(beta2)
        self.eps = float(eps)
        self.m = np.zeros(dim, dtype=float)
        self.v = np.zeros(dim, dtype=float)
        self.t = 0

    def step(self, x: np.ndarray, g: np.ndarray) -> np.ndarray:
        self.t += 1
        self.m = self.beta1 * self.m + (1.0 - self.beta1) * g
        self.v = self.beta2 * self.v + (1.0 - self.beta2) * (g * g)
        m_hat = self.m / (1.0 - self.beta1**self.t)
        v_hat = self.v / (1.0 - self.beta2**self.t)
        return x - self.lr * m_hat / (np.sqrt(v_hat) + self.eps)


@dataclass
class RunResult:
    case: str
    title: str
    exact_energies: np.ndarray
    layerwise_energies: np.ndarray | None
    ground_energy: float
    spectral_gap: float
    num_params: int
    num_terms: int
    num_layer_settings: int


def run_vqe_pair(
    *,
    title: str,
    case: str,
    H: np.ndarray,
    terms: Sequence[PauliTerm],
    n: int,
    steps: int,
    lr: float,
    shots_per_label: int,
    seed: int,
    exact_only: bool,
    sampler: str,
    log_every: int,
) -> RunResult:
    E0, gap = exact_spectrum(H)
    theta0 = random_initial_theta(case, n, seed)

    theta_exact = theta0.copy()
    theta_layer = theta0.copy()
    opt_exact = Adam(theta0.size, lr=lr)
    opt_layer = Adam(theta0.size, lr=lr)
    rng = np.random.default_rng(seed + 1000003)

    exact_energies: List[float] = []
    layerwise_energies: List[float] = []

    specs = layer_specs(case, n)
    num_layer_settings = 1 + 2 * len(specs)
    active_sampler = resolve_sampler(sampler) if not exact_only else "skipped"

    # Sanity check the deterministic layerwise identity at the initial point.
    g_direct = exact_hopf_gradient(theta0, H, case)
    g_layer_exact = exact_layerwise_gradient(theta0, H, case, n)
    identity_residual = float(np.max(np.abs(g_direct - g_layer_exact)))

    print()
    print(f"[{title}]")
    print(f"case = {case}")
    print(f"n = {n}")
    print(f"parameters = {theta0.size}")
    print(f"Pauli terms = {len(terms)}")
    print(f"layerwise settings per Pauli term = {num_layer_settings}")
    print(f"sampler = {active_sampler}")
    print(f"Adam learning rate = {lr}")
    print(f"ground energy = {E0:.12f}")
    print(f"spectral gap = {gap:.12f}")
    print(f"initial layerwise identity residual = {identity_residual:.3e}")
    if not exact_only:
        print(f"shots per label/sign/readout = {shots_per_label}")

    for step in range(steps):
        E_exact = energy(theta_exact, H, case)
        exact_energies.append(E_exact)
        g_exact = exact_hopf_gradient(theta_exact, H, case)
        theta_exact = clip_for_case(opt_exact.step(theta_exact, g_exact), case, n)

        if not exact_only:
            E_layer = energy(theta_layer, H, case)
            layerwise_energies.append(E_layer)
            g_layer = sampled_layerwise_gradient(
                theta_layer,
                terms,
                case=case,
                n=n,
                shots_per_label=shots_per_label,
                rng=rng,
                sampler=sampler,
            )
            theta_layer = clip_for_case(opt_layer.step(theta_layer, g_layer), case, n)

        if log_every > 0 and (step == 0 or (step + 1) % log_every == 0):
            if exact_only:
                print(f"iter {step + 1:4d}: exact E={E_exact:.9f}, gap={E_exact - E0:.3e}")
            else:
                print(
                    f"iter {step + 1:4d}: "
                    f"exact E={E_exact:.9f}, gap={E_exact - E0:.3e}; "
                    f"layerwise E={E_layer:.9f}, gap={E_layer - E0:.3e}"
                )

    exact_energies.append(energy(theta_exact, H, case))
    if not exact_only:
        layerwise_energies.append(energy(theta_layer, H, case))

    return RunResult(
        case=case,
        title=title,
        exact_energies=np.asarray(exact_energies, dtype=float),
        layerwise_energies=None if exact_only else np.asarray(layerwise_energies, dtype=float),
        ground_energy=E0,
        spectral_gap=gap,
        num_params=theta0.size,
        num_terms=len(terms),
        num_layer_settings=num_layer_settings,
    )


# =============================================================================
# Plotting
# =============================================================================


def plot_results(results: Sequence[RunResult], output: str, shots_per_label: int, plot_gap: bool) -> None:
    fig, axes = plt.subplots(1, len(results), figsize=(12.0, 4.8), sharey=False)
    if len(results) == 1:
        axes = [axes]

    for ax, result in zip(axes, results):
        x = np.arange(result.exact_energies.size)
        if plot_gap:
            y_exact = np.maximum(result.exact_energies - result.ground_energy, 1e-14)
            ax.semilogy(x, y_exact, label="Exact Hopf gradient", linewidth=2.0)
            if result.layerwise_energies is not None:
                y_layer = np.maximum(result.layerwise_energies - result.ground_energy, 1e-14)
                ax.semilogy(
                    x,
                    y_layer,
                    label=f"Layerwise circuits, Ns={shots_per_label}",
                    linestyle=":",
                    linewidth=2.0,
                )
            ax.set_ylabel(r"Energy gap $E(\theta)-E_0$")
        else:
            ax.plot(x, result.exact_energies, label="Exact Hopf gradient", linewidth=2.0)
            if result.layerwise_energies is not None:
                ax.plot(
                    x,
                    result.layerwise_energies,
                    label=f"Layerwise circuits, Ns={shots_per_label}",
                    linestyle=":",
                    linewidth=2.0,
                )
            ax.axhline(result.ground_energy, linestyle="--", color='gray', label="Exact ground energy")
            ax.set_ylabel("Energy")

        ax.set_xlabel("Iteration")
        ax.set_title(result.title)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)

    fig.suptitle("n=4 local Hopf VQE: exact gradient vs layerwise sampled estimator")
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print()
    print(f"wrote {output}")


# =============================================================================
# Main
# =============================================================================


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=100, help="Optimizer steps per trajectory.")
    parser.add_argument(
        "--shots",
        "--Ns",
        dest="shots_per_label",
        type=int,
        default=50,
        help="Shots per label/sign/readout for layerwise settings.",
    )
    parser.add_argument("--lr", type=float, default=None, help="Optional fixed Adam learning rate for both toys. Overrides --lr-real and --lr-complex.")
    parser.add_argument("--lr-real", type=float, default=0.15, help="Fixed Adam learning rate for the real-Hopf toy. No line search is used.")
    parser.add_argument("--lr-complex", type=float, default=0.12, help="Fixed Adam learning rate for the complex-Hopf toy. No line search is used.")
    parser.add_argument("--seed", type=int, default=None, help="Optional base seed. Overrides --seed-real and --seed-complex when provided.")
    parser.add_argument("--seed-real", type=int, default=2026, help="Random initialization seed for the real-Hopf toy.")
    parser.add_argument("--seed-complex", type=int, default=2033, help="Random initialization seed for the complex-Hopf toy.")
    parser.add_argument(
        "--sampler",
        choices=("auto", "qibo-explicit", "statevector"),
        default="auto",
        help="Layerwise sampler. 'auto' uses Qibo if available, otherwise statevector fallback.",
    )
    parser.add_argument("--output", default="VQE_qibo.png", help="Output plot path. The file extension controls the Matplotlib output format.")
    parser.add_argument("--exact-only", action="store_true", help="Skip the sampled layerwise trajectories.")
    parser.add_argument("--plot-gap", action="store_true", help="Plot E(theta)-E0 on a logarithmic y-axis instead of raw energies.")
    parser.add_argument("--log-every", type=int, default=10, help="Console logging period. Set to 0 to suppress per-iteration logs.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    lr_real = args.lr_real if args.lr is None else args.lr
    lr_complex = args.lr_complex if args.lr is None else args.lr
    seed_real = args.seed_real if args.seed is None else args.seed
    seed_complex = args.seed_complex if args.seed is None else args.seed + 1

    n = 4
    real_terms = real_local_terms(n)
    complex_terms = complex_local_terms(n)
    H_real = hamiltonian_from_terms(n, real_terms, real_output=True)
    H_complex = hamiltonian_from_terms(n, complex_terms, real_output=False)

    real_result = run_vqe_pair(
        title="Real Hopf: local real chain",
        case="real",
        H=H_real,
        terms=real_terms,
        n=n,
        steps=args.steps,
        lr=lr_real,
        shots_per_label=args.shots_per_label,
        seed=seed_real,
        exact_only=args.exact_only,
        sampler=args.sampler,
        log_every=args.log_every,
    )

    complex_result = run_vqe_pair(
        title="Complex Hopf: local chiral chain",
        case="complex",
        H=H_complex,
        terms=complex_terms,
        n=n,
        steps=args.steps,
        lr=lr_complex,
        shots_per_label=args.shots_per_label,
        seed=seed_complex,
        exact_only=args.exact_only,
        sampler=args.sampler,
        log_every=args.log_every,
    )

    plot_results([real_result, complex_result], args.output, args.shots_per_label, args.plot_gap)


if __name__ == "__main__":
    main()
