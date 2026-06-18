#!/usr/bin/env python3
"""
hopf_complex.py

Complex-Hopf VQE and metrology stress test, diagnostics, and paper-facing plot.

This script merges the complex stress-test generator and its plotting code into
one release-facing entry point.  It repeats the six real-case task definitions
with a complex scrambling circuit, uses ten deterministic complex initial states
per task by default, and compares four tracks:

    1. Hopf-Adam
    2. Hopf-EGT-CG
    3. Hopf-Riemannian-BB
    4. Hopf-Riemannian-LBFGS

No complex Möttönen baseline is implemented.

Default run:
    python hopf_complex.py

Default outputs:
    complex_hopf_stress_data.csv
    complex_stress_summary.pdf

The single figure has the same meaning, mode order, plotting conventions, and
color code as the upper (final-gap) row of hopf_geometric_summary_clean.pdf.
The script depends on hopf_utils.py in the same directory.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics as stats
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.optimize import minimize_scalar
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This script requires scipy. Install it with: pip install scipy") from exc

import hopf_utils as hopf


# -----------------------------------------------------------------------------
# Plot conventions: identical to the upper half of plot_hopf.py, with the
# unavailable complex Möttönen track removed.
# -----------------------------------------------------------------------------

PLOT_MODES = [
    "Hopf-Adam",
    "Hopf-EGT-CG",
    "Hopf-Riemannian-BB",
    "Hopf-Riemannian-LBFGS",
]

MODE_LABEL = {
    "Hopf-Adam": "Hopf Adam",
    "Hopf-EGT-CG": "EGT-CG",
    "Hopf-Riemannian-BB": "R-BB",
    "Hopf-Riemannian-LBFGS": "R-LBFGS",
}

MODE_COLOR = {
    "Hopf-Adam": "#eb6426",
    "Hopf-EGT-CG": "#1f4e79",
    "Hopf-Riemannian-BB": "#4b5fa8",
    "Hopf-Riemannian-LBFGS": "#2a7f8f",
}

MODE_ABBREVIATION = {
    "Hopf-Adam": "HC-Adam",
    "Hopf-EGT-CG": "HC-EGT-CG",
    "Hopf-Riemannian-BB": "HC-R-BB",
    "Hopf-Riemannian-LBFGS": "HC-R-LBFGS",
}

TASK_ORDER = {
    "VQE": ["VQE-1", "VQE-2", "VQE-3"],
    "MET": ["MET-1", "MET-2", "MET-3"],
}

CSV_FIELDS = [
    # Original supplementary CSV compatibility fields.
    "n",
    "app_class",
    "task",
    "seed",
    "optimizer",
    "step",
    "gap",
    # Richer release diagnostics.
    "task_id",
    "problem_seed",
    "seed_index",
    "run_seed",
    "mode",
    "cost",
    "grad_norm",
    "state_grad_norm",
    "state_norm_error",
    "last_step_angle",
    "last_line_evals",
    "wall_time_sec",
]


# -----------------------------------------------------------------------------
# Complex sphere utilities
# -----------------------------------------------------------------------------


def real_inner(x: np.ndarray, y: np.ndarray) -> float:
    """Ambient real pairing Re<x|y> used for the complex state-vector sphere."""
    return float(np.real(np.vdot(np.asarray(x, dtype=complex), np.asarray(y, dtype=complex))))


def unit_normalize(x: np.ndarray, *, eps: float = 1e-15) -> np.ndarray:
    nrm = float(np.linalg.norm(x))
    if nrm < eps:
        raise ValueError("Cannot normalize a near-zero vector.")
    return np.asarray(x, dtype=complex) / nrm


def normalize_with_phase(x: np.ndarray, *, eps: float = 1e-15) -> np.ndarray:
    """Normalize and choose a deterministic global-phase representative."""
    y = unit_normalize(x, eps=eps)
    idx = int(np.argmax(np.abs(y)))
    if abs(y[idx]) > eps:
        y = y * np.exp(-1j * float(np.angle(y[idx])))
    return y


def tangent_project(x: np.ndarray, v: np.ndarray) -> np.ndarray:
    x = unit_normalize(x)
    return np.asarray(v, dtype=complex) - real_inner(x, v) * x


def sphere_exp(x: np.ndarray, direction: np.ndarray, alpha: float) -> np.ndarray:
    """Exact exponential map on the ambient real unit sphere in C^N."""
    x = unit_normalize(x)
    p = tangent_project(x, direction)
    pnorm = float(np.linalg.norm(p))
    if pnorm < 1e-15 or abs(float(alpha)) < 1e-15:
        return x.copy()
    ell = float(alpha) * pnorm
    return unit_normalize(math.cos(ell) * x + math.sin(ell) * (p / pnorm))


def sphere_exp_unit(x: np.ndarray, unit_direction: np.ndarray, angle: float) -> np.ndarray:
    x = unit_normalize(x)
    q = tangent_project(x, unit_direction)
    qnorm = float(np.linalg.norm(q))
    if qnorm < 1e-15 or abs(float(angle)) < 1e-15:
        return x.copy()
    q = q / qnorm
    return unit_normalize(math.cos(float(angle)) * x + math.sin(float(angle)) * q)


def transport_exact(
    x_old: np.ndarray,
    x_new: np.ndarray,
    vector: np.ndarray,
    *,
    eps: float = 1e-12,
) -> np.ndarray:
    """Exact parallel transport along the short great-circle geodesic."""
    x_old = unit_normalize(x_old)
    x_new = unit_normalize(x_new)
    v = tangent_project(x_old, vector)
    denom = 1.0 + real_inner(x_old, x_new)
    if abs(denom) < eps:
        return tangent_project(x_new, v)
    transported = v - (real_inner(x_new, v) / denom) * (x_old + x_new)
    return tangent_project(x_new, transported)


def state_angle(x: np.ndarray, y: np.ndarray) -> float:
    dot = max(-1.0, min(1.0, real_inner(unit_normalize(x), unit_normalize(y))))
    return float(math.acos(dot))


# -----------------------------------------------------------------------------
# Complex scrambling circuit
# -----------------------------------------------------------------------------


@dataclass
class ComplexScrambler:
    """Brickwork R_y/CNOT scrambler augmented by deterministic phase rotations."""

    n: int
    depth: int
    seed: int
    angles: np.ndarray = field(init=False)
    phases: np.ndarray = field(init=False)
    pairs_by_layer: List[List[Tuple[int, int]]] = field(init=False)
    _phase_factors: np.ndarray = field(init=False, repr=False)
    _ry_pairs: List[Tuple[np.ndarray, np.ndarray]] = field(init=False, repr=False)
    _cnot_perms: List[List[np.ndarray]] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.n < 1:
            raise ValueError("n must be >= 1")
        if self.depth < 0:
            raise ValueError("depth must be >= 0")

        rng = np.random.default_rng(self.seed)
        # With the same seed, these R_y angles match the real scrambler before
        # the additional complex phases are drawn.
        self.angles = rng.uniform(-math.pi, math.pi, size=(self.depth, self.n))
        self.phases = rng.uniform(0.0, 2.0 * math.pi, size=(self.depth, self.n))

        self.pairs_by_layer = []
        for layer in range(self.depth):
            pairs: List[Tuple[int, int]] = []
            if self.n >= 2:
                if layer % 2 == 0:
                    for q in range(0, self.n - 1, 2):
                        pairs.append((q, q + 1))
                else:
                    for q in range(1, self.n - 1, 2):
                        pairs.append((q + 1, q))
            self.pairs_by_layer.append(pairs)

        dim = 1 << self.n
        basis = np.arange(dim, dtype=np.int64)
        shifts = np.arange(self.n - 1, -1, -1, dtype=np.int64)
        bits = ((basis[:, None] >> shifts[None, :]) & 1).astype(float)
        self._phase_factors = np.exp(1j * (self.phases @ bits.T))

        self._ry_pairs = []
        for q in range(self.n):
            mask = 1 << (self.n - 1 - q)
            idx0 = basis[(basis & mask) == 0]
            idx1 = idx0 | mask
            self._ry_pairs.append((idx0, idx1))

        self._cnot_perms = []
        for pairs in self.pairs_by_layer:
            layer_perms: List[np.ndarray] = []
            for control, target in pairs:
                cmask = 1 << (self.n - 1 - control)
                tmask = 1 << (self.n - 1 - target)
                perm = basis.copy()
                active = (perm & cmask) != 0
                perm[active] ^= tmask
                layer_perms.append(perm)
            self._cnot_perms.append(layer_perms)

    def apply(self, x: np.ndarray) -> np.ndarray:
        y = np.asarray(x, dtype=complex).copy()
        if y.shape != (1 << self.n,):
            raise ValueError(f"Expected state shape {(1 << self.n,)}, got {y.shape}.")
        for layer in range(self.depth):
            y *= self._phase_factors[layer]
            for q in range(self.n):
                idx0, idx1 = self._ry_pairs[q]
                a = y[idx0].copy()
                b = y[idx1].copy()
                c = math.cos(float(self.angles[layer, q]))
                s = math.sin(float(self.angles[layer, q]))
                y[idx0] = c * a - s * b
                y[idx1] = s * a + c * b
            for perm in self._cnot_perms[layer]:
                y = y[perm]
        return y

    def apply_T(self, x: np.ndarray) -> np.ndarray:
        y = np.asarray(x, dtype=complex).copy()
        if y.shape != (1 << self.n,):
            raise ValueError(f"Expected state shape {(1 << self.n,)}, got {y.shape}.")
        for layer in reversed(range(self.depth)):
            for perm in reversed(self._cnot_perms[layer]):
                y = y[perm]
            for q in reversed(range(self.n)):
                idx0, idx1 = self._ry_pairs[q]
                a = y[idx0].copy()
                b = y[idx1].copy()
                c = math.cos(float(self.angles[layer, q]))
                s = math.sin(float(self.angles[layer, q]))
                y[idx0] = c * a + s * b
                y[idx1] = -s * a + c * b
            y *= np.conjugate(self._phase_factors[layer])
        return y

    def scrambled_basis(self, label: int) -> np.ndarray:
        e = np.zeros(1 << self.n, dtype=complex)
        e[int(label)] = 1.0
        return normalize_with_phase(self.apply(e))


# -----------------------------------------------------------------------------
# Synthetic tasks: exact complex analogues of the real-case stress tests
# -----------------------------------------------------------------------------


@dataclass
class BaseTask:
    app_class: str
    task_id: str
    name: str
    n: int
    problem_seed: int
    scrambler: ComplexScrambler

    def cost(self, psi: np.ndarray) -> float:
        raise NotImplementedError

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        """Gradient under dC = Re<grad|dpsi>."""
        raise NotImplementedError

    def gap(self, psi: np.ndarray) -> float:
        return float(self.cost(psi))


@dataclass
class VQEParentTask(BaseTask):
    target_label: int
    target: np.ndarray

    def cost(self, psi: np.ndarray) -> float:
        overlap = np.vdot(self.target, psi)
        return float(1.0 - abs(overlap) ** 2)

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        overlap = np.vdot(self.target, psi)
        return -2.0 * overlap * self.target


@dataclass
class VQEDiagonalTask(BaseTask):
    target_label: int
    target: np.ndarray
    spectrum: np.ndarray
    distractor_label: Optional[int] = None

    def apply_H(self, psi: np.ndarray) -> np.ndarray:
        raw = self.scrambler.apply_T(psi)
        return self.scrambler.apply(self.spectrum * raw)

    def cost(self, psi: np.ndarray) -> float:
        raw = self.scrambler.apply_T(psi)
        return float(np.dot(self.spectrum, np.abs(raw) ** 2))

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        return 2.0 * self.apply_H(psi)


@dataclass
class METSingleTargetTask(BaseTask):
    target_label: int
    target: np.ndarray

    def expectation_A(self, psi: np.ndarray) -> float:
        return float(abs(np.vdot(self.target, psi)) ** 2)

    def cost(self, psi: np.ndarray) -> float:
        a = self.expectation_A(psi)
        return 1.0 - a * a

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        overlap = np.vdot(self.target, psi)
        a = float(abs(overlap) ** 2)
        return -4.0 * a * overlap * self.target


@dataclass
class METQFITask(BaseTask):
    generator: np.ndarray
    gmin_label: int
    gmax_label: int
    span: float
    tau_plus: np.ndarray
    tau_minus: np.ndarray

    def _Gy(self, psi: np.ndarray, power: int) -> np.ndarray:
        raw = self.scrambler.apply_T(psi)
        return self.scrambler.apply((self.generator ** power) * raw)

    def moments(self, psi: np.ndarray) -> Tuple[float, float]:
        raw = self.scrambler.apply_T(psi)
        prob = np.abs(raw) ** 2
        mu = float(np.dot(self.generator, prob))
        nu = float(np.dot(self.generator ** 2, prob))
        return mu, nu

    def normalized_qfi(self, psi: np.ndarray) -> float:
        mu, nu = self.moments(psi)
        return float(4.0 * max(0.0, nu - mu * mu) / (self.span * self.span))

    def cost(self, psi: np.ndarray) -> float:
        return 1.0 - self.normalized_qfi(psi)

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        mu, _nu = self.moments(psi)
        grad_variance = 2.0 * self._Gy(psi, 2) - 4.0 * mu * self._Gy(psi, 1)
        return -(4.0 / (self.span * self.span)) * grad_variance


@dataclass
class METBalancedTask(BaseTask):
    label1: int
    target1: np.ndarray
    label2: int
    target2: np.ndarray
    beta: float = 20.0

    def expectations(self, psi: np.ndarray) -> Tuple[float, float, complex, complex]:
        o1 = np.vdot(self.target1, psi)
        o2 = np.vdot(self.target2, psi)
        return float(abs(o1) ** 2), float(abs(o2) ** 2), o1, o2

    def fisher_values(self, psi: np.ndarray) -> Tuple[float, float, float]:
        e1, e2, _o1, _o2 = self.expectations(psi)
        f1, f2 = e1 * e1, e2 * e2
        z = -self.beta * np.array([f1, f2], dtype=float)
        zmax = float(np.max(z))
        soft = -(zmax + math.log(0.5 * float(np.sum(np.exp(z - zmax))))) / self.beta
        return f1, f2, soft

    def cost(self, psi: np.ndarray) -> float:
        _f1, _f2, soft = self.fisher_values(psi)
        return 0.25 - soft

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        e1, e2, o1, o2 = self.expectations(psi)
        f1, f2 = e1 * e1, e2 * e2
        z = -self.beta * np.array([f1, f2], dtype=float)
        z -= float(np.max(z))
        weights = np.exp(z)
        weights /= float(np.sum(weights))
        grad_f1 = 4.0 * e1 * o1 * self.target1
        grad_f2 = 4.0 * e2 * o2 * self.target2
        return -(weights[0] * grad_f1 + weights[1] * grad_f2)


def random_distinct_labels(rng: np.random.Generator, dim: int, count: int) -> List[int]:
    return [int(x) for x in rng.choice(dim, size=count, replace=False)]


def make_vqe_tasks(n: int, depth: int) -> List[BaseTask]:
    dim = 1 << n
    seeds = [1101, 1102, 1103]
    tasks: List[BaseTask] = []

    rng = np.random.default_rng(seeds[0])
    scr = ComplexScrambler(n, depth, 100000 + seeds[0])
    x0 = int(rng.integers(dim))
    tau = scr.scrambled_basis(x0)
    tasks.append(VQEParentTask("VQE", "VQE-1", "random_parent", n, seeds[0], scr, x0, tau))

    rng = np.random.default_rng(seeds[1])
    scr = ComplexScrambler(n, depth, 100000 + seeds[1])
    x0 = int(rng.integers(dim))
    tau = scr.scrambled_basis(x0)
    spectrum = np.array(
        [bin(int(x) ^ int(x0)).count("1") / max(1, n) for x in range(dim)],
        dtype=float,
    )
    tasks.append(
        VQEDiagonalTask(
            "VQE",
            "VQE-2",
            "scrambled_hamming_spectrum",
            n,
            seeds[1],
            scr,
            x0,
            tau,
            spectrum,
        )
    )

    rng = np.random.default_rng(seeds[2])
    scr = ComplexScrambler(n, depth, 100000 + seeds[2])
    x0, x1 = random_distinct_labels(rng, dim, 2)
    tau = scr.scrambled_basis(x0)
    spectrum = np.ones(dim, dtype=float)
    spectrum[x0] = 0.0
    spectrum[x1] = 1e-2
    tasks.append(
        VQEDiagonalTask(
            "VQE",
            "VQE-3",
            "small_gap_scrambled_spectrum",
            n,
            seeds[2],
            scr,
            x0,
            tau,
            spectrum,
            x1,
        )
    )
    return tasks


def make_met_tasks(n: int, depth: int) -> List[BaseTask]:
    dim = 1 << n
    seeds = [2101, 2102, 2103]
    tasks: List[BaseTask] = []

    rng = np.random.default_rng(seeds[0])
    scr = ComplexScrambler(n, depth, 200000 + seeds[0])
    x0 = int(rng.integers(dim))
    tau = scr.scrambled_basis(x0)
    tasks.append(
        METSingleTargetTask(
            "MET",
            "MET-1",
            "single_target_fixed_readout_cfi",
            n,
            seeds[0],
            scr,
            x0,
            tau,
        )
    )

    rng = np.random.default_rng(seeds[1])
    scr = ComplexScrambler(n, depth, 200000 + seeds[1])
    gmin, gmax = random_distinct_labels(rng, dim, 2)
    generator = rng.normal(size=dim)
    generator = np.clip(generator, -0.8, 0.8)
    generator[gmin] = -1.0
    generator[gmax] = 1.0
    span = float(generator[gmax] - generator[gmin])
    tau_min = scr.scrambled_basis(gmin)
    tau_max = scr.scrambled_basis(gmax)
    tau_plus = normalize_with_phase(tau_max + tau_min)
    tau_minus = normalize_with_phase(tau_max - tau_min)
    tasks.append(
        METQFITask(
            "MET",
            "MET-2",
            "qfi_extremal_superposition",
            n,
            seeds[1],
            scr,
            generator,
            gmin,
            gmax,
            span,
            tau_plus,
            tau_minus,
        )
    )

    rng = np.random.default_rng(seeds[2])
    scr = ComplexScrambler(n, depth, 200000 + seeds[2])
    x1, x2 = random_distinct_labels(rng, dim, 2)
    tau1 = scr.scrambled_basis(x1)
    tau2 = scr.scrambled_basis(x2)
    tasks.append(
        METBalancedTask(
            "MET",
            "MET-3",
            "two_target_balanced_fisher",
            n,
            seeds[2],
            scr,
            x1,
            tau1,
            x2,
            tau2,
            20.0,
        )
    )
    return tasks


# -----------------------------------------------------------------------------
# Complex Hopf coordinate differential
# -----------------------------------------------------------------------------


def project_complex_theta(theta: np.ndarray, n: int, *, magnitude_eps: float = 1e-9) -> np.ndarray:
    theta = np.asarray(theta, dtype=float).copy()
    dim = 1 << n
    expected = 2 * dim - 1
    if theta.shape != (expected,):
        raise ValueError(f"Complex Hopf theta must have length {expected}, got {theta.size}.")
    theta[: dim - 1] = np.clip(theta[: dim - 1], magnitude_eps, math.pi / 2.0 - magnitude_eps)
    theta[dim - 1 :] = np.mod(theta[dim - 1 :], 2.0 * math.pi)
    return theta


def complex_hopf_coordinate_gradient(theta: np.ndarray, egrad: np.ndarray) -> np.ndarray:
    """Compute Re[J(theta)^† egrad] in O(2^n) time.

    The magnitude part is a bottom-up complex subtree contraction followed by a
    top-down incoming-mass pass.  The phase part follows directly from
    d x_l / d phi_l = i x_l.
    """
    theta = np.asarray(theta, dtype=float)
    egrad = np.asarray(egrad, dtype=complex)
    n = hopf.infer_n_from_theta(theta, case="complex")
    dim = 1 << n
    num_mags = dim - 1
    if egrad.shape != (dim,):
        raise ValueError(f"egrad must have shape {(dim,)}, got {egrad.shape}.")

    mags = theta[:num_mags]
    phases = theta[num_mags:]

    response = np.zeros(2 * dim, dtype=complex)
    response[dim : 2 * dim] = np.exp(-1j * phases) * egrad
    for node in range(num_mags, 0, -1):
        th = float(mags[node - 1])
        response[node] = math.cos(th) * response[2 * node] + math.sin(th) * response[2 * node + 1]

    incoming = np.zeros(num_mags + 1, dtype=float)
    incoming[1] = 1.0
    grad_mags = np.zeros(num_mags, dtype=float)
    for node in range(1, num_mags + 1):
        th = float(mags[node - 1])
        c = math.cos(th)
        s = math.sin(th)
        grad_mags[node - 1] = incoming[node] * float(
            np.real(-s * response[2 * node] + c * response[2 * node + 1])
        )
        left = 2 * node
        right = left + 1
        if left <= num_mags:
            incoming[left] = incoming[node] * c
        if right <= num_mags:
            incoming[right] = incoming[node] * s

    psi = hopf.vector_from_theta(theta, "complex")
    grad_phases = np.imag(np.conjugate(psi) * egrad)
    return np.concatenate([grad_mags, grad_phases])


def complex_metric_diagonal_fast(theta: np.ndarray, psi: Optional[np.ndarray] = None) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    n = hopf.infer_n_from_theta(theta, case="complex")
    dim = 1 << n
    num_mags = dim - 1
    mags = theta[:num_mags]
    g_mag = np.ones(num_mags, dtype=float)
    c2 = np.cos(mags) ** 2
    s2 = np.sin(mags) ** 2
    for node in range(2, num_mags + 1):
        parent = node // 2
        g_mag[node - 1] = g_mag[parent - 1] * (c2[parent - 1] if node % 2 == 0 else s2[parent - 1])
    if psi is None:
        psi = hopf.vector_from_theta(theta, "complex")
    return np.concatenate([g_mag, np.abs(psi) ** 2])


def hopf_params_and_grads(task: BaseTask, psi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return canonical theta, reconstructed state, coordinate gradient, sphere gradient.

    For the full complex Hopf chart, tangent_project(grad_euclidean) is exactly
    J g^{-1} grad_theta on the regular set.  Using the direct projection avoids
    materializing the dense Jacobian at every optimizer step while preserving
    the same state-sphere direction.
    """
    psi = unit_normalize(psi)
    theta = hopf.theta_from_vector(psi, "complex")
    psi_h = unit_normalize(hopf.vector_from_theta(theta, "complex"))
    egrad = task.euclidean_grad(psi_h)
    grad_theta = complex_hopf_coordinate_gradient(theta, egrad)
    state_grad = tangent_project(psi_h, egrad)
    return theta, psi_h, grad_theta, state_grad


# -----------------------------------------------------------------------------
# Optimizer states and line searches
# -----------------------------------------------------------------------------


@dataclass
class GeoCGState:
    prev_x: Optional[np.ndarray] = None
    prev_grad: Optional[np.ndarray] = None
    prev_dir: Optional[np.ndarray] = None


@dataclass
class LBFGSState:
    memory: int = 7
    s_list: List[np.ndarray] = field(default_factory=list)
    y_list: List[np.ndarray] = field(default_factory=list)
    rho_list: List[float] = field(default_factory=list)

    def _refresh(self) -> None:
        kept: List[Tuple[np.ndarray, np.ndarray, float]] = []
        for s, y in zip(self.s_list, self.y_list):
            sy = real_inner(s, y)
            if sy > 1e-12:
                kept.append((s, y, 1.0 / sy))
        kept = kept[-self.memory :]
        self.s_list = [item[0] for item in kept]
        self.y_list = [item[1] for item in kept]
        self.rho_list = [item[2] for item in kept]

    def transport_memory(self, x_old: np.ndarray, x_new: np.ndarray) -> None:
        if not self.s_list:
            return
        self.s_list = [transport_exact(x_old, x_new, s) for s in self.s_list]
        self.y_list = [transport_exact(x_old, x_new, y) for y in self.y_list]
        self._refresh()

    def direction(self, x: np.ndarray, grad: np.ndarray) -> np.ndarray:
        q = tangent_project(x, grad.copy())
        alphas: List[float] = []
        for s, y, rho in reversed(list(zip(self.s_list, self.y_list, self.rho_list))):
            alpha = rho * real_inner(s, q)
            alphas.append(alpha)
            q = q - alpha * y
        if self.s_list:
            sy = real_inner(self.s_list[-1], self.y_list[-1])
            yy = real_inner(self.y_list[-1], self.y_list[-1])
            gamma = sy / yy if yy > 1e-15 else 1.0
        else:
            gamma = 1.0
        r = gamma * q
        for (s, y, rho), alpha in zip(zip(self.s_list, self.y_list, self.rho_list), reversed(alphas)):
            beta = rho * real_inner(y, r)
            r = r + s * (alpha - beta)
        return -tangent_project(x, r)

    def update_memory(
        self,
        x_old: np.ndarray,
        x_new: np.ndarray,
        grad_old: np.ndarray,
        grad_new: np.ndarray,
        step_vector_old: np.ndarray,
    ) -> None:
        self.transport_memory(x_old, x_new)
        s_vec = tangent_project(x_new, transport_exact(x_old, x_new, step_vector_old))
        y_vec = tangent_project(x_new, grad_new - transport_exact(x_old, x_new, grad_old))
        sy = real_inner(s_vec, y_vec)
        if sy > 1e-12:
            self.s_list.append(s_vec)
            self.y_list.append(y_vec)
            self.rho_list.append(1.0 / sy)
            if len(self.s_list) > self.memory:
                self.s_list.pop(0)
                self.y_list.pop(0)
                self.rho_list.pop(0)


@dataclass
class BBState:
    prev_x: Optional[np.ndarray] = None
    prev_grad: Optional[np.ndarray] = None
    prev_step_vec: Optional[np.ndarray] = None
    alpha: float = 1.0
    iteration: int = 0


@dataclass
class AdamState:
    lr: float = 0.03
    beta1: float = 0.9
    beta2: float = 0.999
    eps: float = 1e-8
    m: Optional[np.ndarray] = None
    v: Optional[np.ndarray] = None
    t: int = 0

    def direction(self, grad: np.ndarray, *, normalize_rms: bool = True) -> np.ndarray:
        grad = np.asarray(grad, dtype=float)
        if self.m is None or self.v is None:
            self.m = np.zeros_like(grad)
            self.v = np.zeros_like(grad)
        self.t += 1
        self.m = self.beta1 * self.m + (1.0 - self.beta1) * grad
        self.v = self.beta2 * self.v + (1.0 - self.beta2) * grad * grad
        mhat = self.m / (1.0 - self.beta1 ** self.t)
        vhat = self.v / (1.0 - self.beta2 ** self.t)
        direction = -mhat / (np.sqrt(vhat) + self.eps)
        if float(np.dot(grad, direction)) >= 0.0:
            direction = -grad / (np.sqrt(vhat) + self.eps)
        if float(np.dot(grad, direction)) >= 0.0:
            direction = -grad.copy()
        if normalize_rms:
            rms = float(np.linalg.norm(direction) / math.sqrt(max(1, direction.size)))
            if rms > 1e-15:
                direction = direction / rms
        return direction


def state_gradient(task: BaseTask, psi: np.ndarray) -> np.ndarray:
    return tangent_project(psi, task.euclidean_grad(psi))


def strong_wolfe_geodesic_line_search(
    task: BaseTask,
    x: np.ndarray,
    direction: np.ndarray,
    grad0: np.ndarray,
    args: argparse.Namespace,
    *,
    initial_alpha: float = 1.0,
) -> Tuple[np.ndarray, float, int, np.ndarray, float]:
    """Strong-Wolfe search on a great circle, with bounded cost-only fallback."""
    p = tangent_project(x, direction)
    pnorm = float(np.linalg.norm(p))
    if pnorm < args.grad_tol:
        return x.copy(), 0.0, 0, np.zeros_like(x), 0.0

    grad0 = tangent_project(x, grad0)
    dphi0 = real_inner(grad0, p)
    if dphi0 >= 0.0:
        p = -grad0
        pnorm = float(np.linalg.norm(p))
        if pnorm < args.grad_tol:
            return x.copy(), 0.0, 0, np.zeros_like(x), 0.0
        dphi0 = -real_inner(grad0, grad0)

    f0 = float(task.cost(x))
    alpha_max = max(0.0, float(args.max_line_angle) / pnorm)
    if alpha_max <= 0.0:
        return x.copy(), 0.0, 0, np.zeros_like(x), 0.0

    alpha = min(max(float(initial_alpha), args.line_alpha_min), alpha_max)
    lo = 0.0
    hi: Optional[float] = None
    prev_alpha = 0.0
    prev_f = f0
    best_y = x.copy()
    best_f = f0
    best_alpha = 0.0
    evals = 0

    for iteration in range(int(args.line_maxiter)):
        trial_alpha = alpha
        y = sphere_exp(x, p, trial_alpha)
        f = float(task.cost(y))
        evals += 1
        if math.isfinite(f) and f < best_f:
            best_f = f
            best_y = y
            best_alpha = trial_alpha

        armijo_bad = (not math.isfinite(f)) or f > f0 + args.wolfe_c1 * trial_alpha * dphi0
        nondecrease = iteration > 0 and f >= prev_f
        if armijo_bad or nondecrease:
            hi = trial_alpha
        else:
            grad_y = state_gradient(task, y)
            evals += 1
            transported_p = transport_exact(x, y, p)
            dphi = real_inner(grad_y, transported_p)
            if abs(dphi) <= args.wolfe_c2 * abs(dphi0):
                return y, trial_alpha * pnorm, evals, trial_alpha * p, trial_alpha
            if dphi >= 0.0:
                hi = trial_alpha
            else:
                lo = trial_alpha

        prev_alpha = trial_alpha
        prev_f = f
        if hi is None:
            new_alpha = min(2.0 * trial_alpha, alpha_max)
            if new_alpha <= trial_alpha + 1e-15:
                break
            alpha = new_alpha
        else:
            alpha = 0.5 * (lo + hi)
            if alpha < args.line_alpha_min:
                break

    q = p / pnorm

    def f_scalar(angle: float) -> float:
        return float(task.cost(sphere_exp_unit(x, q, float(angle))))

    try:
        result = minimize_scalar(
            f_scalar,
            bounds=(0.0, float(args.max_line_angle)),
            method="bounded",
            options={"xatol": 1e-8, "maxiter": 80},
        )
        evals += int(getattr(result, "nfev", 0))
        angle = float(result.x)
        y = sphere_exp_unit(x, q, angle)
        f = float(task.cost(y))
        evals += 1
        if math.isfinite(f) and f <= best_f + 1e-15:
            alpha_eff = angle / pnorm
            return y, angle, evals, alpha_eff * p, alpha_eff
    except Exception:
        pass

    if best_alpha > 0.0:
        return best_y, best_alpha * pnorm, evals, best_alpha * p, best_alpha
    return x.copy(), 0.0, evals, np.zeros_like(x), 0.0


def coordinate_cost_line_search(
    task: BaseTask,
    theta: np.ndarray,
    direction: np.ndarray,
    current_cost: float,
    state_from_theta: Callable[[np.ndarray], np.ndarray],
    project_theta: Callable[[np.ndarray], np.ndarray],
    adam: AdamState,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, float, int, float, float]:
    direction = np.asarray(direction, dtype=float)
    if float(np.linalg.norm(direction)) < 1e-15:
        return theta.copy(), 0.0, 0, 0.0, float(current_cost)

    alpha = min(max(float(adam.lr), args.adam_lr_min), args.adam_lr_max)
    best_theta = theta.copy()
    best_cost = float(current_cost)
    best_alpha = 0.0
    evals = 0
    threshold = float(current_cost) - max(
        float(args.adam_accept_atol),
        float(args.adam_accept_rtol) * max(1.0, abs(float(current_cost))),
    )

    for _ in range(int(args.adam_line_evals)):
        trial = project_theta(theta + alpha * direction)
        if float(np.linalg.norm(trial - theta)) < 1e-18:
            alpha *= args.adam_shrink
            continue
        trial_state = state_from_theta(trial)
        cost = float(task.cost(trial_state))
        evals += 1
        if math.isfinite(cost) and cost < best_cost:
            best_cost = cost
            best_theta = trial
            best_alpha = alpha
        if math.isfinite(cost) and cost <= threshold:
            angle = state_angle(state_from_theta(theta), trial_state)
            adam.lr = min(args.adam_lr_max, max(args.adam_lr_min, alpha * args.adam_growth))
            return trial, alpha, evals, angle, cost
        alpha *= args.adam_shrink
        if alpha < args.adam_lr_min:
            break

    if best_alpha > 0.0 and best_cost < current_cost:
        angle = state_angle(state_from_theta(theta), state_from_theta(best_theta))
        adam.lr = min(args.adam_lr_max, max(args.adam_lr_min, best_alpha * args.adam_growth))
        return best_theta, best_alpha, evals, angle, best_cost

    adam.lr = max(args.adam_lr_min, alpha)
    return theta.copy(), 0.0, evals, 0.0, float(current_cost)


def bb_spectral_alpha(state: BBState, x: np.ndarray, grad: np.ndarray, args: argparse.Namespace) -> float:
    if state.prev_x is None or state.prev_grad is None or state.prev_step_vec is None:
        return float(args.bb_initial_alpha)
    s_vec = transport_exact(state.prev_x, x, state.prev_step_vec)
    y_vec = tangent_project(x, grad - transport_exact(state.prev_x, x, state.prev_grad))
    sy = real_inner(s_vec, y_vec)
    yy = real_inner(y_vec, y_vec)
    ss = real_inner(s_vec, s_vec)
    if sy <= 1e-14 or yy <= 1e-14 or ss <= 1e-14:
        return float(state.alpha)
    if args.bb_variant == "bb2" or (args.bb_variant == "alternate" and state.iteration % 2 == 1):
        alpha = sy / yy
    else:
        alpha = ss / sy
    return float(np.clip(alpha, args.bb_min_alpha, args.bb_max_alpha))


# -----------------------------------------------------------------------------
# Run bookkeeping
# -----------------------------------------------------------------------------


def complex_initial_state(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    state = rng.normal(size=1 << n) + 1j * rng.normal(size=1 << n)
    return normalize_with_phase(state)


def run_seeds(task: BaseTask, args: argparse.Namespace) -> List[int]:
    return [int(task.problem_seed) + int(args.seed_offset) + i for i in range(int(args.num_seeds))]


def make_record(
    task: BaseTask,
    mode: str,
    step: int,
    seed_index: int,
    run_seed: int,
    psi: np.ndarray,
    grad_theta: np.ndarray,
    state_grad_vec: np.ndarray,
    elapsed: float,
    *,
    last_step_angle: float | str = "",
    last_line_evals: int | str = "",
) -> Dict[str, Any]:
    cost = float(task.cost(psi))
    gap = float(task.gap(psi))
    return {
        "n": task.n,
        "app_class": task.app_class,
        "task": task.name,
        "seed": seed_index,
        "optimizer": MODE_ABBREVIATION[mode],
        "step": step,
        "gap": gap,
        "task_id": task.task_id,
        "problem_seed": task.problem_seed,
        "seed_index": seed_index,
        "run_seed": run_seed,
        "mode": mode,
        "cost": cost,
        "grad_norm": float(np.linalg.norm(grad_theta)),
        "state_grad_norm": float(np.linalg.norm(state_grad_vec)),
        "state_norm_error": float(abs(np.linalg.norm(psi) - 1.0)),
        "last_step_angle": last_step_angle,
        "last_line_evals": last_line_evals,
        "wall_time_sec": float(elapsed),
    }


def run_hopf_egt_cg(
    task: BaseTask,
    args: argparse.Namespace,
    *,
    run_seed: int,
    seed_index: int,
) -> List[Dict[str, Any]]:
    mode = "Hopf-EGT-CG"
    rows: List[Dict[str, Any]] = []
    start = time.time()
    psi = complex_initial_state(task.n, run_seed)
    state = GeoCGState()
    last_angle: float | str = ""
    last_evals: int | str = ""
    alpha_guess = 1.0

    for step in range(args.steps + 1):
        theta, psi, grad_theta, grad_s = hopf_params_and_grads(task, psi)
        rows.append(
            make_record(
                task,
                mode,
                step,
                seed_index,
                run_seed,
                psi,
                grad_theta,
                grad_s,
                time.time() - start,
                last_step_angle=last_angle,
                last_line_evals=last_evals,
            )
        )
        if step == args.steps:
            break
        if float(np.linalg.norm(grad_s)) < args.grad_tol:
            last_angle, last_evals = 0.0, 0
            continue

        if state.prev_x is None or state.prev_grad is None or state.prev_dir is None:
            direction = -grad_s
        else:
            grad_prev_t = transport_exact(state.prev_x, psi, state.prev_grad)
            dir_prev_t = transport_exact(state.prev_x, psi, state.prev_dir)
            yk = tangent_project(psi, grad_s - grad_prev_t)
            denom = real_inner(dir_prev_t, yk)
            if denom > 1e-14:
                beta_hs = real_inner(grad_s, yk) / denom
                beta_dy = real_inner(grad_s, grad_s) / denom
                beta = max(0.0, min(beta_hs, beta_dy))
            else:
                beta = 0.0
            direction = tangent_project(psi, -grad_s + beta * dir_prev_t)
            if real_inner(direction, grad_s) >= -1e-14 or float(np.linalg.norm(direction)) < 1e-12:
                direction = -grad_s

        old_psi = psi.copy()
        old_grad = grad_s.copy()
        old_direction = direction.copy()
        psi, last_angle, last_evals, _step_vec, alpha_eff = strong_wolfe_geodesic_line_search(
            task,
            psi,
            direction,
            grad_s,
            args,
            initial_alpha=alpha_guess,
        )
        state.prev_x = old_psi
        state.prev_grad = old_grad
        state.prev_dir = old_direction
        alpha_guess = alpha_eff if alpha_eff > 0.0 else 1.0
    return rows


def run_hopf_lbfgs(
    task: BaseTask,
    args: argparse.Namespace,
    *,
    run_seed: int,
    seed_index: int,
) -> List[Dict[str, Any]]:
    mode = "Hopf-Riemannian-LBFGS"
    rows: List[Dict[str, Any]] = []
    start = time.time()
    psi = complex_initial_state(task.n, run_seed)
    optimizer = LBFGSState(memory=args.lbfgs_memory)
    pending: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None
    last_angle: float | str = ""
    last_evals: int | str = ""
    alpha_guess = 1.0

    for step in range(args.steps + 1):
        theta, psi, grad_theta, grad_s = hopf_params_and_grads(task, psi)
        if pending is not None:
            prev_x, prev_grad, prev_step_vec = pending
            optimizer.update_memory(prev_x, psi, prev_grad, grad_s, prev_step_vec)
            pending = None
        rows.append(
            make_record(
                task,
                mode,
                step,
                seed_index,
                run_seed,
                psi,
                grad_theta,
                grad_s,
                time.time() - start,
                last_step_angle=last_angle,
                last_line_evals=last_evals,
            )
        )
        if step == args.steps:
            break
        if float(np.linalg.norm(grad_s)) < args.grad_tol:
            last_angle, last_evals = 0.0, 0
            continue

        direction = optimizer.direction(psi, grad_s)
        if real_inner(direction, grad_s) >= -1e-14 or float(np.linalg.norm(direction)) < 1e-12:
            direction = -grad_s
        old_psi = psi.copy()
        old_grad = grad_s.copy()
        psi, last_angle, last_evals, step_vec, alpha_eff = strong_wolfe_geodesic_line_search(
            task,
            psi,
            direction,
            grad_s,
            args,
            initial_alpha=alpha_guess,
        )
        pending = (old_psi, old_grad, step_vec)
        alpha_guess = alpha_eff if alpha_eff > 0.0 else 1.0
    return rows


def run_hopf_bb(
    task: BaseTask,
    args: argparse.Namespace,
    *,
    run_seed: int,
    seed_index: int,
) -> List[Dict[str, Any]]:
    mode = "Hopf-Riemannian-BB"
    rows: List[Dict[str, Any]] = []
    start = time.time()
    psi = complex_initial_state(task.n, run_seed)
    state = BBState(alpha=float(args.bb_initial_alpha))
    last_angle: float | str = ""
    last_evals: int | str = ""

    for step in range(args.steps + 1):
        theta, psi, grad_theta, grad_s = hopf_params_and_grads(task, psi)
        rows.append(
            make_record(
                task,
                mode,
                step,
                seed_index,
                run_seed,
                psi,
                grad_theta,
                grad_s,
                time.time() - start,
                last_step_angle=last_angle,
                last_line_evals=last_evals,
            )
        )
        if step == args.steps:
            break
        if float(np.linalg.norm(grad_s)) < args.grad_tol:
            last_angle, last_evals = 0.0, 0
            continue

        alpha_guess = bb_spectral_alpha(state, psi, grad_s, args)
        direction = -grad_s
        old_psi = psi.copy()
        old_grad = grad_s.copy()
        psi, last_angle, last_evals, step_vec, alpha_eff = strong_wolfe_geodesic_line_search(
            task,
            psi,
            direction,
            grad_s,
            args,
            initial_alpha=alpha_guess,
        )
        state.prev_x = old_psi
        state.prev_grad = old_grad
        state.prev_step_vec = step_vec
        state.alpha = alpha_eff if alpha_eff > 0.0 else alpha_guess
        state.iteration += 1
    return rows


def run_hopf_adam(
    task: BaseTask,
    args: argparse.Namespace,
    *,
    run_seed: int,
    seed_index: int,
) -> List[Dict[str, Any]]:
    mode = "Hopf-Adam"
    rows: List[Dict[str, Any]] = []
    start = time.time()
    initial = complex_initial_state(task.n, run_seed)
    theta = project_complex_theta(hopf.theta_from_vector(initial, "complex"), task.n)
    adam = AdamState(
        lr=args.adam_lr_init,
        beta1=args.adam_beta1,
        beta2=args.adam_beta2,
        eps=args.adam_eps,
    )
    last_angle: float | str = ""
    last_evals: int | str = ""

    def project(values: np.ndarray) -> np.ndarray:
        return project_complex_theta(values, task.n, magnitude_eps=args.magnitude_eps)

    def state_from(values: np.ndarray) -> np.ndarray:
        return unit_normalize(hopf.vector_from_theta(project(values), "complex"))

    for step in range(args.steps + 1):
        theta = project(theta)
        psi = state_from(theta)
        egrad = task.euclidean_grad(psi)
        grad_theta = complex_hopf_coordinate_gradient(theta, egrad)
        grad_s = tangent_project(psi, egrad)
        rows.append(
            make_record(
                task,
                mode,
                step,
                seed_index,
                run_seed,
                psi,
                grad_theta,
                grad_s,
                time.time() - start,
                last_step_angle=last_angle,
                last_line_evals=last_evals,
            )
        )
        if step == args.steps:
            break
        if float(np.linalg.norm(grad_theta)) < args.grad_tol:
            last_angle, last_evals = 0.0, 0
            continue
        direction = adam.direction(grad_theta, normalize_rms=args.adam_normalize_rms)
        theta, _accepted_lr, last_evals, last_angle, _new_cost = coordinate_cost_line_search(
            task,
            theta,
            direction,
            float(task.cost(psi)),
            state_from,
            project,
            adam,
            args,
        )
    return rows


RUNNERS: Dict[str, Callable[..., List[Dict[str, Any]]]] = {
    "Hopf-Adam": run_hopf_adam,
    "Hopf-EGT-CG": run_hopf_egt_cg,
    "Hopf-Riemannian-BB": run_hopf_bb,
    "Hopf-Riemannian-LBFGS": run_hopf_lbfgs,
}


# -----------------------------------------------------------------------------
# Validation, diagnostics, CSV, and plotting
# -----------------------------------------------------------------------------


def run_self_checks(n: int, depth: int) -> None:
    check_n = min(max(2, n), 4)
    rng = np.random.default_rng(41027)
    scrambler = ComplexScrambler(check_n, min(max(depth, 1), 2), 99173)
    x = normalize_with_phase(rng.normal(size=1 << check_n) + 1j * rng.normal(size=1 << check_n))
    inverse_error = float(np.linalg.norm(scrambler.apply_T(scrambler.apply(x)) - x))

    geometry_n = 3
    dim = 1 << geometry_n
    theta = np.concatenate(
        [
            rng.uniform(0.2, math.pi / 2.0 - 0.2, size=dim - 1),
            rng.uniform(0.0, 2.0 * math.pi, size=dim),
        ]
    )
    h = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    fast = complex_hopf_coordinate_gradient(theta, h)
    jacobian = hopf.jacobian(theta, "complex")
    dense = np.real(jacobian.conjugate().T @ h)
    gradient_error = float(np.max(np.abs(fast - dense)))
    psi = hopf.vector_from_theta(theta, "complex")
    metric = complex_metric_diagonal_fast(theta, psi)
    lifted = jacobian @ (dense / metric)
    lift_error = float(np.linalg.norm(lifted - tangent_project(psi, h)))

    print("Self-checks")
    print(f"  scrambler inverse residual : {inverse_error:.3e}")
    print(f"  coordinate-gradient error  : {gradient_error:.3e}")
    print(f"  metric-lift residual       : {lift_error:.3e}")
    if max(inverse_error, gradient_error, lift_error) > 1e-8:
        raise RuntimeError("Complex Hopf self-check failed.")


def write_csv(path: Path, rows: Iterable[Dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def read_csv(path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: Dict[str, Any] = dict(raw)
            row["n"] = int(float(row["n"]))
            row["seed"] = int(float(row["seed"]))
            row["seed_index"] = int(float(row.get("seed_index") or row["seed"]))
            row["run_seed"] = int(float(row.get("run_seed") or row["seed"]))
            row["step"] = int(float(row["step"]))
            for key in (
                "gap",
                "cost",
                "grad_norm",
                "state_grad_norm",
                "state_norm_error",
                "wall_time_sec",
            ):
                row[key] = float(row[key])
            for key in ("last_step_angle", "last_line_evals"):
                text = str(row.get(key, "")).strip()
                row[key] = float(text) if text else math.nan
            mode = str(row.get("mode", "")).strip()
            if not mode:
                reverse = {value: key for key, value in MODE_ABBREVIATION.items()}
                mode = reverse.get(str(row.get("optimizer", "")).strip(), "")
            row["mode"] = mode
            rows.append(row)
    return rows


def group_rows(rows: Sequence[Dict[str, Any]]) -> Dict[Tuple[str, str, int, str], List[Dict[str, Any]]]:
    groups: Dict[Tuple[str, str, int, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = (
            str(row["app_class"]),
            str(row["task_id"]),
            int(row["run_seed"]),
            str(row["mode"]),
        )
        groups[key].append(row)
    for series in groups.values():
        series.sort(key=lambda item: int(item["step"]))
    return dict(groups)


def plot_summary(
    rows: Sequence[Dict[str, Any]],
    output_path: Path,
    *,
    gap_floor: float,
    threshold: float,
    dpi: int,
) -> None:
    groups = group_rows(rows)
    final_by_app_mode: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    for (app, _task_id, _run_seed, mode), series in groups.items():
        if series:
            final_by_app_mode[(app, mode)].append(max(float(series[-1]["gap"]), gap_floor))

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.6), constrained_layout=True)
    for col, app in enumerate(("VQE", "MET")):
        ax = axes[col]
        modes = [mode for mode in PLOT_MODES if final_by_app_mode.get((app, mode))]
        positions = list(range(1, len(modes) + 1))
        data = [final_by_app_mode[(app, mode)] for mode in modes]
        box = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=True)
        for patch, mode in zip(box.get("boxes", []), modes):
            patch.set_facecolor(MODE_COLOR[mode])
            patch.set_alpha(0.18)
            patch.set_edgecolor(MODE_COLOR[mode])
        for line in box.get("medians", []):
            line.set_linewidth(1.7)
        for position, mode, values in zip(positions, modes, data):
            ax.scatter(
                [position] * len(values),
                values,
                s=18,
                color=MODE_COLOR[mode],
                alpha=0.62,
                zorder=3,
            )
        ax.axhline(threshold, color="0.5", linestyle="--", linewidth=1.0)
        ax.set_yscale("log")
        ax.set_ylim(bottom=gap_floor)
        ax.set_title(f"{app}: final gap")
        ax.set_xticks(positions)
        ax.set_xticklabels([MODE_LABEL[mode] for mode in modes], rotation=25, ha="right")
        ax.set_ylabel("gap")
        ax.grid(True, which="both", alpha=0.18)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def finite(value: float) -> bool:
    return math.isfinite(float(value))


def print_diagnostics(rows: Sequence[Dict[str, Any]], args: argparse.Namespace) -> None:
    groups = group_rows(rows)
    tasks = {(str(row["app_class"]), str(row["task_id"])) for row in rows}
    expected_groups = len(tasks) * int(args.num_seeds) * len(PLOT_MODES)
    expected_steps = set(range(int(args.steps) + 1))

    missing_step_groups = 0
    nonfinite_rows = 0
    negative_gap_rows = 0
    max_norm_error = 0.0
    initial_spreads: List[float] = []

    for series in groups.values():
        steps = {int(row["step"]) for row in series}
        if steps != expected_steps:
            missing_step_groups += 1
        for row in series:
            gap = float(row["gap"])
            if not finite(gap):
                nonfinite_rows += 1
            if gap < -1e-10:
                negative_gap_rows += 1
            max_norm_error = max(max_norm_error, abs(float(row["state_norm_error"])))

    starts: Dict[Tuple[str, str, int], List[float]] = defaultdict(list)
    for (app, task_id, run_seed, _mode), series in groups.items():
        if series:
            starts[(app, task_id, run_seed)].append(float(series[0]["gap"]))
    for values in starts.values():
        if values:
            initial_spreads.append(max(values) - min(values))

    print("\nDiagnostics")
    print(f"  observed traces            : {len(groups)} / expected {expected_groups}")
    print(f"  rows                       : {len(rows)}")
    print(f"  traces with bad step grid  : {missing_step_groups}")
    print(f"  non-finite gap rows        : {nonfinite_rows}")
    print(f"  materially negative gaps  : {negative_gap_rows}")
    print(f"  max state-norm error       : {max_norm_error:.3e}")
    print(f"  max initial spread/modes   : {max(initial_spreads, default=0.0):.3e}")

    header = (
        f"  {'class':<5} {'mode':<25} {'traces':>6} {'mean final':>13} "
        f"{'median':>13} {'min':>13} {'max':>13} {'<= threshold':>12} {'worse':>7}"
    )
    print("\n" + header)
    print("  " + "-" * (len(header) - 2))
    for app in ("VQE", "MET"):
        for mode in PLOT_MODES:
            selected = [series for (a, _t, _s, m), series in groups.items() if a == app and m == mode]
            if not selected:
                continue
            finals = [float(series[-1]["gap"]) for series in selected]
            initial = [float(series[0]["gap"]) for series in selected]
            hits = sum(value <= args.threshold for value in finals)
            worse = sum(final > first + 1e-12 for first, final in zip(initial, finals))
            print(
                f"  {app:<5} {MODE_LABEL[mode]:<25} {len(finals):>6d} "
                f"{stats.fmean(finals):>13.6e} {stats.median(finals):>13.6e} "
                f"{min(finals):>13.6e} {max(finals):>13.6e} "
                f"{hits:>5d}/{len(finals):<6d} {worse:>7d}"
            )

    alerts: List[str] = []
    if len(groups) != expected_groups:
        alerts.append("unexpected number of task/seed/mode traces")
    if missing_step_groups:
        alerts.append("one or more traces do not contain steps 0..steps")
    if nonfinite_rows:
        alerts.append("non-finite gaps detected")
    if max_norm_error > 1e-10:
        alerts.append("state normalization error exceeds 1e-10")
    if max(initial_spreads, default=0.0) > 1e-10:
        alerts.append("optimizers did not start from identical states")
    if alerts:
        print("\n  Alerts: " + "; ".join(alerts) + ".")
    else:
        print("\n  No structural or numerical alerts triggered.")


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------


def run_experiment(args: argparse.Namespace) -> List[Dict[str, Any]]:
    tasks = make_vqe_tasks(args.n, args.scramble_depth) + make_met_tasks(args.n, args.scramble_depth)
    rows: List[Dict[str, Any]] = []
    total_runs = len(tasks) * args.num_seeds * len(PLOT_MODES)
    completed = 0
    experiment_start = time.time()

    print(
        f"Complex Hopf stress test: n={args.n}, tasks={len(tasks)}, "
        f"seeds/task={args.num_seeds}, updates={args.steps}, tracks={len(PLOT_MODES)}"
    )
    for task in tasks:
        seeds = run_seeds(task, args)
        for seed_index, run_seed in enumerate(seeds):
            for mode in PLOT_MODES:
                completed += 1
                run_start = time.time()
                print(
                    f"[{completed:3d}/{total_runs}] {task.app_class} {task.task_id} "
                    f"seed={seed_index:02d} {MODE_LABEL[mode]:<10}",
                    end="",
                    flush=True,
                )
                run_rows = RUNNERS[mode](task, args, run_seed=run_seed, seed_index=seed_index)
                rows.extend(run_rows)
                final_gap = float(run_rows[-1]["gap"])
                print(f"  final={final_gap:.3e}  {time.time() - run_start:.2f}s")

    print(f"Experiment wall time: {time.time() - experiment_start:.1f}s")
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Complex Hopf stress test with four optimizer tracks and one paper-facing final-gap plot."
    )
    parser.add_argument("--n", type=int, default=6, help="Number of qubits (default: 6).")
    parser.add_argument("--steps", type=int, default=200, help="Accepted optimizer updates; step 0 is also recorded.")
    parser.add_argument("--num-seeds", type=int, default=10, help="Deterministic complex initial states per task.")
    parser.add_argument("--seed-offset", type=int, default=77777)
    parser.add_argument("--scramble-depth", type=int, default=4)
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--csv-name", default="complex_hopf_stress_data.csv")
    parser.add_argument("--plot-name", default="hopf_complex.png")
    parser.add_argument("--plot-only", action="store_true", help="Read the existing CSV and regenerate diagnostics/plot.")
    parser.add_argument("--quick", action="store_true", help="Use one seed, two updates, and scramble depth one.")

    parser.add_argument("--gap-floor", type=float, default=1e-16)
    parser.add_argument("--threshold", type=float, default=1e-8)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--grad-tol", type=float, default=1e-12)
    parser.add_argument("--magnitude-eps", type=float, default=1e-9)

    parser.add_argument("--max-line-angle", type=float, default=math.pi / 2.0)
    parser.add_argument("--line-alpha-min", type=float, default=1e-8)
    parser.add_argument("--line-maxiter", type=int, default=25)
    parser.add_argument("--wolfe-c1", type=float, default=0.485)
    parser.add_argument("--wolfe-c2", type=float, default=0.999)

    parser.add_argument("--lbfgs-memory", type=int, default=7)
    parser.add_argument("--bb-initial-alpha", type=float, default=1.0)
    parser.add_argument("--bb-min-alpha", type=float, default=1e-6)
    parser.add_argument("--bb-max-alpha", type=float, default=10.0)
    parser.add_argument("--bb-variant", choices=["bb1", "bb2", "alternate"], default="alternate")

    parser.add_argument("--adam-lr-init", type=float, default=0.03)
    parser.add_argument("--adam-lr-min", type=float, default=1e-8)
    parser.add_argument("--adam-lr-max", type=float, default=0.5)
    parser.add_argument("--adam-beta1", type=float, default=0.9)
    parser.add_argument("--adam-beta2", type=float, default=0.999)
    parser.add_argument("--adam-eps", type=float, default=1e-8)
    parser.add_argument("--adam-line-evals", type=int, default=12)
    parser.add_argument("--adam-shrink", type=float, default=0.5)
    parser.add_argument("--adam-growth", type=float, default=1.2)
    parser.add_argument("--adam-accept-atol", type=float, default=1e-15)
    parser.add_argument("--adam-accept-rtol", type=float, default=1e-12)
    parser.add_argument(
        "--no-adam-rms-normalization",
        dest="adam_normalize_rms",
        action="store_false",
        help="Disable RMS normalization of the Adam direction.",
    )
    parser.set_defaults(adam_normalize_rms=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.quick:
        args.steps = min(args.steps, 2)
        args.num_seeds = 1
        args.scramble_depth = min(args.scramble_depth, 1)
    if args.n < 1:
        raise ValueError("--n must be >= 1")
    if args.steps < 0:
        raise ValueError("--steps must be >= 0")
    if args.num_seeds < 1:
        raise ValueError("--num-seeds must be >= 1")
    if args.gap_floor <= 0.0:
        raise ValueError("--gap-floor must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / args.csv_name
    plot_path = args.outdir / args.plot_name

    if args.plot_only:
        if not csv_path.exists():
            raise FileNotFoundError(f"Cannot use --plot-only: {csv_path} does not exist.")
        rows = read_csv(csv_path)
        print(f"Loaded {len(rows)} rows from {csv_path}")
    else:
        run_self_checks(args.n, args.scramble_depth)
        rows = run_experiment(args)
        count = write_csv(csv_path, rows)
        print(f"Wrote {count} rows to {csv_path}")

    print_diagnostics(rows, args)
    plot_summary(rows, plot_path, gap_floor=args.gap_floor, threshold=args.threshold, dpi=args.dpi)
    print(f"Wrote plot to {plot_path}")


if __name__ == "__main__":
    main()
