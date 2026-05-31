#!/usr/bin/env python3
"""
hopf_data.py

Lightweight synthetic full-real-state stress tests for the Hopf ansatz.

Run:
    python hopf_data.py --n 8

Outputs exactly two CSV files by default:
    vqe_hopf_data_n{n}.csv
    met_hopf_data_n{n}.csv

Each CSV records 3 generated tasks x 10 initial-state seeds x 3 geometry-native optimization modes x (steps + 1) rows by default.
The parameter vector and gradient vector are stored as semicolon-separated strings.

The script assumes that hopf_utils.py is in the same folder or on PYTHONPATH.
It uses the real Hopf chart only.

Modes:
    VQE and MET:
      1. Hopf-EGT-CG
      2. Hopf-Riemannian-LBFGS
      3. Hopf-Riemannian-BB

All three modes use the same quantum-information interface: objective-value
queries and Hopf-coordinate gradient queries. Ritz/subspace oracle updates,
Adam baselines, and Möttönen physical-angle baselines are intentionally omitted
from this cost-and-gradient dataset.

Synthetic tasks:
    VQE-1: random parent Hamiltonian, H = I - |tau><tau|.
    VQE-2: scrambled diagonal Hamming spectrum.
    VQE-3: scrambled small-gap diagonal spectrum.

    MET-1: single-target fixed-readout Fisher, F = <tau|rho|tau>^2.
    MET-2: QFI extremal-superposition task for scrambled diagonal generator.
    MET-3: two-target balanced Fisher soft-min task.

All tasks are scrambled by a fixed real orthogonal circuit per task.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    from scipy.optimize import minimize_scalar
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This script requires scipy. Install it with: pip install scipy") from exc

import hopf_utils as hopf


# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------


def vector_to_string(x: np.ndarray, precision: int = 10) -> str:
    """Serialize a 1D float vector as semicolon-separated values."""
    x = np.asarray(x, dtype=float).ravel()
    fmt = f"%.{precision}g"
    return ";".join(fmt % float(v) for v in x)


def bitstring(label: int, n: int) -> str:
    return format(int(label), f"0{n}b")


def label_from_bitstring(s: str) -> int:
    return int(str(s).strip(), 2)


def popcount_int(value: int) -> int:
    """Python-version/NumPy-safe population count."""
    return bin(int(value)).count("1")


def normalize(x: np.ndarray, *, eps: float = 1e-15) -> np.ndarray:
    nrm = float(np.linalg.norm(x))
    if nrm < eps:
        raise ValueError("Cannot normalize near-zero vector.")
    y = np.asarray(x, dtype=float) / nrm
    # Fix a deterministic global sign for stable target-state generation and CSV traces.
    idx = int(np.argmax(np.abs(y)))
    if y[idx] < 0:
        y = -y
    return y


def unit_normalize(x: np.ndarray, *, eps: float = 1e-15) -> np.ndarray:
    """Normalize without changing the global sign; used for sphere geodesics."""
    nrm = float(np.linalg.norm(x))
    if nrm < eps:
        raise ValueError("Cannot normalize near-zero vector.")
    return np.asarray(x, dtype=float) / nrm


def tangent_project(x: np.ndarray, v: np.ndarray) -> np.ndarray:
    return np.asarray(v, dtype=float) - float(np.dot(x, v)) * np.asarray(x, dtype=float)


def sphere_exp(x: np.ndarray, p: np.ndarray, alpha: float) -> np.ndarray:
    """Exact exponential map on the unit sphere along tangent vector p."""
    x = unit_normalize(x)
    p = tangent_project(x, p)
    nrm = float(np.linalg.norm(p))
    if nrm < 1e-15 or abs(alpha) < 1e-15:
        return x.copy()
    ell = alpha * nrm
    y = math.cos(ell) * x + math.sin(ell) * (p / nrm)
    return unit_normalize(y)


def sphere_exp_unit(x: np.ndarray, q: np.ndarray, ell: float) -> np.ndarray:
    """Sphere geodesic with unit tangent q and angular displacement ell."""
    x = unit_normalize(x)
    q = tangent_project(x, q)
    nrm = float(np.linalg.norm(q))
    if nrm < 1e-15 or abs(ell) < 1e-15:
        return x.copy()
    q = q / nrm
    return unit_normalize(math.cos(ell) * x + math.sin(ell) * q)


def transport_project(x_old: np.ndarray, x_new: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Simple projection vector transport to T_{x_new}S.

    Kept only as a numerical fallback. The Hopf EGT-CG and Riemannian L-BFGS
    modes below use exact sphere geodesic transport.
    """
    return tangent_project(x_new, v)



# -----------------------------------------------------------------------------
# Real orthogonal scrambler
# -----------------------------------------------------------------------------


def _apply_ry_inplace(v: np.ndarray, n: int, qubit: int, theta: float) -> np.ndarray:
    """Apply real R_y-like plane rotation [[c,-s],[s,c]] to a state vector."""
    out = v.copy()
    c, s = math.cos(theta), math.sin(theta)
    mask = 1 << (n - 1 - qubit)
    dim = 1 << n
    for i in range(dim):
        if (i & mask) == 0:
            j = i | mask
            a = v[i]
            b = v[j]
            out[i] = c * a - s * b
            out[j] = s * a + c * b
    return out


def _apply_cnot_inplace(v: np.ndarray, n: int, control: int, target: int) -> np.ndarray:
    """Apply CNOT as a permutation using qubit indices 0..n-1 with q=0 MSB."""
    out = v.copy()
    cmask = 1 << (n - 1 - control)
    tmask = 1 << (n - 1 - target)
    dim = 1 << n
    for i in range(dim):
        if (i & cmask) and ((i & tmask) == 0):
            j = i | tmask
            out[i] = v[j]
            out[j] = v[i]
    return out


@dataclass
class RealScrambler:
    n: int
    depth: int
    seed: int
    angles: np.ndarray = field(init=False)
    pairs_by_layer: List[List[Tuple[int, int]]] = field(init=False)

    def __post_init__(self) -> None:
        rng = np.random.default_rng(self.seed)
        self.angles = rng.uniform(-math.pi, math.pi, size=(self.depth, self.n))
        self.pairs_by_layer = []
        for layer in range(self.depth):
            pairs: List[Tuple[int, int]] = []
            if self.n >= 2:
                if layer % 2 == 0:
                    # Even brickwork: 0->1, 2->3, ...
                    for q in range(0, self.n - 1, 2):
                        pairs.append((q, q + 1))
                else:
                    # Odd brickwork with reversed direction: 2->1, 4->3, ...
                    for q in range(1, self.n - 1, 2):
                        pairs.append((q + 1, q))
            self.pairs_by_layer.append(pairs)

    def apply(self, x: np.ndarray) -> np.ndarray:
        y = np.asarray(x, dtype=float).copy()
        for layer in range(self.depth):
            for q in range(self.n):
                y = _apply_ry_inplace(y, self.n, q, float(self.angles[layer, q]))
            for c, t in self.pairs_by_layer[layer]:
                y = _apply_cnot_inplace(y, self.n, c, t)
        return y

    def apply_T(self, x: np.ndarray) -> np.ndarray:
        y = np.asarray(x, dtype=float).copy()
        for layer in reversed(range(self.depth)):
            for c, t in reversed(self.pairs_by_layer[layer]):
                y = _apply_cnot_inplace(y, self.n, c, t)  # self-inverse
            for q in reversed(range(self.n)):
                y = _apply_ry_inplace(y, self.n, q, -float(self.angles[layer, q]))
        return y

    def basis_state(self, label: int) -> np.ndarray:
        e = np.zeros(1 << self.n, dtype=float)
        e[int(label)] = 1.0
        return e

    def scrambled_basis(self, label: int) -> np.ndarray:
        return normalize(self.apply(self.basis_state(label)))


# -----------------------------------------------------------------------------
# Fast Hopf coordinate gradient: g_theta = J^T grad_euclidean
# -----------------------------------------------------------------------------


def hopf_coordinate_gradient(theta: np.ndarray, euclidean_grad: np.ndarray) -> np.ndarray:
    """Compute real Hopf coordinate gradient without materializing J.

    For node i, support span is [start, end) with halves [start, mid), [mid, end).
    dpsi/dtheta_i = -tan(theta_i) psi on the left half, cot(theta_i) psi on
    the right half.
    """
    theta = np.asarray(theta, dtype=float)
    psi = hopf.vector_from_theta(theta, "real")
    w = psi * np.asarray(euclidean_grad, dtype=float)
    prefix = np.empty(w.size + 1, dtype=float)
    prefix[0] = 0.0
    np.cumsum(w, out=prefix[1:])

    L = theta.size
    n = int(round(math.log2(L + 1)))
    g = np.zeros(L, dtype=float)
    eps = 1e-14
    for j0 in range(L):
        i = j0 + 1
        start, mid, end = hopf._subtree_span(i, n)
        left = prefix[mid] - prefix[start]
        right = prefix[end] - prefix[mid]
        s = math.sin(float(theta[j0]))
        c = math.cos(float(theta[j0]))
        term = 0.0
        if abs(c) > eps:
            term += -(s / c) * left
        if abs(s) > eps:
            term += (c / s) * right
        g[j0] = term
    return g


def theta_from_state_safe(psi: np.ndarray, n: int) -> np.ndarray:
    th = hopf.theta_from_vector(unit_normalize(psi), "real")
    return hopf.clip_theta_hopf_real(th, n=n)


# -----------------------------------------------------------------------------
# Synthetic task classes
# -----------------------------------------------------------------------------


@dataclass
class BaseTask:
    task_type: str
    task_id: str
    name: str
    n: int
    problem_seed: int
    scrambler: RealScrambler
    target_label: int
    target: np.ndarray

    def cost(self, psi: np.ndarray) -> float:
        raise NotImplementedError

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def metrics(self, psi: np.ndarray) -> Dict[str, float | int | str]:
        y = self.scrambler.apply_T(psi)
        probs = y * y
        decoded = int(np.argmax(probs))
        return {
            "target_label": int(self.target_label),
            "target_bitstring": bitstring(self.target_label, self.n),
            "target_overlap": float(np.dot(self.target, psi) ** 2),
            "decoded_argmax": decoded,
            "decoded_bitstring": bitstring(decoded, self.n),
            "decoded_probability": float(probs[decoded]),
            "distractor_label": "",
            "distractor_overlap": "",
            "energy": "",
            "energy_gap": "",
            "energy_per_clause": "",
            "met_objective": "",
            "cfi": "",
            "cfi_gap": "",
            "A_expectation": "",
            "qfi": "",
            "normalized_qfi": "",
            "F1": "",
            "F2": "",
            "F_bal": "",
            "overlap_plus": "",
            "overlap_minus": "",
            "balanced_target_overlap": "",
        }

    def apply_H(self, psi: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Only VQE tasks implement apply_H.")

    def primitive_values(self, psi: np.ndarray) -> np.ndarray:
        """Primitive expectation values used for literal parameter-shift gradients."""
        return np.array([self.cost(psi)], dtype=float)

    def primitive_cost_coefficients(self, psi: np.ndarray) -> np.ndarray:
        """d(cost)/d(primitive_values), evaluated at current psi."""
        return np.array([1.0], dtype=float)


@dataclass
class VQEParentTask(BaseTask):
    def cost(self, psi: np.ndarray) -> float:
        o = float(np.dot(self.target, psi))
        return 1.0 - o * o

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        o = float(np.dot(self.target, psi))
        return -2.0 * o * self.target

    def apply_H(self, psi: np.ndarray) -> np.ndarray:
        return psi - float(np.dot(self.target, psi)) * self.target

    def metrics(self, psi: np.ndarray) -> Dict[str, float | int | str]:
        d = super().metrics(psi)
        energy = self.cost(psi)
        d.update({"energy": energy, "energy_gap": energy, "energy_per_clause": energy})
        return d


@dataclass
class VQEDiagonalTask(BaseTask):
    spectrum: np.ndarray
    distractor_label: Optional[int] = None

    def apply_H(self, psi: np.ndarray) -> np.ndarray:
        y = self.scrambler.apply_T(psi)
        return self.scrambler.apply(self.spectrum * y)

    def cost(self, psi: np.ndarray) -> float:
        y = self.scrambler.apply_T(psi)
        return float(np.dot(self.spectrum, y * y))

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        return 2.0 * self.apply_H(psi)

    def metrics(self, psi: np.ndarray) -> Dict[str, float | int | str]:
        d = super().metrics(psi)
        energy = self.cost(psi)
        d.update({"energy": energy, "energy_gap": energy, "energy_per_clause": energy})
        if self.distractor_label is not None:
            tau1 = self.scrambler.scrambled_basis(self.distractor_label)
            d.update({
                "distractor_label": int(self.distractor_label),
                "distractor_overlap": float(np.dot(tau1, psi) ** 2),
            })
        return d


@dataclass
class METSingleTargetTask(BaseTask):
    def expectation_A(self, psi: np.ndarray) -> float:
        o = float(np.dot(self.target, psi))
        return o * o

    def cost(self, psi: np.ndarray) -> float:
        a = self.expectation_A(psi)
        return 1.0 - a * a

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        o = float(np.dot(self.target, psi))
        # cost = 1 - o^4
        return -4.0 * (o ** 3) * self.target

    def primitive_values(self, psi: np.ndarray) -> np.ndarray:
        return np.array([self.expectation_A(psi)], dtype=float)

    def primitive_cost_coefficients(self, psi: np.ndarray) -> np.ndarray:
        a = self.expectation_A(psi)
        return np.array([-2.0 * a], dtype=float)

    def metrics(self, psi: np.ndarray) -> Dict[str, float | int | str]:
        d = super().metrics(psi)
        a = self.expectation_A(psi)
        cfi = a * a
        d.update({
            "met_objective": cfi,
            "cfi": cfi,
            "cfi_gap": 1.0 - cfi,
            "A_expectation": a,
        })
        return d


@dataclass
class METQFITask(BaseTask):
    generator: np.ndarray
    gmin_label: int
    gmax_label: int
    span: float
    tau_plus: np.ndarray
    tau_minus: np.ndarray

    def _Gy(self, psi: np.ndarray, power: int = 1) -> np.ndarray:
        y = self.scrambler.apply_T(psi)
        if power == 1:
            return self.scrambler.apply(self.generator * y)
        return self.scrambler.apply((self.generator ** power) * y)

    def moments(self, psi: np.ndarray) -> Tuple[float, float]:
        y = self.scrambler.apply_T(psi)
        p = y * y
        mu = float(np.dot(self.generator, p))
        nu = float(np.dot(self.generator ** 2, p))
        return mu, nu

    def normalized_qfi(self, psi: np.ndarray) -> float:
        mu, nu = self.moments(psi)
        return float(4.0 * max(0.0, nu - mu * mu) / (self.span * self.span))

    def cost(self, psi: np.ndarray) -> float:
        return 1.0 - self.normalized_qfi(psi)

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        mu, _nu = self.moments(psi)
        Gpsi = self._Gy(psi, 1)
        G2psi = self._Gy(psi, 2)
        grad_var = 2.0 * G2psi - 4.0 * mu * Gpsi
        return -(4.0 / (self.span * self.span)) * grad_var

    def primitive_values(self, psi: np.ndarray) -> np.ndarray:
        mu, nu = self.moments(psi)
        return np.array([mu, nu], dtype=float)

    def primitive_cost_coefficients(self, psi: np.ndarray) -> np.ndarray:
        mu, _nu = self.moments(psi)
        # cost = 1 - 4(nu - mu^2)/span^2
        return np.array([8.0 * mu / (self.span * self.span), -4.0 / (self.span * self.span)], dtype=float)

    def metrics(self, psi: np.ndarray) -> Dict[str, float | int | str]:
        d = super().metrics(psi)
        qfi = self.normalized_qfi(psi)
        d.update({
            "met_objective": qfi,
            "cfi": "",
            "cfi_gap": 1.0 - qfi,
            "qfi": qfi * self.span * self.span,
            "normalized_qfi": qfi,
            "overlap_plus": float(np.dot(self.tau_plus, psi) ** 2),
            "overlap_minus": float(np.dot(self.tau_minus, psi) ** 2),
            "target_label": f"{self.gmin_label},{self.gmax_label}",
            "target_bitstring": f"{bitstring(self.gmin_label, self.n)}|{bitstring(self.gmax_label, self.n)}",
        })
        return d


@dataclass
class METBalancedTask(BaseTask):
    label2: int
    target2: np.ndarray
    beta: float = 20.0

    def expectations(self, psi: np.ndarray) -> Tuple[float, float]:
        o1 = float(np.dot(self.target, psi))
        o2 = float(np.dot(self.target2, psi))
        return o1 * o1, o2 * o2

    def fisher_values(self, psi: np.ndarray) -> Tuple[float, float, float]:
        e1, e2 = self.expectations(psi)
        F1, F2 = e1 * e1, e2 * e2
        vals = np.array([-self.beta * F1, -self.beta * F2], dtype=float)
        m = float(np.max(vals))
        soft = -(1.0 / self.beta) * (math.log(0.5 * float(np.sum(np.exp(vals - m)))) + m)
        return F1, F2, soft

    def cost(self, psi: np.ndarray) -> float:
        _F1, _F2, soft = self.fisher_values(psi)
        return 0.25 - soft

    def euclidean_grad(self, psi: np.ndarray) -> np.ndarray:
        e1, e2 = self.expectations(psi)
        F1, F2 = e1 * e1, e2 * e2
        weights_raw = np.exp(-self.beta * np.array([F1, F2]))
        weights = weights_raw / float(np.sum(weights_raw))
        o1 = float(np.dot(self.target, psi))
        o2 = float(np.dot(self.target2, psi))
        gradF1 = 4.0 * (o1 ** 3) * self.target
        gradF2 = 4.0 * (o2 ** 3) * self.target2
        # cost = 1/4 - softmin, d softmin/dF_i = weight_i
        return -(weights[0] * gradF1 + weights[1] * gradF2)

    def primitive_values(self, psi: np.ndarray) -> np.ndarray:
        e1, e2 = self.expectations(psi)
        return np.array([e1, e2], dtype=float)

    def primitive_cost_coefficients(self, psi: np.ndarray) -> np.ndarray:
        e1, e2 = self.expectations(psi)
        F1, F2 = e1 * e1, e2 * e2
        weights_raw = np.exp(-self.beta * np.array([F1, F2]))
        weights = weights_raw / float(np.sum(weights_raw))
        # d cost / d e_i = - d soft / dF_i * dF_i/de_i = -weights_i * 2e_i
        return np.array([-2.0 * weights[0] * e1, -2.0 * weights[1] * e2], dtype=float)

    def metrics(self, psi: np.ndarray) -> Dict[str, float | int | str]:
        d = super().metrics(psi)
        F1, F2, soft = self.fisher_values(psi)
        balanced = normalize(self.target + self.target2)
        d.update({
            "met_objective": soft,
            "cfi": soft,
            "cfi_gap": 0.25 - soft,
            "F1": F1,
            "F2": F2,
            "F_bal": soft,
            "balanced_target_overlap": float(np.dot(balanced, psi) ** 2),
            "target_label": f"{self.target_label},{self.label2}",
            "target_bitstring": f"{bitstring(self.target_label, self.n)}|{bitstring(self.label2, self.n)}",
        })
        return d


# -----------------------------------------------------------------------------
# Task generation
# -----------------------------------------------------------------------------


def random_distinct_labels(rng: np.random.Generator, dim: int, count: int) -> List[int]:
    return [int(x) for x in rng.choice(dim, size=count, replace=False)]


def make_vqe_tasks(n: int, depth: int) -> List[BaseTask]:
    dim = 1 << n
    seeds = [1101, 1102, 1103]
    tasks: List[BaseTask] = []

    # VQE-1: random parent Hamiltonian.
    rng = np.random.default_rng(seeds[0])
    scr = RealScrambler(n, depth, 100000 + seeds[0])
    x0 = int(rng.integers(dim))
    tau = scr.scrambled_basis(x0)
    tasks.append(VQEParentTask("VQE", "VQE-1", "random_parent", n, seeds[0], scr, x0, tau))

    # VQE-2: scrambled diagonal Hamming spectrum.
    rng = np.random.default_rng(seeds[1])
    scr = RealScrambler(n, depth, 100000 + seeds[1])
    x0 = int(rng.integers(dim))
    tau = scr.scrambled_basis(x0)
    spectrum = np.zeros(dim, dtype=float)
    for x in range(dim):
        spectrum[x] = popcount_int(int(x) ^ int(x0)) / max(1, n)
    tasks.append(VQEDiagonalTask("VQE", "VQE-2", "scrambled_hamming_spectrum", n, seeds[1], scr, x0, tau, spectrum))

    # VQE-3: small-gap scrambled spectrum.
    rng = np.random.default_rng(seeds[2])
    scr = RealScrambler(n, depth, 100000 + seeds[2])
    x0, x1 = random_distinct_labels(rng, dim, 2)
    tau = scr.scrambled_basis(x0)
    spectrum = np.ones(dim, dtype=float)
    spectrum[x0] = 0.0
    spectrum[x1] = 1e-2
    tasks.append(VQEDiagonalTask("VQE", "VQE-3", "small_gap_scrambled_spectrum", n, seeds[2], scr, x0, tau, spectrum, x1))
    return tasks


def make_met_tasks(n: int, depth: int) -> List[BaseTask]:
    dim = 1 << n
    seeds = [2101, 2102, 2103]
    tasks: List[BaseTask] = []

    # MET-1: single-target fixed-readout Fisher.
    rng = np.random.default_rng(seeds[0])
    scr = RealScrambler(n, depth, 200000 + seeds[0])
    x0 = int(rng.integers(dim))
    tau = scr.scrambled_basis(x0)
    tasks.append(METSingleTargetTask("MET", "MET-1", "single_target_fixed_readout_cfi", n, seeds[0], scr, x0, tau))

    # MET-2: QFI extremal-superposition task.
    rng = np.random.default_rng(seeds[1])
    scr = RealScrambler(n, depth, 200000 + seeds[1])
    gmin, gmax = random_distinct_labels(rng, dim, 2)
    generator = rng.normal(size=dim)
    generator[gmin] = -1.0
    generator[gmax] = 1.0
    # Normalize span exactly to two for a clean optimum.
    generator = np.clip(generator, -0.8, 0.8)
    generator[gmin] = -1.0
    generator[gmax] = 1.0
    span = float(generator[gmax] - generator[gmin])
    tau_min = scr.scrambled_basis(gmin)
    tau_max = scr.scrambled_basis(gmax)
    tau_plus = normalize(tau_max + tau_min)
    tau_minus = normalize(tau_max - tau_min)
    tasks.append(METQFITask(
        "MET", "MET-2", "qfi_extremal_superposition", n, seeds[1], scr,
        gmax, tau_max, generator, gmin, gmax, span, tau_plus, tau_minus
    ))

    # MET-3: two-target balanced Fisher soft-min.
    rng = np.random.default_rng(seeds[2])
    scr = RealScrambler(n, depth, 200000 + seeds[2])
    x1, x2 = random_distinct_labels(rng, dim, 2)
    tau1 = scr.scrambled_basis(x1)
    tau2 = scr.scrambled_basis(x2)
    tasks.append(METBalancedTask("MET", "MET-3", "two_target_balanced_fisher", n, seeds[2], scr, x1, tau1, x2, tau2, 20.0))
    return tasks


# -----------------------------------------------------------------------------
# Optimizer state classes and sphere geometry helpers
# -----------------------------------------------------------------------------


@dataclass
class GeoCGState:
    prev_x: Optional[np.ndarray] = None
    prev_grad: Optional[np.ndarray] = None
    prev_dir: Optional[np.ndarray] = None


@dataclass
class LBFGSState:
    """Limited-memory inverse Hessian approximation in transported tangent spaces."""
    memory: int = 7
    s_list: List[np.ndarray] = field(default_factory=list)
    y_list: List[np.ndarray] = field(default_factory=list)
    rho_list: List[float] = field(default_factory=list)

    def transport_memory(self, x_old: np.ndarray, x_new: np.ndarray) -> None:
        if not self.s_list:
            return
        self.s_list = [transport_exact(x_old, x_new, s) for s in self.s_list]
        self.y_list = [transport_exact(x_old, x_new, y) for y in self.y_list]
        self._refresh_rhos()

    def _refresh_rhos(self) -> None:
        new_s: List[np.ndarray] = []
        new_y: List[np.ndarray] = []
        new_rho: List[float] = []
        for s, y in zip(self.s_list, self.y_list):
            sy = float(np.dot(s, y))
            if sy > 1e-12:
                new_s.append(s)
                new_y.append(y)
                new_rho.append(1.0 / sy)
        self.s_list = new_s[-self.memory:]
        self.y_list = new_y[-self.memory:]
        self.rho_list = new_rho[-self.memory:]

    def direction(self, x: np.ndarray, grad: np.ndarray) -> np.ndarray:
        q = tangent_project(x, grad.copy())
        alpha_vals: List[float] = []
        for s, y, rho in reversed(list(zip(self.s_list, self.y_list, self.rho_list))):
            a = rho * float(np.dot(s, q))
            alpha_vals.append(a)
            q = q - a * y
        if self.s_list:
            sy = float(np.dot(self.s_list[-1], self.y_list[-1]))
            yy = float(np.dot(self.y_list[-1], self.y_list[-1]))
            gamma = sy / yy if yy > 1e-15 else 1.0
        else:
            gamma = 1.0
        r = gamma * q
        for (s, y, rho), a in zip(zip(self.s_list, self.y_list, self.rho_list), reversed(alpha_vals)):
            b = rho * float(np.dot(y, r))
            r = r + s * (a - b)
        return -tangent_project(x, r)

    def update_memory(
        self,
        x_old: np.ndarray,
        x_new: np.ndarray,
        g_old: np.ndarray,
        g_new: np.ndarray,
        step_vec_old: np.ndarray,
    ) -> None:
        # Existing pairs were represented in T_{x_old}S. Move them to T_{x_new}S.
        self.transport_memory(x_old, x_new)
        s_vec = tangent_project(x_new, transport_exact(x_old, x_new, step_vec_old))
        y_vec = tangent_project(x_new, g_new - transport_exact(x_old, x_new, g_old))
        sy = float(np.dot(s_vec, y_vec))
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


def transport_exact(x_old: np.ndarray, x_new: np.ndarray, v: np.ndarray, *, eps: float = 1e-12) -> np.ndarray:
    """Exact parallel transport on the unit sphere along the short geodesic."""
    x_old = unit_normalize(x_old)
    x_new = unit_normalize(x_new)
    v = tangent_project(x_old, v)
    denom = 1.0 + float(np.dot(x_old, x_new))
    if abs(denom) < eps:
        return tangent_project(x_new, v)
    transported = v - (float(np.dot(x_new, v)) / denom) * (x_old + x_new)
    return tangent_project(x_new, transported)


def hopf_metric_diagonal_real_fast(theta: np.ndarray) -> np.ndarray:
    """O(2^n) diagonal Hopf metric for the real chart."""
    theta = np.asarray(theta, dtype=float)
    L = theta.size
    n = int(round(math.log2(L + 1)))
    if (1 << n) - 1 != L:
        raise ValueError("Real Hopf theta length must be 2^n - 1.")
    g = np.ones(L, dtype=float)
    c2 = np.cos(theta) ** 2
    s2 = np.sin(theta) ** 2
    for i in range(1, L + 1):
        left = 2 * i
        right = left + 1
        if left <= L:
            g[left - 1] = g[i - 1] * c2[i - 1]
        if right <= L:
            g[right - 1] = g[i - 1] * s2[i - 1]
    return g


def hopf_state_gradient_from_metric(
    theta: np.ndarray,
    grad_theta: np.ndarray,
    metric_diag: np.ndarray,
    *,
    metric_eps: float = 1e-12,
) -> np.ndarray:
    """Build J g^{-1} grad_theta without materializing the dense Jacobian."""
    theta = np.asarray(theta, dtype=float)
    grad_theta = np.asarray(grad_theta, dtype=float)
    metric_diag = np.asarray(metric_diag, dtype=float)
    L = theta.size
    n = int(round(math.log2(L + 1)))
    psi = hopf.vector_from_theta(theta, "real")
    safe_metric = np.where(metric_diag > metric_eps, metric_diag, metric_eps)
    coeff = grad_theta / safe_metric
    out = np.zeros_like(psi)
    eps = 1e-14

    stack: List[Tuple[int, float]] = [(1, 0.0)]
    while stack:
        node, acc = stack.pop()
        j0 = node - 1
        start, mid, end = hopf._subtree_span(node, n)
        s = math.sin(float(theta[j0]))
        c = math.cos(float(theta[j0]))
        left_acc = acc + (-(s / c) * coeff[j0] if abs(c) > eps else 0.0)
        right_acc = acc + ((c / s) * coeff[j0] if abs(s) > eps else 0.0)
        left_child = 2 * node
        right_child = left_child + 1
        if left_child <= L:
            stack.append((left_child, left_acc))
        else:
            out[start:mid] = psi[start:mid] * left_acc
        if right_child <= L:
            stack.append((right_child, right_acc))
        else:
            out[mid:end] = psi[mid:end] * right_acc
    return tangent_project(psi, out)


def strong_wolfe_geodesic_line_search(
    task: BaseTask,
    x: np.ndarray,
    direction: np.ndarray,
    grad0: np.ndarray,
    args: argparse.Namespace,
    *,
    initial_alpha: float = 1.0,
) -> Tuple[np.ndarray, float, int, np.ndarray, float]:
    """Strong-Wolfe geodesic line search with bounded geodesic fallback."""
    p = tangent_project(x, direction)
    p_norm = float(np.linalg.norm(p))
    if p_norm < 1e-15:
        return x.copy(), 0.0, 0, np.zeros_like(x), 0.0
    grad0 = tangent_project(x, grad0)
    dphi0 = float(np.dot(grad0, p))
    if dphi0 >= 0.0:
        p = -grad0
        p_norm = float(np.linalg.norm(p))
        if p_norm < 1e-15:
            return x.copy(), 0.0, 0, np.zeros_like(x), 0.0
        dphi0 = -float(np.dot(grad0, grad0))

    f0 = float(task.cost(x))
    alpha_max = max(0.0, float(args.max_line_angle) / p_norm)
    if alpha_max <= 0.0:
        return x.copy(), 0.0, 0, np.zeros_like(x), 0.0
    alpha = min(max(float(initial_alpha), args.line_alpha_min), alpha_max)
    c1 = float(args.wolfe_c1)
    c2 = float(args.wolfe_c2)
    nfev = 0
    best_y = x.copy()
    best_f = f0
    best_alpha = 0.0
    lo = 0.0
    hi = alpha_max
    prev_alpha = 0.0
    prev_f = f0

    for _ in range(int(args.line_maxiter)):
        y = sphere_exp(x, p, alpha)
        f = float(task.cost(y))
        nfev += 1
        if np.isfinite(f) and f < best_f:
            best_f = f
            best_y = y
            best_alpha = alpha
        if (not np.isfinite(f)) or f > f0 + c1 * alpha * dphi0 or (alpha > 0 and f >= prev_f and prev_alpha > 0):
            hi = alpha
            alpha = 0.5 * (lo + hi)
            prev_alpha, prev_f = alpha, f
            continue
        _theta_y, _grad_theta_y, grad_y = hopf_params_and_grads(task, y, metric_eps=args.metric_eps)
        nfev += 1
        p_y = transport_exact(x, y, p)
        dphi = float(np.dot(grad_y, p_y))
        if abs(dphi) <= c2 * abs(dphi0):
            return y, alpha * p_norm, nfev, alpha * p, alpha
        if dphi >= 0.0:
            hi = alpha
            alpha = 0.5 * (lo + hi)
        else:
            lo = alpha
            if hi < alpha_max:
                alpha = 0.5 * (lo + hi)
            else:
                new_alpha = min(2.0 * alpha, alpha_max)
                if new_alpha <= alpha + 1e-15:
                    break
                alpha = new_alpha
        prev_alpha, prev_f = alpha, f

    # Safeguard: same exact geodesic, cost-only bounded minimization.
    q = p / p_norm

    def f_scalar(ell: float) -> float:
        return float(task.cost(sphere_exp_unit(x, q, float(ell))))

    try:
        res = minimize_scalar(
            f_scalar,
            bounds=(0.0, float(args.max_line_angle)),
            method="bounded",
            options={"xatol": 1e-8, "maxiter": 80},
        )
        nfev += int(getattr(res, "nfev", 0))
        ell = float(res.x)
        if np.isfinite(ell) and f_scalar(ell) <= best_f + 1e-15:
            y = sphere_exp_unit(x, q, ell)
            alpha_eff = ell / p_norm
            return y, ell, nfev, alpha_eff * p, alpha_eff
    except Exception:
        pass
    if best_alpha > 0.0:
        return best_y, best_alpha * p_norm, nfev, best_alpha * p, best_alpha
    return x.copy(), 0.0, nfev, np.zeros_like(x), 0.0
# -----------------------------------------------------------------------------
# Per-step record construction
# -----------------------------------------------------------------------------


COMMON_FIELDS = [
    "task", "n", "task_id", "task_name", "problem_seed", "seed_index", "run_seed", "mode", "step",
    "scramble_depth", "scramble_seed", "coordinate_type", "num_parameters", "dimension",
    "cost", "grad_norm", "state_grad_norm", "state_norm_error",
    "target_label", "target_bitstring", "target_overlap", "decoded_argmax", "decoded_bitstring", "decoded_probability",
    "distractor_label", "distractor_overlap",
    "energy", "energy_gap", "energy_per_clause", "ground_energy_known",
    "met_objective", "cfi", "cfi_gap", "A_expectation", "qfi", "normalized_qfi", "F1", "F2", "F_bal", "overlap_plus", "overlap_minus", "balanced_target_overlap",
    "last_step_angle", "last_subspace_dim", "last_line_evals", "wall_time_sec",
    "theta", "grad",
]


def make_record(
    task: BaseTask,
    mode: str,
    step: int,
    coordinate_type: str,
    params: np.ndarray,
    grad: np.ndarray,
    psi: np.ndarray,
    state_grad: np.ndarray,
    elapsed: float,
    *,
    run_seed: int | str = "",
    seed_index: int | str = "",
    last_step_angle: float | str = "",
    last_subspace_dim: int | str = "",
    last_line_evals: int | str = "",
    precision: int = 10,
) -> Dict[str, object]:
    metrics = task.metrics(psi)
    row: Dict[str, object] = {k: "" for k in COMMON_FIELDS}
    row.update({
        "task": task.task_type,
        "n": task.n,
        "task_id": task.task_id,
        "task_name": task.name,
        "problem_seed": task.problem_seed,
        "seed_index": seed_index,
        "run_seed": run_seed,
        "mode": mode,
        "step": step,
        "scramble_depth": task.scrambler.depth,
        "scramble_seed": task.scrambler.seed,
        "coordinate_type": coordinate_type,
        "num_parameters": int(params.size),
        "dimension": int(1 << task.n),
        "cost": float(task.cost(psi)),
        "grad_norm": float(np.linalg.norm(grad)),
        "state_grad_norm": float(np.linalg.norm(state_grad)),
        "state_norm_error": float(abs(np.linalg.norm(psi) - 1.0)),
        "ground_energy_known": 0.0 if task.task_type == "VQE" else "",
        "last_step_angle": last_step_angle,
        "last_subspace_dim": last_subspace_dim,
        "last_line_evals": last_line_evals,
        "wall_time_sec": float(elapsed),
        "theta": vector_to_string(params, precision),
        "grad": vector_to_string(grad, precision),
    })
    row.update(metrics)
    return row


# -----------------------------------------------------------------------------
# Runner functions
# -----------------------------------------------------------------------------


def initial_state_for_task(task: BaseTask, run_seed: Optional[int] = None) -> np.ndarray:
    seed = task.problem_seed + 77777 if run_seed is None else int(run_seed)
    return hopf.canonical_initial_state(task.n, seed=seed)


def run_seeds_for_task(task: BaseTask, args: argparse.Namespace) -> List[int]:
    if int(args.num_seeds) < 1:
        raise ValueError("--num-seeds must be >= 1")
    return [int(task.problem_seed) + int(args.seed_offset) + i for i in range(int(args.num_seeds))]


def hopf_params_and_grads(
    task: BaseTask,
    psi: np.ndarray,
    *,
    metric_eps: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta = theta_from_state_safe(psi, task.n)
    psi_h = unit_normalize(hopf.vector_from_theta(theta, "real"))
    grad_e = task.euclidean_grad(psi_h)
    grad_theta = hopf_coordinate_gradient(theta, grad_e)
    metric_diag = hopf_metric_diagonal_real_fast(theta)
    state_grad = hopf_state_gradient_from_metric(theta, grad_theta, metric_diag, metric_eps=metric_eps)
    return theta, grad_theta, state_grad


def run_hopf_egt_cg(task: BaseTask, steps: int, args: argparse.Namespace, *, run_seed: int, seed_index: int) -> List[Dict[str, object]]:
    mode = "Hopf-EGT-CG"
    rows: List[Dict[str, object]] = []
    start_time = time.time()
    psi = normalize(initial_state_for_task(task, run_seed=run_seed))
    state = GeoCGState()
    last_angle: float | str = ""
    last_evals: int | str = ""
    alpha_guess = 1.0
    for step in range(steps + 1):
        theta, grad_theta, grad_s = hopf_params_and_grads(task, psi, metric_eps=args.metric_eps)
        psi = unit_normalize(hopf.vector_from_theta(theta, "real"))
        rows.append(make_record(task, mode, step, "hopf_real", theta, grad_theta, psi, grad_s, time.time() - start_time, last_step_angle=last_angle, last_line_evals=last_evals, run_seed=run_seed, seed_index=seed_index, precision=args.precision))
        if step == steps:
            break
        if state.prev_grad is None or state.prev_x is None or state.prev_dir is None:
            direction = -grad_s
        else:
            g_prev_t = transport_exact(state.prev_x, psi, state.prev_grad)
            p_prev_t = transport_exact(state.prev_x, psi, state.prev_dir)
            yk = tangent_project(psi, grad_s - g_prev_t)
            denom = float(np.dot(p_prev_t, yk))
            if denom > 1e-14:
                beta_hs = float(np.dot(grad_s, yk)) / denom
                beta_dy = float(np.dot(grad_s, grad_s)) / denom
                beta = max(0.0, min(beta_hs, beta_dy))
            else:
                beta = 0.0
            direction = tangent_project(psi, -grad_s + beta * p_prev_t)
            if float(np.dot(direction, grad_s)) >= -1e-14 or np.linalg.norm(direction) < 1e-12:
                direction = -grad_s
        psi_old = psi.copy()
        grad_old = grad_s.copy()
        dir_old = direction.copy()
        psi, last_angle, last_evals, _step_vec, alpha_eff = strong_wolfe_geodesic_line_search(
            task, psi, direction, grad_s, args, initial_alpha=alpha_guess
        )
        state.prev_x = psi_old
        state.prev_grad = grad_old
        state.prev_dir = dir_old
        alpha_guess = alpha_eff if alpha_eff > 0.0 else 1.0
    return rows


def run_hopf_lbfgs(task: BaseTask, steps: int, args: argparse.Namespace, *, run_seed: int, seed_index: int) -> List[Dict[str, object]]:
    mode = "Hopf-Riemannian-LBFGS"
    rows: List[Dict[str, object]] = []
    start_time = time.time()
    psi = normalize(initial_state_for_task(task, run_seed=run_seed))
    opt = LBFGSState(memory=args.lbfgs_memory)
    last_angle: float | str = ""
    last_evals: int | str = ""
    pending: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None
    alpha_guess = 1.0
    for step in range(steps + 1):
        theta, grad_theta, grad_s = hopf_params_and_grads(task, psi, metric_eps=args.metric_eps)
        psi = unit_normalize(hopf.vector_from_theta(theta, "real"))
        if pending is not None:
            prev_x, prev_grad, prev_step_vec = pending
            opt.update_memory(prev_x, psi, prev_grad, grad_s, prev_step_vec)
            pending = None
        rows.append(make_record(task, mode, step, "hopf_real", theta, grad_theta, psi, grad_s, time.time() - start_time, last_step_angle=last_angle, last_line_evals=last_evals, run_seed=run_seed, seed_index=seed_index, precision=args.precision))
        if step == steps:
            break
        direction = opt.direction(psi, grad_s)
        if float(np.dot(direction, grad_s)) >= -1e-14 or np.linalg.norm(direction) < 1e-12:
            direction = -grad_s
        psi_old = psi.copy()
        grad_old = grad_s.copy()
        psi, last_angle, last_evals, step_vec, alpha_eff = strong_wolfe_geodesic_line_search(
            task, psi, direction, grad_s, args, initial_alpha=alpha_guess
        )
        pending = (psi_old, grad_old, step_vec)
        alpha_guess = alpha_eff if alpha_eff > 0.0 else 1.0
    return rows


def bb_spectral_alpha(state: BBState, x: np.ndarray, grad: np.ndarray, args: argparse.Namespace) -> float:
    if state.prev_x is None or state.prev_grad is None or state.prev_step_vec is None:
        return float(args.bb_initial_alpha)
    s_vec = transport_exact(state.prev_x, x, state.prev_step_vec)
    y_vec = tangent_project(x, grad - transport_exact(state.prev_x, x, state.prev_grad))
    sy = float(np.dot(s_vec, y_vec))
    yy = float(np.dot(y_vec, y_vec))
    ss = float(np.dot(s_vec, s_vec))
    if sy <= 1e-14 or yy <= 1e-14 or ss <= 1e-14:
        return float(state.alpha)
    if args.bb_variant == "bb2" or (args.bb_variant == "alternate" and state.iteration % 2 == 1):
        alpha = sy / yy
    else:
        alpha = ss / sy
    return float(np.clip(alpha, args.bb_min_alpha, args.bb_max_alpha))


def run_hopf_bb(task: BaseTask, steps: int, args: argparse.Namespace, *, run_seed: int, seed_index: int) -> List[Dict[str, object]]:
    mode = "Hopf-Riemannian-BB"
    rows: List[Dict[str, object]] = []
    start_time = time.time()
    psi = normalize(initial_state_for_task(task, run_seed=run_seed))
    state = BBState(alpha=float(args.bb_initial_alpha))
    last_angle: float | str = ""
    last_evals: int | str = ""
    for step in range(steps + 1):
        theta, grad_theta, grad_s = hopf_params_and_grads(task, psi, metric_eps=args.metric_eps)
        psi = unit_normalize(hopf.vector_from_theta(theta, "real"))
        rows.append(make_record(task, mode, step, "hopf_real", theta, grad_theta, psi, grad_s, time.time() - start_time, last_step_angle=last_angle, last_line_evals=last_evals, run_seed=run_seed, seed_index=seed_index, precision=args.precision))
        if step == steps:
            break
        alpha_guess = bb_spectral_alpha(state, psi, grad_s, args)
        direction = -grad_s
        psi_old = psi.copy()
        grad_old = grad_s.copy()
        psi, last_angle, last_evals, step_vec, alpha_eff = strong_wolfe_geodesic_line_search(
            task, psi, direction, grad_s, args, initial_alpha=alpha_guess
        )
        state.prev_x = psi_old
        state.prev_grad = grad_old
        state.prev_step_vec = step_vec
        state.alpha = alpha_eff if alpha_eff > 0.0 else alpha_guess
        state.iteration += 1
    return rows


def run_task_modes(task: BaseTask, steps: int, args: argparse.Namespace) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    seeds = run_seeds_for_task(task, args)
    for seed_index, run_seed in enumerate(seeds):
        seed_msg = f"seed {seed_index + 1}/{len(seeds)} run_seed={run_seed}"
        print(f"[{task.task_type}] {task.task_id} {task.name} ({seed_msg}): running Hopf-EGT-CG")
        rows.extend(run_hopf_egt_cg(task, steps, args, run_seed=run_seed, seed_index=seed_index))
        print(f"[{task.task_type}] {task.task_id} {task.name} ({seed_msg}): running Hopf-Riemannian-LBFGS")
        rows.extend(run_hopf_lbfgs(task, steps, args, run_seed=run_seed, seed_index=seed_index))
        print(f"[{task.task_type}] {task.task_id} {task.name} ({seed_msg}): running Hopf-Riemannian-BB")
        rows.extend(run_hopf_bb(task, steps, args, run_seed=run_seed, seed_index=seed_index))
    return rows
# -----------------------------------------------------------------------------
# CSV writing and CLI
# -----------------------------------------------------------------------------


def write_csv(path: str, rows: Iterable[Dict[str, object]]) -> int:
    count = 0
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COMMON_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Synthetic Hopf cost-and-gradient optimizer stress-test data generator.")
    parser.add_argument("--n", type=int, required=True, help="Number of qubits.")
    parser.add_argument("--steps", type=int, default=200, help="Optimization steps; step 0 is also recorded.")
    parser.add_argument("--scramble-depth", type=int, default=4, help="Depth of the real orthogonal scrambler.")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory.")
    parser.add_argument("--precision", type=int, default=10, help="Significant digits for theta/grad vector columns.")
    parser.add_argument("--num-seeds", type=int, default=10, help="Number of deterministic initial-state seeds per synthetic task.")
    parser.add_argument("--seed-offset", type=int, default=77777, help="Initial-state seed offset; run_seed = problem_seed + seed_offset + seed_index.")
    parser.add_argument("--max-line-angle", type=float, default=math.pi / 2.0, help="Maximum geodesic line-search angle.")
    parser.add_argument("--line-alpha-min", type=float, default=1e-8, help="Minimum trial learning-rate scale for geodesic line search.")
    parser.add_argument("--line-maxiter", type=int, default=25, help="Maximum strong-Wolfe backtracking trials before bounded fallback.")
    parser.add_argument("--wolfe-c1", type=float, default=0.485, help="Armijo parameter for the strong-Wolfe line search; default follows the EGT-CG paper.")
    parser.add_argument("--wolfe-c2", type=float, default=0.999, help="Curvature parameter for the strong-Wolfe line search; default follows the EGT-CG paper.")
    parser.add_argument("--metric-eps", type=float, default=1e-12, help="Regularization floor for diagonal Hopf metric divisions.")
    parser.add_argument("--lbfgs-memory", type=int, default=7)
    parser.add_argument("--bb-initial-alpha", type=float, default=1.0, help="Initial BB spectral learning-rate scale.")
    parser.add_argument("--bb-min-alpha", type=float, default=1e-6, help="Minimum BB spectral learning-rate scale.")
    parser.add_argument("--bb-max-alpha", type=float, default=10.0, help="Maximum BB spectral learning-rate scale.")
    parser.add_argument("--bb-variant", choices=["bb1", "bb2", "alternate"], default="alternate", help="Barzilai--Borwein spectral step variant.")
    parser.add_argument("--only", choices=["all", "vqe", "met"], default="all", help="Generate only one task family or both.")
    parser.add_argument("--quick", action="store_true", help="Shortcut: --steps 2 --scramble-depth min(depth, 1).")
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    if args.quick:
        args.steps = min(args.steps, 2)
        args.scramble_depth = min(args.scramble_depth, 1)
    if args.n < 1:
        raise ValueError("--n must be >= 1")
    if args.num_seeds < 1:
        raise ValueError("--num-seeds must be >= 1")
    os.makedirs(args.outdir, exist_ok=True)

    if args.only in ("all", "vqe"):
        vqe_tasks = make_vqe_tasks(args.n, args.scramble_depth)
        all_vqe_rows: List[Dict[str, object]] = []
        for task in vqe_tasks:
            all_vqe_rows.extend(run_task_modes(task, args.steps, args))
        vqe_path = os.path.join(args.outdir, f"vqe_hopf_data_n{args.n}.csv")
        nrows = write_csv(vqe_path, all_vqe_rows)
        print(f"Wrote {nrows} rows to {vqe_path}")

    if args.only in ("all", "met"):
        met_tasks = make_met_tasks(args.n, args.scramble_depth)
        all_met_rows: List[Dict[str, object]] = []
        for task in met_tasks:
            all_met_rows.extend(run_task_modes(task, args.steps, args))
        met_path = os.path.join(args.outdir, f"met_hopf_data_n{args.n}.csv")
        nrows = write_csv(met_path, all_met_rows)
        print(f"Wrote {nrows} rows to {met_path}")


if __name__ == "__main__":
    main()
