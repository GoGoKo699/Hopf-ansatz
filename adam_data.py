#!/usr/bin/env python3
"""
adam_data.py

Adaptive-line-search Adam baseline synthetic full-real-state stress tests for the Hopf ansatz and
Möttönen physical-angle parameterization.

Run:
    python adam_data.py --n 8

Outputs exactly two CSV files by default:
    vqe_adam_data_n{n}.csv
    met_adam_data_n{n}.csv

Each CSV records 3 generated tasks x 10 initial-state seeds x 2 Adam baseline modes x (steps + 1) rows by default.
Each baseline is one adaptive-line-search Adam mode; no learning-rate grid is emitted.
The parameter vector and gradient vector are stored as semicolon-separated strings.

The script assumes that hopf_utils.py is in the same folder or on PYTHONPATH.
It uses the real Hopf chart only.

Modes:
    VQE and MET:
      1. Hopf-Adam
      2. Mottonen-ideal-PS-Adam

The Hopf baseline applies Adam to the Hopf coordinates. The Möttönen baseline
uses the physical post-multiplexing rotation angles and an ideal exact gradient
equivalent to infinite-shot parameter-shift by default. These baselines are
intended for GitHub/reference diagnostics rather than for the paper's
geometry-native cost-and-gradient optimizer section. The Adam step uses cost-only
backtracking, so its quantum interface is still cost evaluations plus gradients.

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
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

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
    # Fix a deterministic global sign for more stable inverse coordinates.
    # This is not physically needed; it only makes CSV traces easier to compare.
    idx = int(np.argmax(np.abs(y)))
    if y[idx] < 0:
        y = -y
    return y


def tangent_project(x: np.ndarray, v: np.ndarray) -> np.ndarray:
    return np.asarray(v, dtype=float) - float(np.dot(x, v)) * np.asarray(x, dtype=float)


def sphere_exp(x: np.ndarray, p: np.ndarray, alpha: float) -> np.ndarray:
    """Exact exponential map on the unit sphere along tangent vector p."""
    p = tangent_project(x, p)
    nrm = float(np.linalg.norm(p))
    if nrm < 1e-15 or abs(alpha) < 1e-15:
        return np.asarray(x, dtype=float).copy()
    ell = alpha * nrm
    y = math.cos(ell) * x + math.sin(ell) * (p / nrm)
    return normalize(y)


def sphere_exp_unit(x: np.ndarray, q: np.ndarray, ell: float) -> np.ndarray:
    """Sphere geodesic with unit tangent q and angular displacement ell."""
    q = tangent_project(x, q)
    nrm = float(np.linalg.norm(q))
    if nrm < 1e-15 or abs(ell) < 1e-15:
        return np.asarray(x, dtype=float).copy()
    q = q / nrm
    return normalize(math.cos(ell) * x + math.sin(ell) * q)


def transport_project(x_old: np.ndarray, x_new: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Simple projection vector transport to T_{x_new}S.

    Kept only as a numerical fallback. The Hopf EGT-CG and Riemannian L-BFGS
    modes below use exact sphere geodesic transport.
    """
    return tangent_project(x_new, v)


def transport_exact_geodesic(
    x: np.ndarray,
    direction: np.ndarray,
    eta: float,
    zeta: np.ndarray,
    *,
    eps: float = 1e-15,
) -> np.ndarray:
    """Exact parallel transport on the unit sphere along Exp_x(eta * direction).

    This is the sphere formula used by exact geodesic transport.  The vector
    ``direction`` and the transported vector ``zeta`` both live in T_x S.  The
    returned vector lives in T_{Exp_x(eta direction)} S.

    T_{eta u}(zeta) = zeta
        - sin(eta ||u||) <u,zeta>/||u|| x
        + (cos(eta ||u||)-1) <u,zeta>/||u||^2 u.
    """
    x = np.asarray(x, dtype=float)
    u = tangent_project(x, direction)
    z = tangent_project(x, zeta)
    un = float(np.linalg.norm(u))
    if un < eps or abs(float(eta)) < eps:
        y = x.copy()
        return tangent_project(y, z)
    inner = float(np.dot(u, z))
    ell = float(eta) * un
    transported = z - math.sin(ell) * (inner / un) * x + (math.cos(ell) - 1.0) * (inner / (un * un)) * u
    y = sphere_exp(x, u, eta)
    return tangent_project(y, transported)


def hopf_metric_diagonal_real(theta: np.ndarray) -> np.ndarray:
    """Efficient real-Hopf diagonal metric.

    The root entry is one.  Each child inherits the parent subtree weight times
    the appropriate branch probability.  This is O(2^n), not O(2^n n).
    """
    theta = np.asarray(theta, dtype=float)
    L = theta.size
    g = np.ones(L, dtype=float)
    if L == 0:
        return g
    c2 = np.cos(theta) ** 2
    s2 = np.sin(theta) ** 2
    for idx in range(2, L + 1):  # one-based tree index
        parent = idx // 2
        if idx % 2 == 0:
            g[idx - 1] = g[parent - 1] * c2[parent - 1]
        else:
            g[idx - 1] = g[parent - 1] * s2[parent - 1]
    return g


def hopf_state_gradient_from_coordinate_gradient(
    theta: np.ndarray,
    grad_theta: np.ndarray,
    metric_diag: np.ndarray,
    *,
    eps: float = 1e-14,
) -> np.ndarray:
    """Lift Hopf coordinate gradients to the state-sphere gradient.

    Computes J g^{-1} grad_theta without materializing J.  For each leaf, only
    the n nodes on its root-to-leaf path contribute, so the cost is O(n 2^n).
    """
    theta = np.asarray(theta, dtype=float)
    grad_theta = np.asarray(grad_theta, dtype=float)
    metric_diag = np.asarray(metric_diag, dtype=float)
    psi = hopf.vector_from_theta(theta, "real")
    L = theta.size
    n = int(round(math.log2(L + 1)))
    dim = 1 << n
    safe_metric = np.maximum(metric_diag, eps)
    nat = grad_theta / safe_metric
    out = np.zeros(dim, dtype=float)
    for leaf in range(dim):
        node = 1
        coeff = 0.0
        for depth in range(n):
            j0 = node - 1
            th = float(theta[j0])
            s = math.sin(th)
            c = math.cos(th)
            bit = (leaf >> (n - 1 - depth)) & 1
            if bit == 0:
                if abs(c) > eps:
                    coeff += -nat[j0] * (s / c)
                node = 2 * node
            else:
                if abs(s) > eps:
                    coeff += nat[j0] * (c / s)
                node = 2 * node + 1
        out[leaf] = psi[leaf] * coeff
    return tangent_project(psi, out)


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
    th = hopf.theta_from_vector(normalize(psi), "real")
    return hopf.clip_theta_hopf_real(th, n=n)


# -----------------------------------------------------------------------------
# Fast Walsh-Hadamard transform for Möttönen physical-angle conversion
# -----------------------------------------------------------------------------


def fwht(x: np.ndarray) -> np.ndarray:
    """Unnormalized Walsh-Hadamard transform in natural order."""
    y = np.asarray(x, dtype=float).copy()
    h = 1
    n = y.size
    while h < n:
        for i in range(0, n, h * 2):
            a = y[i : i + h].copy()
            b = y[i + h : i + 2 * h].copy()
            y[i : i + h] = a + b
            y[i + h : i + 2 * h] = a - b
        h *= 2
    return y


def hopf_theta_to_mottonen_phi(alpha: np.ndarray, n: int) -> np.ndarray:
    """Convert logical tree/Hopf split angles alpha to physical multiplexed Ry angles."""
    alpha = np.asarray(alpha, dtype=float)
    out = np.empty_like(alpha)
    base = 0
    for level in range(n):
        k = 1 << level
        block = alpha[base : base + k]
        if level == 0:
            out[base] = 2.0 * block[0]
        else:
            out[base : base + k] = (2.0 ** (1 - level)) * fwht(block)
        base += k
    return out


def mottonen_phi_to_hopf_theta(phi: np.ndarray, n: int) -> np.ndarray:
    """Convert physical multiplexed Ry angles to logical tree/Hopf split angles."""
    phi = np.asarray(phi, dtype=float)
    out = np.empty_like(phi)
    base = 0
    for level in range(n):
        k = 1 << level
        block = phi[base : base + k]
        if level == 0:
            out[base] = 0.5 * block[0]
        else:
            out[base : base + k] = 0.5 * fwht(block)
        base += k
    return out


def mottonen_state_from_phi(phi: np.ndarray, n: int) -> np.ndarray:
    alpha = mottonen_phi_to_hopf_theta(phi, n)
    return normalize(hopf.vector_from_theta(alpha, "real"))


def mottonen_phi_from_state(psi: np.ndarray, n: int) -> np.ndarray:
    alpha = theta_from_state_safe(psi, n)
    return hopf_theta_to_mottonen_phi(alpha, n)


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
# Optimizer state
# -----------------------------------------------------------------------------


@dataclass
class AdamState:
    """Adam moments plus an online cost-only backtracking step length.

    The moment update is Adam. The scalar step length is adapted by evaluating
    the objective along the Adam direction. This avoids a hand-tuned global
    learning rate while keeping the same quantum information interface as the
    geometry-native optimizers: cost evaluations and gradients only.
    """

    lr: float = 0.03
    beta1: float = 0.9
    beta2: float = 0.999
    eps: float = 1e-8
    m: Optional[np.ndarray] = None
    v: Optional[np.ndarray] = None
    t: int = 0

    def direction(self, grad: np.ndarray, *, normalize_rms: bool = True) -> np.ndarray:
        grad = np.asarray(grad, dtype=float)
        if self.m is None:
            self.m = np.zeros_like(grad)
            self.v = np.zeros_like(grad)
        self.t += 1
        self.m = self.beta1 * self.m + (1.0 - self.beta1) * grad
        self.v = self.beta2 * self.v + (1.0 - self.beta2) * (grad * grad)
        mhat = self.m / (1.0 - self.beta1 ** self.t)
        vhat = self.v / (1.0 - self.beta2 ** self.t)
        d = -mhat / (np.sqrt(vhat) + self.eps)

        # Momentum can occasionally turn the Adam direction non-descent. Restart
        # the first moment for this direction only; keep the second-moment scale.
        if float(np.dot(grad, d)) >= 0.0:
            d = -grad / (np.sqrt(vhat) + self.eps)
        if float(np.dot(grad, d)) >= 0.0:
            d = -grad.copy()

        if normalize_rms:
            rms = float(np.linalg.norm(d) / math.sqrt(max(1, d.size)))
            if rms > 1e-15:
                d = d / rms
        return d

# -----------------------------------------------------------------------------


def state_angle(x: np.ndarray, y: np.ndarray) -> float:
    dot = float(np.dot(normalize(x), normalize(y)))
    dot = max(-1.0, min(1.0, dot))
    return float(math.acos(dot))


def cost_line_search_update(
    params: np.ndarray,
    direction: np.ndarray,
    current_cost: float,
    cost_from_params,
    project_params,
    state_from_params,
    adam: AdamState,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, float, int, float, float]:
    """Backtracking line search along an Adam coordinate direction.

    Returns (new_params, accepted_lr, number_of_cost_evaluations, state_angle,
    new_cost). The line search is deliberately cost-only: it does not request a
    Hessian, quantum metric tensor, projected Hamiltonian matrix, or any new
    quantum primitive beyond trial cost evaluations.
    """
    params = np.asarray(params, dtype=float)
    direction = np.asarray(direction, dtype=float)
    dnorm = float(np.linalg.norm(direction))
    if dnorm < 1e-15 or not math.isfinite(dnorm):
        return params.copy(), 0.0, 0, 0.0, float(current_cost)

    alpha = min(max(float(adam.lr), args.adam_lr_min), args.adam_lr_max)
    best_params = params.copy()
    best_cost = float(current_cost)
    best_alpha = 0.0
    evals = 0

    # A tiny absolute margin prevents accepting pure roundoff as progress.
    atol = float(args.adam_accept_atol)
    rtol = float(args.adam_accept_rtol)
    threshold = float(current_cost) - max(atol, rtol * max(1.0, abs(float(current_cost))))

    for _ in range(args.adam_line_evals):
        trial = project_params(params + alpha * direction)
        if float(np.linalg.norm(trial - params)) < 1e-18:
            alpha *= args.adam_shrink
            continue
        c = float(cost_from_params(trial))
        evals += 1
        if math.isfinite(c) and c < best_cost:
            best_cost = c
            best_params = trial
            best_alpha = alpha
        if math.isfinite(c) and c <= threshold:
            x0 = state_from_params(params)
            x1 = state_from_params(trial)
            angle = state_angle(x0, x1)
            adam.lr = min(args.adam_lr_max, max(args.adam_lr_min, alpha * args.adam_growth))
            return trial, alpha, evals, angle, c
        alpha *= args.adam_shrink
        if alpha < args.adam_lr_min:
            break

    # If no step meets the strict threshold, accept the best strictly improving
    # trial. Otherwise stay put and reduce the next trial length.
    if best_alpha > 0.0 and best_cost < float(current_cost):
        x0 = state_from_params(params)
        x1 = state_from_params(best_params)
        angle = state_angle(x0, x1)
        adam.lr = min(args.adam_lr_max, max(args.adam_lr_min, best_alpha * args.adam_growth))
        return best_params, best_alpha, evals, angle, best_cost

    adam.lr = max(args.adam_lr_min, alpha)
    return params.copy(), 0.0, evals, 0.0, float(current_cost)


def mottonen_adjoint_ps_equivalent_gradient(task: BaseTask, phi: np.ndarray, n: int) -> np.ndarray:
    """Fast exact gradient w.r.t. physical Möttönen Ry angles.

    This is mathematically equal to the infinite-shot ideal parameter-shift
    gradient, but uses the linear Walsh-Hadamard relation between physical
    multiplexed Ry angles and logical tree/Hopf split angles.
    """
    alpha = mottonen_phi_to_hopf_theta(phi, n)
    psi = normalize(hopf.vector_from_theta(alpha, "real"))
    grad_e = task.euclidean_grad(psi)
    g_alpha = hopf_coordinate_gradient(alpha, grad_e)
    g_phi = np.empty_like(g_alpha)
    base = 0
    for level in range(n):
        k = 1 << level
        block = g_alpha[base : base + k]
        if level == 0:
            g_phi[base] = 0.5 * block[0]
        else:
            g_phi[base : base + k] = 0.5 * fwht(block)
        base += k
    return g_phi

def mottonen_literal_ps_gradient(task: BaseTask, phi: np.ndarray, n: int) -> np.ndarray:
    phi = np.asarray(phi, dtype=float)
    P = phi.size
    grad = np.zeros(P, dtype=float)

    if task.task_type == "VQE":
        for j in range(P):
            pp = phi.copy(); pm = phi.copy()
            pp[j] += math.pi / 2.0
            pm[j] -= math.pi / 2.0
            fp = task.cost(mottonen_state_from_phi(pp, n))
            fm = task.cost(mottonen_state_from_phi(pm, n))
            grad[j] = 0.5 * (fp - fm)
        return grad

    # MET: parameter shift the primitive expectation values, then apply the
    # exact chain rule for the nonlinear Fisher objective at the current state.
    psi0 = mottonen_state_from_phi(phi, n)
    coeff = task.primitive_cost_coefficients(psi0)
    for j in range(P):
        pp = phi.copy(); pm = phi.copy()
        pp[j] += math.pi / 2.0
        pm[j] -= math.pi / 2.0
        prim_p = task.primitive_values(mottonen_state_from_phi(pp, n))
        prim_m = task.primitive_values(mottonen_state_from_phi(pm, n))
        dprim = 0.5 * (prim_p - prim_m)
        grad[j] = float(np.dot(coeff, dprim))
    return grad


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
    "adam_lr", "accepted_lr",
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
    adam_lr: float | str = "",
    accepted_lr: float | str = "",
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
        "adam_lr": adam_lr,
        "accepted_lr": accepted_lr,
        "theta": vector_to_string(params, precision),
        "grad": vector_to_string(grad, precision),
    })
    row.update(metrics)
    return row


def initial_state_for_task(task: BaseTask, run_seed: Optional[int] = None) -> np.ndarray:
    seed = task.problem_seed + 77777 if run_seed is None else int(run_seed)
    return hopf.canonical_initial_state(task.n, seed=seed)


def run_seeds_for_task(task: BaseTask, args: argparse.Namespace) -> List[int]:
    if int(args.num_seeds) < 1:
        raise ValueError("--num-seeds must be >= 1")
    return [int(task.problem_seed) + int(args.seed_offset) + i for i in range(int(args.num_seeds))]


def run_hopf_adam(task: BaseTask, steps: int, args: argparse.Namespace, *, run_seed: int, seed_index: int) -> List[Dict[str, object]]:
    mode = "Hopf-Adam"
    rows: List[Dict[str, object]] = []
    start_time = time.time()
    psi = initial_state_for_task(task, run_seed=run_seed)
    theta = theta_from_state_safe(psi, task.n)
    adam = AdamState(lr=args.adam_lr_init, beta1=args.adam_beta1, beta2=args.adam_beta2, eps=args.adam_eps)
    last_angle: float | str = ""
    last_evals: int | str = ""
    last_accepted_lr: float | str = ""

    def project(th: np.ndarray) -> np.ndarray:
        return hopf.clip_theta_hopf_real(th, n=task.n)

    def state_from(th: np.ndarray) -> np.ndarray:
        return normalize(hopf.vector_from_theta(project(th), "real"))

    def cost_from(th: np.ndarray) -> float:
        return float(task.cost(state_from(th)))

    for step in range(steps + 1):
        theta = project(theta)
        psi = state_from(theta)
        grad_e = task.euclidean_grad(psi)
        grad_theta = hopf_coordinate_gradient(theta, grad_e)
        state_grad = tangent_project(psi, grad_e)
        rows.append(make_record(
            task, mode, step, "hopf_real", theta, grad_theta, psi, state_grad,
            time.time() - start_time,
            last_step_angle=last_angle,
            last_line_evals=last_evals,
            run_seed=run_seed,
            seed_index=seed_index,
            adam_lr=adam.lr,
            accepted_lr=last_accepted_lr,
            precision=args.precision,
        ))
        if step == steps:
            break

        current_cost = float(task.cost(psi))
        direction = adam.direction(grad_theta, normalize_rms=args.adam_normalize_rms)
        theta, acc_lr, evals, angle, _new_cost = cost_line_search_update(
            theta, direction, current_cost, cost_from, project, state_from, adam, args
        )
        last_angle = angle
        last_evals = evals
        last_accepted_lr = acc_lr if acc_lr > 0 else ""
    return rows


def run_mottonen_ps_adam(task: BaseTask, steps: int, args: argparse.Namespace, *, run_seed: int, seed_index: int) -> List[Dict[str, object]]:
    mode = "Mottonen-ideal-PS-Adam"
    rows: List[Dict[str, object]] = []
    start_time = time.time()
    psi0 = initial_state_for_task(task, run_seed=run_seed)
    phi = mottonen_phi_from_state(psi0, task.n)
    if args.wrap_mottonen:
        phi = ((phi + math.pi) % (2.0 * math.pi)) - math.pi
    adam = AdamState(lr=args.adam_lr_init, beta1=args.adam_beta1, beta2=args.adam_beta2, eps=args.adam_eps)
    last_angle: float | str = ""
    last_evals: int | str = ""
    last_accepted_lr: float | str = ""

    def project(ph: np.ndarray) -> np.ndarray:
        if args.wrap_mottonen:
            return ((np.asarray(ph, dtype=float) + math.pi) % (2.0 * math.pi)) - math.pi
        return np.asarray(ph, dtype=float)

    def state_from(ph: np.ndarray) -> np.ndarray:
        return mottonen_state_from_phi(project(ph), task.n)

    def cost_from(ph: np.ndarray) -> float:
        return float(task.cost(state_from(ph)))

    for step in range(steps + 1):
        phi = project(phi)
        psi = state_from(phi)
        if args.mottonen_gradient == "literal":
            grad_phi = mottonen_literal_ps_gradient(task, phi, task.n)
        else:
            grad_phi = mottonen_adjoint_ps_equivalent_gradient(task, phi, task.n)
        grad_e = task.euclidean_grad(psi)
        state_grad = tangent_project(psi, grad_e)
        rows.append(make_record(
            task, mode, step, "mottonen_physical_ry", phi, grad_phi, psi, state_grad,
            time.time() - start_time,
            last_step_angle=last_angle,
            last_line_evals=last_evals,
            run_seed=run_seed,
            seed_index=seed_index,
            adam_lr=adam.lr,
            accepted_lr=last_accepted_lr,
            precision=args.precision,
        ))
        if step == steps:
            break

        current_cost = float(task.cost(psi))
        direction = adam.direction(grad_phi, normalize_rms=args.adam_normalize_rms)
        phi, acc_lr, evals, angle, _new_cost = cost_line_search_update(
            phi, direction, current_cost, cost_from, project, state_from, adam, args
        )
        last_angle = angle
        last_evals = evals
        last_accepted_lr = acc_lr if acc_lr > 0 else ""
    return rows


def run_task_modes(task: BaseTask, steps: int, args: argparse.Namespace) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    seeds = run_seeds_for_task(task, args)
    for seed_index, run_seed in enumerate(seeds):
        seed_msg = f"seed {seed_index + 1}/{len(seeds)} run_seed={run_seed}"
        print(f"[{task.task_type}] {task.task_id} {task.name} ({seed_msg}): running Hopf-Adam")
        rows.extend(run_hopf_adam(task, steps, args, run_seed=run_seed, seed_index=seed_index))
        print(f"[{task.task_type}] {task.task_id} {task.name} ({seed_msg}): running Mottonen-ideal-PS-Adam")
        rows.extend(run_mottonen_ps_adam(task, steps, args, run_seed=run_seed, seed_index=seed_index))
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
    parser = argparse.ArgumentParser(description="Adaptive-line-search Adam-baseline full-real-state stress-test data generator.")
    parser.add_argument("--n", type=int, required=True, help="Number of qubits.")
    parser.add_argument("--steps", type=int, default=200, help="Optimization steps; step 0 is also recorded.")
    parser.add_argument("--scramble-depth", type=int, default=4, help="Depth of the real orthogonal scrambler.")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory.")
    parser.add_argument("--precision", type=int, default=10, help="Significant digits for theta/grad vector columns.")
    parser.add_argument("--num-seeds", type=int, default=10, help="Number of deterministic initial-state seeds per synthetic task.")
    parser.add_argument("--seed-offset", type=int, default=77777, help="Initial-state seed offset; run_seed = problem_seed + seed_offset + seed_index.")
    parser.add_argument("--adam-lr-init", type=float, default=0.03, help="Initial scalar step length for adaptive-line-search Adam.")
    parser.add_argument("--adam-lr-min", type=float, default=1e-8, help="Minimum scalar step length used by Adam backtracking.")
    parser.add_argument("--adam-lr-max", type=float, default=0.5, help="Maximum scalar step length used by Adam backtracking.")
    parser.add_argument("--adam-line-evals", type=int, default=25, help="Maximum trial cost evaluations per Adam step.")
    parser.add_argument("--adam-shrink", type=float, default=0.5, help="Backtracking shrink factor.")
    parser.add_argument("--adam-growth", type=float, default=1.25, help="Growth factor for the next initial step after an accepted move.")
    parser.add_argument("--adam-accept-atol", type=float, default=1e-15, help="Absolute improvement required for immediate line-search acceptance.")
    parser.add_argument("--adam-accept-rtol", type=float, default=0.0, help="Relative improvement required for immediate line-search acceptance.")
    parser.add_argument("--adam-normalize-rms", dest="adam_normalize_rms", action="store_true", default=True, help="Normalize Adam direction to unit RMS before scalar line search.")
    parser.add_argument("--no-adam-normalize-rms", dest="adam_normalize_rms", action="store_false", help="Use the raw Adam direction without RMS normalization.")
    parser.add_argument("--adam-beta1", type=float, default=0.9)
    parser.add_argument("--adam-beta2", type=float, default=0.999)
    parser.add_argument("--adam-eps", type=float, default=1e-8)
    parser.add_argument("--mottonen-gradient", choices=["adjoint", "literal"], default="adjoint",
                        help="Use fast exact adjoint gradient equivalent to infinite-shot parameter shift, or literal two-shift parameter shift.")
    parser.add_argument("--wrap-mottonen", action="store_true", help="Wrap Möttönen physical angles to [-pi,pi) after each accepted step.")
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
        vqe_path = os.path.join(args.outdir, f"vqe_adam_data_n{args.n}.csv")
        nrows = write_csv(vqe_path, all_vqe_rows)
        print(f"Wrote {nrows} rows to {vqe_path}")

    if args.only in ("all", "met"):
        met_tasks = make_met_tasks(args.n, args.scramble_depth)
        all_met_rows: List[Dict[str, object]] = []
        for task in met_tasks:
            all_met_rows.extend(run_task_modes(task, args.steps, args))
        met_path = os.path.join(args.outdir, f"met_adam_data_n{args.n}.csv")
        nrows = write_csv(met_path, all_met_rows)
        print(f"Wrote {nrows} rows to {met_path}")


if __name__ == "__main__":
    main()
