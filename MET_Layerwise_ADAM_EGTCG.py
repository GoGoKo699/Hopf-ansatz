r"""MET_Layerwise_ADAM_EGTCG.py
================================================================================
Numerical fixed-readout Ramsey-metrology panel for the real Hopf ansatz.

This is the metrological analogue of VQE_Layerwise_ADAM_EGTCG.py. The Hopf
state is used as a probe for common-phase Ramsey sensing with

    U_phi = exp[-i phi G],       G = (1/2) sum_j Z_j,

followed by H^{\otimes n} and computational-basis measurement. The full
bitstring measurement is coarse grained to even/odd parity, giving a specified
binary readout. The optimized objective is the regularized inverse classical
Fisher information of this fixed readout,

    C(theta) = 1/(F_C(theta) + eta),
    F_C(theta) = q(theta)^2/[p(theta)(1-p(theta))],

where p is the even-parity probability at phi0 and q = dp/dphi at phi0.

The exact gradient uses the chain rule for p and q. The finite-shot estimator
uses the same Hopf tangent-state and signed-branch structure as the VQE script,
but replaces the Hamiltonian by the Ramsey chain-rule observable. The local
slope q is estimated with centered finite differences at phi0 +/- delta, and
stochastic inverse-Fisher prefactors are stabilized as described in the paper.

For n=6 and phi0=pi/(2n), the ideal noiseless fixed-readout GHZ-type reference
has F_C=36 and cost 1/(36+eta). This is a local reference value, not a noisy or
global phase-estimation guarantee.

Output:
    * MET_Layerwise_ADAM_EGTCG.pdf.

Dependencies: numpy, matplotlib, hopf_utils.
================================================================================
"""

from __future__ import annotations
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
import math
import numpy as np
import matplotlib.pyplot as plt
import hopf_utils as hopf
try:
    from threadpoolctl import threadpool_limits
except Exception:  # pragma: no cover
    threadpool_limits = None

# ============================================================
# 1. Baseline Utilities & Ramsey Metrology Problem
# ============================================================

COST_EPS = 1e-9
PROB_EPS = 1e-9
MIN_FISHER_FOR_MC_PREFACTOR = 1e-6
FINITE_DIFF_DELTA = 0.12
COEFF_CLIP = 5.0e3

# Each Ramsey moment in the finite-shot estimator is measured at three
# physical phase settings: phi0, phi0 + delta, and phi0 - delta.
RAMSEY_PHASE_SETTINGS_PER_MOMENT = 3

# The current-probe moment budget is used twice in the MC gradient estimator:
# once to estimate p and q for the chain-rule prefactors, and once to estimate
# the current-state weighted moment used in the signed-branch control variate.
CURRENT_PROBE_MOMENT_BLOCKS_PER_GRADIENT = 2


class Baseline:
    """Fixed-readout Ramsey baseline tools for the metrology numerical panel."""

    @staticmethod
    def hadamard_matrix(n):
        """Dense H^{otimes n} readout matrix used for exact Ramsey calculations."""
        H1 = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / math.sqrt(2.0)
        H = np.array([[1.0 + 0.0j]])
        for _ in range(n):
            H = np.kron(H, H1)
        return H

    @staticmethod
    def make_common_phase_problem(n=6):
        """
        Construct the fixed-readout Ramsey problem.

            U_phi = exp[-i phi G],   G = (1/2) sum_j Z_j,
            readout = H^{otimes n} followed by computational-basis measurement.

        Computational-basis convention:
            Z|0> = +|0>, Z|1> = -|1>.
        """
        dim = 1 << n
        G_diag = np.zeros(dim, dtype=float)
        parity_even = np.zeros(dim, dtype=bool)
        complement = np.zeros(dim, dtype=int)
        full_mask = dim - 1
        for basis in range(dim):
            ones = int(basis.bit_count())
            G_diag[basis] = 0.5 * (n - 2 * ones)
            parity_even[basis] = (ones % 2 == 0)
            complement[basis] = full_mask ^ basis

        phi0 = math.pi / (2.0 * n)
        H_readout = Baseline.hadamard_matrix(n)

        phase0 = np.exp(-1j * phi0 * G_diag)
        phase_plus = np.exp(-1j * (phi0 + FINITE_DIFF_DELTA) * G_diag)
        phase_minus = np.exp(-1j * (phi0 - FINITE_DIFF_DELTA) * G_diag)

        exact_fisher = float((G_diag.max() - G_diag.min()) ** 2)
        exact_cost = 1.0 / (exact_fisher + COST_EPS)
        return {
            "n": n,
            "dim": dim,
            "G_diag": G_diag,
            "parity_even": parity_even,
            "complement": complement,
            "phi0": phi0,
            "H_readout": H_readout,
            "phase0": phase0,
            "phase_plus": phase_plus,
            "phase_minus": phase_minus,
            "exact_fisher": exact_fisher,
            "exact_cost": exact_cost,
        }

    @staticmethod
    def ramsey_amplitudes_from_phase(psi, problem, phase):
        """State amplitudes after a cached U_phi and H^{otimes n} readout."""
        psi = np.asarray(psi, dtype=complex)
        return problem["H_readout"] @ (phase * psi)

    @staticmethod
    def ramsey_probabilities_from_phase(psi, problem, phase):
        """Exact computational-basis probabilities after the full Ramsey circuit."""
        amp = Baseline.ramsey_amplitudes_from_phase(psi, problem, phase)
        probs = np.real(np.conj(amp) * amp)
        probs = np.maximum(probs, 0.0)
        s = float(np.sum(probs))
        if s <= 0.0:
            raise ValueError("Invalid Ramsey probabilities: total probability is non-positive.")
        return probs / s

    @staticmethod
    def ramsey_even_probability_from_phase(psi, problem, phase):
        """
        Even-parity probability after the full Ramsey circuit.

        This evaluates the exact parity coarse graining of

            U_phi -> H^{otimes n} -> computational-basis measurement

        using the equivalent parity identity

            p_even(phi) = [1 + Re <psi|U_phi^\\dagger X^{\\otimes n} U_phi|psi>]/2.

        It is statistically identical to sampling full Ramsey bitstrings and
        then coarse graining by parity, but avoids repeated dense readout
        simulations inside the Monte Carlo gradient loops.
        """
        psi = np.asarray(psi, dtype=complex)
        comp = problem["complement"]
        # phase = exp(-i phi G), hence U_phi^\dagger X_all U_phi contributes
        # exp(+2 i phi G_b) on the complement-pair matrix element.
        parity_phase = np.conj(phase) ** 2
        corr = np.vdot(psi, parity_phase * psi[comp])
        p_even = 0.5 * (1.0 + float(np.real(corr)))
        return float(np.clip(p_even, 0.0, 1.0))

    @staticmethod
    def ramsey_even_probability_and_slope(psi, problem):
        """Exact p_even(phi0) and q_even(phi0)=d p_even / dphi."""
        psi = np.asarray(psi, dtype=complex)
        G = problem["G_diag"]
        phase = problem["phase0"]
        comp = problem["complement"]
        parity_phase = np.conj(phase) ** 2
        pair = np.conj(psi) * psi[comp] * parity_phase
        corr = np.sum(pair)
        dcorr = np.sum((1j * 2.0 * G) * pair)
        p = 0.5 * (1.0 + float(np.real(corr)))
        q = 0.5 * float(np.real(dcorr))
        return float(np.clip(p, 0.0, 1.0)), q

    @staticmethod
    def classical_fisher_from_pq(p, q):
        p_safe = float(np.clip(p, PROB_EPS, 1.0 - PROB_EPS))
        denom = max(p_safe * (1.0 - p_safe), PROB_EPS)
        return float((q * q) / denom)

    @staticmethod
    def fisher_from_psi(psi, problem):
        p, q = Baseline.ramsey_even_probability_and_slope(psi, problem)
        return Baseline.classical_fisher_from_pq(p, q)

    @staticmethod
    def cost_from_psi(psi, problem):
        fisher = Baseline.fisher_from_psi(psi, problem)
        return 1.0 / (fisher + COST_EPS)

    @staticmethod
    def fisher_cost_coefficients(p, q, clip=True, fisher_floor=MIN_FISHER_FOR_MC_PREFACTOR):
        """
        Chain-rule coefficients for the minimized cost

            C = 1/(F_C + eps),     F_C = q^2/[p(1-p)].

        Returns coeff_p, coeff_q, and the Fisher estimate used for the prefactor,
        so that

            dC = coeff_p * dp + coeff_q * dq.
        """
        p_safe = float(np.clip(p, PROB_EPS, 1.0 - PROB_EPS))
        denom = max(p_safe * (1.0 - p_safe), PROB_EPS)
        fisher = float((q * q) / denom)
        fisher_for_prefactor = max(fisher, fisher_floor)

        dF_dp = -q * q * (1.0 - 2.0 * p_safe) / (denom * denom)
        dF_dq = 2.0 * q / denom
        scale = -1.0 / ((fisher_for_prefactor + COST_EPS) ** 2)
        coeff_p = scale * dF_dp
        coeff_q = scale * dF_dq
        if clip:
            coeff_p = np.clip(coeff_p, -COEFF_CLIP, COEFF_CLIP)
            coeff_q = np.clip(coeff_q, -COEFF_CLIP, COEFF_CLIP)
        return float(coeff_p), float(coeff_q), fisher

    @staticmethod
    def state_gradient_cost(psi, problem):
        """Euclidean gradient dC/dpsi for a real state-vector line search."""
        psi = np.asarray(psi, dtype=float)
        G = problem["G_diag"]
        H = problem["H_readout"]
        phase = problem["phase0"]
        even = problem["parity_even"]

        amp = H @ (phase * psi)
        amp_phi = H @ ((-1j * G * phase) * psi)
        p = float(np.sum(np.abs(amp[even]) ** 2).real)
        q = float(2.0 * np.real(np.vdot(amp[even], amp_phi[even])))
        coeff_p, coeff_q, fisher = Baseline.fisher_cost_coefficients(p, q, clip=False, fisher_floor=0.0)

        mask = even.astype(float)
        # dp = 2 Re <B dpsi, Pi_even B psi>
        tmp_p = mask * amp
        back_p = np.conj(phase) * (H.conj().T @ tmp_p)
        grad_p = 2.0 * np.real(back_p)

        # dq = 2 Re[ <B dpsi, Pi_even B_phi psi> + <B psi, Pi_even B_phi dpsi> ]
        tmp_q1 = mask * amp_phi
        tmp_q2 = mask * amp
        back_q1 = np.conj(phase) * (H.conj().T @ tmp_q1)
        back_q2 = (1j * G * np.conj(phase)) * (H.conj().T @ tmp_q2)
        grad_q = 2.0 * np.real(back_q1 + back_q2)

        return coeff_p * grad_p + coeff_q * grad_q

    @staticmethod
    def expected_cost(theta, problem):
        psi = hopf.vector_from_theta(theta, "real")
        return Baseline.cost_from_psi(psi, problem)

    @staticmethod
    def expected_fisher(theta, problem):
        psi = hopf.vector_from_theta(theta, "real")
        return Baseline.fisher_from_psi(psi, problem)

    @staticmethod
    def exact_cost_gradient(theta, problem):
        psi = hopf.vector_from_theta(theta, "real")
        J = hopf.jacobian(theta, "real")
        grad_psi = Baseline.state_gradient_cost(psi, problem)
        return np.real(grad_psi @ J).flatten()

    @staticmethod
    def init_theta(n, seed=42):
        # Ramsey-standard product probe |+>^{otimes n}.  This avoids starting
        # from a random state with nearly zero readout slope, where 1/F_C is
        # numerically singular and not informative for this fixed-readout benchmark.
        dim = 1 << n
        psi0 = np.ones(dim, dtype=float) / math.sqrt(dim)
        theta0 = hopf.theta_from_vector(psi0, "real")
        theta0 = hopf.clip_theta_hopf_real(theta0, n=n)
        return theta0

    @staticmethod
    def adam_update(theta, grad, state, lr=0.01, b1=0.9, b2=0.999, eps=1e-8):
        m, v, t = state
        t += 1
        m = b1 * m + (1 - b1) * grad
        v = b2 * v + (1 - b2) * (grad**2)
        m_hat = m / (1 - b1**t)
        v_hat = v / (1 - b2**t)
        theta_new = theta - lr * m_hat / (np.sqrt(v_hat) + eps)
        theta_new = hopf.clip_theta_hopf_real(theta_new)
        return theta_new, (m, v, t)


# ============================================================
# 2. Geometry & EGT-CG Optimizer
# ============================================================

def expmap_sphere(x, v, t):
    vn = float(np.linalg.norm(v))
    if vn < 1e-15: return x.copy()
    a = t * vn
    return math.cos(a) * x + math.sin(a) * (v / vn)


def transport_sphere(x, y, w):
    c = float(np.dot(x, y))
    denom = 1.0 + c
    if denom < 1e-12:
        wp = w - float(np.dot(w, y)) * y
        return wp
    wp = w - (float(np.dot(w, y)) / denom) * (x + y)
    wp = wp - float(np.dot(wp, y)) * y
    return wp


class EGT_CG:
    def __init__(self, dim, t0=5.0, c1=1e-4, c2=0.99, max_ls=25, eps=1e-12):
        self.dim = dim
        self.t0, self.c1, self.c2, self.max_ls, self.eps = t0, c1, c2, max_ls, eps
        self.psi_prev = None
        self.v_prev = None
        self.u_prev = None

    def _build_v(self, theta, grad_theta):
        psi = hopf.vector_from_theta(theta, "real")
        v = np.zeros_like(psi)
        g_diag = hopf.metric_diagonal(list(theta), "real")

        for i in range(1, self.dim + 1):
            dpi = hopf.partial_derivative_column_real(theta, i)
            gii = g_diag[i-1]
            if gii < self.eps: continue
            nat_i = float(grad_theta[i - 1]) / gii
            v -= nat_i * dpi

        v -= float(np.dot(psi, v)) * psi
        return psi, v

    def _C_psi(self, psi, problem):
        return Baseline.cost_from_psi(psi, problem)

    def _dphi(self, psi, u, problem):
        grad_psi = Baseline.state_gradient_cost(psi, problem)
        return float(np.dot(grad_psi, u))

    def _line_search(self, psi0, u0, problem):
        u0n = float(np.linalg.norm(u0))
        if u0n < self.eps: return 0.0, psi0
        phi0 = self._C_psi(psi0, problem)
        dphi0 = self._dphi(psi0, u0, problem)
        if dphi0 >= 0.0: return 0.0, psi0

        t = float(self.t0)
        for _ in range(self.max_ls):
            psi_t = expmap_sphere(psi0, u0, t)
            phi_t = self._C_psi(psi_t, problem)
            if phi_t > phi0 + self.c1 * t * dphi0:
                t *= 0.5
                continue
            u_t = transport_sphere(psi0, psi_t, u0)
            dphi_t = self._dphi(psi_t, u_t, problem)
            if abs(dphi_t) <= self.c2 * abs(dphi0): return t, psi_t
            t *= 0.5
        return t, expmap_sphere(psi0, u0, t)

    def step(self, theta, grad_theta, problem):
        theta = np.asarray(theta, dtype=float)
        grad_theta = np.asarray(grad_theta, dtype=float)
        psi, v = self._build_v(theta, grad_theta)

        if self.psi_prev is None:
            u = v.copy()
            self.psi_prev = psi.copy()
            self.v_prev = v.copy()
            self.u_prev = u.copy()
        else:
            v_prev_t = transport_sphere(self.psi_prev, psi, self.v_prev)
            u_prev_t = transport_sphere(self.psi_prev, psi, self.u_prev)
            num = float(np.dot(v, (v - v_prev_t)))
            den = float(np.dot(v_prev_t, v_prev_t)) + 1e-15
            beta = max(0.0, num / den)
            u = v + beta * u_prev_t
            u -= float(np.dot(psi, u)) * psi
            self.psi_prev = psi.copy()
            self.v_prev = v.copy()
            self.u_prev = u.copy()

        # Numerical safeguard: restart if the transported CG direction is not descent.
        if self._dphi(psi, u, problem) >= 0.0:
            u = v.copy()
            self.u_prev = u.copy()

        _, psi_new = self._line_search(psi, u, problem)
        theta_new = hopf.theta_from_vector(psi_new, "real")
        theta_new = hopf.clip_theta_hopf_real(theta_new)
        return theta_new


# ============================================================
# 3. Estimators & Budget
# ============================================================

def get_layer_budgets(plan: str, n: int, total_grad_shots: int):
    """Allocate tangent-state and branch-state shots by Hopf tree layer."""
    budgets = []
    budget_per_type = total_grad_shots / 2.0
    if plan == "param":
        M = (1 << n) - 1
        shots_per_param = budget_per_type / float(M)
        for d in range(n):
            layer_size = 1 << d
            val = max(1, int(round(shots_per_param * layer_size)))
            budgets.append(val)
    else: raise ValueError(f"Unknown plan: {plan}")
    return budgets, budgets


def sample_even_probability_from_phase(psi, problem, phase, shots, rng):
    """
    Sample the coarse-grained even-parity outcome after the full Ramsey circuit.

    This is equivalent to sampling full Ramsey bitstrings and then checking
    parity, but it uses the binomial coarse-grained distribution for speed.
    """
    shots = int(shots)
    if shots <= 0:
        raise ValueError("shots must be positive")
    p_even = Baseline.ramsey_even_probability_from_phase(psi, problem, phase)
    p_even = float(np.clip(p_even, 0.0, 1.0))
    return float(rng.binomial(shots, p_even) / float(shots))


def estimate_ramsey_pq_mc(psi, problem, shots, rng):
    """
    Shot-estimated p_even and finite-difference q_even for the Ramsey readout.

    p is estimated at phi0.  q = dp/dphi is estimated by a central
    finite-difference pair of full Ramsey measurement circuits at
    phi0 +/- FINITE_DIFF_DELTA.
    """
    p0 = sample_even_probability_from_phase(psi, problem, problem["phase0"], shots, rng)
    pp = sample_even_probability_from_phase(psi, problem, problem["phase_plus"], shots, rng)
    pm = sample_even_probability_from_phase(psi, problem, problem["phase_minus"], shots, rng)
    q = (pp - pm) / (2.0 * FINITE_DIFF_DELTA)
    return p0, float(q)


def estimate_weighted_ramsey_moment_mc(psi, problem, coeff_p, coeff_q, shots, rng):
    """
    Estimate the weighted moment of the Ramsey chain-rule observable

        O = coeff_p M_even + coeff_q D_even,

    where M_even is the full Ramsey POVM coarse-grained to even parity and
    D_even = dM_even/dphi.  The M part is sampled at phi0; the D part is sampled
    by the central finite-difference Ramsey circuits at phi0 +/- delta.
    """
    if int(shots) <= 0:
        return 0.0
    p0 = sample_even_probability_from_phase(psi, problem, problem["phase0"], shots, rng)
    pp = sample_even_probability_from_phase(psi, problem, problem["phase_plus"], shots, rng)
    pm = sample_even_probability_from_phase(psi, problem, problem["phase_minus"], shots, rng)
    q = (pp - pm) / (2.0 * FINITE_DIFF_DELTA)
    return float(coeff_p * p0 + coeff_q * q)


def estimate_gradient_layerwise(theta, n, problem, Npsi, b_part, b_phi, rng):
    """
    Layerwise signed-branch estimator for the Ramsey inverse-CFI gradient.

    First estimate p and q for the current state.  These define the scalar
    chain-rule observable for C = 1/(F_C+eta).  Then estimate the Hopf transition
    elements of that scalar observable using the same signed-branch layerwise
    structure as the VQE script.
    """
    M = (1 << n) - 1
    psi = hopf.vector_from_theta(theta, "real")

    p_hat, q_hat = estimate_ramsey_pq_mc(psi, problem, Npsi, rng)
    coeff_p, coeff_q, _ = Baseline.fisher_cost_coefficients(p_hat, q_hat)
    Epsi = estimate_weighted_ramsey_moment_mc(psi, problem, coeff_p, coeff_q, Npsi, rng)

    g_diag = hopf.metric_diagonal(list(theta), "real")
    Epartial = np.zeros(M + 1)
    grad = np.zeros(M)
    layer_start = 1

    for d in range(n):
        layer_size = 1 << d
        nd_p, nd_v = b_part[d], b_phi[d]

        tangent_cache = []
        valid_cache = []
        for offset in range(layer_size):
            i = layer_start + offset
            gi = g_diag[i-1]
            if gi > 1e-12:
                dpsi = hopf.partial_derivative_column_real(theta, i)
                tangent_cache.append(dpsi / np.sqrt(gi))
                valid_cache.append(True)
            else:
                tangent_cache.append(None)
                valid_cache.append(False)

        # Phase 1: tangent-state weighted moments.
        offsets_p = rng.integers(0, layer_size, size=nd_p)
        counts_p = np.bincount(offsets_p, minlength=layer_size)
        for offset, cnt in enumerate(counts_p):
            if cnt <= 0 or not valid_cache[offset]:
                continue
            i = layer_start + offset
            Epartial[i] = estimate_weighted_ramsey_moment_mc(
                tangent_cache[offset], problem, coeff_p, coeff_q, int(cnt), rng
            )

        # Phase 2: signed-branch transition overlaps.
        offsets_v = rng.integers(0, layer_size, size=nd_v)
        signs_v = np.where(rng.random(nd_v) < 0.5, 1.0, -1.0)
        acc_grad_sum, acc_grad_cnt = np.zeros(layer_size), np.zeros(layer_size)
        for offset in range(layer_size):
            if not valid_cache[offset]:
                continue
            i = layer_start + offset
            mask_offset = offsets_v == offset
            if not np.any(mask_offset):
                continue
            gi = g_diag[i-1]
            control_variate = 0.5 * (Epsi + Epartial[i])
            for s in (1.0, -1.0):
                cnt = int(np.count_nonzero(mask_offset & (signs_v == s)))
                if cnt <= 0:
                    continue
                v_vec = (psi + s * tangent_cache[offset]) / math.sqrt(2.0)
                val = estimate_weighted_ramsey_moment_mc(v_vec, problem, coeff_p, coeff_q, cnt, rng)
                overlap = s * (val - control_variate)
                acc_grad_sum[offset] += cnt * (2.0 * np.sqrt(gi) * overlap)
                acc_grad_cnt[offset] += cnt

        for offset in range(layer_size):
            i = layer_start + offset
            if acc_grad_cnt[offset] > 0:
                grad[i-1] = acc_grad_sum[offset] / acc_grad_cnt[offset]
            else:
                grad[i-1] = 0.0
        layer_start += layer_size
    return grad


# ============================================================
# 4. Experiments
# ============================================================

def run_adam_experiment(plan, total_grad_shots, n, problem, theta_init, steps, lr, Npsi, seed):
    rng = np.random.default_rng(seed)
    theta = theta_init.copy()
    state = (np.zeros_like(theta), np.zeros_like(theta), 0)
    b_part, b_phi = get_layer_budgets(plan, n, total_grad_shots)
    costs = []
    for t in range(steps):
        costs.append(Baseline.expected_cost(theta, problem))
        grad = estimate_gradient_layerwise(theta, n, problem, Npsi, b_part, b_phi, rng)
        theta, state = Baseline.adam_update(theta, grad, state, lr)
    return np.array(costs)


def run_egtcg_experiment(plan, total_grad_shots, n, problem, theta_init, steps, t0, Npsi, seed):
    rng = np.random.default_rng(seed)
    theta = theta_init.copy()
    opt = EGT_CG(dim=theta.size, t0=t0)
    b_part, b_phi = get_layer_budgets(plan, n, total_grad_shots)
    costs = []
    for t in range(steps):
        costs.append(Baseline.expected_cost(theta, problem))
        grad = estimate_gradient_layerwise(theta, n, problem, Npsi, b_part, b_phi, rng)
        theta = opt.step(theta, grad, problem)
    return np.array(costs)


def run_exact_baseline_adam(theta_init, problem, steps, lr):
    theta = theta_init.copy()
    state = (np.zeros_like(theta), np.zeros_like(theta), 0)
    costs = []
    for t in range(steps):
        costs.append(Baseline.expected_cost(theta, problem))
        grad = Baseline.exact_cost_gradient(theta, problem)
        theta, state = Baseline.adam_update(theta, grad, state, lr)
    return np.array(costs)


def run_exact_baseline_egtcg(theta_init, problem, steps, t0):
    theta = theta_init.copy()
    opt = EGT_CG(dim=theta.size, t0=t0)
    costs = []
    for t in range(steps):
        costs.append(Baseline.expected_cost(theta, problem))
        grad = Baseline.exact_cost_gradient(theta, problem)
        theta = opt.step(theta, grad, problem)
    return np.array(costs)


# ============================================================
# 5. Plotting
# ============================================================

def plot_results_comparison(data, steps, total_grad_shots, n_psi_shots, outfile="MET_Layerwise_ADAM_EGTCG.pdf"):
    """Plot using physical Ramsey hardware shot counts in the title."""
    fig, ax = plt.subplots(figsize=(5.5, 6))
    x = np.arange(steps)
    ax.axhline(data["exact_cost"], color='gray', linestyle='--')
    ax.plot(x, data["adam_exact"], '-',  color="tab:blue", label='Adam (exact)', linewidth=1.5)
    ax.plot(x, data["egt_exact"], '-',  color="tab:green", label='EGT-CG (exact)', linewidth=2.0)
    ax.plot(x, data["adam_param"], '--', color="tab:blue", label='Adam (finite-shot)', linewidth=2.0)
    ax.plot(x, data["egt_param"], '--', color="tab:green", label='EGT-CG (finite-shot)', linewidth=2.0)
    ax.set_xlabel("Step")
    ax.set_ylabel(r"Ramsey inverse-CFI cost $1/(F_C+\eta)$")
    ax.grid(True, alpha=0.25)
    fig.suptitle(f"Common-phase Ramsey metrology (n=6)\nGradient Shots: {total_grad_shots}, Moment Shots: {n_psi_shots}", fontsize=12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    if outfile: plt.savefig(outfile, dpi=300)


# ============================================================
# Main
# ============================================================

def main():
    if threadpool_limits is not None:
        threadpool_limits(1)

    n = 6
    steps = 200
    lr = 0.02

    # EGT-CG line-search scale used for the inverse-CFI objective.
    t0 = 5.0

    # Internal budgets are per Ramsey phase setting.  The plotted shot counts
    # are converted to total physical Ramsey circuit executions by accounting
    # for the three phase settings and the current-probe moment blocks.
    N_MOMENT_PER_PHASE = 150 * 3
    TOTAL_GRAD_PER_PHASE = 3150 * 3

    physical_grad_shots = RAMSEY_PHASE_SETTINGS_PER_MOMENT * TOTAL_GRAD_PER_PHASE
    physical_moment_shots = (
        CURRENT_PROBE_MOMENT_BLOCKS_PER_GRADIENT
        * RAMSEY_PHASE_SETTINGS_PER_MOMENT
        * N_MOMENT_PER_PHASE
    )
    physical_total_shots_per_mc_step = physical_grad_shots + physical_moment_shots

    problem = Baseline.make_common_phase_problem(n)
    theta_init = Baseline.init_theta(n, seed=123)

    data = {"exact_fisher": problem["exact_fisher"], "exact_cost": problem["exact_cost"]}
    data["adam_exact"] = run_exact_baseline_adam(theta_init, problem, steps, lr)
    data["egt_exact"]  = run_exact_baseline_egtcg(theta_init, problem, steps, t0)
    data["adam_param"] = run_adam_experiment("param", TOTAL_GRAD_PER_PHASE, n, problem, theta_init, steps, lr, N_MOMENT_PER_PHASE, seed=202)
    data["egt_param"] = run_egtcg_experiment("param", TOTAL_GRAD_PER_PHASE, n, problem, theta_init, steps, t0, N_MOMENT_PER_PHASE, seed=302)

    plot_results_comparison(data, steps, physical_grad_shots, physical_moment_shots)

    print("6-qubit common-phase Ramsey metrology with full measurement circuit")
    print(f"Operating point phi0 = pi/(2n) = {problem['phi0']:.9f}")
    print(f"Theoretical optimum: F_C = {problem['exact_fisher']:.6f}, cost = 1/F_C = {problem['exact_cost']:.9f}")
    print(f"Per-phase gradient budget used internally: {TOTAL_GRAD_PER_PHASE}")
    print(f"Per-phase current-probe moment budget used internally: {N_MOMENT_PER_PHASE}")
    print(f"Physical gradient shots shown in plot: {physical_grad_shots}")
    print(f"Physical moment shots shown in plot: {physical_moment_shots}")
    print(f"Total physical Ramsey shots per MC gradient step: {physical_total_shots_per_mc_step}")
    print(f"Final Adam Exact cost: {data['adam_exact'][-1]:.9f}")
    print(f"Final EGT-CG Exact cost: {data['egt_exact'][-1]:.9f}")
    print(f"Final Adam MC cost:    {data['adam_param'][-1]:.9f}")
    print(f"Final EGT-CG MC cost:  {data['egt_param'][-1]:.9f}")


if __name__ == "__main__":
    main()
