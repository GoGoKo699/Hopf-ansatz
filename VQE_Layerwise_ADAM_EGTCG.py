"""VQE_Layerwise_ADAM_EGTCG.py
================================================================================
Numerical VQE panel for the real Hopf ansatz.

The problem is the six-qubit open-boundary transverse-field Ising Hamiltonian
used in the paper,

    H_TFIM = -J sum_j Z_j Z_{j+1} - h sum_j X_j,
    J = 1, h = 0.6.

The script compares Adam and EGT-CG under two gradient-access models:

    * exact state-vector gradients, and
    * finite-shot layerwise Hopf-gradient estimates.

The finite-shot estimator implements the paper's signed-branch transition-moment
identity. For each magnitude layer of the Hopf tree, it samples normalized
coordinate tangent states and signed branch states, then forms the gradient
component using the analytic diagonal metric. Shots are allocated across layers
proportionally to layer size so that each coordinate receives approximately the
same expected number of samples.

Shot accounting follows the numerical section of the paper. Internal budgets
are per Pauli term; the plotted physical shot counts multiply those budgets by
the number of TFIM Pauli measurement streams.

Output:
    * VQE_Layerwise_ADAM_EGTCG.pdf.

Dependencies: numpy, matplotlib, hopf_utils.
================================================================================
"""

from __future__ import annotations
import os

# Keep the small dense linear algebra single-threaded; for n=6 the BLAS
# thread-launch overhead is larger than the matrix work itself.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import math
import numpy as np
try:
    from threadpoolctl import threadpool_limits
    _THREADPOOL_LIMIT = threadpool_limits(1)
    _THREADPOOL_LIMIT.__enter__()
except Exception:
    _THREADPOOL_LIMIT = None
import matplotlib.pyplot as plt
import hopf_utils as hopf

# ============================================================
# 1. Baseline Utilities & Hamiltonian
# ============================================================

class Baseline:
    """VQE baseline tools for the TFIM numerical panel."""
    
    @staticmethod
    def pauli_matrix(name):
        if name == 'I': return np.eye(2, dtype=complex)
        if name == 'X': return np.array([[0, 1], [1, 0]], dtype=complex)
        if name == 'Y': return np.array([[0, -1j], [1j, 0]], dtype=complex)
        if name == 'Z': return np.array([[1, 0], [0, -1]], dtype=complex)
        raise ValueError(f"Unknown Pauli {name}")

    @staticmethod
    def kron_n(ops):
        res = ops[0]
        for op in ops[1:]:
            res = np.kron(res, op)
        return res

    @staticmethod
    def make_tfim_problem(n=6, J=1.0, h=0.6):
        """Construct the open-boundary one-dimensional TFIM Hamiltonian.

        ``terms`` keeps the dense Pauli matrix for compatibility, but also
        stores a lightweight descriptor used by the MC estimator:

            (coefficient, dense_matrix, kind, data)

        where kind is ``"ZZ"`` with a precomputed sign vector, or ``"X"``
        with a precomputed basis permutation.
        """
        dim = 1 << n
        H = np.zeros((dim, dim), dtype=complex)
        terms = []
        basis = np.arange(dim)

        for i in range(n - 1):
            ops = ['I'] * n
            ops[i], ops[i+1] = 'Z', 'Z'
            mat = Baseline.kron_n([Baseline.pauli_matrix(o) for o in ops])
            H -= J * mat

            mask_i = 1 << (n - 1 - i)
            mask_j = 1 << (n - 1 - (i + 1))
            z_i = np.where((basis & mask_i) == 0, 1.0, -1.0)
            z_j = np.where((basis & mask_j) == 0, 1.0, -1.0)
            zz_signs = z_i * z_j
            terms.append((-J, mat, "ZZ", zz_signs))

        for i in range(n):
            ops = ['I'] * n
            ops[i] = 'X'
            mat = Baseline.kron_n([Baseline.pauli_matrix(o) for o in ops])
            H -= h * mat

            mask_i = 1 << (n - 1 - i)
            x_perm = basis ^ mask_i
            terms.append((-h, mat, "X", x_perm))

        eigvals = np.linalg.eigvalsh(H)
        exact_E = eigvals[0]
        return H, exact_E, terms

    @staticmethod
    def expected_value(theta, H):
        psi = hopf.vector_from_theta(theta, "real")
        return float(np.real(psi.conj().T @ H @ psi))

    @staticmethod
    def init_theta(n, seed=42):
        psi0 = hopf.canonical_initial_state(n, seed)
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
    def __init__(self, dim, t0=1.0, c1=1e-4, c2=0.9, max_ls=25, eps=1e-12):
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

    def _E_psi(self, psi, H): return float(np.real(psi @ (H @ psi)))
    def _dphi(self, psi, u, H): return float(np.real(np.dot(2.0 * (H @ psi), u)))

    def _line_search(self, psi0, u0, H):
        u0n = float(np.linalg.norm(u0))
        if u0n < self.eps: return 0.0, psi0
        phi0 = self._E_psi(psi0, H)
        dphi0 = self._dphi(psi0, u0, H)
        if dphi0 >= 0.0: return 0.0, psi0

        t = float(self.t0)
        for _ in range(self.max_ls):
            psi_t = expmap_sphere(psi0, u0, t)
            phi_t = self._E_psi(psi_t, H)
            if phi_t > phi0 + self.c1 * t * dphi0:
                t *= 0.5
                continue
            u_t = transport_sphere(psi0, psi_t, u0)
            dphi_t = self._dphi(psi_t, u_t, H)
            if abs(dphi_t) <= self.c2 * abs(dphi0): return t, psi_t
            t *= 0.5
        return t, expmap_sphere(psi0, u0, t)

    def step(self, theta, grad_theta, H):
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

        _, psi_new = self._line_search(psi, u, H)
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

def pauli_expectation_batch_fast(term, states):
    """Fast vectorized <P> for one TFIM Pauli descriptor.

    ``states`` has shape (B, dim). Terms are stored as
    ``(coefficient, dense_matrix, kind, data)`` where ``kind`` is:
      - ``"ZZ"``: ``data`` is a computational-basis sign vector;
      - ``"X"``: ``data`` is the bit-flip permutation for that qubit.
    """
    c, P, kind, data = term
    states = np.asarray(states, dtype=complex)
    if states.ndim == 1:
        states = states[None, :]

    if kind == "ZZ":
        probs = np.abs(states) ** 2
        return np.real(probs @ data)
    if kind == "X":
        return np.real(np.sum(np.conj(states) * states[:, data], axis=1))

    # Fallback for unexpected term formats.
    return np.real(np.einsum("bi,ij,bj->b", states.conj(), P, states, optimize=True))


_TERM_PACK_CACHE = {}


def _pack_terms(terms):
    """Return vectorized TFIM term data, cached by object identity."""
    key = id(terms)
    cached = _TERM_PACK_CACHE.get(key)
    if cached is not None:
        return cached

    zz_diags, zz_coeffs = [], []
    x_perms, x_coeffs = [], []
    for term in terms:
        # Current make_tfim_problem format: (coefficient, dense_matrix, kind, data).
        if len(term) == 4:
            c, _P, kind, data = term
        # Fallback lightweight format: (coefficient, kind, data).
        elif len(term) == 3:
            c, kind, data = term
        else:
            raise ValueError("Expected TFIM term descriptor with 3 or 4 entries")

        if kind in ("ZZ", "diag"):
            zz_coeffs.append(float(c))
            zz_diags.append(np.asarray(data, dtype=float))
        elif kind in ("X", "perm"):
            x_coeffs.append(float(c))
            x_perms.append(np.asarray(data, dtype=np.int64))
        else:
            raise ValueError(f"Unknown term kind: {kind}")

    dim = len(zz_diags[0]) if zz_diags else len(x_perms[0])
    pack = {
        "zz_diags": np.asarray(zz_diags, dtype=float).reshape(len(zz_diags), dim),
        "zz_coeffs": np.asarray(zz_coeffs, dtype=float),
        "x_perms": np.asarray(x_perms, dtype=np.int64).reshape(len(x_perms), dim),
        "x_coeffs": np.asarray(x_coeffs, dtype=float),
        "coeffs": np.asarray(zz_coeffs + x_coeffs, dtype=float),
    }
    _TERM_PACK_CACHE[key] = pack
    return pack


def pauli_expectation_fast(term, psi):
    """Fast expectation value for one TFIM Pauli descriptor."""
    if len(term) == 4:
        _c, _P, kind, data = term
    elif len(term) == 3:
        _c, kind, data = term
    else:
        raise ValueError("Expected TFIM term descriptor with 3 or 4 entries")

    if kind in ("ZZ", "diag"):
        return float(np.dot(np.abs(psi) ** 2, data))
    if kind in ("X", "perm"):
        return float(np.real(np.vdot(psi, psi[data])))
    raise ValueError(f"Unknown term kind: {kind}")


def estimate_energy_batch_mc(terms, states, shots, rng):
    """Vectorized block/binomial Hamiltonian estimator.

    The estimator uses the same Pauli-term binomial statistics as direct
    term-wise sampling, but batches identical Hopf tangent/branch configurations
    inside each layer.
    """
    states = np.asarray(states, dtype=complex)
    single = states.ndim == 1
    if single:
        states = states[None, :]

    B = states.shape[0]
    shots_arr = np.asarray(shots, dtype=int)
    if shots_arr.ndim == 0:
        shots_arr = np.full(B, int(shots_arr), dtype=int)
    else:
        shots_arr = shots_arr.reshape(-1).astype(int)
        if shots_arr.size != B:
            raise ValueError("shots must be a scalar or have one entry per state")

    out = np.zeros(B, dtype=float)
    active = shots_arr > 0
    if not np.any(active):
        return float(out[0]) if single else out

    pack = _pack_terms(terms)
    S = states[active]
    active_shots = shots_arr[active]

    pieces = []
    if pack["zz_diags"].size:
        probs = np.abs(S) ** 2
        pieces.append(probs @ pack["zz_diags"].T)
    if pack["x_perms"].size:
        x_states = S[:, pack["x_perms"]]
        pieces.append(np.real(np.sum(np.conj(S)[:, None, :] * x_states, axis=2)))

    if pieces:
        mu = np.concatenate(pieces, axis=1)
    else:
        mu = np.zeros((S.shape[0], 0), dtype=float)

    p_plus = np.clip(0.5 * (1.0 + mu), 0.0, 1.0)
    n_plus = rng.binomial(active_shots[:, None], p_plus)
    means = 2.0 * n_plus / active_shots[:, None] - 1.0
    out[active] = means @ pack["coeffs"]
    return float(out[0]) if single else out


def estimate_energy_block_mc(terms, psi, shots, rng):
    """Scalar compatibility wrapper around the vectorized estimator."""
    return estimate_energy_batch_mc(terms, psi, shots, rng)


# Alias used by the estimator entry points.
def estimate_energy_mc(terms, psi, shots, rng):
    return estimate_energy_block_mc(terms, psi, shots, rng)


def estimate_gradient_layerwise(theta, n, terms, Npsi, b_part, b_phi, rng):
    """Finite-shot layerwise Hopf-gradient estimator using grouped label sampling."""
    M = (1 << n) - 1
    psi = hopf.vector_from_theta(theta, "real")
    Epsi = estimate_energy_batch_mc(terms, psi, Npsi, rng)
    g_diag = hopf.metric_diagonal(list(theta), "real")
    Epartial = np.zeros(M + 1, dtype=float)
    grad = np.zeros(M, dtype=float)
    layer_start = 1
    inv_sqrt2 = 1.0 / math.sqrt(2.0)

    for d in range(n):
        layer_size = 1 << d
        nd_p, nd_v = int(b_part[d]), int(b_phi[d])

        # Cache normalized tangent states for this layer exactly once.
        tangent_cache = [None] * layer_size
        active = np.zeros(layer_size, dtype=bool)
        for offset in range(layer_size):
            i = layer_start + offset
            gi = float(g_diag[i - 1])
            if gi > 1e-12:
                dpsi = hopf.partial_derivative_column_real(theta, i)
                tangent_cache[offset] = dpsi / math.sqrt(gi)
                active[offset] = True

        # Phase 1: E_partial[i] = <u_i|H|u_i>.
        if nd_p > 0:
            offsets_p = rng.integers(0, layer_size, size=nd_p)
            counts_p = np.bincount(offsets_p, minlength=layer_size).astype(int)
        else:
            counts_p = np.zeros(layer_size, dtype=int)

        idx_p = [offset for offset in range(layer_size) if counts_p[offset] > 0 and active[offset]]
        if idx_p:
            states_p = np.asarray([tangent_cache[offset] for offset in idx_p])
            shots_p = counts_p[idx_p]
            vals_p = np.asarray(estimate_energy_batch_mc(terms, states_p, shots_p, rng), dtype=float)
            for offset, val in zip(idx_p, vals_p):
                Epartial[layer_start + offset] = val

        # Phase 2: transition overlaps from signed branch states.
        if nd_v > 0:
            offsets_v = rng.integers(0, layer_size, size=nd_v)
            signs_bool = rng.random(nd_v) < 0.5
            counts_plus = np.bincount(offsets_v[signs_bool], minlength=layer_size).astype(int)
            counts_minus = np.bincount(offsets_v[~signs_bool], minlength=layer_size).astype(int)
        else:
            counts_plus = np.zeros(layer_size, dtype=int)
            counts_minus = np.zeros(layer_size, dtype=int)

        acc_grad_sum = np.zeros(layer_size, dtype=float)
        acc_grad_cnt = counts_plus.astype(float) + counts_minus.astype(float)

        idx_plus = [offset for offset in range(layer_size) if counts_plus[offset] > 0 and active[offset]]
        if idx_plus:
            states_plus = np.asarray([(psi + tangent_cache[offset]) * inv_sqrt2 for offset in idx_plus])
            shots_plus = counts_plus[idx_plus]
            vals_plus = np.asarray(estimate_energy_batch_mc(terms, states_plus, shots_plus, rng), dtype=float)
            for offset, val in zip(idx_plus, vals_plus):
                cv = 0.5 * (Epsi + Epartial[layer_start + offset])
                acc_grad_sum[offset] += counts_plus[offset] * (val - cv)

        idx_minus = [offset for offset in range(layer_size) if counts_minus[offset] > 0 and active[offset]]
        if idx_minus:
            states_minus = np.asarray([(psi - tangent_cache[offset]) * inv_sqrt2 for offset in idx_minus])
            shots_minus = counts_minus[idx_minus]
            vals_minus = np.asarray(estimate_energy_batch_mc(terms, states_minus, shots_minus, rng), dtype=float)
            for offset, val in zip(idx_minus, vals_minus):
                cv = 0.5 * (Epsi + Epartial[layer_start + offset])
                acc_grad_sum[offset] -= counts_minus[offset] * (val - cv)

        for offset in range(layer_size):
            i = layer_start + offset
            gi = float(g_diag[i - 1])
            if gi > 1e-12 and acc_grad_cnt[offset] > 0 and active[offset]:
                mean_overlap = acc_grad_sum[offset] / acc_grad_cnt[offset]
                grad[i - 1] = 2.0 * math.sqrt(gi) * mean_overlap
            else:
                grad[i - 1] = 0.0

        layer_start += layer_size

    return grad

# ============================================================
# Main
# ============================================================

def run_adam_experiment(plan, total_grad_shots, n, H, terms, theta_init, steps, lr, Npsi, seed):
    rng = np.random.default_rng(seed)
    theta = theta_init.copy()
    state = (np.zeros_like(theta), np.zeros_like(theta), 0)
    b_part, b_phi = get_layer_budgets(plan, n, total_grad_shots)
    energies = []
    for t in range(steps):
        energies.append(Baseline.expected_value(theta, H))
        grad = estimate_gradient_layerwise(theta, n, terms, Npsi, b_part, b_phi, rng)
        theta, state = Baseline.adam_update(theta, grad, state, lr)
    return np.array(energies)

def run_egtcg_experiment(plan, total_grad_shots, n, H, terms, theta_init, steps, t0, Npsi, seed):
    rng = np.random.default_rng(seed)
    theta = theta_init.copy()
    opt = EGT_CG(dim=theta.size, t0=t0)
    b_part, b_phi = get_layer_budgets(plan, n, total_grad_shots)
    energies = []
    for t in range(steps):
        energies.append(Baseline.expected_value(theta, H))
        grad = estimate_gradient_layerwise(theta, n, terms, Npsi, b_part, b_phi, rng)
        theta = opt.step(theta, grad, H)
    return np.array(energies)

def run_exact_baseline_adam(theta_init, H, steps, lr):
    theta = theta_init.copy()
    state = (np.zeros_like(theta), np.zeros_like(theta), 0)
    energies = []
    for t in range(steps):
        energies.append(Baseline.expected_value(theta, H))
        psi = hopf.vector_from_theta(theta, "real")
        J = hopf.jacobian(theta, "real")
        grad = 2.0 * np.real(psi.T @ H @ J).flatten()
        theta, state = Baseline.adam_update(theta, grad, state, lr)
    return np.array(energies)

def run_exact_baseline_egtcg(theta_init, H, steps, t0):
    theta = theta_init.copy()
    opt = EGT_CG(dim=theta.size, t0=t0)
    energies = []
    for t in range(steps):
        energies.append(Baseline.expected_value(theta, H))
        psi = hopf.vector_from_theta(theta, "real")
        J = hopf.jacobian(theta, "real")
        grad = 2.0 * np.real(psi.T @ H @ J).flatten()
        theta = opt.step(theta, grad, H)
    return np.array(energies)

def plot_results_comparison(data, steps, total_grad_shots, n_psi_shots, outfile="VQE_Layerwise_ADAM_EGTCG.pdf"):
    fig, ax = plt.subplots(figsize=(5.5, 6))
    x = np.arange(steps)
    ax.axhline(data["exact_E"], color='gray', linestyle='--')
    ax.plot(x, data["adam_exact"], '-',  color="tab:blue", label='Adam (exact)', linewidth=2.0)
    ax.plot(x, data["egt_exact"], '-',  color="tab:green", label='EGT-CG (exact)', linewidth=2.0)
    ax.plot(x, data["adam_param"], '--', color="tab:blue", label='Adam (finite-shot)', linewidth=2.0)
    ax.plot(x, data["egt_param"], '--', color="tab:green", label='EGT-CG (finite-shot)', linewidth=2.0)
    ax.set_xlabel("Step")
    ax.set_ylabel("Energy")
    ax.grid(True, alpha=0.25)
    fig.suptitle(f"VQE for TFIM (n=6)\nGradient Shots: {total_grad_shots}, Energy Shots: {n_psi_shots}", fontsize=12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    if outfile: plt.savefig(outfile, dpi=300)

def main():
    n = 6
    steps = 200
    lr = 0.02
    t0 = 1.0

    # Internal estimator budgets are shots per Pauli term.  The TFIM
    # Hamiltonian has one measurement stream for each Pauli term in ``terms``,
    # so the physical circuit-execution count is larger by ``len(terms)``.
    # The plotted and printed shot counts use this physical accounting.
    N_psi_per_term = 150 * 2
    TOTAL_GRAD_per_term = 3150 * 2

    H, exact_E, terms = Baseline.make_tfim_problem(n)
    n_pauli_terms = len(terms)
    N_psi_physical = n_pauli_terms * N_psi_per_term
    TOTAL_GRAD_physical = n_pauli_terms * TOTAL_GRAD_per_term

    theta_init = Baseline.init_theta(n, seed=123)
    data = {"exact_E": exact_E}
    data["adam_exact"] = run_exact_baseline_adam(theta_init, H, steps, lr)
    data["egt_exact"]  = run_exact_baseline_egtcg(theta_init, H, steps, t0)
    data["adam_param"] = run_adam_experiment(
        "param", TOTAL_GRAD_per_term, n, H, terms, theta_init, steps, lr,
        N_psi_per_term, seed=202
    )
    data["egt_param"] = run_egtcg_experiment(
        "param", TOTAL_GRAD_per_term, n, H, terms, theta_init, steps, t0,
        N_psi_per_term, seed=302
    )

    print("6-qubit TFIM VQE with layerwise Hopf-gradient estimation")
    print(f"Number of Pauli measurement terms: {n_pauli_terms}")
    print(f"Per-term gradient budget used internally: {TOTAL_GRAD_per_term}")
    print(f"Per-term current-state energy budget used internally: {N_psi_per_term}")
    print(f"Physical gradient shots shown in plot: {TOTAL_GRAD_physical}")
    print(f"Physical energy shots shown in plot: {N_psi_physical}")
    print(
        "Total physical Pauli-measurement shots per MC gradient step: "
        f"{TOTAL_GRAD_physical + N_psi_physical}"
    )
    print(f"Exact ground-state energy: {exact_E:.9f}")
    print(f"Final Adam Exact energy: {data['adam_exact'][-1]:.9f}")
    print(f"Final EGT-CG Exact energy: {data['egt_exact'][-1]:.9f}")
    print(f"Final Adam MC energy:    {data['adam_param'][-1]:.9f}")
    print(f"Final EGT-CG MC energy:  {data['egt_param'][-1]:.9f}")

    plot_results_comparison(data, steps, TOTAL_GRAD_physical, N_psi_physical)

if __name__ == "__main__":
    main()
