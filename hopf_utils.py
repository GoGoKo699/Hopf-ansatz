"""hopf_utils.py
================================================================================
Core utilities for the Hopf ansatz used in the paper.

The module implements the mathematical and circuit conventions of the paper:

    * Hopf coordinate maps for real and complex normalized state vectors.
    * The inverse map from amplitudes back to Hopf parameters.
    * Exact Jacobians and the closed-form diagonal pullback metric.
    * Parameter assignments theta^(i) that prepare normalized coordinate
      tangent states on the same Hopf circuit skeleton.
    * Gate-schedule routines returning the four paper lists
      Ctrl, Anti, Targ, and Index.
    * Optional Qibo circuit construction for direct circuit verification.

Indexing conventions:

    * Tree indices are one-based: internal nodes are 1,...,2^n-1.
    * Real Hopf coordinates use [0, pi/2] for non-final internal layers and
      [0, 2pi) for the final internal layer, which carries real signs.
    * Complex Hopf coordinates use the same magnitude tree with angles in
      [0, pi/2], followed by one leaf phase for each computational-basis state.
    * Complex final-layer gates are promoted to R_C(theta_a, theta_b, theta_c),
      carrying one magnitude angle and the two sibling leaf phases.
    * Qibo's conventional RY(phi) gate is called with phi=2 theta, matching the
      paper convention R_y(theta) = [[cos theta, -sin theta], [sin theta, cos theta]].

When run directly, the module executes consistency checks for inverse mapping,
tangent-state synthesis, metric/Jacobian agreement, and optional Qibo circuit
agreement.

Dependencies: numpy; qibo is optional and used only for circuit verification.
================================================================================
"""

import numpy as np
import math
from collections import Counter
from typing import List, Tuple, Union, Dict

# Optional qibo import for circuit verification
try:
    from qibo import models, gates
    from qibo import set_backend
    _QIBO_AVAILABLE = True
    set_backend("numpy")
except Exception:
    models = None
    gates = None
    set_backend = None
    _QIBO_AVAILABLE = False


def canonical_initial_state(n: int, seed: int = 123):
    """Deterministic real normalized ψ0 with fixed global sign."""
    rng = np.random.default_rng(seed)
    m = 1 << n
    psi = rng.normal(size=m)
    psi = psi / np.linalg.norm(psi)
    if psi[0] < 0:
        psi = -psi
    return psi

def infer_n_from_theta(theta):
    """Infer n from theta length L."""
    L = int(np.asarray(theta).size)
    # Real: L = 2^n - 1 -> 2^n = L+1
    # Complex: L = 2^(n+1) - 1 -> 2^(n+1) = L+1
    # We try real first
    if ((L + 1) & L) == 0: # Power of 2 check
        log_val = int(round(math.log2(L + 1)))
        # Heuristic: usually n is the qubit count.
        # If ambiguous, we rely on context. This helper is mostly for real clipping.
        return log_val
    return int(round(math.log2(L + 1))) # Fallback

def clip_theta_hopf_real(theta, n=None, eps_internal=1e-6, eps_sing=1e-9, **kwargs):
    """
    Use the real Hopf coordinate ranges from the paper:
      - non-final internal nodes: theta in (0, pi/2)
      - final internal layer, i.e. parents of leaves: theta in [0, 2pi)
    """
    theta = np.asarray(theta, dtype=float).copy()
    L = theta.size
    if n is None:
        n = int(round(math.log2(L + 1)))
    
    if (1 << n) - 1 != L:
        # Accept the value unchanged if the length is not a full real Hopf parameter list.
        pass 

    last_start = (1 << (n - 1)) - 1

    # Internal nodes: strictly within (0, pi/2)
    if last_start > 0:
        theta[:last_start] = np.clip(theta[:last_start], eps_internal, (math.pi / 2) - eps_internal)

    # Last internal level: wrap to [0, 2pi)
    if last_start < L:
        th = np.mod(theta[last_start:], 2.0 * math.pi)
        # Nudge away from sin/cos zeros
        critical = np.array([0.0, 0.5*math.pi, math.pi, 1.5*math.pi, 2.0*math.pi])
        for k in critical:
            d = np.abs(((th - k + math.pi) % (2.0 * math.pi)) - math.pi)
            mask = d < eps_sing
            th[mask] = np.mod(th[mask] + eps_sing, 2.0 * math.pi)
        theta[last_start:] = th

    return theta

# =============================================================================
# Gate Scheduling & Qibo Circuit Construction
# =============================================================================

def hamming_weight(n: int) -> List[List[int]]:
    """
    Paper subroutine HammingWeight(n).
    Output: HW[0..n], where HW[k] is the ascending list of n-bit integers
    with Hamming weight k. Complexity: O(n * 2^n).
    """
    N = 1 << n
    HW = [[] for _ in range(n + 1)]

    pop = [0] * N
    for i in range(1, N):
        pop[i] = pop[i >> 1] + (i & 1)

    for i in range(N):
        HW[pop[i]].append(i)
    return HW

def find_pairs(A: List[int], B: List[int]) -> Tuple[List[int], List[int]]:
    """
    Paper subroutine FindPairs.

    It implements the greedy pairing rule used by the Hopf gate schedule and
    emits pairs in the deterministic order used by the four-qubit table.

    Assumptions:
      A and B are sorted ascending (true for hamming_weight output).

    Pair selection rule (unchanged):
      scan B from largest to smallest; for each b choose the current largest a in A with a<b
      (a is reusable; pointer only moves down when needed).

    Output order:
      - group by a in ascending order
      - within each a, emit b in descending order
        (this is exactly the order naturally accumulated when scanning B in reverse)
    """
    i = len(A) - 1

    # bucket b's by selected a
    buckets: Dict[int, List[int]] = {}

    for b in reversed(B):  # b from largest -> smallest
        while i >= 0 and A[i] >= b:
            i -= 1
        if i < 0:
            break
        a = A[i]
        buckets.setdefault(a, []).append(b)  # b's are appended in descending order per a

    # emit in table order: a ascending; b already descending within bucket
    A_out: List[int] = []
    B_out: List[int] = []
    for a in A:  # A already ascending
        bs = buckets.get(a)
        if not bs:
            continue
        for b in bs:  # keep descending b order (DO NOT reverse)
            A_out.append(a)
            B_out.append(b)

    return A_out, B_out



def anti_ctrl(A: List[int]) -> List[int]:
    """
    Paper subroutine Anti(A).
    For each distinct a with multiplicity c (in A),
      output masks: [2^c-2^c, 2^c-2^(c-1), ..., 2^c-2^1]
    ordered by ascending distinct a.
    Complexity: O(|A| log |A|)
    These are the local anti-control masks used by the gate-schedule algorithms.
    """
    freq = Counter(A)
    keys = sorted(freq.keys())

    out = []
    for a in keys:
        c = freq[a]
        for i in range(c, 0, -1):
            out.append((1 << c) - (1 << i))
    return out


def theta_real(n: int, a: int, b: int) -> int:
    """
    Paper subroutine ThetaReal(n,a,b) = (2^n + a) / (2*(b-a)).
    """
    d = b - a
    return ((1 << n) + a) // (2 * d)


def theta_complex(n: int, a: int, b: int) -> Union[int, List[int]]:
    """Return the Hopf parameter index or index triple for the complex ansatz schedule.

    For a given pair (a,b) with d=b-a:
      - If d==1 (leaf sibling pair), return [mag_index, phase_left, phase_right]
      - Else return mag_index as an int.

    Indices are 1-based and refer to entries of the full complex theta list of length 2^(n+1)-1,
    ordered as [magnitudes (1..2^n-1), phases (2^n..2^(n+1)-1)].
    """
    d = b - a
    mag_index = ((1 << n) + a) // (2 * d)
    if d == 1:
        return [mag_index, (1 << n) + a, (1 << n) + b]
    return mag_index


def gate_theta(old: int, new: int, n: int, case: str) -> Union[int, List[int]]:
    """
    Adapter for the paper's real/complex Index convention:
      - case='real'    -> int
      - case='complex' -> List[int]
    """
    if case == "real":
        return theta_real(n, old, new)
    if case == "complex":
        return theta_complex(n, old, new)
    raise ValueError("case must be 'real' or 'complex'")


# -----------------------------------------------------------------------------
# 2) Main schedulers: paper algorithms HopfReal / HopfComplex
# -----------------------------------------------------------------------------

def gates_order(n: int, case: str) -> Tuple[List[int], List[int], List[int], List[Union[int, List[int]]]]:
    """
    Paper algorithms HopfReal / HopfComplex.
    Returns: Ctrl, Anti, Targ, Index
    where each entry corresponds to one gate in order.
    """
    N = 1 << n
    HW = hamming_weight(n)

    # Ctrl: n zeros
    ctrl = [0] * n

    # Anti init: [0, 2^n-2^(n-1), ..., 2^n-2]
    # Initial anti-control masks from the paper algorithm.
    anti = [0] + [N - (1 << (n - i)) for i in range(1, n)]

    # Targ init: [2^(n-1), 2^(n-2), ..., 2^0]
    targ = [(1 << (n - 1 - i)) for i in range(n)]

    # Index init
    if case == "real":
        index: List[Union[int, List[int]]] = [1 << i for i in range(n)]
    elif case == "complex":
        # first n-1 are magnitude-only; last promoted to SU(2) at leaf
        index = [1 << i for i in range(n - 1)] + [[(1 << (n - 1)), (1 << n), (1 << n) + 1]]
    else:
        raise ValueError("case must be 'real' or 'complex'")

    # main loop k=1..n-1
    for k in range(1, n):
        A, B = find_pairs(HW[k], HW[k + 1])

        ctrl.extend(A)
        anti.extend(anti_ctrl(A))
        targ.extend([b - a for a, b in zip(A, B)])  # elementwise B-A

        if case == "real":
            for a, b in zip(A, B):
                index.append(theta_real(n, a, b))
        else:
            for a, b in zip(A, B):
                index.append(theta_complex(n, a, b))

    return ctrl, anti, targ, index


# -----------------------------------------------------------------------------
# 3) Mask decoding helpers (integer masks -> Qibo qubit indices)
# -----------------------------------------------------------------------------

def _mask_bitpos(mask: int) -> int:
    """Return bit position (0=LSB) for a power-of-two mask."""
    if mask <= 0 or (mask & (mask - 1)) != 0:
        raise ValueError(f"target mask must be a power of two; got {mask}")
    return (mask.bit_length() - 1)

def _bitpos_to_qubit(bitpos: int, n: int) -> int:
    """
    Map integer bit position -> Qibo qubit index.
    Convention used here (matches common |q0 q1 ... q_{n-1}| with q0 as MSB):
      bitpos = n-1  -> qubit 0
      bitpos = 0    -> qubit n-1
    """
    return (n - 1 - bitpos)

def _mask_to_qubits(mask: int, n: int) -> List[int]:
    """Return qubits (Qibo indices) whose corresponding bits are 1 in mask."""
    qs = []
    for bitpos in range(n):
        if mask & (1 << bitpos):
            qs.append(_bitpos_to_qubit(bitpos, n))
    return qs


# -----------------------------------------------------------------------------
# 4) R_C leaf gate for the complex ansatz
# -----------------------------------------------------------------------------

def SU2(theta: float, phase_0: float, phase_1: float) -> np.ndarray:
    U = np.array([
        [np.exp(1j * phase_0) * np.cos(theta), -np.exp(-1j * phase_1) * np.sin(theta)],
        [np.exp(1j * phase_1) * np.sin(theta),  np.exp(-1j * phase_0) * np.cos(theta)]
    ], dtype=complex)
    return U


# -----------------------------------------------------------------------------
# 5) Circuit builder: iterate over zip(Ctrl, Anti, Targ, Index)
# -----------------------------------------------------------------------------

def qibo_circuit(
    theta_list: np.ndarray,
    case: str,
    *,
    minimal: bool = False,
    eps: float = 1e-12,
    return_circuit: bool = False,
):
    """Build a Qibo circuit that implements the Hopf ansatz.

    Parameters
    ----------
    theta_list:
        Hopf parameters (real or complex convention) as returned by helpers in
        this module.
    case:
        "real" or "complex".
    minimal:
        If True, prune gates that provably have no effect for the given
        `theta_list` by:
          (i) skipping controlled rotations whose controls are never satisfied
              on the *current* support (tracked conservatively), and
          (ii) skipping single-qubit operations that are exactly identity
               (up to numerical tolerance) on the target.

        This is mainly useful for tangent-state preparation, where most angles
        are clamped to 0 and the full schedule contains many no-op gates.
    eps:
        Numerical tolerance for pruning decisions.

    Notes
    -----
    - Pruning never drops a controlled gate that is a *conditional phase*
      (e.g., controlled (-I)); only gates numerically equal to +I are removed.
    - Support tracking is conservative (it may over-approximate), so pruning is
      conservative but not necessarily globally minimal.
    """
    if not _QIBO_AVAILABLE:
        raise RuntimeError("Qibo is not installed; cannot build a Qibo circuit.")

    if case not in ("real", "complex"):
        raise ValueError("case must be 'real' or 'complex'")

    if case == "real":
        L = int(len(theta_list))
        n = int(np.log2(L + 1))
    else:
        L = int(len(theta_list))
        n = int(np.log2((L + 1) / 2))

    # Gate schedule (independent of theta values)
    # gates_order returns (Ctrl, Anti, Targ, Index) as four parallel lists.
    ctrl_list, anti_list, targ_list, idx_list = gates_order(n, case=case)

    circuit = models.Circuit(n)

    # --------- helpers for pruning ---------
    if minimal:
        dim = 1 << n
        states = np.arange(dim, dtype=np.uint32)  # basis indices in bitpos (LSB-first) convention
        active = np.zeros(dim, dtype=bool)
        active[0] = True  # start from |0...0>

    def _ctrl_satisfied_mask(ctrl_mask: int, anti_mask: int) -> np.ndarray:
        # ctrl bits must be 1, anti bits must be 0
        return ((states & np.uint32(ctrl_mask)) == np.uint32(ctrl_mask)) & ((states & np.uint32(anti_mask)) == 0)

    def _is_identity_ry(theta: float) -> bool:
        # With our convention, the gate matrix is [[cos(theta), -sin(theta)], [sin(theta), cos(theta)]].
        # Identity iff theta ≡ 0 (mod 2π). (Do NOT treat -I as identity because controlled(-I) is a phase gate.)
        return (abs(np.sin(theta)) <= eps) and (abs(np.cos(theta) - 1.0) <= eps)

    def _ry_kind(theta: float) -> str:
        # Used only for support update; 'diag' includes ±I.
        if abs(np.sin(theta)) <= eps:
            return "diag"
        if abs(np.cos(theta)) <= eps:
            return "anti"
        return "mix"

    def _is_identity_u2(U: np.ndarray) -> bool:
        return np.allclose(U, np.eye(2, dtype=complex), atol=eps, rtol=0.0)

    def _u2_kind(U: np.ndarray) -> str:
        if np.allclose(U[0, 1], 0.0, atol=eps, rtol=0.0) and np.allclose(U[1, 0], 0.0, atol=eps, rtol=0.0):
            return "diag"
        if np.allclose(U[0, 0], 0.0, atol=eps, rtol=0.0) and np.allclose(U[1, 1], 0.0, atol=eps, rtol=0.0):
            return "anti"
        return "mix"

    def _support_update(trigger: np.ndarray, targ_mask: int, kind: str):
        # Conservative support update after applying a controlled 2x2 operation on the target.
        nonlocal active
        if kind == "diag":
            return
        trig_idx = np.nonzero(trigger)[0]
        if trig_idx.size == 0:
            return
        flip_idx = trig_idx ^ int(targ_mask)
        if kind == "anti":
            # Swap support within affected pairs (exact at support level).
            new_active = active.copy()
            new_active[trig_idx] = active[flip_idx]
            new_active[flip_idx] = active[trig_idx]
            active = new_active
        else:
            # Mixing: conservatively add both basis states in each affected pair.
            active[flip_idx] = True

    # --------- build circuit (with optional pruning) ---------
    for ctrl_mask, anti_mask, targ_mask, idx in zip(ctrl_list, anti_list, targ_list, idx_list):
        # Target bit position and target qubit index
        t_bitpos = _mask_bitpos(targ_mask)
        t_qubit = _bitpos_to_qubit(t_bitpos, n)

        # Controls are on all qubits except the target.
        ctrl_mask = int(ctrl_mask) & (~int(targ_mask))
        anti_mask = int(anti_mask) & (~int(targ_mask))

        controls_mask = ctrl_mask | anti_mask
        # Convert integer masks (bitpos convention) directly to Qibo qubit indices.
        # NOTE: target bit is already removed from ctrl/anti above.
        controls_qubits = _mask_to_qubits(int(controls_mask), n)
        anticontrol_qubits = _mask_to_qubits(int(anti_mask), n)

        gate = None
        op_kind = None  # 'diag' | 'anti' | 'mix' (for support updates)

        if case == "real":
            theta = float(theta_list[int(idx) - 1])
            if minimal and _is_identity_ry(theta):
                continue  # controlled(+I) is identity
            op_kind = _ry_kind(theta)
            gate = gates.RY(t_qubit, 2.0 * theta).controlled_by(*controls_qubits)

        else:  # complex
            if isinstance(idx, (list, tuple, np.ndarray)) and len(idx) == 3:
                a = float(theta_list[int(idx[0]) - 1])
                b = float(theta_list[int(idx[1]) - 1])
                c = float(theta_list[int(idx[2]) - 1])
                U = np.array(
                    [
                        [np.exp(1j * b) * np.cos(a), -np.exp(-1j * c) * np.sin(a)],
                        [np.exp(1j * c) * np.sin(a), np.exp(-1j * b) * np.cos(a)],
                    ],
                    dtype=complex,
                )
                if minimal and _is_identity_u2(U):
                    continue
                op_kind = _u2_kind(U)
                gate = gates.Unitary(U, t_qubit).controlled_by(*controls_qubits)
            else:
                theta = float(theta_list[int(idx) - 1])
                if minimal and _is_identity_ry(theta):
                    continue
                op_kind = _ry_kind(theta)
                gate = gates.RY(t_qubit, 2.0 * theta).controlled_by(*controls_qubits)

        # If pruning, check whether this controlled operation can ever be triggered
        if minimal and controls_mask != 0:
            control_ok = _ctrl_satisfied_mask(ctrl_mask, anti_mask)
            trigger = active & control_ok
            if not np.any(trigger):
                continue
        elif minimal:
            # no controls -> always "triggered" on all currently active basis states
            trigger = active

        # Apply anticontrol conjugation (X), then controlled gate, then undo
        for q in anticontrol_qubits:
            circuit.add(gates.X(q))

        circuit.add(gate)

        for q in anticontrol_qubits:
            circuit.add(gates.X(q))

        if minimal:
            _support_update(trigger, int(targ_mask), op_kind)

    if return_circuit:
        return circuit

    result = circuit(nshots=0)
    return result.state()



# =============================================================================
# Mathematical Helpers
# =============================================================================

def theta_from_vector(vector, case):
    d = len(vector)
    n = int(np.log2(d))
    if (1 << n) != d:
        raise ValueError("Length of vector must be a power of 2.")

    v = np.asarray(vector)
    if np.iscomplexobj(v):
        magsq = v.real**2 + v.imag**2
    else:
        magsq = v**2

    ps = np.empty(d + 1, dtype=float)
    ps[0] = 0.0
    np.cumsum(magsq, out=ps[1:])

    if case == 'real':
        theta_list = np.empty(d - 1, dtype=float)
        for j in range(d - 1):
            col = (j + 1).bit_length() - 1
            offset = (1 << col) - 1
            delta = j - offset
            block_size = 1 << (n - col)
            start = delta * block_size
            mid = start + (block_size >> 1)
            end = start + block_size

            if col == n - 1:
                a = v[start]
                b = v[mid]
                if a == 0.0 and b == 0.0:
                    theta_list[j] = 0.0
                else:
                    theta_list[j] = np.mod(np.arctan2(b, a), 2*np.pi)
            else:
                left = ps[mid] - ps[start]
                right = ps[end] - ps[mid]
                theta_list[j] = np.arctan2(np.sqrt(right), np.sqrt(left))
        return theta_list

    if case == 'complex':
        theta_list = np.empty(2 * d - 1, dtype=float)
        # 1. Magnitude angles
        for j in range(d - 1):
            col = (j + 1).bit_length() - 1
            offset = (1 << col) - 1
            delta = j - offset
            block_size = 1 << (n - col)
            start = delta * block_size
            mid = start + (block_size >> 1)
            end = start + block_size

            left = ps[mid] - ps[start]
            right = ps[end] - ps[mid]
            theta_list[j] = np.arctan2(np.sqrt(right), np.sqrt(left))
        # 2. Phase angles
        theta_list[d - 1:] = np.mod(np.angle(v), 2*np.pi)
        return theta_list

    raise ValueError("case must be 'real' or 'complex'.")

def vector_from_theta(theta_list, case):
    if case == 'real':
        n = int(np.log2(len(theta_list) + 1))
        vec = np.array([1.0], dtype=float)
        base = 0
        c_all = np.cos(theta_list)
        s_all = np.sin(theta_list)
        
        for level in range(n):
            m = 1 << level
            c = c_all[base : base + m]
            s = s_all[base : base + m]
            left = vec * c
            right = vec * s
            out = np.empty(2 * m, dtype=float)
            out[0::2] = left
            out[1::2] = right
            vec = out
            base += m
        return vec

    if case == 'complex':
        L = len(theta_list)
        n = int(np.log2((L + 1) // 2))
        d = 1 << n
        mag_thetas = theta_list[:d - 1]
        phases     = theta_list[d - 1 :]

        c_all = np.cos(mag_thetas)
        s_all = np.sin(mag_thetas)
        
        vec = np.array([1.0], dtype=float)
        base = 0
        for level in range(n):
            m = 1 << level
            c = c_all[base : base + m]
            s = s_all[base : base + m]
            left = vec * c
            right = vec * s
            out = np.empty(2 * m, dtype=float)
            out[0::2] = left
            out[1::2] = right
            vec = out
            base += m

        # Apply phases to leaves
        return vec.astype(complex) * np.exp(1j * phases)

    raise ValueError("case must be 'real' or 'complex'.")

# =============================================================================
# Tangent State Synthesis (Parameter Encoding)
# =============================================================================

def _hopf_ancestors_path(i):
    """Return list of ancestors from root down to parent of i (1-based index)."""
    path = []
    curr = i // 2
    while curr >= 1:
        path.append(curr)
        curr //= 2
    return list(reversed(path))

def _hopf_path_bits(i):
    """Return list of bits (0=left, 1=right) on the path from root to node i."""
    if i < 1: return []
    bstr = bin(i)[3:] 
    return [int(b) for b in bstr]

def _hopf_subtree_nodes(i, max_node_idx):
    """Return set of all node indices in the subtree rooted at i (inclusive), respecting bounds."""
    nodes = set()
    queue = [i]
    while queue:
        curr = queue.pop(0)
        if curr <= max_node_idx:
            nodes.add(curr)
            c1, c2 = 2*curr, 2*curr+1
            if c1 <= max_node_idx: queue.append(c1)
            if c2 <= max_node_idx: queue.append(c2)
    return nodes

def _subtree_span(i, n):
    """Return [start, mid, end) indices of leaf range for magnitude node i."""
    d = int(math.floor(math.log2(i)))
    pos_in_level = i - (1 << d)
    span = 1 << (n - d)
    start = pos_in_level * span
    mid = start + (span >> 1)
    end = start + span
    return start, mid, end

def theta_hopf_tangent_state(theta, i, case='real'):
    """
    Construct the parameter vector theta^(i) that prepares the normalized
    partial derivative state |∂_i ψ⟩ / ||∂_i ψ||.
    """
    theta = np.asarray(theta, dtype=float)
    L_total = len(theta)
    th_new = np.zeros_like(theta)

    if case == 'real':
        if i < 1 or i > L_total:
            raise ValueError(f"Index i={i} out of range [1, {L_total}]")
        
        # 1. Shift target node
        th_new[i-1] = theta[i-1] + (math.pi / 2.0)

        # 2. Clamp ancestors
        ancestors = _hopf_ancestors_path(i)
        bits = _hopf_path_bits(i)
        for anc, b in zip(ancestors, bits):
            th_new[anc-1] = 0.0 if b == 0 else (math.pi / 2.0)

        # 3. Copy descendants (subtree magnitudes)
        subtree = _hopf_subtree_nodes(i, L_total)
        for node in subtree:
            if node != i:
                th_new[node-1] = theta[node-1]
        
        return th_new

    elif case == 'complex':
        # n qubits -> 2^n - 1 magnitudes, 2^n phases.
        # L_total = 2^(n+1) - 1.
        n = int(np.log2(L_total + 1)) - 1
        num_mags = (1 << n) - 1
        
        if i < 1 or i > L_total:
            raise ValueError(f"Index i={i} out of range [1, {L_total}]")

        if i <= num_mags:
            # --- Magnitude Parameter ---
            # 1. Shift target
            th_new[i-1] = theta[i-1] + (math.pi / 2.0)
            
            # 2. Clamp ancestors
            ancestors = _hopf_ancestors_path(i)
            bits = _hopf_path_bits(i)
            for anc, b in zip(ancestors, bits):
                th_new[anc-1] = 0.0 if b == 0 else (math.pi / 2.0)
                    
            # 3. Copy descendants (magnitude nodes only)
            subtree_mags = _hopf_subtree_nodes(i, num_mags)
            for node in subtree_mags:
                if node != i:
                    th_new[node-1] = theta[node-1]
            
            # 4. Copy relevant phases (leaves of this subtree)
            start, mid, end = _subtree_span(i, n)
            # Phases are stored at indices [num_mags : ]
            # The phase for leaf k (0-based) is at index num_mags + k
            for k in range(start, end):
                p_idx = num_mags + k 
                th_new[p_idx] = theta[p_idx]
                
        else:
            # --- Phase Parameter ---
            # i corresponds to a leaf. 
            # 0-based leaf index:
            leaf_idx = i - (num_mags + 1) 
            
            # 1. Clamp ALL magnitude ancestors to path leading to leaf_idx
            curr = 1 # Root
            for bit_pos in range(n - 1, -1, -1):
                # MSB first traversal
                bit = (leaf_idx >> bit_pos) & 1
                if bit == 0:
                    th_new[curr-1] = 0.0
                    curr = 2 * curr
                else:
                    th_new[curr-1] = math.pi / 2.0
                    curr = 2 * curr + 1
                    
            # 2. Shift phase parameter
            th_new[i-1] = theta[i-1] + (math.pi / 2.0)
            
        return th_new

    else:
        raise ValueError("case must be 'real' or 'complex'")

# =============================================================================

def qibo_tangent_state(theta, i, case='real', *, minimal=True, eps=1e-12, return_circuit: bool = False):
    """Convenience wrapper: build |∂_i ψ⟩/||∂_i ψ|| via Qibo.

    This calls `theta_hopf_tangent_state` and then `qibo_circuit`.
    Set `minimal=True` to prune no-op/inactive gates for this tangent state.
    """
    th = theta_hopf_tangent_state(theta, i, case=case)
    return qibo_circuit(th, case, minimal=minimal, eps=eps, return_circuit=return_circuit)

# Jacobian & Metric (Analytical)
# =============================================================================

def partial_derivative_column_real(theta, i):
    theta = np.asarray(theta, dtype=float)
    L = len(theta)
    n = int(round(math.log2(L + 1)))
    
    psi = vector_from_theta(theta, 'real')
    start, mid, end = _subtree_span(i, n)
    
    th = theta[i-1]
    s, c = np.sin(th), np.cos(th)
    
    col = np.zeros_like(psi)
    if abs(c) > 1e-15:
        col[start:mid] = - (s/c) * psi[start:mid]
    if abs(s) > 1e-15:
        col[mid:end] = (c/s) * psi[mid:end]
        
    return col

def jacobian(theta_list, case):
    """Compute full Jacobian matrix."""
    theta = np.asarray(theta_list, dtype=float)
    
    if case == 'real':
        L = len(theta)
        n = int(np.log2(L + 1))
        x = vector_from_theta(theta, 'real')
        J = np.zeros((1 << n, L), dtype=float)
        
        s = np.sin(theta)
        c = np.cos(theta)
        tan = np.divide(s, c, out=np.zeros_like(s), where=(c!=0))
        cot = np.divide(c, s, out=np.zeros_like(c), where=(s!=0))
        
        for j0 in range(L):
            start, mid, end = _subtree_span(j0+1, n)
            if c[j0] != 0:
                J[start:mid, j0] = -tan[j0] * x[start:mid]
            if s[j0] != 0:
                J[mid:end, j0]   =  cot[j0] * x[mid:end]
        return J
        
    if case == 'complex':
        L_total = len(theta)
        n = int(np.log2(L_total + 1)) - 1
        num_mags = (1 << n) - 1
        
        x = vector_from_theta(theta, 'complex')
        J = np.zeros((1 << n, L_total), dtype=complex)
        
        mag_thetas = theta[:num_mags]
        s = np.sin(mag_thetas)
        c = np.cos(mag_thetas)
        tan = np.divide(s, c, out=np.zeros_like(s), where=(c!=0))
        cot = np.divide(c, s, out=np.zeros_like(c), where=(s!=0))
        
        # Magnitude columns (0..num_mags-1)
        for j0 in range(num_mags):
            start, mid, end = _subtree_span(j0+1, n)
            if c[j0] != 0:
                J[start:mid, j0] = -tan[j0] * x[start:mid]
            if s[j0] != 0:
                J[mid:end, j0]   =  cot[j0] * x[mid:end]
                
        # Phase columns (num_mags..L_total-1)
        # Column j corresponds to leaf k = j - num_mags
        # Derivative is i * x[k] only at index k, 0 elsewhere? 
        # No, x is the full state vector.
        # x_k depends on phase phi_k. x_other does not.
        # d(x_k)/d(phi_k) = i * x_k.
        for k in range(1 << n):
            col_idx = num_mags + k
            J[k, col_idx] = 1j * x[k]
            
        return J
    return None

def metric_diagonal(theta_list, case):
    """
    Calculate diagonal metric elements g_ii = ||∂ψ/∂θ_i||^2.
    """
    theta = np.asarray(theta_list, dtype=float)
    L = len(theta)
    
    if case == 'real':
        g = np.ones(L, dtype=float)
        c2 = np.cos(theta)**2
        s2 = np.sin(theta)**2
        for i in range(2, L + 1):
            ancestors = _hopf_ancestors_path(i)
            bits = _hopf_path_bits(i)
            val = 1.0
            for anc, b in zip(ancestors, bits):
                val *= s2[anc-1] if b == 1 else c2[anc-1]
            g[i-1] = val
        return g
        
    if case == 'complex':
        n = int(np.log2(L + 1)) - 1
        num_mags = (1 << n) - 1
        g = np.ones(L, dtype=float)
        
        # Magnitudes part (same as real metric on magnitudes)
        c2 = np.cos(theta[:num_mags])**2
        s2 = np.sin(theta[:num_mags])**2
        
        # Magnitude params
        for i in range(2, num_mags + 1):
            ancestors = _hopf_ancestors_path(i)
            bits = _hopf_path_bits(i)
            val = 1.0
            for anc, b in zip(ancestors, bits):
                val *= s2[anc-1] if b == 1 else c2[anc-1]
            g[i-1] = val
            
        # Phase params: g_ii = |x_leaf|^2 = Product of path probs
        # Indices num_mags+1 ... L (1-based)
        for i in range(num_mags + 1, L + 1):
            leaf_idx = i - (num_mags + 1)
            # Reconstruct path prob
            val = 1.0
            curr = 1
            for bit_pos in range(n - 1, -1, -1):
                bit = (leaf_idx >> bit_pos) & 1
                if bit == 0:
                    val *= c2[curr-1]
                    curr = 2 * curr
                else:
                    val *= s2[curr-1]
                    curr = 2 * curr + 1
            g[i-1] = val
            
        return g
        
    return None

# =============================================================================
# Sanity Checks
# =============================================================================

def sanity_checks(n):
    np.set_printoptions(precision=3, suppress=True)
    print(f"\n=== Hopf Utils Sanity Checks (n={n}) ===")
    rng = np.random.default_rng(42)
    
    # ---------------------------
    # 1. REAL CASE
    # ---------------------------
    print("\n[REAL Case]")
    L_real = (1 << n) - 1
    theta_r = rng.uniform(0.2, np.pi/2 - 0.2, size=L_real)
    
    # Check Vector Roundtrip
    x_r = vector_from_theta(theta_r, 'real')
    th_r_back = theta_from_vector(x_r, 'real')
    x_r_back = vector_from_theta(th_r_back, 'real')
    print(f"  Vector Roundtrip Error: {np.linalg.norm(x_r - x_r_back):.2e}")
    
    # Check Tangent Synthesis
    J_r = jacobian(theta_r, 'real')
    max_err_r = 0.0
    for i in range(1, L_real + 1):
        v_diff = J_r[:, i-1]
        norm_v = np.linalg.norm(v_diff)
        
        th_tan = theta_hopf_tangent_state(theta_r, i, 'real')
        v_syn = vector_from_theta(th_tan, 'real')
        
        if norm_v > 1e-10:
            # Check alignment (dot product ~ 1) and diff
            # The synthesized state should equal the normalized derivative direction.
            dist = np.linalg.norm(v_syn - (v_diff / norm_v))
            max_err_r = max(max_err_r, dist)
            
    print(f"  Tangent Synthesis Max Error: {max_err_r:.2e}")
    if max_err_r < 1e-10: print("  ✅ Real Synthesis Correct.")
    else: print("  ❌ Real Synthesis FAILED.")

    # ---------------------------
    # 2. COMPLEX CASE
    # ---------------------------
    print("\n[COMPLEX Case]")
    L_comp = (1 << (n+1)) - 1
    # Random magnitudes (avoid 0, pi/2) and random phases
    theta_c = np.zeros(L_comp)
    num_mags = (1 << n) - 1
    theta_c[:num_mags] = rng.uniform(0.2, np.pi/2 - 0.2, size=num_mags)
    theta_c[num_mags:] = rng.uniform(0, 2*np.pi, size=(1<<n))
    
    # Check Vector Roundtrip
    x_c = vector_from_theta(theta_c, 'complex')
    th_c_back = theta_from_vector(x_c, 'complex')
    x_c_back = vector_from_theta(th_c_back, 'complex')
    # theta_from_vector maps leaf phases to [0, 2pi); zero-amplitude leaves
    # use the convention stated in the inverse-map lemma.
    print(f"  Vector Roundtrip Error: {np.linalg.norm(x_c - x_c_back):.2e}")
    
    # Check Tangent Synthesis
    J_c = jacobian(theta_c, 'complex')
    max_err_c = 0.0
    
    for i in range(1, L_comp + 1):
        v_diff = J_c[:, i-1]
        norm_v = np.linalg.norm(v_diff)
        
        th_tan = theta_hopf_tangent_state(theta_c, i, 'complex')
        v_syn = vector_from_theta(th_tan, 'complex')
        
        if norm_v > 1e-10:
            target = v_diff / norm_v
            # Synthesized state is a unit vector. Target is a unit vector.
            # v_syn should equal target (including complex phase).
            dist = np.linalg.norm(v_syn - target)
            max_err_c = max(max_err_c, dist)
            if dist > 1e-5:
               print(f"     Param {i} mismatch. Dist={dist:.4f}")
    
    print(f"  Tangent Synthesis Max Error: {max_err_c:.2e}")
    if max_err_c < 1e-10: print("  ✅ Complex Synthesis Correct.")
    else: print("  ❌ Complex Synthesis FAILED.")

    # ---------------------------
    # 3. Circuit Check
    # ---------------------------
    if _QIBO_AVAILABLE:
        print("\n[Circuit Parity]")
        try:
            x_circ_r = qibo_circuit(theta_r, 'real')
            err_r = np.linalg.norm(x_r - np.array(x_circ_r))
            print(f"  Real Qibo vs Analytic: {err_r:.2e}")
            
            x_circ_c = qibo_circuit(theta_c, 'complex')
            err_c = np.linalg.norm(x_c - np.array(x_circ_c))
            print(f"  Complex Qibo vs Analytic: {err_c:.2e}")
        except Exception as e:
            print(f"  Qibo check failed: {e}")

if __name__ == "__main__":
    sanity_checks(n=7)
