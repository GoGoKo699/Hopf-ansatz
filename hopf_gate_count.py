#!/usr/bin/env python3

"""hopf_gate_count.py

CNOT-count safeguard for the Hopf ansatz gate schedules.

The script uses hopf_utils.gates_order to generate the real and complex Hopf
schedules, then evaluates the no-clean-ancilla CNOT model stated in the paper.
It prints a table and saves a normalized plot of #CNOT/(n 2^n).

Control semantics:
    * Ctrl bit 1  -> positive control on that physical qubit.
    * Anti bit 1  -> negative control on that physical qubit.
    * m = popcount(Ctrl | Anti) is the number of controls.
    * An Index entry of length three denotes a promoted final-layer R_C gate;
      an integer Index entry denotes a controlled R_y gate.

CNOT model:
    * m=0: 0.
    * 1 <= m <= 4: 2^(m+1)-2.
    * m >= 5:
        controlled R_y: 16*(m+1)-40,
        controlled R_C: 20*(m+1)-38 if (m+1) is odd, otherwise 20*(m+1)-42.

Closed-form counts checked against the generated schedules:
    * Real:
        G_R(n) = sum_{m=0}^{n-1} binom(n,m) c_Ry(m).
    * Complex:
        G_C(n) = sum_{m=0}^{n-1}
                 [binom(n-1,m)c_Ry(m) + binom(n-1,m-1)c_RC(m)]
                 + c_RC(n-1) - c_Ry(n-1),
      where the final correction accounts for the top-control promoted R_C gate
      structure of HopfComplex.

Usage:
    python3 hopf_gate_count.py
    python3 hopf_gate_count.py --nmin 4 --nmax 20 --out hopf_gate_count.pdf
    python3 hopf_gate_count.py --show
"""

from __future__ import annotations

from typing import Any, Dict, List
import argparse
import math
import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import hopf_utils as hu


# ------------------ helpers ------------------

def popcount(x: int) -> int:
    x = int(x)
    try:
        return x.bit_count()          # Python >= 3.8
    except AttributeError:
        return bin(x).count("1")      # Python <= 3.7


# ------------------ CNOT cost model ------------------

def cnot_cost_small(m: int) -> int:
    if m <= 0:
        return 0
    return (1 << (m + 1)) - 2


def cnot_cost_ry(m: int) -> int:
    if m <= 4:
        return cnot_cost_small(m)
    ng = m + 1
    return 16 * ng - 40


def cnot_cost_rc(m: int) -> int:
    if m <= 4:
        return cnot_cost_small(m)
    ng = m + 1
    return (20 * ng - 38) if (ng % 2 == 1) else (20 * ng - 42)


def is_rc_gate(idx: Any) -> bool:
    return isinstance(idx, (list, tuple)) and len(idx) == 3


# ------------------ numerical counting from gates_order ------------------

def count_for_n(n: int) -> Dict[str, int]:
    """Numerical CNOT count from hopf_utils.gates_order under the paper's Ctrl/Anti semantics."""
    ctrl_r, anti_r, targ_r, idx_r = hu.gates_order(n, "real")
    ctrl_c, anti_c, targ_c, idx_c = hu.gates_order(n, "complex")

    total_r = 0
    for cm, am in zip(ctrl_r, anti_r):
        m = popcount(cm | am)
        total_r += cnot_cost_ry(m)

    total_c = 0
    for cm, am, ix in zip(ctrl_c, anti_c, idx_c):
        m = popcount(cm | am)
        total_c += cnot_cost_rc(m) if is_rc_gate(ix) else cnot_cost_ry(m)

    return {"n": n, "HopfReal_CNOT": total_r, "HopfComplex_CNOT": total_c}


# ------------------ theory (binomial sums) ------------------

def C(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def G_real_theory(n: int) -> int:
    return sum(C(n, m) * cnot_cost_ry(m) for m in range(0, n))


def G_complex_theory(n: int) -> int:
    total = 0
    for m in range(0, n):
        total += C(n - 1, m) * cnot_cost_ry(m)
        total += C(n - 1, m - 1) * cnot_cost_rc(m)

    # Top-control correction:
    # the naive binomial sum counts one (n-1)-controlled Ry and only n-1
    # promoted R_C gates; the actual HopfComplex schedule has no such Ry
    # and has n promoted R_C gates at m=n-1.
    total += cnot_cost_rc(n - 1) - cnot_cost_ry(n - 1)
    return total


def norm_factor(n: int) -> int:
    return n * (1 << n)


# ------------------ plotting ------------------

def make_plot(nmin: int, nmax: int, out: str, show: bool) -> None:
    import matplotlib.pyplot as plt

    ns: List[int] = list(range(nmin, nmax + 1))

    y_r_num: List[float] = []
    y_c_num: List[float] = []
    y_r_th: List[float] = []
    y_c_th: List[float] = []

    for n in ns:
        denom = float(norm_factor(n))

        r = count_for_n(n)
        y_r_num.append(r["HopfReal_CNOT"] / denom)
        y_c_num.append(r["HopfComplex_CNOT"] / denom)

        y_r_th.append(G_real_theory(n) / denom)
        y_c_th.append(G_complex_theory(n) / denom)

    plt.figure()
    plt.plot(ns, y_r_th, label="Hopf real (theory)", color="tab:blue", linewidth=2)
    plt.plot(ns, y_c_th, label="Hopf complex (theory)",  color="tab:green", linewidth=2)
    plt.scatter(ns, y_r_num, label="Hopf real (numerical)", color="tab:blue", marker="o")
    plt.scatter(ns, y_c_num, label="Hopf complex (numerical)", color="tab:green", marker="o")

    plt.axhline(8.0, linestyle="--", color="tab:blue", alpha=0.6)
    plt.axhline(9.0, linestyle="--", color="tab:green", alpha=0.6)

    plt.xlabel("number of qubits n")
    plt.ylabel(r"#CNOT / (n * 2^n)")
    plt.xticks(ns)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    print(f"Saved plot: {out}")
    if show:
        plt.show()


# ------------------ CLI ------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=20)
    ap.add_argument("--out", type=str, default="hopf_gate_count.pdf",
                    help="Output plot filename (pdf/png/svg).")
    ap.add_argument("--show", action="store_true", help="Show plot window.")
    args = ap.parse_args()

    if args.nmin < 2 or args.nmax < args.nmin:
        raise ValueError("Invalid n range.")

    # 1) print table
    print("n  HopfReal_CNOT  HopfComplex_CNOT")
    for n in range(args.nmin, args.nmax + 1):
        r = count_for_n(n)
        print(f"{r['n']:>2} {r['HopfReal_CNOT']:>13} {r['HopfComplex_CNOT']:>16}")

    # 2) save plot
    make_plot(args.nmin, args.nmax, args.out, args.show)


if __name__ == "__main__":
    main()
