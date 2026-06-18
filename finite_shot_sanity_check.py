#!/usr/bin/env python3
"""
finite_shot_sanity_check.py

Finite-shot sanity check for the symmetric signed-branch Hopf gradient
estimator used in the paper.

The experiment fixes one real n-qubit Hopf state and one diagonal normalized
Hamming-spectrum Hamiltonian. For every Hopf coordinate it prepares the exact
normalized tangent state, forms the + and - signed branch distributions, and
estimates the two branch energies by finite-shot sampling. The resulting full
coordinate gradient is compared with an independent exact tree-gradient
reference over repeated trials.

The command-line value ``--shots S`` follows the convention used for the paper
figure: each of the + and - branch-energy estimates for each gradient component
receives S samples. Thus each component uses 2S branch-state measurements in
this symmetric estimator.

Default run:
    python finite_shot_sanity_check.py

Default outputs:
    finite_shot_sanity_data.csv
    finite_shot_sanity.png

Pass ``--formats png pdf`` to generate both raster and vector summaries.

The script depends only on NumPy, Matplotlib, and ``hopf_utils.py`` in the same
repository. It is a fixed-state estimator check, not a complete hardware or
full-gradient sampling-complexity benchmark.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import hopf_utils as hopf


DEFAULT_SHOTS = (100, 1_000, 10_000)
DEFAULT_TRIALS = 50


def _normalized_probabilities(values: np.ndarray) -> Tuple[np.ndarray, float]:
    """Return a normalized probability vector and its pre-normalization error."""
    probs = np.asarray(values, dtype=float)
    if probs.ndim != 1 or np.any(~np.isfinite(probs)):
        raise ValueError("Probabilities must be a finite one-dimensional array.")
    if np.any(probs < -1e-14):
        raise ValueError("Encountered a materially negative probability.")
    probs = np.maximum(probs, 0.0)
    total = float(np.sum(probs))
    if total <= 0.0:
        raise ValueError("Probability vector has zero total mass.")
    return probs / total, abs(total - 1.0)


def _hopf_coordinate_gradient_real(theta: np.ndarray, euclidean_grad: np.ndarray) -> np.ndarray:
    """Independent O(2^n) real-Hopf tree gradient used as the exact reference."""
    theta = np.asarray(theta, dtype=float)
    euclidean_grad = np.asarray(euclidean_grad, dtype=float)
    length = theta.size
    n = int(round(math.log2(length + 1)))
    dim = 1 << n
    if length != dim - 1:
        raise ValueError("Real Hopf theta length must be 2^n - 1.")
    if euclidean_grad.shape != (dim,):
        raise ValueError(f"euclidean_grad must have shape ({dim},).")

    subtree_response = np.zeros(length + 1, dtype=float)
    for node in range(length, 0, -1):
        start, mid, _ = hopf._subtree_span(node, n)
        left_child = 2 * node
        right_child = left_child + 1
        left = subtree_response[left_child] if left_child <= length else euclidean_grad[start]
        right = subtree_response[right_child] if right_child <= length else euclidean_grad[mid]
        angle = float(theta[node - 1])
        subtree_response[node] = math.cos(angle) * left + math.sin(angle) * right

    gradient = np.zeros(length, dtype=float)
    prefix = np.zeros(length + 1, dtype=float)
    prefix[1] = 1.0
    for node in range(1, length + 1):
        start, mid, _ = hopf._subtree_span(node, n)
        left_child = 2 * node
        right_child = left_child + 1
        left = subtree_response[left_child] if left_child <= length else euclidean_grad[start]
        right = subtree_response[right_child] if right_child <= length else euclidean_grad[mid]
        angle = float(theta[node - 1])
        sine = math.sin(angle)
        cosine = math.cos(angle)
        gradient[node - 1] = prefix[node] * (-sine * left + cosine * right)
        if left_child <= length:
            prefix[left_child] = prefix[node] * cosine
        if right_child <= length:
            prefix[right_child] = prefix[node] * sine
    return gradient


def prepare_fixed_experiment(
    n: int,
    state_seed: int,
) -> Dict[str, np.ndarray | float]:
    """Prepare the fixed state, exact gradient, and all signed-branch distributions."""
    if n < 1:
        raise ValueError("n must be at least 1.")

    dim = 1 << n
    length = dim - 1

    # Normalized diagonal Hamming spectrum in the computational basis.
    diagonal = np.fromiter((int(i).bit_count() for i in range(dim)), dtype=float, count=dim)
    diagonal /= max(1.0, float(np.max(diagonal)))

    # Preserve the state convention used for the paper's finite-shot figure.
    rng_state = np.random.default_rng(state_seed)
    theta = rng_state.uniform(0.2, math.pi / 2.0 - 0.2, size=length)
    psi = np.asarray(hopf.vector_from_theta(theta, "real"), dtype=float)
    psi /= np.linalg.norm(psi)

    metric = np.asarray(hopf.metric_diagonal(theta, "real"), dtype=float)
    exact_gradient = _hopf_coordinate_gradient_real(theta, 2.0 * diagonal * psi)

    plus_probabilities = np.empty((length, dim), dtype=float)
    minus_probabilities = np.empty((length, dim), dtype=float)
    exact_branch_gradient = np.empty(length, dtype=float)

    max_tangent_norm_error = 0.0
    max_tangent_overlap = 0.0
    max_probability_norm_error = 0.0

    for index in range(1, length + 1):
        tangent_theta = hopf.theta_hopf_tangent_state(theta, index, case="real")
        tangent = np.asarray(hopf.vector_from_theta(tangent_theta, "real"), dtype=float)
        tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm <= 1e-15:
            raise RuntimeError(f"Tangent state {index} has near-zero norm.")
        tangent /= tangent_norm

        max_tangent_norm_error = max(max_tangent_norm_error, abs(float(np.linalg.norm(tangent)) - 1.0))
        max_tangent_overlap = max(max_tangent_overlap, abs(float(np.dot(psi, tangent))))

        plus = (psi + tangent) / math.sqrt(2.0)
        minus = (psi - tangent) / math.sqrt(2.0)
        p_plus, plus_error = _normalized_probabilities(plus * plus)
        p_minus, minus_error = _normalized_probabilities(minus * minus)
        max_probability_norm_error = max(max_probability_norm_error, plus_error, minus_error)

        plus_probabilities[index - 1] = p_plus
        minus_probabilities[index - 1] = p_minus

        e_plus = float(np.dot(p_plus, diagonal))
        e_minus = float(np.dot(p_minus, diagonal))
        exact_branch_gradient[index - 1] = math.sqrt(max(metric[index - 1], 0.0)) * (e_plus - e_minus)

    exact_norm = float(np.linalg.norm(exact_gradient))
    if exact_norm <= 1e-15:
        raise RuntimeError("The selected exact gradient is too small for a relative-error test.")

    jacobian_gradient = 2.0 * np.asarray(hopf.jacobian(theta, "real"), dtype=float).T @ (diagonal * psi)

    return {
        "theta": theta,
        "psi": psi,
        "diagonal": diagonal,
        "metric": metric,
        "plus_probabilities": plus_probabilities,
        "minus_probabilities": minus_probabilities,
        "exact_gradient": exact_gradient,
        "exact_gradient_norm": exact_norm,
        "objective": float(np.dot(diagonal, psi * psi)),
        "exact_branch_residual": float(np.linalg.norm(exact_branch_gradient - exact_gradient)),
        "jacobian_reference_residual": float(np.linalg.norm(jacobian_gradient - exact_gradient)),
        "max_tangent_norm_error": max_tangent_norm_error,
        "max_tangent_overlap": max_tangent_overlap,
        "max_probability_norm_error": max_probability_norm_error,
        "state_norm_error": abs(float(np.linalg.norm(psi)) - 1.0),
    }


def estimate_gradient(
    diagonal: np.ndarray,
    metric: np.ndarray,
    plus_probabilities: np.ndarray,
    minus_probabilities: np.ndarray,
    shots_per_branch: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample the symmetric signed-branch estimator for the full coordinate gradient."""
    if shots_per_branch < 1:
        raise ValueError("shots_per_branch must be positive.")

    length = metric.size
    estimate = np.empty(length, dtype=float)
    for index in range(length):
        plus_counts = rng.multinomial(shots_per_branch, plus_probabilities[index])
        minus_counts = rng.multinomial(shots_per_branch, minus_probabilities[index])
        e_plus = float(np.dot(plus_counts, diagonal)) / shots_per_branch
        e_minus = float(np.dot(minus_counts, diagonal)) / shots_per_branch
        estimate[index] = math.sqrt(max(float(metric[index]), 0.0)) * (e_plus - e_minus)
    return estimate


def summarize_errors(rows: Iterable[Dict[str, float | int]]) -> Dict[int, Dict[str, float]]:
    grouped: Dict[int, List[float]] = {}
    for row in rows:
        shots = int(row["shots"])
        grouped.setdefault(shots, []).append(float(row["rel_error"]))

    summary: Dict[int, Dict[str, float]] = {}
    for shots, values in sorted(grouped.items()):
        array = np.asarray(values, dtype=float)
        sample_std = float(np.std(array, ddof=1)) if array.size > 1 else 0.0
        summary[shots] = {
            "mean": float(np.mean(array)),
            "std": sample_std,
            "sem": sample_std / math.sqrt(array.size) if array.size else math.nan,
            "median": float(np.median(array)),
            "trials": float(array.size),
        }
    return summary


def fitted_log_slope(summary: Dict[int, Dict[str, float]]) -> float:
    shots = np.asarray(sorted(summary), dtype=float)
    means = np.asarray([summary[int(value)]["mean"] for value in shots], dtype=float)
    if shots.size < 2 or np.any(means <= 0.0):
        return math.nan
    return float(np.polyfit(np.log10(shots), np.log10(means), 1)[0])


def save_csv(rows: Sequence[Dict[str, float | int]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "shots",
        "shots_per_branch",
        "total_branch_shots_per_component",
        "trial",
        "rel_error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_csv(path: Path) -> List[Dict[str, float | int]]:
    rows: List[Dict[str, float | int]] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            shots = int(raw.get("shots") or raw.get("shots_per_branch") or 0)
            rows.append(
                {
                    "shots": shots,
                    "shots_per_branch": int(raw.get("shots_per_branch") or shots),
                    "total_branch_shots_per_component": int(
                        raw.get("total_branch_shots_per_component") or 2 * shots
                    ),
                    "trial": int(raw["trial"]),
                    "rel_error": float(raw["rel_error"]),
                }
            )
    if not rows:
        raise ValueError(f"No data rows found in {path}.")
    return rows


def plot_result(
    summary: Dict[int, Dict[str, float]],
    output_dir: Path,
    formats: Sequence[str],
    dpi: int,
) -> List[Path]:
    shot_values = np.asarray(sorted(summary), dtype=float)
    means = np.asarray([summary[int(value)]["mean"] for value in shot_values], dtype=float)
    sems = np.asarray([summary[int(value)]["sem"] for value in shot_values], dtype=float)

    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    ax.errorbar(
        shot_values,
        means,
        yerr=sems,
        fmt="o-",
        color="#1f4e79",
        capsize=4,
        markersize=5,
        linewidth=1.5,
        label="finite-shot estimate",
    )

    reference = means[0] * math.sqrt(shot_values[0]) / np.sqrt(shot_values)
    ax.plot(
        shot_values,
        reference,
        linestyle="--",
        color="0.45",
        linewidth=1.1,
        label=r"$\propto S^{-1/2}$",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Shots per signed branch and gradient component")
    ax.set_ylabel("Relative gradient error")
    ax.grid(True, which="both", alpha=0.24)
    ax.legend(frameon=False)
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for extension in formats:
        path = output_dir / f"finite_shot_sanity.{extension}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)
    plt.close(fig)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finite-shot Hopf gradient estimator sanity check")
    parser.add_argument("--n", type=int, default=6, help="Number of qubits (default: 6).")
    parser.add_argument(
        "--shots",
        type=int,
        nargs="+",
        default=list(DEFAULT_SHOTS),
        help="Shots allocated to each signed branch per component (default: 100 1000 10000).",
    )
    parser.add_argument("--trials", type=int, default=DEFAULT_TRIALS, help="Trials per shot count (default: 50).")
    parser.add_argument("--state-seed", type=int, default=123, help="Fixed-state seed (default: 123).")
    parser.add_argument("--estimator-seed", type=int, default=999, help="Sampling seed (default: 999).")
    parser.add_argument("--output-dir", type=Path, default=HERE, help="Output directory (default: script directory).")
    parser.add_argument("--formats", nargs="+", choices=("png", "pdf"), default=["png"])
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate the plot and summary from finite_shot_sanity_data.csv without rerunning samples.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    csv_path = output_dir / "finite_shot_sanity_data.csv"

    if args.plot_only:
        rows = load_csv(csv_path)
        summary = summarize_errors(rows)
        slope = fitted_log_slope(summary)
        written = plot_result(summary, output_dir, args.formats, args.dpi)
        print(f"Loaded {len(rows)} rows from {csv_path}")
        for shots, values in summary.items():
            print(
                f"  S={shots:>6d}: mean={values['mean']:.6f}, "
                f"std={values['std']:.6f}, sem={values['sem']:.6f}"
            )
        print(f"Log-log slope of mean relative error: {slope:.6f}")
        for path in written:
            print(f"Plot saved to {path}")
        return

    shots_list = sorted(set(int(value) for value in args.shots))
    if not shots_list or any(value < 1 for value in shots_list):
        raise ValueError("Every --shots value must be a positive integer.")
    if args.trials < 1:
        raise ValueError("--trials must be positive.")

    fixed = prepare_fixed_experiment(args.n, args.state_seed)
    exact_gradient = np.asarray(fixed["exact_gradient"], dtype=float)
    exact_norm = float(fixed["exact_gradient_norm"])
    diagonal = np.asarray(fixed["diagonal"], dtype=float)
    metric = np.asarray(fixed["metric"], dtype=float)
    plus_probabilities = np.asarray(fixed["plus_probabilities"], dtype=float)
    minus_probabilities = np.asarray(fixed["minus_probabilities"], dtype=float)

    print(f"Finite-shot signed-branch sanity check: n={args.n}, parameters={exact_gradient.size}")
    print(f"Shot counts per signed branch: {shots_list}")
    print(f"Trials per shot count: {args.trials}")
    print(f"Objective value: {float(fixed['objective']):.9f}")
    print(f"Exact gradient norm: {exact_norm:.9e}")
    print(f"Exact branch-identity residual: {float(fixed['exact_branch_residual']):.3e}")
    print(f"Tree/Jacobian reference residual: {float(fixed['jacobian_reference_residual']):.3e}")
    print(f"Maximum state/tangent norm error: {max(float(fixed['state_norm_error']), float(fixed['max_tangent_norm_error'])):.3e}")
    print(f"Maximum |<psi|tangent>|: {float(fixed['max_tangent_overlap']):.3e}")
    print(f"Maximum branch-probability normalization error: {float(fixed['max_probability_norm_error']):.3e}")

    rows: List[Dict[str, float | int]] = []
    rng = np.random.default_rng(args.estimator_seed)

    for shots in shots_list:
        estimates = np.empty((args.trials, exact_gradient.size), dtype=float)
        errors = np.empty(args.trials, dtype=float)
        for trial in range(args.trials):
            estimate = estimate_gradient(
                diagonal,
                metric,
                plus_probabilities,
                minus_probabilities,
                shots,
                rng,
            )
            estimates[trial] = estimate
            errors[trial] = float(np.linalg.norm(estimate - exact_gradient) / exact_norm)
            rows.append(
                {
                    "shots": shots,
                    "shots_per_branch": shots,
                    "total_branch_shots_per_component": 2 * shots,
                    "trial": trial,
                    "rel_error": errors[trial],
                }
            )

        sample_std = float(np.std(errors, ddof=1)) if args.trials > 1 else 0.0
        sem = sample_std / math.sqrt(args.trials)
        relative_bias = float(np.linalg.norm(np.mean(estimates, axis=0) - exact_gradient) / exact_norm)
        print(
            f"S={shots:>6d}: mean rel. error={float(np.mean(errors)):.6f}, "
            f"std={sample_std:.6f}, sem={sem:.6f}, rel. mean-estimator bias={relative_bias:.6f}"
        )

    save_csv(rows, csv_path)
    summary = summarize_errors(rows)
    slope = fitted_log_slope(summary)
    written = plot_result(summary, output_dir, args.formats, args.dpi)

    print(f"Log-log slope of mean relative error: {slope:.6f} (ideal: -0.5)")
    print(f"Data saved to {csv_path}")
    for path in written:
        print(f"Plot saved to {path}")


if __name__ == "__main__":
    main()
