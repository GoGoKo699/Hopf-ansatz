#!/usr/bin/env python3
"""
plot_all_github_clean.py

Github-facing plain-scale cost plots for one selected n.

Expected input files:
    vqe_hopf_data_n{n}.csv, met_hopf_data_n{n}.csv
    vqe_adam_data_n{n}.csv, met_adam_data_n{n}.csv

Default outputs:
    all_vqe_n{n}_clean.pdf/.png
    all_met_n{n}_clean.pdf/.png
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics as stats
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt

try:
    csv.field_size_limit(min(sys.maxsize, 2_147_483_647))
except Exception:
    pass

TASK_ORDER = {
    "VQE": ["VQE-1", "VQE-2", "VQE-3"],
    "MET": ["MET-1", "MET-2", "MET-3"],
}

TASK_LABEL = {
    "VQE-1": "Parent Hamiltonian",
    "VQE-2": "Hamming spectrum",
    "VQE-3": "Small-gap spectrum",
    "MET-1": "Single-target Fisher",
    "MET-2": "QFI superposition",
    "MET-3": "Balanced Fisher",
    "random_parent": "Parent Hamiltonian",
    "scrambled_hamming_spectrum": "Hamming spectrum",
    "small_gap_scrambled_spectrum": "Small-gap spectrum",
    "single_target_fixed_readout_cfi": "Single-target Fisher",
    "qfi_extremal_superposition": "QFI superposition",
    "two_target_balanced_fisher": "Balanced Fisher",
}

MODE_ORDER = [
    "Hopf-EGT-CG",
    "Hopf-Riemannian-BB",
    "Hopf-Riemannian-LBFGS",
    "Hopf-Adam",
    "Mottonen-ideal-PS-Adam",
    "Hopf-geodesic-CG",
    "Hopf-Ritz3",
    "Hopf-Subspace3",
]

MODE_LABEL = {
    "Hopf-EGT-CG": "Hopf EGT-CG",
    "Hopf-Riemannian-BB": "Hopf R-BB",
    "Hopf-Riemannian-LBFGS": "Hopf R-LBFGS",
    "Hopf-Adam": "Hopf Adam",
    "Mottonen-ideal-PS-Adam": "Möttönen Adam",
    "Hopf-geodesic-CG": "Hopf geodesic-CG",
    "Hopf-Ritz3": "Hopf Ritz3",
    "Hopf-Subspace3": "Hopf Subspace3",
}

MODE_COLOR = {
    "Hopf-EGT-CG": "#1f4e79",
    "Hopf-Riemannian-BB": "#4b5fa8",
    "Hopf-Riemannian-LBFGS": "#2a7f8f",
    "Hopf-Adam": "#eb6426",
    "Mottonen-ideal-PS-Adam": "#f01d1d",
    "Hopf-geodesic-CG": "#376092",
    "Hopf-Ritz3": "#756bb1",
    "Hopf-Subspace3": "#2ca25f",
}

MODE_STYLE = {
    "Hopf-EGT-CG": "-",
    "Hopf-Riemannian-BB": "-.",
    "Hopf-Riemannian-LBFGS": "--",
    "Hopf-Adam": ":",
    "Mottonen-ideal-PS-Adam": (0, (5, 2)),
    "Hopf-geodesic-CG": "-",
    "Hopf-Ritz3": "-.",
    "Hopf-Subspace3": "--",
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def existing(path: Path) -> Path | None:
    if path.exists():
        return path
    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gz
    return None


def sf(x: Any, default: float = math.nan) -> float:
    try:
        s = str(x).strip()
        if s == "" or s.lower() in {"nan", "none", "null"}:
            return default
        return float(s)
    except Exception:
        return default


def si(x: Any, default: int = -1) -> int:
    try:
        return int(float(str(x).strip()))
    except Exception:
        return default


def read_csv(path: Path, task: str, n: int) -> List[Dict[str, Any]]:
    found = existing(path)
    if found is None:
        return []
    out: List[Dict[str, Any]] = []
    with open_text(found) as f:
        reader = csv.DictReader(f)
        for raw in reader:
            row = dict(raw)
            row["source_file"] = found.name
            row["task"] = str(row.get("task") or task).upper()
            row["n"] = si(row.get("n"), n)
            row["step"] = si(row.get("step"), -1)
            row["task_id"] = str(row.get("task_id") or "").strip()
            row["task_name"] = str(row.get("task_name") or row.get("task_id") or "").strip()
            row["mode"] = str(row.get("mode") or "").strip()
            row["run_seed"] = str(row.get("run_seed") or row.get("initial_seed") or "legacy_single_seed").strip()
            row["cost"] = sf(row.get("cost"))
            if row["task"] == task and row["n"] == n and row["step"] >= 0 and math.isfinite(row["cost"]):
                out.append(row)
    return out


def load_rows(task: str, n: int, hopf_dir: Path, adam_dir: Path) -> List[Dict[str, Any]]:
    prefix = task.lower()
    rows: List[Dict[str, Any]] = []
    rows.extend(read_csv(hopf_dir / f"{prefix}_hopf_data_n{n}.csv", task, n))
    rows.extend(read_csv(adam_dir / f"{prefix}_adam_data_n{n}.csv", task, n))
    if not rows:
        raise FileNotFoundError(f"No {prefix} rows found for n={n}")
    seen = set()
    deduped: List[Dict[str, Any]] = []
    order = {m: i for i, m in enumerate(MODE_ORDER)}
    for row in sorted(rows, key=lambda r: (r["task_id"], order.get(r["mode"], 999), r["source_file"], r["step"])):
        key = (row["task"], row["n"], row["task_id"], row["run_seed"], row["mode"], row["step"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def rows_for(rows: Iterable[Dict[str, Any]], task_id: str, mode: str) -> List[Dict[str, Any]]:
    out = [r for r in rows if r.get("task_id") == task_id and r.get("mode") == mode]
    out.sort(key=lambda r: (r["step"], str(r.get("run_seed", ""))))
    return out


def mean_and_error(values: Sequence[float], error_band: str) -> Tuple[float, float]:
    vals = [float(v) for v in values if math.isfinite(float(v))]
    if not vals:
        return math.nan, math.nan
    mu = sum(vals) / len(vals)
    if error_band == "none" or len(vals) < 2:
        return mu, 0.0
    sd = stats.stdev(vals)
    if error_band == "sem":
        return mu, sd / math.sqrt(len(vals))
    if error_band == "std":
        return mu, sd
    raise ValueError("error_band must be one of: none, sem, std")


def aggregate_step_band(series: Sequence[Dict[str, Any]], *, error_band: str) -> Tuple[List[int], List[float], List[float], List[float]]:
    by_step: Dict[int, List[float]] = {}
    for row in series:
        y = float(row["cost"])
        if math.isfinite(y):
            by_step.setdefault(int(row["step"]), []).append(y)
    xs = sorted(by_step)
    means: List[float] = []
    lows: List[float] = []
    highs: List[float] = []
    for step in xs:
        mu, err = mean_and_error(by_step[step], error_band)
        means.append(mu)
        lows.append(mu - err)
        highs.append(mu + err)
    return xs, means, lows, highs


def task_title(rows: Sequence[Dict[str, Any]], task_id: str) -> str:
    if task_id in TASK_LABEL:
        return TASK_LABEL[task_id]
    name = next((str(r.get("task_name")) for r in rows if r.get("task_id") == task_id), task_id)
    return TASK_LABEL.get(name, name.replace("_", " "))


def save_formats(fig: plt.Figure, outbase: Path, formats: Sequence[str], dpi: int) -> None:
    outbase.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(outbase.with_suffix(f".{fmt}"), bbox_inches="tight", dpi=dpi)


def plot_task(task: str, n: int, rows: List[Dict[str, Any]], outdir: Path, formats: Sequence[str], dpi: int, error_band: str, band_alpha: float) -> None:
    task_ids = [tid for tid in TASK_ORDER[task] if any(r.get("task_id") == tid for r in rows)]
    if len(task_ids) < 3:
        task_ids = sorted({str(r.get("task_id")) for r in rows})[:3]
    modes_present = []
    for mode in MODE_ORDER:
        if any(r.get("mode") == mode for r in rows):
            modes_present.append(mode)
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.0), sharex=True, constrained_layout=True)
    for ax, task_id in zip(axes, task_ids[:3]):
        for mode in modes_present:
            series = rows_for(rows, task_id, mode)
            if not series:
                continue
            xs, ys, lo, hi = aggregate_step_band(series, error_band=error_band)
            color = MODE_COLOR.get(mode, "0.25")
            if error_band != "none" and len(xs) >= 2:
                ax.fill_between(xs, lo, hi, color=color, alpha=band_alpha, linewidth=0)
            ax.plot(xs, ys, color=color, linestyle=MODE_STYLE.get(mode, "-"), linewidth=1.65, label=MODE_LABEL.get(mode, mode))
        ax.axhline(0.0, color="0.55", linestyle="--", linewidth=1.0)
        ax.set_title(task_title(rows, task_id))
        ax.set_xlabel("step")
        ax.set_ylabel("cost")
        ax.grid(True, alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{task} cost traces, n={n}")
    save_formats(fig, outdir / f"all_{task.lower()}_n{n}_clean", formats, dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hopf-dir", type=Path, default=Path("."))
    parser.add_argument("--adam-dir", type=Path, default=Path("."))
    parser.add_argument("--outdir", type=Path, default=Path("figures_all_clean"))
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--formats", nargs="+", default=["pdf", "png"])
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--error-band", choices=["sem", "std", "none"], default="sem", help="Shaded band around each mean trace across run_seed values.")
    parser.add_argument("--band-alpha", type=float, default=0.18, help="Transparency for shaded error bands.")
    args = parser.parse_args()

    for task in ["VQE", "MET"]:
        rows = load_rows(task, args.n, args.hopf_dir, args.adam_dir)
        plot_task(task, args.n, rows, args.outdir, args.formats, args.dpi, args.error_band, args.band_alpha)
    print(f"Wrote combined figures to {args.outdir}")


if __name__ == "__main__":
    main()
