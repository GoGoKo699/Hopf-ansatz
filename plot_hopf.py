#!/usr/bin/env python3
"""
plot_hopf_geometric_clean.py

Paper-facing plots for the geometry-native Hopf optimizer data.

Expected input files:
    vqe_hopf_data_n{n}.csv
    met_hopf_data_n{n}.csv

Default outputs:
    hopf_geometric_summary_clean.pdf/.png
    hopf_geometric_n{detail_n}_convergence_clean.pdf/.png
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics as stats
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt

try:
    csv.field_size_limit(min(sys.maxsize, 2_147_483_647))
except Exception:
    pass

GEOM_MODES = [
    "Hopf-EGT-CG",
    "Hopf-Riemannian-BB",
    "Hopf-Riemannian-LBFGS",
]

MODE_LABEL = {
    "Hopf-EGT-CG": "EGT-CG",
    "Hopf-Riemannian-BB": "R-BB",
    "Hopf-Riemannian-LBFGS": "R-LBFGS",
}

MODE_COLOR = {
    "Hopf-EGT-CG": "#1f4e79",
    "Hopf-Riemannian-BB": "#4b5fa8",
    "Hopf-Riemannian-LBFGS": "#2a7f8f",
}

MODE_STYLE = {
    "Hopf-EGT-CG": "-",
    "Hopf-Riemannian-BB": "-.",
    "Hopf-Riemannian-LBFGS": "--",
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

TASK_ORDER = {
    "VQE": ["VQE-1", "VQE-2", "VQE-3"],
    "MET": ["MET-1", "MET-2", "MET-3"],
}


def parse_ns(text: str) -> List[int]:
    out: List[int] = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


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


def finite(x: float) -> bool:
    return math.isfinite(float(x))


def gap_value(row: Dict[str, Any]) -> float:
    task = str(row.get("task", "")).upper()
    if task == "VQE":
        for key in ("energy_gap", "cost", "energy"):
            v = sf(row.get(key))
            if finite(v):
                return v
    if task == "MET":
        for key in ("cfi_gap", "cost"):
            v = sf(row.get(key))
            if finite(v):
                return v
    return sf(row.get("cost"))


def read_rows(path: Path, task: str, n: int) -> List[Dict[str, Any]]:
    found = existing(path)
    if found is None:
        return []
    rows: List[Dict[str, Any]] = []
    with open_text(found) as f:
        reader = csv.DictReader(f)
        for raw in reader:
            mode = str(raw.get("mode", "")).strip()
            if mode == "Hopf-geodesic-CG":
                mode = "Hopf-EGT-CG"
            if mode not in GEOM_MODES:
                continue
            row = dict(raw)
            row["mode"] = mode
            row["task"] = str(row.get("task") or task).upper()
            row["n"] = si(row.get("n"), n)
            row["step"] = si(row.get("step"), -1)
            row["task_id"] = str(row.get("task_id") or "").strip()
            row["task_name"] = str(row.get("task_name") or row.get("task_id") or "").strip()
            row["run_seed"] = str(row.get("run_seed") or row.get("initial_seed") or "legacy_single_seed").strip()
            row["gap"] = gap_value(row)
            row["cost"] = sf(row.get("cost"))
            if row["task"] == task and row["n"] == n and row["step"] >= 0 and finite(row["gap"]):
                rows.append(row)
    return rows


def load_all(indir: Path, ns: Sequence[int]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for n in ns:
        rows.extend(read_rows(indir / f"vqe_hopf_data_n{n}.csv", "VQE", n))
        rows.extend(read_rows(indir / f"met_hopf_data_n{n}.csv", "MET", n))
    if not rows:
        raise FileNotFoundError(f"No Hopf data files found in {indir}")
    return rows


def group_rows(rows: Iterable[Dict[str, Any]], keys: Sequence[str]) -> Dict[Tuple[Any, ...], List[Dict[str, Any]]]:
    out: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = defaultdict(list)
    for row in rows:
        out[tuple(row.get(k) for k in keys)].append(row)
    for vals in out.values():
        vals.sort(key=lambda r: int(r.get("step", -1)))
    return dict(out)


def median(vals: Sequence[float]) -> float:
    clean = [float(v) for v in vals if finite(float(v))]
    return stats.median(clean) if clean else math.nan


def first_step_below(series: Sequence[Dict[str, Any]], threshold: float) -> float:
    for row in sorted(series, key=lambda r: int(r.get("step", -1))):
        if row["gap"] <= threshold:
            return float(row["step"])
    return math.nan


def clipped_gap(v: float, floor: float) -> float:
    if not finite(v):
        return math.nan
    return max(float(v), floor)


def mean_and_error(values: Sequence[float], error_band: str) -> Tuple[float, float]:
    vals = [float(v) for v in values if finite(float(v))]
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


def aggregate_step_band(series: Sequence[Dict[str, Any]], *, gap_floor: float, error_band: str) -> Tuple[List[int], List[float], List[float], List[float]]:
    by_step: Dict[int, List[float]] = defaultdict(list)
    for row in series:
        y = clipped_gap(float(row["gap"]), gap_floor)
        if finite(y):
            by_step[int(row["step"])].append(y)
    xs = sorted(by_step)
    means: List[float] = []
    lows: List[float] = []
    highs: List[float] = []
    for step in xs:
        mu, err = mean_and_error(by_step[step], error_band)
        means.append(clipped_gap(mu, gap_floor))
        lows.append(clipped_gap(mu - err, gap_floor))
        highs.append(clipped_gap(mu + err, gap_floor))
    return xs, means, lows, highs


def save_formats(fig: plt.Figure, outbase: Path, formats: Sequence[str], dpi: int) -> None:
    outbase.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(outbase.with_suffix(f".{fmt}"), bbox_inches="tight", dpi=dpi)


def plot_summary(rows: List[Dict[str, Any]], outdir: Path, formats: Sequence[str], dpi: int, gap_floor: float, threshold: float) -> None:
    groups = group_rows(rows, ["task", "n", "task_id", "run_seed", "mode"])
    final_by_task_mode: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    step_by_task_mode: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    for (task, _n, _tid, _run_seed, mode), series in groups.items():
        final_by_task_mode[(task, mode)].append(clipped_gap(series[-1]["gap"], gap_floor))
        st = first_step_below(series, threshold)
        if finite(st):
            step_by_task_mode[(task, mode)].append(st)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2), constrained_layout=True)
    for col, task in enumerate(["VQE", "MET"]):
        ax = axes[0][col]
        positions = list(range(1, len(GEOM_MODES) + 1))
        data = [final_by_task_mode.get((task, mode), []) for mode in GEOM_MODES]
        box = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=True)
        for patch, mode in zip(box.get("boxes", []), GEOM_MODES):
            patch.set_facecolor(MODE_COLOR[mode])
            patch.set_alpha(0.18)
            patch.set_edgecolor(MODE_COLOR[mode])
        for line in box.get("medians", []):
            line.set_linewidth(1.7)
        for i, mode in enumerate(GEOM_MODES, start=1):
            vals = data[i - 1]
            ax.scatter([i] * len(vals), vals, s=18, color=MODE_COLOR[mode], alpha=0.62, zorder=3)
        ax.axhline(threshold, color="0.5", linestyle="--", linewidth=1.0)
        ax.set_yscale("log")
        ax.set_ylim(bottom=gap_floor)
        ax.set_title(f"{task}: final gap")
        ax.set_xticks(positions)
        ax.set_xticklabels([MODE_LABEL[m] for m in GEOM_MODES], rotation=25, ha="right")
        ax.set_ylabel("gap")
        ax.grid(True, which="both", alpha=0.18)

        ax = axes[1][col]
        data2 = [step_by_task_mode.get((task, mode), []) for mode in GEOM_MODES]
        box2 = ax.boxplot(data2, positions=positions, widths=0.55, patch_artist=True, showfliers=True)
        for patch, mode in zip(box2.get("boxes", []), GEOM_MODES):
            patch.set_facecolor(MODE_COLOR[mode])
            patch.set_alpha(0.18)
            patch.set_edgecolor(MODE_COLOR[mode])
        for line in box2.get("medians", []):
            line.set_linewidth(1.7)
        for i, mode in enumerate(GEOM_MODES, start=1):
            vals = data2[i - 1]
            ax.scatter([i] * len(vals), vals, s=18, color=MODE_COLOR[mode], alpha=0.62, zorder=3)
        ax.set_title(f"{task}: first step below {threshold:g}")
        ax.set_xticks(positions)
        ax.set_xticklabels([MODE_LABEL[m] for m in GEOM_MODES], rotation=25, ha="right")
        ax.set_ylabel("step")
        ax.grid(True, alpha=0.22)
    save_formats(fig, outdir / "hopf_geometric_summary_clean", formats, dpi)
    plt.close(fig)


def plot_detail(rows: List[Dict[str, Any]], outdir: Path, formats: Sequence[str], dpi: int, detail_n: int, gap_floor: float, threshold: float, error_band: str, band_alpha: float) -> None:
    rows_n = [r for r in rows if int(r.get("n", -1)) == detail_n]
    if not rows_n:
        return
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 6.9), sharex=True, constrained_layout=True)
    for row_idx, task in enumerate(["VQE", "MET"]):
        task_rows = [r for r in rows_n if r.get("task") == task]
        task_ids = [tid for tid in TASK_ORDER[task] if any(r.get("task_id") == tid for r in task_rows)]
        if len(task_ids) < 3:
            task_ids = sorted({str(r.get("task_id")) for r in task_rows})[:3]
        for col, task_id in enumerate(task_ids[:3]):
            ax = axes[row_idx][col]
            for mode in GEOM_MODES:
                series = [r for r in task_rows if r.get("task_id") == task_id and r.get("mode") == mode]
                series.sort(key=lambda r: int(r.get("step", -1)))
                if not series:
                    continue
                xs, ys, lo, hi = aggregate_step_band(series, gap_floor=gap_floor, error_band=error_band)
                if error_band != "none" and len(xs) >= 2:
                    ax.fill_between(xs, lo, hi, color=MODE_COLOR[mode], alpha=band_alpha, linewidth=0)
                ax.plot(xs, ys, color=MODE_COLOR[mode], linestyle=MODE_STYLE[mode], linewidth=1.75, label=MODE_LABEL[mode])
            title_name = TASK_LABEL.get(task_id)
            if title_name is None:
                any_row = next((r for r in task_rows if r.get("task_id") == task_id), {})
                title_name = TASK_LABEL.get(str(any_row.get("task_name")), task_id)
            ax.axhline(threshold, color="0.5", linestyle="--", linewidth=1.0)
            ax.set_yscale("log")
            ax.set_ylim(bottom=gap_floor)
            ax.set_title(f"{task}: {title_name}")
            ax.grid(True, which="both", alpha=0.18)
            if row_idx == 1:
                ax.set_xlabel("step")
            if col == 0:
                ax.set_ylabel("gap")
            if row_idx == 0 and col == 0:
                ax.legend(frameon=False, fontsize=8)
    save_formats(fig, outdir / f"hopf_geometric_n{detail_n}_convergence_clean", formats, dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=Path, default=Path("."))
    parser.add_argument("--outdir", type=Path, default=Path("figures_hopf_clean"))
    parser.add_argument("--ns", type=parse_ns, default=parse_ns("6-10"))
    parser.add_argument("--detail-n", type=int, default=10)
    parser.add_argument("--formats", nargs="+", default=["pdf", "png"])
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--gap-floor", type=float, default=1e-16)
    parser.add_argument("--threshold", type=float, default=1e-8)
    parser.add_argument("--error-band", choices=["sem", "std", "none"], default="sem", help="Shaded band around mean convergence traces across run_seed values.")
    parser.add_argument("--band-alpha", type=float, default=0.18, help="Transparency for shaded error bands.")
    args = parser.parse_args()

    rows = load_all(args.indir, args.ns)
    plot_summary(rows, args.outdir, args.formats, args.dpi, args.gap_floor, args.threshold)
    plot_detail(rows, args.outdir, args.formats, args.dpi, args.detail_n, args.gap_floor, args.threshold, args.error_band, args.band_alpha)
    print(f"Wrote Hopf figures to {args.outdir}")


if __name__ == "__main__":
    main()
