#!/usr/bin/env python3
"""
diagnose_hopf.py

Streaming diagnostics for synthetic Hopf cost-and-gradient stress-test CSV files produced by hopf_data.py.

Default expected input files in --indir:
    vqe_hopf_data_n6.csv,  met_hopf_data_n6.csv,
    ...
    vqe_hopf_data_n10.csv, met_hopf_data_n10.csv

Default output:
    hopf_data_diagnostics.txt

The script intentionally avoids pandas and does not parse the full theta/grad vectors by default,
because the vector columns can dominate file size. It expects only the three Hopf
cost-and-gradient optimizer modes: EGT-CG, Riemannian L-BFGS, and Riemannian
Barzilai--Borwein. It computes convergence, ranking, completeness, and pathology
summaries from scalar columns while optionally checking sampled vector lengths by
counting semicolon-separated entries.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import math
import os
from pathlib import Path
import statistics as stats
import sys
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Any

# Large theta/grad CSV fields can exceed Python's default CSV field-size limit.
try:
    csv.field_size_limit(min(sys.maxsize, 2_147_483_647))
except Exception:
    pass

VQE_EXPECTED_MODES = [
    "Hopf-EGT-CG",
    "Hopf-Riemannian-LBFGS",
    "Hopf-Riemannian-BB",
]

MET_EXPECTED_MODES = [
    "Hopf-EGT-CG",
    "Hopf-Riemannian-LBFGS",
    "Hopf-Riemannian-BB",
]

THRESHOLDS = [1e-2, 1e-4, 1e-6, 1e-8]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def safe_float(x: Any, default: float = math.nan) -> float:
    if x is None:
        return default
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in {"none", "null"}:
        return default
    try:
        return float(s)
    except Exception:
        return default


def safe_int(x: Any, default: Optional[int] = None) -> Optional[int]:
    if x is None:
        return default
    if isinstance(x, int):
        return x
    s = str(x).strip()
    if s == "":
        return default
    try:
        return int(float(s))
    except Exception:
        return default


def finite(x: float) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def fmt_float(x: float, digits: int = 6) -> str:
    if not finite(x):
        return "nan"
    ax = abs(x)
    if ax == 0:
        return "0"
    if ax < 1e-4 or ax >= 1e5:
        return f"{x:.{digits}e}"
    return f"{x:.{digits}g}"


def fmt_int(x: Optional[int]) -> str:
    return "" if x is None else str(x)


def mean(values: Sequence[float]) -> float:
    vals = [v for v in values if finite(v)]
    return sum(vals) / len(vals) if vals else math.nan


def median(values: Sequence[float]) -> float:
    vals = [v for v in values if finite(v)]
    return stats.median(vals) if vals else math.nan


def stdev(values: Sequence[float]) -> float:
    vals = [v for v in values if finite(v)]
    return stats.stdev(vals) if len(vals) >= 2 else 0.0 if len(vals) == 1 else math.nan


def quantile(values: Sequence[float], q: float) -> float:
    vals = sorted(v for v in values if finite(v))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return (1 - frac) * vals[lo] + frac * vals[hi]


def vector_length_from_string(s: str) -> int:
    s = (s or "").strip()
    if not s:
        return 0
    return s.count(";") + 1


def metric_value(row: Dict[str, str]) -> float:
    """Lower is better for both tasks."""
    task = (row.get("task") or "").strip().upper()
    if task == "VQE":
        for col in ("energy_gap", "cost", "energy"):
            v = safe_float(row.get(col))
            if finite(v):
                return v
        return math.nan
    if task == "MET":
        for col in ("cfi_gap", "cost"):
            v = safe_float(row.get(col))
            if finite(v):
                return v
        return math.nan
    # Fallback: lower cost is better.
    return safe_float(row.get("cost"))


def task_score_label(task: str) -> str:
    if task.upper() == "VQE":
        return "gap = energy_gap (lower is better)"
    if task.upper() == "MET":
        return "gap = cfi_gap or cost (lower is better)"
    return "cost/gap (lower is better)"


@dataclass
class RowPoint:
    step: int
    metric: float
    cost: float
    grad_norm: float
    state_grad_norm: float
    target_overlap: float
    decoded_probability: float
    wall_time_sec: float
    last_step_angle: float
    last_line_evals: float
    last_subspace_dim: float
    raw: Dict[str, str] = field(default_factory=dict)


@dataclass
class GroupStats:
    key: Tuple[str, int, str, str, str, str]  # task, n, task_id, task_name, run_seed, mode
    rows: int = 0
    steps: set = field(default_factory=set)
    points: List[RowPoint] = field(default_factory=list)
    nonfinite_metric_rows: int = 0
    nonfinite_grad_rows: int = 0
    max_state_norm_error: float = 0.0
    max_grad_norm: float = math.nan
    max_state_grad_norm: float = math.nan
    vector_mismatch_count: int = 0
    vector_checks: int = 0
    coordinate_type: str = ""
    num_parameters: Optional[int] = None
    dimension: Optional[int] = None
    problem_seed: str = ""
    scramble_depth: str = ""
    scramble_seed: str = ""
    seed_index: str = ""
    run_seed: str = ""
    source_file: str = ""

    def add_point(self, p: RowPoint) -> None:
        self.rows += 1
        self.steps.add(p.step)
        self.points.append(p)
        if not finite(p.metric):
            self.nonfinite_metric_rows += 1
        if not finite(p.grad_norm):
            self.nonfinite_grad_rows += 1
        sne = safe_float(p.raw.get("state_norm_error"), 0.0)
        if finite(sne):
            self.max_state_norm_error = max(self.max_state_norm_error, abs(sne))
        if finite(p.grad_norm):
            self.max_grad_norm = p.grad_norm if not finite(self.max_grad_norm) else max(self.max_grad_norm, p.grad_norm)
        if finite(p.state_grad_norm):
            self.max_state_grad_norm = p.state_grad_norm if not finite(self.max_state_grad_norm) else max(self.max_state_grad_norm, p.state_grad_norm)

    @property
    def task(self) -> str:
        return self.key[0]

    @property
    def n(self) -> int:
        return self.key[1]

    @property
    def task_id(self) -> str:
        return self.key[2]

    @property
    def task_name(self) -> str:
        return self.key[3]

    @property
    def run_seed_value(self) -> str:
        return self.key[4]

    @property
    def mode(self) -> str:
        return self.key[5]

    def sorted_points(self) -> List[RowPoint]:
        return sorted(self.points, key=lambda p: p.step)

    def initial(self) -> Optional[RowPoint]:
        pts = self.sorted_points()
        return pts[0] if pts else None

    def final(self) -> Optional[RowPoint]:
        pts = self.sorted_points()
        return pts[-1] if pts else None

    def best(self) -> Optional[RowPoint]:
        pts = [p for p in self.points if finite(p.metric)]
        return min(pts, key=lambda p: p.metric) if pts else None

    def metric_series(self) -> List[float]:
        return [p.metric for p in self.sorted_points()]

    def final_metric(self) -> float:
        p = self.final()
        return p.metric if p else math.nan

    def initial_metric(self) -> float:
        p = self.initial()
        return p.metric if p else math.nan

    def best_metric(self) -> float:
        p = self.best()
        return p.metric if p else math.nan

    def best_step(self) -> Optional[int]:
        p = self.best()
        return p.step if p else None

    def final_target_overlap(self) -> float:
        p = self.final()
        return p.target_overlap if p else math.nan

    def final_decoded_probability(self) -> float:
        p = self.final()
        return p.decoded_probability if p else math.nan

    def final_grad_norm(self) -> float:
        p = self.final()
        return p.grad_norm if p else math.nan

    def final_state_grad_norm(self) -> float:
        p = self.final()
        return p.state_grad_norm if p else math.nan

    def total_line_evals(self) -> float:
        vals = [p.last_line_evals for p in self.points if finite(p.last_line_evals)]
        return sum(vals) if vals else 0.0

    def mean_line_evals(self) -> float:
        vals = [p.last_line_evals for p in self.points if finite(p.last_line_evals)]
        return mean(vals) if vals else 0.0

    def mean_step_angle(self) -> float:
        vals = [p.last_step_angle for p in self.points if finite(p.last_step_angle)]
        return mean(vals) if vals else math.nan

    def monotone_fraction(self) -> float:
        vals = [v for v in self.metric_series() if finite(v)]
        if len(vals) <= 1:
            return math.nan
        nonincreasing = sum(1 for a, b in zip(vals[:-1], vals[1:]) if b <= a + 1e-15)
        return nonincreasing / (len(vals) - 1)

    def improvement(self) -> float:
        a = self.initial_metric()
        b = self.final_metric()
        if not finite(a) or not finite(b):
            return math.nan
        return a - b

    def relative_improvement(self) -> float:
        a = self.initial_metric()
        imp = self.improvement()
        if not finite(a) or abs(a) < 1e-15 or not finite(imp):
            return math.nan
        return imp / abs(a)

    def first_step_below(self, threshold: float) -> Optional[int]:
        for p in self.sorted_points():
            if finite(p.metric) and p.metric <= threshold:
                return p.step
        return None

    def terminal_slope(self, window: int = 20) -> float:
        pts = [p for p in self.sorted_points() if finite(p.metric)]
        if len(pts) < 3:
            return math.nan
        vals = [max(p.metric, 1e-300) for p in pts[-window:]]
        xs = list(range(len(vals)))
        if len(vals) < 2:
            return math.nan
        ys = [math.log10(v) for v in vals]
        xbar = mean(xs)
        ybar = mean(ys)
        den = sum((x - xbar) ** 2 for x in xs)
        if den == 0:
            return math.nan
        return sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / den


@dataclass
class FileStats:
    path: Path
    task_kind: str
    n: int
    rows: int = 0
    columns: List[str] = field(default_factory=list)
    groups: Dict[Tuple[str, int, str, str, str, str], GroupStats] = field(default_factory=dict)
    parse_errors: int = 0
    vector_mismatch_count: int = 0
    vector_checks: int = 0


def infer_task_and_n(path: Path) -> Tuple[str, Optional[int]]:
    name = path.name.lower()
    task = "VQE" if name.startswith("vqe") else "MET" if name.startswith("met") else "UNKNOWN"
    nval = None
    # Expected suffix *_n{n}.csv or *_n{n}.csv.gz
    stem = path.name
    for token in stem.replace(".", "_").split("_"):
        if token.startswith("n") and token[1:].isdigit():
            nval = int(token[1:])
    return task, nval


def find_input_files(indir: Path, ns: Sequence[int]) -> List[Path]:
    files: List[Path] = []
    for n in ns:
        for base in (f"vqe_hopf_data_n{n}.csv", f"met_hopf_data_n{n}.csv"):
            p = indir / base
            pgz = indir / (base + ".gz")
            if p.exists():
                files.append(p)
            elif pgz.exists():
                files.append(pgz)
            else:
                # Missing files are reported later from expected list.
                pass
    return files


def expected_file_list(indir: Path, ns: Sequence[int]) -> List[Path]:
    out = []
    for n in ns:
        out.append(indir / f"vqe_hopf_data_n{n}.csv")
        out.append(indir / f"met_hopf_data_n{n}.csv")
    return out


def row_to_group_key(row: Dict[str, str], fallback_task: str, fallback_n: int) -> Tuple[str, int, str, str, str, str]:
    task = (row.get("task") or fallback_task or "UNKNOWN").strip()
    n = safe_int(row.get("n"), fallback_n)
    if n is None:
        n = fallback_n
    task_id = str(row.get("task_id") or row.get("task_name") or "").strip()
    task_name = str(row.get("task_name") or task_id or "UNKNOWN_TASK").strip()
    run_seed = str(row.get("run_seed") or row.get("initial_seed") or "legacy_single_seed").strip()
    mode = str(row.get("mode") or "UNKNOWN_MODE").strip()
    return (task, int(n), task_id, task_name, run_seed, mode)


def process_file(path: Path, expected_steps: int, vector_check: str) -> FileStats:
    fallback_task, n_guess = infer_task_and_n(path)
    if n_guess is None:
        n_guess = -1
    fs = FileStats(path=path, task_kind=fallback_task, n=n_guess)

    with open_text(path) as f:
        reader = csv.DictReader(f)
        fs.columns = reader.fieldnames or []
        for row in reader:
            fs.rows += 1
            try:
                key = row_to_group_key(row, fallback_task, n_guess)
                if key not in fs.groups:
                    gs = GroupStats(key=key, source_file=str(path))
                    gs.coordinate_type = str(row.get("coordinate_type") or "")
                    gs.num_parameters = safe_int(row.get("num_parameters"), None)
                    gs.dimension = safe_int(row.get("dimension"), None)
                    gs.problem_seed = str(row.get("problem_seed") or "")
                    gs.scramble_depth = str(row.get("scramble_depth") or "")
                    gs.scramble_seed = str(row.get("scramble_seed") or "")
                    gs.seed_index = str(row.get("seed_index") or "")
                    gs.run_seed = str(row.get("run_seed") or row.get("initial_seed") or "legacy_single_seed")
                    fs.groups[key] = gs
                gs = fs.groups[key]
                step = safe_int(row.get("step"), -1)
                if step is None:
                    step = -1
                metric = metric_value(row)
                p = RowPoint(
                    step=int(step),
                    metric=metric,
                    cost=safe_float(row.get("cost")),
                    grad_norm=safe_float(row.get("grad_norm")),
                    state_grad_norm=safe_float(row.get("state_grad_norm")),
                    target_overlap=safe_float(row.get("target_overlap")),
                    decoded_probability=safe_float(row.get("decoded_probability")),
                    wall_time_sec=safe_float(row.get("wall_time_sec")),
                    last_step_angle=safe_float(row.get("last_step_angle")),
                    last_line_evals=safe_float(row.get("last_line_evals")),
                    last_subspace_dim=safe_float(row.get("last_subspace_dim")),
                    raw=row,
                )
                gs.add_point(p)

                # Optional vector-length checking. This does not parse floats, only counts entries.
                do_vec = False
                if vector_check == "all":
                    do_vec = True
                elif vector_check == "sample" and step in {0, expected_steps}:
                    do_vec = True
                if do_vec:
                    expected_p = safe_int(row.get("num_parameters"), None)
                    if expected_p is not None:
                        theta_len = vector_length_from_string(row.get("theta") or "")
                        grad_len = vector_length_from_string(row.get("grad") or "")
                        gs.vector_checks += 1
                        fs.vector_checks += 1
                        if theta_len != expected_p or grad_len != expected_p:
                            gs.vector_mismatch_count += 1
                            fs.vector_mismatch_count += 1
            except Exception:
                fs.parse_errors += 1
    return fs


def completeness_summary(fs: FileStats, expected_steps: int, expected_num_seeds: int) -> Dict[str, Any]:
    expected_modes = VQE_EXPECTED_MODES if fs.task_kind == "VQE" else MET_EXPECTED_MODES if fs.task_kind == "MET" else []
    task_ids = sorted({g.task_id for g in fs.groups.values()})
    modes = sorted({g.mode for g in fs.groups.values()})
    missing_modes = [m for m in expected_modes if m not in modes]
    expected_rows = len(task_ids) * len(expected_modes) * int(expected_num_seeds) * (expected_steps + 1) if expected_modes else None
    missing_groups = []
    seed_count_problems = []
    for tid in task_ids:
        task_names = sorted({g.task_name for g in fs.groups.values() if g.task_id == tid})
        tname = task_names[0] if task_names else tid
        observed_seeds = sorted({g.run_seed_value for g in fs.groups.values() if g.task_id == tid})
        if len(observed_seeds) != int(expected_num_seeds):
            seed_count_problems.append((tid, tname, len(observed_seeds), int(expected_num_seeds), observed_seeds[:10]))
        for seed in observed_seeds:
            for m in expected_modes:
                if not any(g.task_id == tid and g.run_seed_value == seed and g.mode == m for g in fs.groups.values()):
                    missing_groups.append((tid, tname, seed, m))
    incomplete_step_groups = []
    expected_step_set = set(range(0, expected_steps + 1))
    for g in fs.groups.values():
        missing_steps = sorted(expected_step_set - g.steps)
        extra_steps = sorted(s for s in g.steps if s not in expected_step_set)
        if missing_steps or extra_steps:
            incomplete_step_groups.append((g, missing_steps[:10], len(missing_steps), extra_steps[:10], len(extra_steps)))
    return {
        "task_ids": task_ids,
        "modes": modes,
        "missing_modes": missing_modes,
        "expected_rows": expected_rows,
        "missing_groups": missing_groups,
        "seed_count_problems": seed_count_problems,
        "incomplete_step_groups": incomplete_step_groups,
    }


def rank_groups_for_task(groups: Sequence[GroupStats]) -> List[GroupStats]:
    return sorted(groups, key=lambda g: (math.inf if not finite(g.final_metric()) else g.final_metric(), g.mode))


def aggregate_by_mode(groups: Sequence[GroupStats]) -> Dict[str, Dict[str, Any]]:
    by_mode: Dict[str, List[GroupStats]] = defaultdict(list)
    for g in groups:
        by_mode[g.mode].append(g)
    out: Dict[str, Dict[str, Any]] = {}
    for mode, gs in by_mode.items():
        finals = [g.final_metric() for g in gs]
        bests = [g.best_metric() for g in gs]
        imps = [g.improvement() for g in gs]
        rel_imps = [g.relative_improvement() for g in gs]
        overlaps = [g.final_target_overlap() for g in gs]
        probs = [g.final_decoded_probability() for g in gs]
        grad_final = [g.final_grad_norm() for g in gs]
        out[mode] = {
            "count": len(gs),
            "mean_final": mean(finals),
            "median_final": median(finals),
            "std_final": stdev(finals),
            "p25_final": quantile(finals, 0.25),
            "p75_final": quantile(finals, 0.75),
            "mean_best": mean(bests),
            "mean_improvement": mean(imps),
            "mean_relative_improvement": mean(rel_imps),
            "mean_overlap": mean(overlaps),
            "mean_decoded_probability": mean(probs),
            "mean_final_grad_norm": mean(grad_final),
        }
    return out


def win_counts(groups: Sequence[GroupStats], by: str = "final") -> Counter:
    """Count wins by mode per (task, n, task_id, run_seed)."""
    bucket: Dict[Tuple[str, int, str, str], List[GroupStats]] = defaultdict(list)
    for g in groups:
        bucket[(g.task, g.n, g.task_id, g.run_seed_value)].append(g)
    c = Counter()
    for _key, gs in bucket.items():
        if by == "best":
            sorted_g = sorted(gs, key=lambda g: (math.inf if not finite(g.best_metric()) else g.best_metric(), g.mode))
        else:
            sorted_g = sorted(gs, key=lambda g: (math.inf if not finite(g.final_metric()) else g.final_metric(), g.mode))
        if sorted_g:
            c[sorted_g[0].mode] += 1
    return c


def mode_rank_distribution(groups: Sequence[GroupStats], by: str = "final") -> Dict[str, List[int]]:
    bucket: Dict[Tuple[str, int, str, str], List[GroupStats]] = defaultdict(list)
    for g in groups:
        bucket[(g.task, g.n, g.task_id, g.run_seed_value)].append(g)
    out: Dict[str, List[int]] = defaultdict(list)
    for _key, gs in bucket.items():
        if by == "best":
            sorted_g = sorted(gs, key=lambda g: (math.inf if not finite(g.best_metric()) else g.best_metric(), g.mode))
        else:
            sorted_g = sorted(gs, key=lambda g: (math.inf if not finite(g.final_metric()) else g.final_metric(), g.mode))
        for rank, g in enumerate(sorted_g, start=1):
            out[g.mode].append(rank)
    return out


def line(width: int = 120) -> str:
    return "-" * width


def table(rows: List[Sequence[Any]], headers: Sequence[str], max_width: int = 120) -> str:
    """Simple text table."""
    str_rows = [[str(x) for x in r] for r in rows]
    cols = len(headers)
    widths = [len(str(h)) for h in headers]
    for r in str_rows:
        for i in range(cols):
            widths[i] = max(widths[i], len(r[i]) if i < len(r) else 0)
    # Keep wide tables readable by capping some long columns.
    widths = [min(w, 38) for w in widths]

    def clip(s: str, w: int) -> str:
        return s if len(s) <= w else s[: max(0, w - 1)] + "…"

    def fmt_row(r: Sequence[str]) -> str:
        return "  ".join(clip(str(r[i]) if i < len(r) else "", widths[i]).ljust(widths[i]) for i in range(cols))

    out = [fmt_row([str(h) for h in headers]), fmt_row(["-" * w for w in widths])]
    out.extend(fmt_row(r) for r in str_rows)
    return "\n".join(out)


def build_report(
    file_stats: Sequence[FileStats],
    expected_files: Sequence[Path],
    args: argparse.Namespace,
) -> str:
    all_groups: List[GroupStats] = []
    for fs in file_stats:
        all_groups.extend(fs.groups.values())

    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: List[str] = []
    lines.append("HOPF SYNTHETIC STRESS-TEST CSV DIAGNOSTICS")
    lines.append(line())
    lines.append(f"Generated: {now}")
    lines.append(f"Input directory: {Path(args.indir).resolve()}")
    lines.append(f"Requested n values: {', '.join(map(str, args.ns))}")
    lines.append(f"Expected steps: 0..{args.steps} ({args.steps + 1} records per seed/mode run)")
    lines.append(f"Expected initial-state seeds per task: {args.num_seeds}")
    lines.append(f"Vector check mode: {args.vector_check}")
    lines.append(f"Metric convention: lower is better. VQE uses energy_gap when available; MET uses cfi_gap when available.")
    lines.append("")

    # File inventory.
    lines.append("1. FILE INVENTORY AND COMPLETENESS")
    lines.append(line())
    existing = {fs.path.resolve() for fs in file_stats}
    inv_rows = []
    for p in expected_files:
        candidates = [p, Path(str(p) + ".gz")]
        found = next((c for c in candidates if c.resolve() in existing), None)
        if found is None:
            inv_rows.append([p.name, "MISSING", "", "", "", ""])
        else:
            fs = next(fs0 for fs0 in file_stats if fs0.path.resolve() == found.resolve())
            comp = completeness_summary(fs, args.steps, args.num_seeds)
            exp = comp["expected_rows"]
            exp_s = "" if exp is None else str(exp)
            inv_rows.append([
                found.name,
                "OK",
                fs.task_kind,
                fs.n,
                fs.rows,
                exp_s,
            ])
    lines.append(table(inv_rows, ["file", "status", "task", "n", "rows", "expected_rows"]))
    lines.append("")

    # Per-file details.
    detail_rows = []
    for fs in file_stats:
        comp = completeness_summary(fs, args.steps, args.num_seeds)
        detail_rows.append([
            fs.path.name,
            len(comp["task_ids"]),
            len(comp["modes"]),
            len(fs.groups),
            fs.parse_errors,
            fs.vector_checks,
            fs.vector_mismatch_count,
            len(comp["missing_groups"]),
            len(comp["seed_count_problems"]),
            len(comp["incomplete_step_groups"]),
        ])
    if detail_rows:
        lines.append("Per-file structure summary:")
        lines.append(table(detail_rows, ["file", "tasks", "modes", "groups", "parse_err", "vec_checks", "vec_bad", "missing_groups", "bad_seed_counts", "bad_steps"]))
        lines.append("")

    # Alerts.
    alerts: List[str] = []
    for p in expected_files:
        if not p.exists() and not Path(str(p) + ".gz").exists():
            alerts.append(f"Missing expected file: {p.name}")
    for fs in file_stats:
        comp = completeness_summary(fs, args.steps, args.num_seeds)
        if fs.parse_errors:
            alerts.append(f"{fs.path.name}: {fs.parse_errors} row parse errors.")
        if fs.vector_mismatch_count:
            alerts.append(f"{fs.path.name}: {fs.vector_mismatch_count}/{fs.vector_checks} sampled theta/grad vector length checks failed.")
        if comp["missing_modes"]:
            alerts.append(f"{fs.path.name}: missing expected modes: {', '.join(comp['missing_modes'])}")
        if comp["seed_count_problems"]:
            tid, _tname, observed, expected, seed_preview = comp["seed_count_problems"][0]
            alerts.append(f"{fs.path.name}: {len(comp['seed_count_problems'])} tasks have an unexpected number of run seeds. Example task {tid}: observed {observed}, expected {expected}, preview={seed_preview}.")
        if comp["missing_groups"]:
            preview = "; ".join(f"task {tid}/seed {seed}/{mode}" for tid, _tn, seed, mode in comp["missing_groups"][:5])
            alerts.append(f"{fs.path.name}: {len(comp['missing_groups'])} missing task/seed/mode groups. First: {preview}")
        if comp["incomplete_step_groups"]:
            g, missing_preview, nmiss, extra_preview, nextra = comp["incomplete_step_groups"][0]
            alerts.append(f"{fs.path.name}: {len(comp['incomplete_step_groups'])} groups with missing/extra steps. Example {g.task_name}/{g.mode}: missing={nmiss}, extra={nextra}.")
    for g in all_groups:
        if g.nonfinite_metric_rows:
            alerts.append(f"{Path(g.source_file).name}: {g.task_name}/{g.mode} has {g.nonfinite_metric_rows} nonfinite metric rows.")
        if g.nonfinite_grad_rows:
            alerts.append(f"{Path(g.source_file).name}: {g.task_name}/{g.mode} has {g.nonfinite_grad_rows} nonfinite grad_norm rows.")
        if g.max_state_norm_error > args.state_norm_tol:
            alerts.append(f"{Path(g.source_file).name}: {g.task_name}/{g.mode} max state_norm_error={fmt_float(g.max_state_norm_error)} > {args.state_norm_tol}.")
        if finite(g.final_metric()) and finite(g.initial_metric()) and g.final_metric() > g.initial_metric() + args.worse_tol:
            alerts.append(f"{Path(g.source_file).name}: {g.task_name}/{g.mode} final metric worse than initial: {fmt_float(g.initial_metric())} -> {fmt_float(g.final_metric())}.")
        if finite(g.max_grad_norm) and g.max_grad_norm > args.grad_explosion_threshold:
            alerts.append(f"{Path(g.source_file).name}: {g.task_name}/{g.mode} max grad_norm={fmt_float(g.max_grad_norm)} exceeds threshold {args.grad_explosion_threshold}.")

    lines.append("2. ALERTS / DATA-QUALITY FLAGS")
    lines.append(line())
    if alerts:
        for i, a in enumerate(alerts, start=1):
            lines.append(f"[{i}] {a}")
    else:
        lines.append("No structural or numerical alerts triggered by the configured thresholds.")
    lines.append("")

    # Overall aggregate rankings by task kind.
    lines.append("3. OVERALL MODE AGGREGATES")
    lines.append(line())
    for task_kind in ("VQE", "MET"):
        gs_task = [g for g in all_groups if g.task.upper() == task_kind]
        if not gs_task:
            continue
        lines.append(f"{task_kind}: {task_score_label(task_kind)}")
        agg = aggregate_by_mode(gs_task)
        wins_final = win_counts(gs_task, by="final")
        wins_best = win_counts(gs_task, by="best")
        ranks = mode_rank_distribution(gs_task, by="final")
        rows = []
        for mode, d in sorted(agg.items(), key=lambda kv: (math.inf if not finite(kv[1]["mean_final"]) else kv[1]["mean_final"], kv[0])):
            rank_list = ranks.get(mode, [])
            rows.append([
                mode,
                d["count"],
                fmt_float(d["mean_final"]),
                fmt_float(d["median_final"]),
                fmt_float(d["std_final"]),
                fmt_float(d["mean_best"]),
                fmt_float(d["mean_improvement"]),
                fmt_float(d["mean_overlap"]),
                fmt_float(d["mean_decoded_probability"]),
                wins_final.get(mode, 0),
                wins_best.get(mode, 0),
                fmt_float(mean(rank_list)),
            ])
        lines.append(table(rows, ["mode", "groups", "mean_final", "median_final", "std_final", "mean_best", "mean_improve", "mean_overlap", "mean_dec_prob", "final_wins", "best_wins", "mean_rank"]))
        lines.append("")

    # Aggregate by n and task kind.
    lines.append("4. AGGREGATES BY n")
    lines.append(line())
    for task_kind in ("VQE", "MET"):
        for n in sorted({g.n for g in all_groups if g.task.upper() == task_kind}):
            gs_nt = [g for g in all_groups if g.task.upper() == task_kind and g.n == n]
            if not gs_nt:
                continue
            lines.append(f"{task_kind}, n={n}: {task_score_label(task_kind)}")
            agg = aggregate_by_mode(gs_nt)
            wins_final = win_counts(gs_nt, by="final")
            rows = []
            for mode, d in sorted(agg.items(), key=lambda kv: (math.inf if not finite(kv[1]["mean_final"]) else kv[1]["mean_final"], kv[0])):
                rows.append([
                    mode,
                    d["count"],
                    fmt_float(d["mean_final"]),
                    fmt_float(d["median_final"]),
                    fmt_float(d["p25_final"]),
                    fmt_float(d["p75_final"]),
                    fmt_float(d["mean_best"]),
                    wins_final.get(mode, 0),
                ])
            lines.append(table(rows, ["mode", "groups", "mean_final", "median_final", "p25", "p75", "mean_best", "wins"]))
            lines.append("")

    # Per task final rankings.
    lines.append("5. PER-TASK FINAL RANKINGS")
    lines.append(line())
    bucket: Dict[Tuple[str, int, str, str, str], List[GroupStats]] = defaultdict(list)
    for g in all_groups:
        bucket[(g.task, g.n, g.task_id, g.task_name, g.run_seed_value)].append(g)
    for key in sorted(bucket.keys(), key=lambda x: (x[0], x[1], x[2], x[3], x[4])):
        task_kind, n, tid, tname, run_seed = key
        gs = rank_groups_for_task(bucket[key])
        lines.append(f"{task_kind} n={n} task_id={tid} task_name={tname} run_seed={run_seed} | {task_score_label(task_kind)}")
        rows = []
        for rank, g in enumerate(gs, start=1):
            rows.append([
                rank,
                g.mode,
                fmt_float(g.initial_metric()),
                fmt_float(g.final_metric()),
                fmt_float(g.best_metric()),
                fmt_int(g.best_step()),
                fmt_float(g.improvement()),
                fmt_float(g.relative_improvement()),
                fmt_float(g.final_target_overlap()),
                fmt_float(g.final_decoded_probability()),
                fmt_float(g.final_grad_norm()),
                fmt_float(g.mean_step_angle()),
                fmt_float(g.total_line_evals()),
            ])
        lines.append(table(rows, ["rank", "mode", "initial", "final", "best", "best_step", "improve", "rel_improve", "final_overlap", "final_dec_prob", "final_grad", "mean_angle", "line_evals_total"]))
        lines.append("")

    # Threshold diagnostics.
    lines.append("6. THRESHOLD-HIT DIAGNOSTICS")
    lines.append(line())
    for task_kind in ("VQE", "MET"):
        gs_task = [g for g in all_groups if g.task.upper() == task_kind]
        if not gs_task:
            continue
        lines.append(f"{task_kind}: first step with metric <= threshold; blank means never reached.")
        rows = []
        for g in sorted(gs_task, key=lambda x: (x.n, x.task_id, x.mode)):
            rows.append([
                g.n,
                g.task_id,
                g.task_name,
                g.run_seed_value,
                g.mode,
                *[fmt_int(g.first_step_below(t)) for t in THRESHOLDS],
            ])
        lines.append(table(rows, ["n", "task_id", "task_name", "run_seed", "mode", "<=1e-2", "<=1e-4", "<=1e-6", "<=1e-8"]))
        lines.append("")

    # Convergence quality / pathology metrics per group.
    lines.append("7. CONVERGENCE-SHAPE AND NUMERICAL DIAGNOSTICS")
    lines.append(line())
    rows = []
    for g in sorted(all_groups, key=lambda x: (x.task, x.n, x.task_id, x.mode)):
        rows.append([
            g.task,
            g.n,
            g.task_id,
            g.run_seed_value,
            g.mode,
            g.rows,
            len(g.steps),
            fmt_float(g.monotone_fraction()),
            fmt_float(g.terminal_slope()),
            fmt_float(g.max_state_norm_error),
            fmt_float(g.max_grad_norm),
            fmt_float(g.max_state_grad_norm),
            g.vector_checks,
            g.vector_mismatch_count,
        ])
    lines.append(table(rows, ["task", "n", "task_id", "run_seed", "mode", "rows", "steps", "monotone_frac", "tail_log10_slope", "max_norm_err", "max_grad", "max_state_grad", "vec_checks", "vec_bad"]))
    lines.append("")

    # Recommended plot candidates based on diagnostics.
    lines.append("8. AUTOMATIC PLOT-SCOUTING NOTES")
    lines.append(line())
    for task_kind in ("VQE", "MET"):
        gs_task = [g for g in all_groups if g.task.upper() == task_kind]
        if not gs_task:
            continue
        wins = win_counts(gs_task, by="final")
        agg = aggregate_by_mode(gs_task)
        if agg:
            best_mean_mode = min(agg, key=lambda m: math.inf if not finite(agg[m]["mean_final"]) else agg[m]["mean_final"])
            lines.append(f"{task_kind}: best mean-final mode = {best_mean_mode} with mean_final={fmt_float(agg[best_mean_mode]['mean_final'])}; final win counts = {dict(wins)}")
        # Most diagnostic task: largest spread between best and worst final mode.
        spreads = []
        for key, gs in bucket.items():
            if key[0].upper() != task_kind:
                continue
            finals = [g.final_metric() for g in gs if finite(g.final_metric())]
            if len(finals) >= 2:
                spreads.append((max(finals) - min(finals), key, gs))
        if spreads:
            spread, key, _gs = max(spreads, key=lambda x: x[0])
            lines.append(f"{task_kind}: most visually separating task appears to be n={key[1]}, task_id={key[2]}, task_name={key[3]}, final spread={fmt_float(spread)}.")
    lines.append("")

    lines.append("9. INTERPRETATION CAUTIONS")
    lines.append(line())
    lines.append("- This report ranks modes using final scalar gaps. It does not inspect the full theta/grad vectors except optional length checks.")
    lines.append("- All rows are Hopf-coordinate rows. The expected modes are Hopf-EGT-CG, Hopf-Riemannian-LBFGS, and Hopf-Riemannian-BB only.")
    lines.append("- The modes use cost and Hopf-coordinate-gradient information. Nonzero last_line_evals are line-search work calls, not shot counts in the infinite-shot synthetic data; for strong-Wolfe steps this counter includes trial cost evaluations and curvature-check gradient-oracle calls.")
    lines.append("- CSV vector columns are intentionally left unparsed for speed; rerun with --vector-check all only for small n/files.")
    lines.append("")

    return "\n".join(lines)


def parse_ns(s: str) -> List[int]:
    out: List[int] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Diagnose synthetic Hopf stress-test CSV files.")
    p.add_argument("--indir", default=".", help="Directory containing vqe_hopf_data_n*.csv and met_hopf_data_n*.csv files.")
    p.add_argument("--ns", default="6-10", help="Comma/range list of n values, e.g. '6-10' or '6,7,8'.")
    p.add_argument("--steps", type=int, default=200, help="Expected number of optimization updates. Step 0 is expected too.")
    p.add_argument("--num-seeds", type=int, default=10, help="Expected number of initial-state seeds per task.")
    p.add_argument("--out", default="hopf_data_diagnostics.txt", help="Output diagnostic text file.")
    p.add_argument("--vector-check", choices=["none", "sample", "all"], default="sample", help="Check theta/grad vector lengths. sample checks step 0 and final step only.")
    p.add_argument("--state-norm-tol", type=float, default=1e-8, help="Alert if max state_norm_error exceeds this.")
    p.add_argument("--worse-tol", type=float, default=1e-12, help="Alert if final metric exceeds initial metric by more than this.")
    p.add_argument("--grad-explosion-threshold", type=float, default=1e8, help="Alert if grad_norm exceeds this.")
    p.add_argument("--no-print", action="store_true", help="Write the report to file but do not print it to stdout.")
    args = p.parse_args()
    args.ns = parse_ns(args.ns)
    if args.num_seeds < 1:
        raise ValueError("--num-seeds must be >= 1")
    return args


def main() -> None:
    args = parse_args()
    indir = Path(args.indir)
    expected = expected_file_list(indir, args.ns)
    files = find_input_files(indir, args.ns)

    file_stats: List[FileStats] = []
    for path in files:
        print(f"Reading {path} ...", file=sys.stderr)
        file_stats.append(process_file(path, args.steps, args.vector_check))

    report = build_report(file_stats, expected, args)
    outpath = Path(args.out)
    if not outpath.is_absolute():
        outpath = indir / outpath
    outpath.write_text(report, encoding="utf-8")
    print(f"Wrote diagnostic report: {outpath}", file=sys.stderr)
    if not args.no_print:
        print(report)


if __name__ == "__main__":
    main()
