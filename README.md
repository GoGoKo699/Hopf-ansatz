# Hopf Ansatz Stress-Test Code

This repository contains reproducible synthetic stress-test code for the Hopf ansatz.

The code generates synthetic VQE and metrology-inspired optimization traces, runs seed-aware diagnostics, and produces clean GitHub-facing plots. Generated datasets, diagnostic text files, and paper-facing figure outputs are intentionally not committed because they can be regenerated from the scripts.

## Repository contents

| File | Purpose |
|---|---|
| `hopf_utils.py` | Core Hopf ansatz utilities: real/complex coordinate maps, inverse maps, Jacobians, diagonal metrics, tangent-state assignments, gate schedules, and optional Qibo circuit checks. |
| `hopf_data.py` | Generates geometry-native Hopf optimizer data for synthetic VQE and metrology-inspired tasks. |
| `adam_data.py` | Generates adaptive-line-search Adam baseline data for Hopf coordinates and Möttönen physical-angle coordinates. |
| `diagnose_hopf_data.py` | Streaming diagnostics for Hopf optimizer CSV files. |
| `diagnose_adam.py` | Streaming diagnostics for Adam baseline CSV files. |
| `plot_hopf.py` | Generates geometry-native Hopf plots from Hopf CSV files. |
| `plot_all.py` | Generates plain-scale GitHub-facing comparison plots including Hopf optimizers and Adam baselines. |

## Preview plots

The following plots are GitHub-facing comparison plots at `n=10`. They are not the paper-facing figures.

![VQE cost traces, n=10](all_vqe_n10_clean.png)

![MET cost traces, n=10](all_met_n10_clean.png)

## Requirements

Use Python 3.10 or newer.

Install the required packages:

```bash
pip install numpy scipy matplotlib
```

Optional package:

```bash
pip install qibo
```

`qibo` is only needed for optional circuit verification in `hopf_utils.py`. The data-generation, diagnostic, and plotting scripts do not require it.

## Quick smoke test

Run a small two-qubit test before generating the full datasets.

```bash
mkdir -p data_smoke figures_smoke

python hopf_data.py --n 2 --steps 2 --num-seeds 2 --scramble-depth 1 --outdir data_smoke
python adam_data.py --n 2 --steps 2 --num-seeds 2 --scramble-depth 1 --outdir data_smoke

python diagnose_hopf_data.py --indir data_smoke --ns 2 --steps 2 --num-seeds 2
python diagnose_adam.py --indir data_smoke --ns 2 --steps 2 --num-seeds 2

python plot_all.py --hopf-dir data_smoke --adam-dir data_smoke --n 2 --outdir figures_smoke
```

Expected generated files include:

```text
data_smoke/vqe_hopf_data_n2.csv
data_smoke/met_hopf_data_n2.csv
data_smoke/vqe_adam_data_n2.csv
data_smoke/met_adam_data_n2.csv
hopf_data_diagnostics.txt
adam_data_diagnostics.txt
figures_smoke/all_vqe_n2_clean.png
figures_smoke/all_met_n2_clean.png
```

## Full data generation

The default experiment uses:

- system sizes `n = 6, 7, 8, 9, 10`;
- `200` optimization updates, with step `0` also recorded;
- `10` deterministic initial-state seeds per synthetic task;
- three VQE tasks and three metrology-inspired tasks per `n`;
- three geometry-native Hopf optimizer modes;
- two Adam baseline modes.

Generate the full Hopf and Adam datasets:

```bash
mkdir -p data

for n in 6 7 8 9 10; do
    python hopf_data.py --n "$n" --outdir data
    python adam_data.py --n "$n" --outdir data
done
```

Default Hopf outputs:

```text
data/vqe_hopf_data_n6.csv
data/met_hopf_data_n6.csv
...
data/vqe_hopf_data_n10.csv
data/met_hopf_data_n10.csv
```

Default Adam outputs:

```text
data/vqe_adam_data_n6.csv
data/met_adam_data_n6.csv
...
data/vqe_adam_data_n10.csv
data/met_adam_data_n10.csv
```

Per system size and task family, the default row counts are:

```text
Hopf CSV: 3 tasks × 10 seeds × 3 modes × 201 steps = 18090 rows
Adam CSV: 3 tasks × 10 seeds × 2 modes × 201 steps = 12060 rows
```

The CSV files can be large because each row stores scalar diagnostics as well as semicolon-separated parameter and gradient vectors.

The diagnostic and plotting scripts also accept `.csv.gz` files, so generated CSVs may be compressed after creation:

```bash
gzip data/*.csv
```

## Synthetic tasks

Each task is scrambled by a fixed real orthogonal circuit, so the optimum is generic in the computational basis.

### VQE tasks

| Task ID | Name | Description |
|---|---|---|
| `VQE-1` | Parent Hamiltonian | `H = I - |tau><tau|` |
| `VQE-2` | Scrambled Hamming spectrum | A diagonal Hamming-distance spectrum conjugated by the real scrambler. |
| `VQE-3` | Small-gap spectrum | A scrambled diagonal Hamiltonian with one ground state, one nearby excited state at gap `1e-2`, and all other levels at unit energy. |

### Metrology-inspired tasks

| Task ID | Name | Description |
|---|---|---|
| `MET-1` | Single-target Fisher | Fixed-readout Fisher proxy, minimized as `1 - F`. |
| `MET-2` | QFI superposition | Normalized-QFI objective whose optimum is an extremal-generator superposition. |
| `MET-3` | Balanced Fisher | Nonlinear two-target soft-min Fisher objective. |

Lower gap or cost is better in all diagnostic summaries.

## Optimizer modes

### Geometry-native Hopf optimizers

Generated by `hopf_data.py`:

```text
Hopf-EGT-CG
Hopf-Riemannian-LBFGS
Hopf-Riemannian-BB
```

These optimizers use the same cost-plus-Hopf-coordinate-gradient interface. The Hopf diagonal metric is used to lift coordinate gradients to state-sphere gradients, and trial states are moved on the normalized real state sphere.

### Adam baselines

Generated by `adam_data.py`:

```text
Hopf-Adam
Mottonen-ideal-PS-Adam
```

`Hopf-Adam` applies Adam to Hopf coordinates.

`Mottonen-ideal-PS-Adam` applies Adam to physical post-multiplexing Möttönen rotation angles with an exact gradient equivalent to infinite-shot parameter shift.

Both Adam baselines use adaptive cost-only backtracking. Trial evaluations are additional objective calls, not additional gradient, metric-estimation, Hessian, or Hamiltonian-action primitives.

## Seed convention

By default, every fixed synthetic task is run from `10` deterministic initial-state seeds.

The default seed rule is:

```text
run_seed = problem_seed + seed_offset + seed_index
```

with:

```text
seed_offset = 77777
seed_index = 0, 1, ..., num_seeds - 1
```

The default `num_seeds` is `10`.

To change the number of initial states, pass the same value to both the data-generation and diagnostic scripts:

```bash
python hopf_data.py --n 8 --num-seeds 5 --outdir data
python adam_data.py --n 8 --num-seeds 5 --outdir data

python diagnose_hopf_data.py --indir data --ns 8 --num-seeds 5
python diagnose_adam.py --indir data --ns 8 --num-seeds 5
```

## Diagnostics

Run diagnostics after generating the datasets.

```bash
python diagnose_hopf_data.py \
    --indir data \
    --ns 6-10 \
    --steps 200 \
    --num-seeds 10 \
    --out hopf_data_diagnostics.txt

python diagnose_adam.py \
    --indir data \
    --ns 6-10 \
    --steps 200 \
    --num-seeds 10 \
    --out adam_data_diagnostics.txt
```

The reports check:

- expected files and row counts;
- expected task, mode, seed, and step completeness;
- CSV parse errors;
- sampled parameter and gradient vector lengths;
- nonfinite metrics or gradients;
- state-norm errors;
- aggregate final gaps;
- win counts and final rankings.

The diagnostic scripts are streaming and avoid loading the full vector columns into pandas.

## Plotting

### Geometry-native Hopf plots

Use `plot_hopf.py` for plots involving only the three geometry-native Hopf optimizers.

```bash
python plot_hopf.py \
    --indir data \
    --outdir figures_hopf_clean \
    --ns 6-10 \
    --detail-n 10 \
    --formats pdf png \
    --error-band sem \
    --band-alpha 0.18
```

Default outputs:

```text
figures_hopf_clean/hopf_geometric_summary_clean.pdf
figures_hopf_clean/hopf_geometric_summary_clean.png
figures_hopf_clean/hopf_geometric_n10_convergence_clean.pdf
figures_hopf_clean/hopf_geometric_n10_convergence_clean.png
```

These are paper-facing-style Hopf-only plots. The generated files are reproducible and are not committed by default.

### Hopf-plus-Adam GitHub comparison plots

Use `plot_all.py` for plain-scale comparison plots including both Hopf optimizers and Adam baselines.

```bash
python plot_all.py \
    --hopf-dir data \
    --adam-dir data \
    --outdir figures_all_clean \
    --n 10 \
    --formats pdf png \
    --error-band sem \
    --band-alpha 0.18
```

Default outputs:

```text
figures_all_clean/all_vqe_n10_clean.pdf
figures_all_clean/all_vqe_n10_clean.png
figures_all_clean/all_met_n10_clean.pdf
figures_all_clean/all_met_n10_clean.png
```

To use the generated PNGs as the README preview images, copy them to the repository root:

```bash
cp figures_all_clean/all_vqe_n10_clean.png .
cp figures_all_clean/all_met_n10_clean.png .
```

The plotting scripts support the following error-band options:

```text
sem
std
none
```

The default is `sem`. The `--band-alpha` option controls the transparency of the shaded bands.

## Optional Hopf utility checks

Run:

```bash
python hopf_utils.py
```

This executes consistency checks for inverse mapping, tangent-state synthesis, metric/Jacobian agreement, and optional Qibo circuit agreement when Qibo is installed.

## Reproducibility notes

The scripts use deterministic random seeds by default. Re-running the same commands with the same Python, NumPy, and SciPy versions should reproduce the same synthetic task definitions and initial-state seeds.

The experiments are exact state-vector simulations. They do not model finite-shot sampling noise.

The generated stress-test datasets use the real Hopf chart. `hopf_utils.py` also contains real and complex Hopf utilities used by the paper.

## Suggested `.gitignore`

Generated data, diagnostics, and paper-facing figures are reproducible and should usually stay out of version control.

```gitignore
# Python
__pycache__/
*.py[cod]
.venv/
venv/

# Generated CSV data
data/
data_smoke/
*_hopf_data_n*.csv
*_adam_data_n*.csv
*_hopf_data_n*.csv.gz
*_adam_data_n*.csv.gz

# Diagnostics
*_diagnostics.txt
hopf_data_diagnostics.txt
adam_data_diagnostics.txt

# Generated figures
figures_smoke/
figures_hopf_clean/
figures_all_clean/

# Local notebook/editor files
.ipynb_checkpoints/
.DS_Store
```

The two GitHub preview images may be committed at the repository root:

```text
all_vqe_n10_clean.png
all_met_n10_clean.png
```

## Citation

If you use this code, cite the accompanying Hopf ansatz paper.

Add the final arXiv, DOI, or BibTeX entry here when available.
