# Reproducibility checklist

This repository intentionally does not track the full generated CSV datasets. They are large derived artifacts and can be regenerated from deterministic scripts in this repository.

The commands below assume they are run from the repository root.

## Environment

Use Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Optional explicit Qibo circuit checks require:

```bash
python -m pip install -r requirements-optional.txt
```

## Minimal smoke test

These commands check the core utilities and small safeguard paths without generating the full synthetic CSV archive.

```bash
python hopf_utils.py

mkdir -p figures_safeguards
python hopf_gate_count.py \
    --nmin 2 \
    --nmax 6 \
    --out figures_safeguards/hopf_gate_count_smoke.pdf

MPLBACKEND=Agg python VQE_qibo.py \
    --steps 2 \
    --shots 20 \
    --sampler statevector \
    --output figures_safeguards/VQE_qibo_smoke.png \
    --log-every 0
```

The `VQE_qibo.py --sampler statevector` path is the dependency-light fallback. Use `--sampler qibo-explicit` after installing Qibo to force the explicit Qibo circuit path.

## Generate the full synthetic CSV data

```bash
mkdir -p data

for n in 6 7 8 9 10; do
    python hopf_data.py --n "$n" --steps 200 --num-seeds 10 --outdir data
    python adam_data.py --n "$n" --steps 200 --num-seeds 10 --outdir data
done
```

For a smaller development run, reduce `--num-seeds`, reduce `--steps`, or use `--quick`. Pass the same `--steps` and `--num-seeds` values to the diagnostic scripts.

## Diagnostics

```bash
repo_root="$(pwd)"
mkdir -p diagnostics

python diagnose_hopf.py \
    --indir data \
    --ns 6-10 \
    --steps 200 \
    --num-seeds 10 \
    --out "$repo_root/diagnostics/hopf_data_diagnostics.txt"

python diagnose_adam.py \
    --indir data \
    --ns 6-10 \
    --steps 200 \
    --num-seeds 10 \
    --out "$repo_root/diagnostics/adam_data_diagnostics.txt"
```

The diagnostic scripts check file completeness, seed coverage, step coverage, scalar convergence summaries, sampled vector lengths, nonfinite values, state-norm errors, final rankings, and threshold hits. They avoid loading the full vector columns into pandas.

## Plots

```bash
mkdir -p figures
python plot_hopf.py \
    --indir data \
    --adam-dir data \
    --outdir figures \
    --ns 6-10 \
    --detail-n 10 \
    --formats pdf png
```

## Standalone safeguards

```bash
python hopf_gate_count.py
MPLBACKEND=Agg python VQE_qibo.py
```

The CNOT script checks generated schedules against the closed-form count formulas. The VQE/Qibo script is a small real/complex layerwise gradient-access realizability check; it is separate from the large synthetic optimizer traces.
