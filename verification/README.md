# Verification and reproduction

The deterministic scripts below regenerate the CSV datasets, diagnostics, and
figures from the repository sources.

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

python hopf_complex.py \
    --quick \
    --outdir figures_safeguards \
    --csv-name complex_hopf_smoke.csv \
    --plot-name hopf_complex_smoke.png

python finite_shot_sanity_check.py \
    --shots 100 1000 \
    --trials 3 \
    --output-dir figures_safeguards
```

The `VQE_qibo.py --sampler statevector` path is the dependency-light fallback. Use `--sampler qibo-explicit` after installing Qibo to force the explicit Qibo circuit path.

The [engineering guide](../docs/ENGINEERING_GUIDE.md#numerical-clipping-for-optimizer-traces)
specifies the chart-clipping tolerances. The default diagonal-metric division
floor is `1e-12`, set by `--metric-eps`.


## What to inspect in the output

The [reviewer route](../review/README.md) connects these checks to the complete
[scientific account](../science/README.md). Run small identity checks before
regenerating the full numerical study. All-size statements are established by
the linked proofs; finite runs check their implementation.

| Check | Expected result | Interpretation |
|---|---|---|
| `python hopf_utils.py` | Real and complex state round-trip residuals near floating-point precision; normalized tangent synthesis error below `1e-10`; both synthesis messages report correct. | The script prints diagnostics, so inspect them as well as its exit status. Optional Qibo residuals compare the prepared state column. |
| `hopf_gate_count.py --nmin 2 --nmax 6` | Real assigned counts `4, 24, 100, 360, 1168`; complex counts `4, 24, 100, 360, 1300`. | Integer schedule counts and closed forms use the declared no-clean-ancilla charges. |
| `finite_shot_sanity_check.py` | Exact branch and tree/Jacobian residuals near roundoff; all state/probability normalization diagnostics small. | Sampling error is separate. The default error-versus-shots slope is approximately `-0.5`; a three-trial smoke run is not a precise slope measurement. |
| `hopf_complex.py --quick` | Scrambler, gradient and metric-lift residuals pass the built-in `1e-8` check; the small CSV and plot are produced. | Inspect the printed completeness and finite-value diagnostics. The full study's normalization and matched-start alert threshold is `1e-10`. |
| `VQE_qibo.py --sampler qibo-explicit` | The script reports its exact layerwise energy-identity check and generates real/complex trajectories. | This is a four-qubit circuit and sampling check, with its shot convention in [Experiments](../docs/EXPERIMENTS.md#10-local-qibo-layerwise-safeguard). |

The default fixed-state finite-shot run reports exact branch-identity and
tree/Jacobian residuals around `2.5e-16`, and mean relative errors around
`0.1860`, `0.05634`, and `0.01860` at `100`, `1000`, and `10000` shots per
signed branch. These are reference diagnostics for the fixed seed, not
universal tolerances for statistical data. See [experiment definitions and
reported results](../docs/EXPERIMENTS.md).

## Generate the full synthetic CSV data

```bash
mkdir -p data

for n in 6 7 8 9 10; do
    python hopf_data.py --n "$n" --steps 200 --num-seeds 10 --outdir data
    python adam_data.py --n "$n" --steps 200 --num-seeds 10 --outdir data
done
```

For a smaller run, reduce `--num-seeds`, reduce `--steps`, or use `--quick`. Pass the same `--steps` and `--num-seeds` values to the diagnostic scripts.

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

## Standalone experiments and safeguards

```bash
python hopf_gate_count.py
MPLBACKEND=Agg python VQE_qibo.py
python finite_shot_sanity_check.py
python hopf_complex.py
```

The CNOT script checks generated schedules against the closed-form count formulas. The VQE/Qibo script is a small real/complex layerwise gradient-access realizability check. The finite-shot script checks the signed-branch estimator at a fixed state, and `hopf_complex.py` runs the focused `n=6` complex-Hopf stress test. These scripts are separate from the multi-size real-Hopf CSV pipeline.

## Paper figures and numerical output

| Paper item | Source in the scientific account | Regeneration |
|---|---|---|
| Figure 1: aggregate final precision and threshold hitting | [Numerical studies](../science/numerical-studies.md#S4.F1) | Full real-Hopf and Adam CSV pipelines above, then `plot_hopf.py`; output `hopf_geometric_summary_clean.png`. |
| Figure 2: per-task convergence at `n = 10` | [Numerical studies](../science/numerical-studies.md#S4.F2) | Same CSVs with `--detail-n 10`; output `hopf_geometric_n10_convergence_clean.png`. |
| Figures 3–6: tree, ansatz circuit, indexed tangent and branch circuits | [Four-qubit construction](../science/four-qubit-example.md) | Original paper assets and captions accompany the explicit construction; inspect the gate/parameter assignments against `hopf_utils.py` and `VQE_qibo.py`. These are construction diagrams rather than sampled data. |
| Tables I–III: gate schedule, metric and tangent assignments | [Four-qubit construction](../science/four-qubit-example.md) | Compare `gates_order`, `metric_diagonal`, `jacobian`, and `theta_hopf_tangent_state` with the displayed entries. |

The two paper data figures are distinct from the additional committed complex,
finite-shot, and local-circuit summaries. The [complete paper map](../review/PAPER_MAP.md)
links every numbered item to its Markdown home.

For the full real study, each size produces `18,090` rows per geometry-method
application-class CSV and `12,060` rows per Adam-baseline application-class CSV.
The expected aggregate contains `150` traces per optimizer per application
class. Use the [reported aggregate targets](../docs/EXPERIMENTS.md#paper-reported-aggregate-reproduction-targets)
only after the file, seed and step diagnostics pass. Tiny medians are sensitive
to floating-point rounding; the declared task, seed, size, update and optimizer
settings define the comparison.
