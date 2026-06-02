# Hopf Ansatz Stress-Test Code

This repository contains reproducible synthetic stress-test code for the Hopf ansatz.

The code generates synthetic VQE and metrology-inspired optimization traces, runs seed-aware diagnostics, produces GitHub-facing comparison plots, and includes standalone safeguard scripts for the Hopf CNOT-count formulas and the circuit-level real-Hopf gradient construction.

Full generated CSV datasets and local safeguard figure outputs are not included in this archive. Included diagnostic text files and the paper PDF are release artifacts; they can be regenerated from the scripts once the full CSV datasets are present.

## Repository contents

| File | Role |
|---|---|
| `hopf_utils.py` | Core Hopf utilities: real and complex coordinate maps, inverse maps, Jacobians, diagonal metrics, tangent-state assignments, gate schedules, and optional Qibo circuit checks. |
| `hopf_data.py` | Generates synthetic stress-test CSVs for the three geometry-native Hopf optimizers. |
| `adam_data.py` | Generates synthetic stress-test CSVs for the two adaptive Adam baselines. |
| `diagnose_hopf.py` | Checks completeness and numerical quality of Hopf optimizer CSVs. |
| `diagnose_adam.py` | Checks completeness and numerical quality of Adam baseline CSVs. |
| `plot_hopf.py` | Generates summary and convergence plots from Hopf CSVs. |
| `hopf_gate_count.py` | Safeguard script for the CNOT-count formulas. It builds real and complex Hopf gate schedules and checks numerical counts against the closed-form binomial formulas. |
| `VQE_qibo.py` | Circuit-level safeguard script for the real-Hopf gradient formula. It builds explicit Qibo ancilla circuits for transition-moment estimation in a small VQE example. |

The VQE test contains three scrambled Hamiltonian tests:

**Parent Hamiltonian:** recover a scrambled target state by minimizing the parent-Hamiltonian gap. The task cost is the missing target overlap,

```math
\mathcal{C}_{\mathrm{parent}}(\psi)
=
1 - |\langle \tau|\psi\rangle|^2 .
```

Here τ is the scrambled target state. The minimum is zero, reached when the variational state ψ equals τ up to the real global sign convention used in the code.

**Hamming spectrum:** minimize a structured diagonal Hamming-distance spectrum hidden by a real orthogonal scrambling circuit. The task cost is the expected normalized Hamming distance in the unscrambled basis,

```math
\mathcal{C}_{\mathrm{Ham}}(\psi)
=
\sum_x
\frac{d_{\mathrm{H}}(x,x_0)}{n}
\,
|(S^\top\psi)_x|^2 .
```

Here S is the fixed real scrambler and x₀ is the hidden target bit string. Physically, this tests whether the optimizer can find the scrambled computational-basis ground state of a simple but hidden diagonal spectrum.

**Small-gap spectrum:** minimize a scrambled Hamiltonian with a deliberately small spectral gap near the ground state. The task cost is the expected value of a scrambled diagonal spectrum,

```math
\mathcal{C}_{\mathrm{gap}}(\psi)
=
\sum_x
E_x
|(S^\top\psi)_x|^2 ,
```

with

```math
E_{x_0}=0,\qquad
E_{x_1}=10^{-2},\qquad
E_x=1\ \text{for all other }x .
```

The physical stress feature is the nearby distractor state x₁: the optimizer must distinguish the true ground state from a low-lying excited state.

The metrology test contains three nonlinear sensing-inspired tests:

**Single-target Fisher:** maximize a fixed-readout Fisher proxy for one scrambled target state. The task cost is

```math
\mathcal{C}_{\mathrm{single}}(\psi)
=
1 - F,
\qquad
F =
|\langle \tau|\psi\rangle|^4 .
```

Equivalently, the code first computes the target probability A = |⟨τ|ψ⟩|² and then uses F = A². Physically, this rewards concentration of the probe state on one scrambled readout pattern.

**QFI superposition:** maximize normalized pure-state QFI for a scrambled diagonal generator. In the unscrambled basis, the code computes

```math
\mu =
\sum_x g_x |(S^\top\psi)_x|^2,
\qquad
\nu =
\sum_x g_x^2 |(S^\top\psi)_x|^2,
```

and then uses the normalized QFI objective

```math
F_Q^{\mathrm{norm}}(\psi)
=
\frac{4(\nu-\mu^2)}{\mathrm{span}(G)^2}.
```

The task cost is

```math
\mathcal{C}_{\mathrm{QFI}}(\psi)
=
1 - F_Q^{\mathrm{norm}}(\psi).
```

Physically, this rewards a probe state with large generator variance. The optimum is an equal superposition of the scrambled minimum- and maximum-generator eigenstates.

**Balanced Fisher:** optimize a two-target Fisher objective where the optimum requires balanced overlap with two scrambled target states. The code computes

```math
e_1 = |\langle \tau_1|\psi\rangle|^2,
\qquad
e_2 = |\langle \tau_2|\psi\rangle|^2,
```

then

```math
F_1=e_1^2,
\qquad
F_2=e_2^2.
```

The balanced Fisher score is a soft minimum,

```math
F_{\mathrm{bal}}(\psi)
=
-\frac{1}{\beta}
\log\!\left(
\frac{
e^{-\beta F_1}
+
e^{-\beta F_2}
}{2}
\right),
\qquad
\beta=20 .
```

The task cost is

```math
\mathcal{C}_{\mathrm{bal}}(\psi)
=
\frac{1}{4}
-
F_{\mathrm{bal}}(\psi).
```

The optimum has balanced probability on the two scrambled targets, giving approximately e₁ = e₂ = 1/2, hence F₁ = F₂ = 1/4 and zero cost.

To regenerate the plots locally, generate the datasets and run `plot_hopf.py`; see the sections below.

## Requirements

Use Python 3.10 or newer.

Install the required packages for data generation, diagnostics, CNOT counting, and plotting:

```bash
pip install numpy scipy matplotlib
```

Optional package for circuit-level checks:

```bash
pip install qibo
```

`qibo` is only needed for optional circuit checks in `hopf_utils.py` and for the standalone circuit-safeguard script `VQE_qibo.py`. The synthetic data-generation, diagnostic, CNOT-count, and plotting scripts do not require Qibo.


## Synthetic tasks

Each task is scrambled by a fixed real orthogonal circuit. This prevents the optimizers from exploiting an obvious computational-basis target.
### VQE tasks

**`VQE-1`: Parent Hamiltonian**

Minimizes the parent-Hamiltonian gap for a scrambled target state. This is the direct target-state recovery task.

**Objective:** minimize one minus the squared target overlap.

**Stress feature:** direct scrambled target-state recovery.

**`VQE-2`: Scrambled Hamming spectrum**

Minimizes a diagonal Hamming-distance spectrum after conjugation by the real scrambler.

**Objective:** minimize the scrambled Hamming-spectrum energy.

**Stress feature:** a structured diagonal spectrum hidden by a real orthogonal circuit.

**`VQE-3`: Small-gap spectrum**

Minimizes a scrambled diagonal Hamiltonian with one ground state, one nearby excited state, and all other levels at unit energy.

**Objective:** minimize the scrambled small-gap Hamiltonian energy.

**Stress feature:** a small spectral gap of `1e-2` near the optimum.

### Metrology-inspired tasks

**`MET-1`: Single-target Fisher**

Optimizes a nonlinear fixed-readout Fisher proxy for a scrambled target state.

**Objective:** minimize one minus the squared Fisher proxy.

**Stress feature:** nonlinear single-target fixed-readout objective.

**`MET-2`: QFI superposition**

Optimizes normalized pure-state QFI for a scrambled diagonal generator.

**Objective:** maximize normalized QFI, equivalently minimize the normalized-QFI gap.

**Stress feature:** the optimum is an equal superposition of the extremal generator eigenstates.

**`MET-3`: Balanced Fisher**

Optimizes a nonlinear two-target Fisher soft-min objective.

**Objective:** minimize the balanced two-target Fisher gap.

**Stress feature:** the optimum requires balanced overlap with two scrambled target states.

Lower cost or gap is better in all diagnostic summaries.

## Optimizer modes

### Geometry-native Hopf optimizers

Generated by `hopf_data.py`:

```text
Hopf-EGT-CG
Hopf-Riemannian-LBFGS
Hopf-Riemannian-BB
```

These modes use the same cost-plus-Hopf-coordinate-gradient interface. The Hopf diagonal metric is used to lift coordinate gradients to state-sphere gradients, and accepted states are mapped back to Hopf coordinates.

### Adam baselines

Generated by `adam_data.py`:

```text
Hopf-Adam
Mottonen-ideal-PS-Adam
```

`Hopf-Adam` applies Adam directly to Hopf coordinates.

`Mottonen-ideal-PS-Adam` applies Adam to the physical post-multiplexing Möttönen rotation angles with an exact gradient equivalent to infinite-shot parameter shift.

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

python diagnose_hopf.py --indir data --ns 8 --num-seeds 5
python diagnose_adam.py --indir data --ns 8 --num-seeds 5
```

## Diagnostics

Run diagnostics after generating the datasets.

```bash
mkdir -p diagnostics

python diagnose_hopf.py \
    --indir data \
    --ns 6-10 \
    --steps 200 \
    --num-seeds 10 \
    --out diagnostics/hopf_data_diagnostics.txt

python diagnose_adam.py \
    --indir data \
    --ns 6-10 \
    --steps 200 \
    --num-seeds 10 \
    --out diagnostics/adam_data_diagnostics.txt
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

## Optional Hopf utility checks

Run:

```bash
python hopf_utils.py
```

This executes consistency checks for inverse mapping, tangent-state synthesis, metric/Jacobian agreement, and optional Qibo circuit agreement when Qibo is installed.

## Safeguard scripts

The repository includes two standalone safeguard scripts. They are not part of the main synthetic-data pipeline, and their generated figure files are not committed.

### CNOT-count safeguard

Run:

```bash
python hopf_gate_count.py
```

This prints a table of real and complex Hopf CNOT counts for the default range `n = 4, ..., 20`, using the no-clean-ancilla CNOT model stated in the paper. It also saves a local plot named:

```text
hopf_gate_count.pdf
```

To choose a smaller range and write the local output into a generated-figure folder:

```bash
mkdir -p figures_safeguards

python hopf_gate_count.py \
    --nmin 2 \
    --nmax 10 \
    --out figures_safeguards/hopf_gate_count.pdf
```

The script compares counts obtained directly from the generated Hopf gate schedules against the closed-form binomial count formulas. The saved plot is only a local check and is not shown in this README.

### Qibo circuit-level gradient safeguard

Run:

```bash
MPLBACKEND=Agg python VQE_qibo.py
```

This script uses Qibo to build explicit ancilla Hadamard-test circuits for real-Hopf tangent-state transition moments of the form:

```text
Re <partial_i psi(theta)| P_alpha |psi(theta)>
```

inside a small TFIM VQE example. It compares a sampled circuit-gradient trajectory with an exact state-vector-gradient trajectory and writes:

```text
VQE_ADAM_MC_qibo.pdf
```

This output is a local circuit-realizability check and is not shown in this README.

## Reproducibility notes

The scripts use deterministic random seeds by default. Re-running the same commands with the same Python, NumPy, and SciPy versions should reproduce the same synthetic task definitions and initial-state seeds.

In Hopf CSVs, `last_line_evals` is a line-search work counter. It includes trial objective evaluations and, when strong-Wolfe curvature checks are reached, gradient-oracle calls used only by the line search. In Adam CSVs, `last_line_evals` counts adaptive cost-only trial evaluations.

The generated stress-test datasets use the real Hopf chart. `hopf_utils.py` also contains real and complex Hopf utilities used by the paper.


## Citation

If you use this code, cite the accompanying Hopf ansatz paper.
