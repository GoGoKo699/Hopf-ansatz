# The Hopf Ansatz: reference implementation and reproducibility

Reference code, engineering conventions, validation safeguards, and reproducible
experiments for:

**A Compass on the Quantum State Sphere: The Hopf Ansatz for Arbitrary
Pure-State Optimization**  
Ruge Lin and Guangxi Li (2026)

[Read the paper on arXiv](https://arxiv.org/abs/2607.14231)

The Hopf ansatz is a fixed binary-tree circuit family for arbitrary normalized
real and complex pure states. In addition to forward state preparation, it
provides an explicit inverse map, a diagonal pullback metric, exactly preparable
normalized coordinate tangents, structured gradient access, and
geometry-aware optimization on the state sphere.

## Choose a route

| Goal | Start here |
|---|---|
| Check whether the repository supports a paper claim | [Claim support and evidence map](docs/CLAIM_SUPPORT.md) |
| Implement the chart, inverse map, metric, tangents, or gate schedule | [Engineering guide](docs/ENGINEERING_GUIDE.md) |
| Understand the numerical studies and safeguards | [Experiments and evidence](docs/EXPERIMENTS.md) |
| Reproduce the scripts and generated data | [Reproducibility](REPRODUCIBILITY.md) |
| Inspect the companion reverse-gradient constructions | [Hopf-QBP repository](https://github.com/GoGoKo699/Hopf-QBP) |

## Implementation

For `n` qubits and `N = 2**n`, the real chart uses `N - 1` angles; the
complex chart adds `N` leaf phases. The core API in `hopf_utils.py` provides
forward and inverse maps, exact Jacobians, diagonal metrics, normalized tangent
assignments, and deterministic real/complex gate schedules with optional Qibo
circuit checks.

The optimization scripts combine boundary-safe coordinate gradients with
state-sphere geometry. The experiments cover real and complex synthetic tasks,
a fixed-state finite-shot estimator, and a local `n = 4` layerwise circuit
check. The resource ledger assigns logical CNOT costs under a declared
no-clean-ancilla model.

Coordinate conventions, singular-boundary handling, and implementation examples
are in the [engineering guide](docs/ENGINEERING_GUIDE.md). The
[claim map](docs/CLAIM_SUPPORT.md) connects each result to its checks and scope;
the paper supplies the formal proofs.

## Quick start

Use Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Install Qibo only for explicit circuit checks:

```bash
python -m pip install -r requirements-optional.txt
```

Run the core map, inverse, metric, tangent, and optional Qibo checks:

```bash
python hopf_utils.py
```

For the small checks, full experiments, diagnostics, and figure generation, use
[REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Evidence snapshots

The committed figures summarize the complex optimization study, fixed-state
finite-shot estimator, and local layerwise circuit check.

<p align="center">
  <img src="hopf_complex.png" width="31%" alt="Complex Hopf optimization stress-test summary">
  <img src="finite_shot_sanity.png" width="31%" alt="Finite-shot signed-branch gradient error">
  <img src="VQE_qibo.png" width="31%" alt="Local Qibo layerwise-gradient VQE safeguard">
</p>

- `hopf_complex.png` summarizes a focused complex-chart optimization study.
- `finite_shot_sanity.png` isolates the statistical convergence of one signed-branch estimator.
- `VQE_qibo.png` checks local real and complex layerwise circuit realizability at `n = 4`.

The exact meaning, settings, and limitations of each panel are documented in
[docs/EXPERIMENTS.md](docs/EXPERIMENTS.md).

## Repository map

| Path | Role |
|---|---|
| `hopf_utils.py` | Coordinate maps, inverse maps, Jacobians, metrics, tangent assignments, native gate schedules, and optional Qibo circuit checks. |
| `hopf_gate_count.py` | Generated-schedule versus closed-form assigned CNOT safeguard. |
| `hopf_data.py` | Multi-size real-Hopf geometry-native optimization data generator. |
| `adam_data.py` | Real Hopf-Adam and ideal Möttönen-parameter-shift-equivalent baselines. |
| `diagnose_hopf.py` | Streaming completeness and numerical diagnostics for geometry-native CSV data. |
| `diagnose_adam.py` | Streaming diagnostics for the Adam baseline CSV data. |
| `plot_hopf.py` | Aggregate and per-task plots from regenerated real-Hopf and Adam datasets. |
| `hopf_complex.py` | Focused complex-Hopf stress test, diagnostics, and committed summary figure. |
| `finite_shot_sanity_check.py` | Fixed-state signed-branch estimator check and committed statistical figure. |
| `VQE_qibo.py` | Local real/complex layerwise-gradient circuit safeguard and committed VQE figure. |
| `docs/CLAIM_SUPPORT.md` | Claim-to-code and claim-to-evidence map. |
| `docs/ENGINEERING_GUIDE.md` | Self-contained implementation and adaptation guide. |
| `docs/EXPERIMENTS.md` | Experimental designs, outputs, reported results, and interpretation limits. |
| `REPRODUCIBILITY.md` | Clean-environment commands and complete data-regeneration workflow. |

## Related work

The companion [Hopf-QBP repository](https://github.com/GoGoKo699/Hopf-QBP)
provides exact-logical global-frame, direct-phase, and checkpointed
reverse-gradient constructions. The two repositories have no runtime dependency
on one another.

## Citation

For the Hopf ansatz and the scientific results in this repository, cite:

**Ruge Lin and Guangxi Li, “A Compass on the Quantum State Sphere: The Hopf
Ansatz for Arbitrary Pure-State Optimization,” arXiv:2607.14231 (2026).**

```bibtex
@article{lin2026hopf,
  title   = {A Compass on the Quantum State Sphere:
             The Hopf Ansatz for Arbitrary Pure-State Optimization},
  author  = {Lin, Ruge and Li, Guangxi},
  journal = {arXiv preprint arXiv:2607.14231},
  year    = {2026},
  url     = {https://arxiv.org/abs/2607.14231}
}
```

## License

This software is released under the [MIT License](LICENSE).
