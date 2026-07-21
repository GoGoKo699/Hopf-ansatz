# Hopf Ansatz

This repository accompanies two papers built on the same Hopf-coordinate construction.
The numerical optimization studies of the first paper and the exact-logical validation
of the second paper are kept in one repository but are separated by scope.

| Paper | Scope | Repository material |
|---|---|---|
| **[A Compass on the Quantum State Sphere: The Hopf Ansatz for Arbitrary Pure-State Optimization](https://arxiv.org/abs/2607.14231)** | Hopf coordinates, inverse maps, geometry, optimization algorithms, synthetic VQE and metrology-inspired stress tests | Existing root-level scripts, diagnostics, plots, and [`PAPER1_NUMERICS.md`](PAPER1_NUMERICS.md) |
| **The Compass in Reverse: Quantum Backpropagation with the Hopf Ansatz** | Exact-logical complete-gradient estimators, global frames, complex phase decoding, polyspherical frames, and checkpointed adjoints | Deterministic Qibo circuit and decoder checks in [`qbp_validation/`](qbp_validation/) |

## Scope separation

The first paper contains numerical optimization studies. Its scripts generate optimizer
traces, diagnostics, and figures.

The second paper is theoretical and makes no numerical-performance claim. Its repository
material consists only of deterministic sanity checks of the exact-logical state, frame,
measurement, decoder, singular-chart, checkpoint, and compiler-ledger identities. It does
not contain optimization comparisons, empirical scaling plots, hardware experiments, or
noise studies.

## Repository map

- [`PAPER1_NUMERICS.md`](PAPER1_NUMERICS.md): guide to the first paper's numerical and reproducibility material.
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md): commands for regenerating the first paper's datasets, diagnostics, safeguards, and figures.
- [`hopf_utils.py`](hopf_utils.py): shared real and complex Hopf maps, inverse maps, Jacobians, metrics, tangent-state assignments, and gate schedules.
- [`VQE_qibo.py`](VQE_qibo.py): the first paper's local real/complex VQE circuit-realizability demonstration.
- [`qbp_validation/`](qbp_validation/): exact-distribution Qibo validation for the second manuscript.

## Quick start

Use Python 3.10 or newer.

### Paper I: numerical safeguards

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-optional.txt
python hopf_utils.py
MPLBACKEND=Agg python VQE_qibo.py --steps 2 --shots 20 --sampler statevector --log-every 0
```

See [`PAPER1_NUMERICS.md`](PAPER1_NUMERICS.md) and
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the full data-generation and plotting workflow.

### Paper II: exact-logical Qibo validation

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-optional.txt
python -m unittest discover -s qbp_validation/tests -v
```

The Paper-II suite evaluates exact Qibo statevectors and complete measurement
probability distributions for small systems. It creates no plots or benchmark datasets.
See [`qbp_validation/README.md`](qbp_validation/README.md) for the validation matrix and
scope of each test group.
