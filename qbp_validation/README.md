# Exact-logical Qibo validation for Hopf quantum backpropagation

This directory contains deterministic Qibo circuit-level validation for
*The Compass in Reverse: Quantum Backpropagation with the Hopf Ansatz*.
The manuscript is theoretical and makes no numerical-performance claim. The suite checks
exact-logical state, frame, measurement, and decoder identities using small statevector
circuits and independent analytic references.

It covers the real global frame, complex magnitude and phase families, real and complex
magnitude checkpoints, singular chart cases, representative polyspherical trees, the
complete four-qubit construction, and the assigned compiler-ledger arithmetic.

It does **not** benchmark optimization, finite-shot convergence, execution time, memory,
routing, approximate synthesis, device noise, or hardware performance. Compiler-ledger
arithmetic is checked separately and is not inferred from Qibo-transpiled gate counts.

## What counts as circuit-level validation

A circuit check in this directory satisfies both conditions:

1. the tested state or measurement distribution is produced by a
   `qibo.models.Circuit` executed with Qibo's NumPy statevector backend;
2. the expected result is computed independently by the recursive NumPy formulas in
   `reference.py`, which contains no Qibo imports.

The tests do not inject an expected statevector into Qibo. They build the Hopf frame and
breadth-first checkpoint circuits from controlled `RY(2 theta)` gates with explicit open
controls. The complex diagonal phase layer and phase-calibrated controlled observable are
represented as exact logical `Unitary` gates, matching the manuscript's access model.
Representative polyspherical pair rotations are also exact two-level `Unitary` gates; those
checks validate the finite-dimensional frame algebra, not an arbitrary-tree synthesis claim.

All estimator checks use the full output statevector after the stated basis rotations and
sum the complete probability distribution exactly. No Monte Carlo shots are used.

## Layout

| File | Role |
|---|---|
| `conventions.py` | Big-endian basis labels, node/marker maps, gate specifications, polyspherical topology, and assigned ledger formulas. |
| `reference.py` | Independent recursive states, frame columns, metric factors, coordinate derivatives, and exact gradients. |
| `circuits.py` | Qibo builders for real frames, complex phase layers, global estimators, checkpoints, and polyspherical frames. |
| `decoders.py` | Walsh, signed-histogram, signed one-hot, and checkpoint decoders. |
| `cases.py` | Deterministic interior, singular, observable, and tree cases. |
| `tests/` | Circuit, decoder, singular-case, four-qubit, polyspherical, and ledger checks. |

## Validation matrix

| Test group | Manuscript statement checked | Qibo calculation | Independent reference |
|---|---|---|---|
| `test_conventions.py` | Basis ordering, marker labels, open controls, and the corrected `S^dagger`-then-`H` Y readout | Small deterministic circuits | Explicit bits and matrices |
| `test_real_frame.py` | `W_R |0> = |psi>` and `W_R |lambda(j)> = |e_j>` including singular completions | Qibo circuit unitary and statevector | Recursive state and local complements |
| `test_real_global_estimator.py` | One all-X distribution returns every real magnitude derivative | Exact Qibo probabilities | Analytic coordinate gradient |
| `test_complex_magnitude.py` | Phase-dressed frame columns and the asymmetric magnitude decoder | Qibo frame/phase circuits and probabilities | Complex magnitude derivatives |
| `test_complex_phase.py` | Ancilla-Y/system-Z signed one-hot phase estimator and common-phase cancellation | Exact Qibo probabilities | Analytic leaf-phase derivatives |
| `test_checkpoints.py` | Real and complex signed one-hot magnitude checkpoints at every depth | Qibo Y/Z/Y distributions | Analytic layer gradients |
| `test_singular_cases.py` | Zero metric factors give zero coordinate derivatives without breaking the frame | Boundary-angle circuits | Recursive zero derivatives |
| `test_polyspherical.py` | Local-complement frame, diagonal metric, translation, and parity estimator | Qibo two-level circuits | Independent ordered-tree recursion |
| `test_decoders.py` | FWHT, signed leaf histogram, suffix marginalization, and one-hot decoding | Deterministic probability arrays | Direct sums |
| `test_four_qubit_example.py` | All appendix markers, node-5 formulas, and all four checkpoint settings | Fixed `n=4` circuits | Displayed manuscript formulas |
| `test_resource_ledger.py` | Assigned controlled-rotation, frame, depth, and checkpoint charges | No transpiler count is used | Exact integer formulas |

## Coverage

The balanced Hopf checks use `n = 1, 2, 3, 4`. They include:

- deterministic generic interior parameters;
- real final-layer angles beyond `pi/2`, exercising the sign convention;
- zero and `pi/2` upstream angles, producing zero metric factors;
- complex states with zero-amplitude leaves;
- Pauli, diagonal-reflection, and fixed-seed Householder observables satisfying
  `O = O^dagger` and `O^2 = I`.

The polyspherical checks use two ordered full binary trees on eight encoded leaves, including
an unbalanced tree with a nonzero root anchor so that the theorem's Pauli-X translation is
exercised explicitly.

## Run

From the repository root:

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-optional.txt
python -m unittest discover -s qbp_validation/tests -v
```

The suite was developed and checked with Qibo 0.3.4 and its NumPy backend. A separate pinned
environment, validation CLI, and GitHub Actions workflow are intentionally deferred.

A successful run executes the complete deterministic test suite and creates no figures,
CSVs, or benchmark summaries.

## Interpretation

The tests validate finite-size instances of the exact algebraic premises used in the
manuscript. They do not numerically prove concentration inequalities or asymptotic resource
statements. In particular:

- Hoeffding and fixed-norm concentration results remain mathematical theorems;
- decoder complexity and storage claims are not timed empirically;
- the polyspherical resource extension remains conditional on the stated efficient-frame
  realization assumptions;
- Qibo gate counts are not compared with the manuscript's assigned compiler ledger.
