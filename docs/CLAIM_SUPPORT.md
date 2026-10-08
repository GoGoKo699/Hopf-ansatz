# Claim support and evidence map

This page connects the results in *A Compass on the Quantum State Sphere: The
Hopf Ansatz for Arbitrary Pure-State Optimization* to the implementation and
numerical evidence. The paper provides the all-size statements and proofs;
the repository tests finite realizations.

Use the [engineering guide](ENGINEERING_GUIDE.md) for formulas, API contracts,
coordinate conventions, and compiler assumptions. Use
[Experiments](EXPERIMENTS.md) for study designs and reported results, and
[Reproducibility](../REPRODUCIBILITY.md) to run the checks.

## Claim-to-evidence map

| Paper-level statement | Repository support | Main files | Scope |
|---|---|---|---|
| The real and complex Hopf recursions produce normalized states with the stated coordinate order | Forward-map and round-trip checks | `hopf_utils.py` | Finite numerical instances; the all-size recursion is analytic. |
| Any normalized real or complex state has an explicit Hopf representative | Forward-after-inverse round-trip | `hopf_utils.py` | Zero subtrees and zero leaves use a fixed convention; coordinates are not globally unique. |
| The inverse map has linear classical time and memory | Tree-based implementation and paper algorithm | `hopf_utils.py` | Algorithmic time and memory complexity. |
| The pullback metric is diagonal | Dense Jacobian Gram matrix versus analytic diagonal entries | `hopf_utils.py`, `hopf_complex.py` | Checked at finite dimensions; complex geometry uses the round-sphere real pairing. |
| Normalized coordinate tangents can be prepared on the same Hopf skeleton | Tangent assignment versus normalized Jacobian columns | `hopf_utils.py` | Only regular coordinates have a normalized tangent; zero metric gives a vanishing raw differential. |
| The native gate schedules reproduce the analytic states | Generated schedules and optional Qibo parity | `hopf_utils.py` | Prepared state from the initialized input; full-unitary completions may differ. |
| Real and complex assigned CNOT formulas match the generated schedules | Schedule-by-schedule counting versus closed forms | `hopf_gate_count.py` | Assigned logical charges in the declared no-clean-ancilla model. |
| Hopf coordinate gradients can be evaluated by tree contractions | Fast routines compared with dense Jacobian references or exact identities | `hopf_data.py`, `hopf_complex.py`, `finite_shot_sanity_check.py` | Classical reference calculations; quantum access uses the paper’s stated assumptions. |
| Signed branch states recover coordinate gradients | Exact branch identity and finite-shot sampling | `finite_shot_sanity_check.py`, `VQE_qibo.py` | Fixed observable decompositions and stated shot convention. |
| Layerwise preparations organize a full gradient into logarithmically many compiled families | Explicit `n = 4` layerwise circuit safeguard and general layer specification | `VQE_qibo.py` | Local `n = 4` circuit validation; the all-size count is analytic. |
| Geometry-native optimizers use only cost and Hopf-coordinate-gradient queries | Executable EGT-CG, R-LBFGS, and R-BB pipelines | `hopf_data.py`, `hopf_complex.py` | Classical statevector studies with the stated cost/gradient interface. |
| The real-chart optimization study is stable across the stated task family | Multi-size deterministic generation and streaming diagnostics | `hopf_data.py`, `diagnose_hopf.py`, `plot_hopf.py` | Synthetic tasks at the specified sizes and seeds. |
| The complex extension is operational beyond isolated identities | Focused `n = 6` complex stress test and self-checks | `hopf_complex.py` | One size and six synthetic task families, with four Hopf optimizer modes. |
| Finite-shot error follows the expected inverse-square-root trend in the fixed test | Repeated sampling at three shot counts | `finite_shot_sanity_check.py` | Shots are allocated per branch and per component. |
| A local real/complex layerwise estimator tracks exact-gradient VQE trajectories | Explicit Qibo or equivalent ideal statevector simulation | `VQE_qibo.py` | Local `n = 4` Hamiltonians and fixed-rate Adam trajectories. |

## Reading the evidence

- **Core identities:** `hopf_utils.py` checks forward-after-inverse state
  reconstruction and normalized Jacobian columns against tangent preparations.
  State-level comparison accommodates nonunique boundary coordinates. Optional
  Qibo checks compare prepared statevectors with the analytic maps.
- **Gradient and geometry checks:** `finite_shot_sanity_check.py` compares an
  independent real tree gradient with a dense Jacobian and exact signed-branch
  energies. `hopf_complex.py` compares the fast complex gradient with
  `Re(J† h)` and the metric lift with the projected state-space gradient.
- **Resource accounting:** `hopf_gate_count.py` compares generated gate-list
  charges with closed-form schedule counts. The
  [assigned CNOT ledger](ENGINEERING_GUIDE.md#15-assigned-cnot-ledger) specifies
  the control-dependent model, including its readout and compilation boundary.
- **Sampling:** the fixed-state finite-shot experiment assigns `S` shots to
  each sign of each coordinate, hence `2*S` branch-state samples per coordinate.
  The local VQE script instead uses shots per label, sign, and Pauli readout.
  [Experiments](EXPERIMENTS.md) states each convention alongside its results.
- **Optimization:** generators record task, seed, mode, and step grids.
  Streaming diagnostics check completeness, nonfinite values, normalization,
  threshold hits, and final rankings; plotting scripts read regenerated CSVs.
  The metric division floor is a numerical safeguard near singular chart
  boundaries, while the analytic metric remains unchanged.

The companion [Hopf-QBP repository](https://github.com/GoGoKo699/Hopf-QBP)
contains exact-logical reverse-gradient constructions and their validation.
