# Reviewer route

This route connects **A Compass on the Quantum State Sphere: The Hopf Ansatz
for Arbitrary Pure-State Optimization**, [arXiv:2607.14231v1](https://arxiv.org/abs/2607.14231v1),
to its complete Markdown account and implementation. The [science entry](../science/README.md)
gives a connected reading route. The [complete paper map](PAPER_MAP.md) resolves
original section, equation, statement, algorithm, figure, table and reference
identifiers to their authoritative repository locations.

## A direct path through the argument

1. Read the [chart and metric conventions](../science/chart-and-geometry.md):
   normalized real and complex state vectors, tree indexing, angle ranges,
   coordinate singularities, and the retained complex global phase.
2. Check the [gradient construction](../science/gradient-access.md) and its
   [proofs](../science/proofs.md): metric factors, normalized tangents, signed
   branches, and indexed preparation. Keep the number of compiled circuit
   families separate from the precision-dependent measurement budget.
3. Inspect the [inverse and native schedule](../science/algorithms.md), then the
   [four-qubit worked example](../science/four-qubit-example.md). The
   [engineering guide](../docs/ENGINEERING_GUIDE.md) translates the notation
   into array, bit and physical-angle conventions.
4. Follow [geometry-native optimization](../science/optimizers.md) and
   [metrological chain-rule observables](../science/metrological-objectives.md)
   to the [numerical study](../science/numerical-studies.md). The
   [experiment specification](../docs/EXPERIMENTS.md) fixes tasks, seeds,
   outputs and the matched baseline comparisons.
5. Run the [identity checks and reproduction commands](../verification/README.md).
   Use the numerical output as evidence about the implementation and stated
   experiments; use the written arguments for all-size statements.

## Results, proofs and implementation

| Result and paper statement | Complete proof or construction | Implementation / check | Scope of the evidence |
|---|---|---|---|
| Normalized real and complex charts, Definitions 1–3 | [Chart definitions](../science/chart-and-geometry.md#Thmdefinition1); [product formula, Lemma 3](../science/proofs.md#Thmlemma3) | `vector_from_theta`, `gates_order`, optional Qibo parity in `hopf_utils.py` | Finite state reconstruction checks supplement the all-size recursion. Native schedules agree on the initialized state column; full-unitary completions may differ. |
| Inverse coordinates and linear classical cost, Lemma 1 | [Inverse map](../science/algorithms.md#Thmlemma1); [proof](../science/proofs.md#A5.SSx6) | `theta_from_vector` followed by `vector_from_theta`; `python hopf_utils.py` | The algorithm has `O(N)` time and memory. Zero subtrees and zero leaves use fixed coordinate conventions; the inverse is not globally unique. |
| Diagonal metric, Definition 4 and Theorem 1 | [Metric convention and statement](../science/chart-and-geometry.md#Thmdefinition4); [proof](../science/proofs.md#A5.SSx2) | [Executable Gram-matrix check](../verification/README.md#check-the-metric-against-the-jacobian) of `metric_diagonal` against `jacobian`; `hopf_complex.py` self-checks | Real Fubini–Study/sphere pairing and complex ambient round-sphere real pairing are distinct conventions. Numerical checks are finite-dimensional. |
| VQE gradient and signed-branch identity, Theorem 2 | [Statement](../science/gradient-access.md#Thmtheorem2); [proof](../science/proofs.md#A5.SSx3) | `prepare_fixed_experiment` in `finite_shot_sanity_check.py`; independent tree gradient, dense Jacobian and exact branch energies | Exact identities and finite-shot sampling are separate checks. Quantum access uses the stated objective and measurement assumptions. |
| Preparable normalized tangents, Definition 5 and Theorem 3 | [Definition and construction](../science/gradient-access.md#Thmdefinition5); [proof](../science/proofs.md#A5.SSx4) | `theta_hopf_tangent_state` against normalized `jacobian` columns; `python hopf_utils.py` | A normalized coordinate tangent requires nonzero metric weight; a zero metric gives a vanishing raw differential. |
| Layerwise families and uniform per-setting cost, Theorem 4 | [Statement](../science/gradient-access.md#Thmtheorem4); [proof](../science/proofs.md#A5.SSx5) | `VQE_qibo.py` implements the real/complex indexed-layer circuit check | Circuit validation is at `n = 4`; the all-size family count and cost follow from the proof. Compiled families and precision-dependent shot budgets remain separate resources. |
| Assigned native CNOT counts, Lemma 2 and Theorem 5 | [Mixed-sign control model](../science/algorithms.md#Thmlemma2); [count formula](../science/algorithms.md#Thmtheorem5) | `count_for_n` against `G_real_theory` / `G_complex_theory` in `hopf_gate_count.py`; [expected integers](../verification/README.md#what-to-inspect-in-the-output) | Exact schedule counts under the declared logical no-clean-ancilla model, with the readout and compilation boundary in the [engineering ledger](../docs/ENGINEERING_GUIDE.md#15-assigned-cnot-ledger). |
| Jacobian and Gram forms, Theorem 6 | [Statement](../science/optimizers.md#Thmtheorem6); [proof](../science/proofs.md#A5.SSx7) | `jacobian`, `metric_diagonal`; complex fast gradient versus `Re(J† h)` and metric lift versus projected state-space gradient | Classical reference calculations verify the finite formulas and preserve the chosen sphere metric. |
| Geometry-native optimization, Appendix C | [Optimizer construction](../science/optimizers.md) | EGT-CG, R-LBFGS and R-BB in `hopf_data.py` and `hopf_complex.py` | The optimizer uses cost and Hopf-coordinate-gradient queries plus classical geometric processing. Experiments are statevector studies with the stated interface. |
| QFI/CFI chain-rule observables, Appendix D | [Metrological derivations](../science/metrological-objectives.md) | Compare the displayed chain rule and local-estimation assumptions with the objective definitions | The synthetic metrology-inspired tasks below have their own explicit definitions; their numerical evidence is not a replacement for the QFI/CFI derivations. |

## Numerical and worked-example evidence

| Material | Where to inspect | Distinct evidence and reproduction route |
|---|---|---|
| Algorithms 1–2 and Subroutines 1–5 | [Algorithms](../science/algorithms.md) | Compare the complete pseudocode with `gates_order`, `hamming_weight`, `find_pairs`, `anti_ctrl`, `theta_real`, and `theta_complex` in `hopf_utils.py`. |
| Four-qubit construction; Figures 3–6 and Tables I–III | [Worked example](../science/four-qubit-example.md) | Inspect explicit gate, metric and tangent rows and the indexed circuits against the core functions and `VQE_qibo.py`. |
| Real optimization study; Section IV and Figures 1–2 | [Numerical studies](../science/numerical-studies.md) | Matched real-Hopf and Adam pipelines use `n = 6,…,10`, six tasks and ten seeds. Use the [archived numerical implementation](../verification/README.md#source-revision-for-the-reported-numerical-study), generate CSVs, run streaming diagnostics, and reproduce the two plots with the [reported aggregate targets](../docs/EXPERIMENTS.md#paper-reported-aggregate-reproduction-targets). |
| Complex extension | [Focused study](../docs/EXPERIMENTS.md#8-focused-complex-hopf-study) | `hopf_complex.py` checks one size, `n = 6`, six synthetic task families and four optimizer modes, with completeness and finite-value diagnostics. |
| Fixed-state finite-shot estimator | [Sampling experiment](../docs/EXPERIMENTS.md#9-finite-shot-signed-branch-check) | `finite_shot_sanity_check.py` assigns `S` samples to each sign of each coordinate: `2*S` branch-state samples per coordinate. The default three-shot-count study checks the inverse-square-root error trend at that fixed state. |
| Local layerwise VQE | [Circuit experiment](../docs/EXPERIMENTS.md#10-local-qibo-layerwise-safeguard) | Explicit Qibo or ideal statevector simulation at `n = 4`, two Hamiltonians and fixed-rate Adam trajectories. Shots are per label, sign and Pauli readout. |

[Verification and reproduction](../verification/README.md) owns the commands,
expected outputs and tolerances. The streaming diagnostics check task, seed,
mode and step completeness, finite values, state normalization, threshold hits
and final rankings. Metric division floors are numerical safeguards near
singular coordinates; they do not change the analytic metric.

The companion [Hopf-QBP repository](https://github.com/GoGoKo699/Hopf-QBP)
contains the complete reverse-gradient construction and its validation.
