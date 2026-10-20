# Complete paper map

This index maps [2607.14231v1](https://arxiv.org/abs/2607.14231v1) to its complete Markdown account. Labels follow the printed paper; links retain stable source identifiers. The [reviewer route](README.md) connects the principal claims to proofs and executable checks. The same inventory is available as [structured data](paper-map.json).

| Item type | Count |
|---|---:|
| Section | 48 |
| Statement | 14 |
| Equation | 121 |
| Figure | 6 |
| Table | 3 |
| Algorithm | 7 |
| Reference | 56 |

## Sections

| Paper item | Markdown location |
|---|---|
| Abstract | [README.md](../science/README.md#abstract1) |
| I Introduction | [motivation.md](../science/motivation.md#S1) |
| II A universal Hopf variational chart | [chart-and-geometry.md](../science/chart-and-geometry.md#S2) |
| Hopf geometry on the state sphere | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.SSx1) |
| III Hopf gradient access and scaling | [gradient-access.md](../science/gradient-access.md#S3) |
| Gradient estimation from signed branch states | [gradient-access.md](../science/gradient-access.md#S3.SSx1) |
| Magnitude-layer batching | [gradient-access.md](../science/gradient-access.md#S3.SSx2) |
| Single phase-layer batching for the complex ansatz | [gradient-access.md](../science/gradient-access.md#S3.SSx3) |
| Complete estimator and scaling analysis | [gradient-access.md](../science/gradient-access.md#S3.SSx4) |
| IV Numerical simulations | [numerical-studies.md](../science/numerical-studies.md#S4) |
| Stress-test design | [numerical-studies.md](../science/numerical-studies.md#S4.SSx1) |
| Aggregate results | [numerical-studies.md](../science/numerical-studies.md#S4.SSx2) |
| Per-task convergence | [numerical-studies.md](../science/numerical-studies.md#S4.SSx3) |
| V Conclusion | [conclusion.md](../science/conclusion.md#S5) |
| Data and code availability | [conclusion.md](../science/conclusion.md#Sx1) |
| References | [paper-bibliography.md](../references/paper-bibliography.md#bib) |
| Appendix A Algorithms and subroutines | [algorithms.md](../science/algorithms.md#A1) |
| Inverse Hopf map | [algorithms.md](../science/algorithms.md#A1.SSx1) |
| Gate-schedule routines | [algorithms.md](../science/algorithms.md#A1.SSx2) |
| CNOT counting model | [algorithms.md](../science/algorithms.md#A1.SSx3) |
| Appendix B A four-qubit example of the Hopf ansatz | [four-qubit-example.md](../science/four-qubit-example.md#A2) |
| Hopf tree and leaf-phase indexing | [four-qubit-example.md](../science/four-qubit-example.md#A2.SSx1) |
| Ansatz circuit | [four-qubit-example.md](../science/four-qubit-example.md#A2.SSx2) |
| Metric and Jacobian entries | [four-qubit-example.md](../science/four-qubit-example.md#A2.SSx3) |
| Representative normalized coordinate tangents | [four-qubit-example.md](../science/four-qubit-example.md#A2.SSx4) |
| Layerwise derivative and branch circuits | [four-qubit-example.md](../science/four-qubit-example.md#A2.SSx5) |
| Appendix C Geometry-native optimizers for the Hopf ansatz | [optimizers.md](../science/optimizers.md#A3) |
| Jacobian formulas and inverse differential map | [optimizers.md](../science/optimizers.md#A3.SSx1) |
| State-sphere conventions | [optimizers.md](../science/optimizers.md#A3.SSx2) |
| Lifting Hopf coordinate gradients to the state sphere | [optimizers.md](../science/optimizers.md#A3.SSx3) |
| Hopf exact geodesic transport with conjugate gradients | [optimizers.md](../science/optimizers.md#A3.SSx4) |
| Hopf Riemannian L-BFGS | [optimizers.md](../science/optimizers.md#A3.SSx5) |
| Hopf Riemannian Barzilai–Borwein | [optimizers.md](../science/optimizers.md#A3.SSx6) |
| Summary of the optimizer family | [optimizers.md](../science/optimizers.md#A3.SSx7) |
| Appendix D Metrological objectives as chain-rule observables | [metrological-objectives.md](../science/metrological-objectives.md#A4) |
| A general chain-rule form | [metrological-objectives.md](../science/metrological-objectives.md#A4.SSx1) |
| Pure-state QFI objective | [metrological-objectives.md](../science/metrological-objectives.md#A4.SSx2) |
| Fixed-readout Ramsey CFI objective | [metrological-objectives.md](../science/metrological-objectives.md#A4.SSx3) |
| Compiled-setting scaling | [metrological-objectives.md](../science/metrological-objectives.md#A4.SSx4) |
| Appendix E Proofs of the main results | [proofs.md](../science/proofs.md#A5) |
| Closed-form product formula for Hopf amplitudes | [proofs.md](../science/proofs.md#A5.SSx1) |
| Proof of Theorem 1 | [proofs.md](../science/proofs.md#A5.SSx2) |
| Proof of Theorem 2 | [proofs.md](../science/proofs.md#A5.SSx3) |
| Proof of Theorem 3 | [proofs.md](../science/proofs.md#A5.SSx4) |
| Proof of Theorem 4 | [proofs.md](../science/proofs.md#A5.SSx5) |
| Proof of Lemma 1 | [proofs.md](../science/proofs.md#A5.SSx6) |
| Proof of Theorem 6 | [proofs.md](../science/proofs.md#A5.SSx7) |
| Acknowledgements | [conclusion.md](../science/conclusion.md#S5.acknowledgements1) |

## Statements

| Paper item | Markdown location |
|---|---|
| Definition 1 (Hopf binary tree) . | [chart-and-geometry.md](../science/chart-and-geometry.md#Thmdefinition1) |
| Definition 2 (Hopf coordinate system) . | [chart-and-geometry.md](../science/chart-and-geometry.md#Thmdefinition2) |
| Definition 3 (Hopf ansatz) . | [chart-and-geometry.md](../science/chart-and-geometry.md#Thmdefinition3) |
| Definition 4 (Metrics in parameter space) . | [chart-and-geometry.md](../science/chart-and-geometry.md#Thmdefinition4) |
| Theorem 1 (Diagonal metric of the Hopf ansatz) . | [chart-and-geometry.md](../science/chart-and-geometry.md#Thmtheorem1) |
| Definition 5 (Normalized coordinate tangent) . | [gradient-access.md](../science/gradient-access.md#Thmdefinition5) |
| Theorem 2 (Gradient estimation for VQE) . | [gradient-access.md](../science/gradient-access.md#Thmtheorem2) |
| Theorem 3 (Normalized coordinate tangents as quantum states (Hopf)) . | [gradient-access.md](../science/gradient-access.md#Thmtheorem3) |
| Theorem 4 (Uniform per-setting cost of Hopf gradient configurations) . | [gradient-access.md](../science/gradient-access.md#Thmtheorem4) |
| Lemma 1 (Inverse Hopf map and classical cost) . | [algorithms.md](../science/algorithms.md#Thmlemma1) |
| Lemma 2 (CNOT counting model for mixed-sign multi-controls) . | [algorithms.md](../science/algorithms.md#Thmlemma2) |
| Theorem 5 (CNOT counts for Hopf ansatze under the piecewise no-ancilla model) . | [algorithms.md](../science/algorithms.md#Thmtheorem5) |
| Theorem 6 (Jacobian formulas and Gram forms) . | [optimizers.md](../science/optimizers.md#Thmtheorem6) |
| Lemma 3 (Closed-form product formula for Hopf amplitudes) . | [proofs.md](../science/proofs.md#Thmlemma3) |

## Algorithms

| Paper item | Markdown location |
|---|---|
| Algorithm 1 | [algorithms.md](../science/algorithms.md#A1.SSx2.p3.3) |
| Subroutine 1 | [algorithms.md](../science/algorithms.md#A1.SSx2.p5.1) |
| Subroutine 2 | [algorithms.md](../science/algorithms.md#A1.SSx2.p6.3) |
| Subroutine 3 | [algorithms.md](../science/algorithms.md#A1.SSx2.p7.3) |
| Subroutine 4 | [algorithms.md](../science/algorithms.md#A1.SSx2.p8.3) |
| Subroutine 5 | [algorithms.md](../science/algorithms.md#A1.SSx2.p9.3) |
| Algorithm 2 | [algorithms.md](../science/algorithms.md#A1.SSx2.p10.3) |

## Equations

| Paper item | Markdown location |
|---|---|
| Equation (1) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E1) |
| Equation (2) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E2) |
| Equation (3) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E3) |
| Equation (4) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E4) |
| Equation (5) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E5) |
| Equation (6) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E6) |
| Equation (7) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E7) |
| Equation (8) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E8) |
| Equation (9) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E9) |
| Equation (10) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E10) |
| Equation (11) | [chart-and-geometry.md](../science/chart-and-geometry.md#S2.E11) |
| Equation (12) | [gradient-access.md](../science/gradient-access.md#S3.E12) |
| Equation (13) | [gradient-access.md](../science/gradient-access.md#S3.E13) |
| Equation (14) | [gradient-access.md](../science/gradient-access.md#S3.E14) |
| Equation (15) | [gradient-access.md](../science/gradient-access.md#S3.E15) |
| Equation (16) | [gradient-access.md](../science/gradient-access.md#S3.E16) |
| Equation (17) | [gradient-access.md](../science/gradient-access.md#S3.E17) |
| Equation (18) | [gradient-access.md](../science/gradient-access.md#S3.E18) |
| Equation (19) | [gradient-access.md](../science/gradient-access.md#S3.E19) |
| Equation (20) | [gradient-access.md](../science/gradient-access.md#S3.E20) |
| Equation (21) | [gradient-access.md](../science/gradient-access.md#S3.E21) |
| Equation (22) | [gradient-access.md](../science/gradient-access.md#S3.E22) |
| Equation (23) | [gradient-access.md](../science/gradient-access.md#S3.E23) |
| Equation (24) | [gradient-access.md](../science/gradient-access.md#S3.E24) |
| Equation (25) | [gradient-access.md](../science/gradient-access.md#S3.E25) |
| Equation (26) | [gradient-access.md](../science/gradient-access.md#S3.E26) |
| Equation (27) | [gradient-access.md](../science/gradient-access.md#S3.E27) |
| Equation (28) | [gradient-access.md](../science/gradient-access.md#S3.E28) |
| Equation (29) | [gradient-access.md](../science/gradient-access.md#S3.E29) |
| Equation (30) | [gradient-access.md](../science/gradient-access.md#S3.E30) |
| Equation (31) | [gradient-access.md](../science/gradient-access.md#S3.E31) |
| Equation (32) | [gradient-access.md](../science/gradient-access.md#S3.E32) |
| Equation (33) | [gradient-access.md](../science/gradient-access.md#S3.E33) |
| Equation (34) | [gradient-access.md](../science/gradient-access.md#S3.E34) |
| Equation (35) | [gradient-access.md](../science/gradient-access.md#S3.E35) |
| Equation (36) | [gradient-access.md](../science/gradient-access.md#S3.E36) |
| Equation (37) | [gradient-access.md](../science/gradient-access.md#S3.E37) |
| Equation (38) | [gradient-access.md](../science/gradient-access.md#S3.E38) |
| Equation (39) | [gradient-access.md](../science/gradient-access.md#S3.E39) |
| Equation (40) | [gradient-access.md](../science/gradient-access.md#S3.E40) |
| Equation (41) | [gradient-access.md](../science/gradient-access.md#S3.E41) |
| Equation (42) | [gradient-access.md](../science/gradient-access.md#S3.E42) |
| Equation (A1) | [algorithms.md](../science/algorithms.md#A1.E43) |
| Equation (A2) | [algorithms.md](../science/algorithms.md#A1.E44) |
| Equation (A3) | [algorithms.md](../science/algorithms.md#A1.E45) |
| Equation (A4) | [algorithms.md](../science/algorithms.md#A1.E46) |
| Equation (A5) | [algorithms.md](../science/algorithms.md#A1.E47) |
| Equation (B1) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E48) |
| Equation (B2) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E49) |
| Equation (B3) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E50) |
| Equation (B4) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E51) |
| Equation (B5) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E52) |
| Equation (B6) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E53) |
| Equation (B7) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E54) |
| Equation (B8) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E55) |
| Equation (B9) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E56) |
| Equation (B10) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E57) |
| Equation (B11) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E58) |
| Equation (B12) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E59) |
| Equation (B13) | [four-qubit-example.md](../science/four-qubit-example.md#A2.E60) |
| Equation (C1) | [optimizers.md](../science/optimizers.md#A3.E61) |
| Equation (C2) | [optimizers.md](../science/optimizers.md#A3.E62) |
| Equation (C3) | [optimizers.md](../science/optimizers.md#A3.E63) |
| Equation (C4) | [optimizers.md](../science/optimizers.md#A3.E64) |
| Equation (C5) | [optimizers.md](../science/optimizers.md#A3.E65) |
| Equation (C6) | [optimizers.md](../science/optimizers.md#A3.E66) |
| Equation (C7) | [optimizers.md](../science/optimizers.md#A3.E67) |
| Equation (C8) | [optimizers.md](../science/optimizers.md#A3.E68) |
| Equation (C9) | [optimizers.md](../science/optimizers.md#A3.E69) |
| Equation (C10) | [optimizers.md](../science/optimizers.md#A3.E70) |
| Equation (C11) | [optimizers.md](../science/optimizers.md#A3.E71) |
| Equation (C12) | [optimizers.md](../science/optimizers.md#A3.E72) |
| Equation (C13) | [optimizers.md](../science/optimizers.md#A3.E73) |
| Equation (C14) | [optimizers.md](../science/optimizers.md#A3.E74) |
| Equation (C15) | [optimizers.md](../science/optimizers.md#A3.E75) |
| Equation (C16) | [optimizers.md](../science/optimizers.md#A3.E76) |
| Equation (C17) | [optimizers.md](../science/optimizers.md#A3.E77) |
| Equation (C18) | [optimizers.md](../science/optimizers.md#A3.E78) |
| Equation (C19) | [optimizers.md](../science/optimizers.md#A3.E79) |
| Equation (C20) | [optimizers.md](../science/optimizers.md#A3.E80) |
| Equation (C21) | [optimizers.md](../science/optimizers.md#A3.E81) |
| Equation (C22) | [optimizers.md](../science/optimizers.md#A3.E82) |
| Equation (C23) | [optimizers.md](../science/optimizers.md#A3.E83) |
| Equation (C24) | [optimizers.md](../science/optimizers.md#A3.E84) |
| Equation (C25) | [optimizers.md](../science/optimizers.md#A3.E85) |
| Equation (C26) | [optimizers.md](../science/optimizers.md#A3.E86) |
| Equation (C27) | [optimizers.md](../science/optimizers.md#A3.E87) |
| Equation (C28) | [optimizers.md](../science/optimizers.md#A3.E88) |
| Equation (C29) | [optimizers.md](../science/optimizers.md#A3.E89) |
| Equation (C30) | [optimizers.md](../science/optimizers.md#A3.E90) |
| Equation (C31) | [optimizers.md](../science/optimizers.md#A3.E91) |
| Equation (C32) | [optimizers.md](../science/optimizers.md#A3.E92) |
| Equation (C33) | [optimizers.md](../science/optimizers.md#A3.E93) |
| Equation (C34) | [optimizers.md](../science/optimizers.md#A3.E94) |
| Equation (C35) | [optimizers.md](../science/optimizers.md#A3.E95) |
| Equation (C36) | [optimizers.md](../science/optimizers.md#A3.E96) |
| Equation (C37) | [optimizers.md](../science/optimizers.md#A3.E97) |
| Equation (C38) | [optimizers.md](../science/optimizers.md#A3.E98) |
| Equation (C39) | [optimizers.md](../science/optimizers.md#A3.E99) |
| Equation (C40) | [optimizers.md](../science/optimizers.md#A3.E100) |
| Equation (C41) | [optimizers.md](../science/optimizers.md#A3.E101) |
| Equation (D1) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E102) |
| Equation (D2) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E103) |
| Equation (D3) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E104) |
| Equation (D4) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E105) |
| Equation (D5) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E106) |
| Equation (D6) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E107) |
| Equation (D7) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E108) |
| Equation (D8) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E109) |
| Equation (D9) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E110) |
| Equation (D10) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E111) |
| Equation (D11) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E112) |
| Equation (D12) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E113) |
| Equation (D13) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E114) |
| Equation (D14) | [metrological-objectives.md](../science/metrological-objectives.md#A4.E115) |
| Equation (E1) | [proofs.md](../science/proofs.md#A5.E116) |
| Equation (E2) | [proofs.md](../science/proofs.md#A5.E117) |
| Equation (E3) | [proofs.md](../science/proofs.md#A5.E118) |
| Equation (E4) | [proofs.md](../science/proofs.md#A5.E119) |
| Equation (E5) | [proofs.md](../science/proofs.md#A5.E120) |
| Equation (E6) | [proofs.md](../science/proofs.md#A5.E121) |

## Figures

| Paper item | Markdown location |
|---|---|
| Figure 1 | [numerical-studies.md](../science/numerical-studies.md#S4.F1) |
| Figure 2 | [numerical-studies.md](../science/numerical-studies.md#S4.F2) |
| Figure 3 | [four-qubit-example.md](../science/four-qubit-example.md#A2.F3) |
| Figure 4 | [four-qubit-example.md](../science/four-qubit-example.md#A2.F4) |
| Figure 5 | [four-qubit-example.md](../science/four-qubit-example.md#A2.F5) |
| Figure 6 | [four-qubit-example.md](../science/four-qubit-example.md#A2.F6) |

## Tables

| Paper item | Markdown location |
|---|---|
| Table I | [four-qubit-example.md](../science/four-qubit-example.md#A2.T1) |
| Table II | [four-qubit-example.md](../science/four-qubit-example.md#A2.T2) |
| Table III | [four-qubit-example.md](../science/four-qubit-example.md#A2.T3) |

## References

| Paper item | Markdown location |
|---|---|
| Reference [1] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib16) |
| Reference [2] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib17) |
| Reference [3] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib18) |
| Reference [4] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib39) |
| Reference [5] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib31) |
| Reference [6] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib25) |
| Reference [7] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib26) |
| Reference [8] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib27) |
| Reference [9] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib28) |
| Reference [10] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib29) |
| Reference [11] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib14) |
| Reference [12] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib12) |
| Reference [13] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib11) |
| Reference [14] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib13) |
| Reference [15] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib15) |
| Reference [16] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib5) |
| Reference [17] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib32) |
| Reference [18] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib55) |
| Reference [19] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib56) |
| Reference [20] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib1) |
| Reference [21] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib4) |
| Reference [22] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib3) |
| Reference [23] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib6) |
| Reference [24] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib2) |
| Reference [25] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib19) |
| Reference [26] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib20) |
| Reference [27] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib21) |
| Reference [28] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib7) |
| Reference [29] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib22) |
| Reference [30] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib24) |
| Reference [31] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib23) |
| Reference [32] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib33) |
| Reference [33] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib34) |
| Reference [34] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib35) |
| Reference [35] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib36) |
| Reference [36] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib37) |
| Reference [37] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib38) |
| Reference [38] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib51) |
| Reference [39] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib49) |
| Reference [40] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib50) |
| Reference [41] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib54) |
| Reference [42] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib42) |
| Reference [43] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib30) |
| Reference [44] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib41) |
| Reference [45] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib40) |
| Reference [46] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib44) |
| Reference [47] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib45) |
| Reference [48] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib43) |
| Reference [49] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib46) |
| Reference [50] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib53) |
| Reference [51] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib10) |
| Reference [52] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib9) |
| Reference [53] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib8) |
| Reference [54] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib47) |
| Reference [55] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib48) |
| Reference [56] | [paper-bibliography.md](../references/paper-bibliography.md#bib.bib52) |
