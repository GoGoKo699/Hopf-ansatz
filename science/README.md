# A Compass on the Quantum State Sphere: The Hopf Ansatz for Arbitrary Pure-State Optimization

Ruge Lin  ([ORCID](https://orcid.org/0000-0002-9297-7622))

Corresponding author: [gogoko699@gmail.com](mailto:gogoko699@gmail.com)

Affiliation: Thrust of Artificial Intelligence, Information Hub, The Hong Kong University of Science and Technology (Guangzhou), Guangzhou 511453, China

Guangxi Li  ([ORCID](https://orcid.org/0009-0009-5481-2411))

Affiliation: Quantum Science Center of Guangdong-Hong Kong-Macao Greater Bay Area, Shenzhen 518045, China

<a id="abstract1"></a>

## Abstract

Optimizing arbitrary quantum state vectors is like navigating the unit sphere in Hilbert space: beyond target reachability, optimization asks for coordinates, local distance information, and measurable directions. We introduce the Hopf ansatz, a binary-tree circuit for arbitrary normalized real and complex state vectors. Internal angles steer probability between subtrees, leaf phases carry the complex degrees of freedom, and the same tree gives state preparation and an explicit inverse map from amplitudes to physical angles. Together these structures act as a compass for the search: the inverse map gives coordinates, the diagonal induced metric gives local distance information, and nonzero coordinate tangents become preparable normalized states.

For Hamiltonian objectives, and for objectives with the same local transition-moment form, each gradient component is a known scale factor times a transition moment between the current state and a tangent state. A branch-state construction expresses these moments through expectation-value measurements, while the tree organizes the compiled gradient settings by magnitude layer and leaf-indexed phase family. Thus the number of distinct compiled gradient-access circuit families grows only logarithmically with Hilbert-space dimension, while the measurement budget required for a chosen precision remains a separate statistical cost. In deterministic real-state benchmarks with exact costs, exact gradients, and known global optima, metric-aware Hopf optimizers reach numerical-precision median gaps, with the clearest baseline gains in smaller VQE mean final gaps and stronger concentration of metrology-inspired traces at numerical precision. The Hopf ansatz turns universal state preparation into a navigable framework for arbitrary pure-state optimization.

## Read the paper

This is the complete Markdown text of [arXiv:2607.14231](https://arxiv.org/abs/2607.14231), [version 1](https://arxiv.org/abs/2607.14231v1). Numbering follows the paper PDF.

[Reviewer route](../review/README.md) · [Learning route](../learning/MILESTONES.md) · [Repository](../README.md)

| Paper part | Read |
| --- | --- |
| I Introduction | [motivation.md](motivation.md#S1) |
| II A universal Hopf variational chart | [chart-and-geometry.md](chart-and-geometry.md#S2) |
| III Hopf gradient access and scaling | [gradient-access.md](gradient-access.md#S3) |
| IV Numerical simulations | [numerical-studies.md](numerical-studies.md#S4) |
| V Conclusion | [conclusion.md](conclusion.md#S5) |
| Data and code availability | [conclusion.md](conclusion.md#Sx1) |
| References | [paper-bibliography.md](../references/paper-bibliography.md#bib) |
| Appendix A Algorithms and subroutines | [algorithms.md](algorithms.md#A1) |
| Appendix B A four-qubit example of the Hopf ansatz | [four-qubit-example.md](four-qubit-example.md#A2) |
| Appendix C Geometry-native optimizers for the Hopf ansatz | [optimizers.md](optimizers.md#A3) |
| Appendix D Metrological objectives as chain-rule observables | [metrological-objectives.md](metrological-objectives.md#A4) |
| Appendix E Proofs of the main results | [proofs.md](proofs.md#A5) |
