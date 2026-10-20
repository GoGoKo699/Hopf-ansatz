[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="S5"></a>

# V Conclusion

The Hopf ansatz gives arbitrary pure-state optimization an executable geometry. A single binary tree prepares normalized real and complex state vectors, recovers circuit angles from amplitudes, and equips those angles with local geometric data for metric-aware updates. In deterministic real-state benchmarks with exact cost and gradient evaluations, metric-aware Hopf optimizers reach numerical-precision median gaps; relative to coordinate-update baselines, the clearest gains are smaller VQE mean final gaps and stronger concentration of metrology-inspired traces at numerical precision.

At the optimizer interface, local calculus becomes circuit operations. Coordinate derivatives appear as known metric factors times tangent-state transition moments, and branch-state constructions express them through measurable expectation values. The tree organizes the compiled gradient-access circuit families: magnitude tangents are grouped by layer, and complex phase tangents by a leaf-indexed construction. Consequently, the number of distinct compiled gradient-access circuit families grows only logarithmically with Hilbert-space dimension, while the measurement budget required for a chosen precision remains a separate statistical cost.

<a id="S5.p3.1"></a>

The framework applies most directly when the search object is an arbitrary normalized quantum state vector and the quantity being optimized admits a usable local derivative formula. The same geometric role can be explored in state discrimination [[42](../references/paper-bibliography.md#bib.bib42)], variational quantum hypothesis testing [[43](../references/paper-bibliography.md#bib.bib30)], variational sensing and metrology [[4](../references/paper-bibliography.md#bib.bib39), [44](../references/paper-bibliography.md#bib.bib41), [5](../references/paper-bibliography.md#bib.bib31), [6](../references/paper-bibliography.md#bib.bib25), [9](../references/paper-bibliography.md#bib.bib28), [45](../references/paper-bibliography.md#bib.bib40), [10](../references/paper-bibliography.md#bib.bib29)], quantum control [[46](../references/paper-bibliography.md#bib.bib44), [47](../references/paper-bibliography.md#bib.bib45)], and variational compilation [[48](../references/paper-bibliography.md#bib.bib43)]. Across these settings, objectives, readouts, and performance gains remain problem-specific, but the need for an executable bridge between state-space geometry and circuit parameters is shared.

The main message is structural. When arbitrary-state freedom is part of the model, the circuit can carry its own coordinate map, metric, and differential probes. The Hopf ansatz turns universal state preparation into a navigable circuit framework, a compass for moving on the quantum state sphere.

<a id="S5.acknowledgements1"></a>

### Acknowledgements.

 The authors are sincerely grateful to Allan Tosta, André J. Ferreira-Martins, Thiago O. Maciel, and Tobias Haug for thoughtful discussions and perspectives that helped shape this work, and to Leandro Aolita and the Technology Innovation Institute (TII) for their valued support of the project. Special thanks are due to Giancarlo Camilo and Renato M. S. Farias for sustained guidance, encouragement, and perspective, and to Ranyiliu Chen and Xin Wang for help with the final revision. Guangxi Li was supported by National Natural Science Foundation of China (Grant No. 62402206). This work was partially supported by the Guangdong Provincial Quantum Science Strategic Initiative (Grant Nos. GDZX2403008, GDZX2503001, and GDZX2403001).

<a id="Sx1"></a>

# Data and code availability

<a id="Sx1.p1.1"></a>

The manuscript uses analytic derivations and deterministic state-vector simulations; no experimental data were used. The accompanying public repository [[49](../references/paper-bibliography.md#bib.bib46)] provides documented commands to regenerate the reported real-Hopf stress-test data, diagnostics, and figures. It also includes auxiliary checks: an $`n=6`$ complex-Hopf stress test using the same three VQE and three metrology-inspired task families, finite-shot signed-branch gradient tests, CNOT-count verification, and Qibo [[50](../references/paper-bibliography.md#bib.bib53)] circuit-level gradient-access examples for local four-qubit real and complex VQE toys.

Generated CSV files are omitted from source control because they are derived artifacts and can be large, but they can be regenerated from the documented commands. Selected figures and diagnostic panels are included as release artifacts for direct comparison with the manuscript results.
