# Claim support and evidence map

This page is for reviewers and readers checking whether the repository supports
the claims associated with *A Compass on the Quantum State Sphere: The Hopf
Ansatz for Arbitrary Pure-State Optimization*.

The repository contains several kinds of evidence. They should not be
conflated.

| Support type | Meaning |
|---|---|
| **Analytic identity implemented and checked** | Two independent finite-dimensional calculations or constructions are compared numerically. |
| **Explicit circuit/statevector check** | A Qibo circuit or an equivalent ideal statevector path is compared with the analytic Hopf construction. |
| **Generated-schedule accounting** | A gate list is generated and its assigned resource count is compared with a closed-form formula. |
| **Finite-shot evidence** | A fixed estimator is sampled repeatedly and compared with its exact reference. |
| **Optimization evidence** | Deterministic synthetic tasks are run across specified sizes, seeds, and optimizer modes. |
| **Mathematical statement exposed by code** | The code implements finite premises or formulas, while the all-size conclusion is proved in the paper rather than by numerical testing. |

Numerical experiments do not prove universality, diagonal geometry, or
asymptotic complexity. The paper provides the formal statements and proofs; the
repository makes the constructions inspectable and tests their finite
realizations.

## Audit summary

| Paper-level statement | Repository support | Main files | Evidence boundary |
|---|---|---|---|
| The real and complex Hopf recursions produce normalized states with the stated coordinate order | Forward-map and round-trip checks | `hopf_utils.py` | Finite numerical instances; the all-size recursion is analytic. |
| Any normalized real or complex state has an explicit Hopf representative | Forward-after-inverse round-trip | `hopf_utils.py` | Zero subtrees and zero leaves use a fixed convention; coordinates are not globally unique. |
| The inverse map has linear classical time and memory | Tree-based implementation and paper algorithm | `hopf_utils.py` | Complexity is an algorithmic deduction, not a timing benchmark. |
| The pullback metric is diagonal | Dense Jacobian Gram matrix versus analytic diagonal entries | `hopf_utils.py`, `hopf_complex.py` | Checked at finite dimensions; complex geometry uses the round-sphere real pairing. |
| Normalized coordinate tangents can be prepared on the same Hopf skeleton | Tangent assignment versus normalized Jacobian columns | `hopf_utils.py` | Only regular coordinates have a normalized tangent; zero metric gives a vanishing raw differential. |
| The native gate schedules reproduce the analytic states | Generated schedules and optional Qibo parity | `hopf_utils.py` | State preparation from the initialized input; not a statement that all unitary completions are unique. |
| Real and complex assigned CNOT formulas match the generated schedules | Schedule-by-schedule counting versus closed forms | `hopf_gate_count.py` | Compiler-relative, no-clean-ancilla model; not transpiler or hardware counts. |
| Hopf coordinate gradients can be evaluated by tree contractions | Fast routines compared with dense Jacobian references or exact identities | `hopf_data.py`, `hopf_complex.py`, `finite_shot_sanity_check.py` | Classical reference calculations; they do not replace the quantum-access assumptions of the paper. |
| Signed branch states recover coordinate gradients | Exact branch identity and finite-shot sampling | `finite_shot_sanity_check.py`, `VQE_qibo.py` | Fixed observable decompositions and stated shot convention. |
| Layerwise preparations organize a full gradient into logarithmically many compiled families | Explicit `n = 4` layerwise circuit safeguard and general layer specification | `VQE_qibo.py` | The local demo is not a general all-`n` compiled benchmark; the all-size count is analytic. |
| Geometry-native optimizers use only cost and Hopf-coordinate-gradient queries | Executable EGT-CG, R-LBFGS, and R-BB pipelines | `hopf_data.py`, `hopf_complex.py` | Classical statevector studies, not hardware optimizer benchmarks. |
| The real-chart optimization study is stable across the stated task family | Multi-size deterministic generation and streaming diagnostics | `hopf_data.py`, `diagnose_hopf.py`, `plot_hopf.py` | Synthetic tasks at the specified sizes and seeds. |
| The complex extension is operational beyond isolated identities | Focused `n = 6` complex stress test and self-checks | `hopf_complex.py` | One size and six synthetic task families; no complex Möttönen baseline. |
| Finite-shot error follows the expected inverse-square-root trend in the fixed test | Repeated sampling at three shot counts | `finite_shot_sanity_check.py` | Shots are allocated per branch and per component; not a global hardware budget. |
| A local real/complex layerwise estimator tracks exact-gradient VQE trajectories | Explicit Qibo or equivalent ideal statevector simulation | `VQE_qibo.py` | Local `n = 4` toy Hamiltonians; not performance or asymptotic evidence. |

## 1. Forward coordinate maps

### Operational statement

For `n` qubits, let `N = 2**n`. The real chart uses `N - 1` internal-node
angles. The complex chart uses the same magnitude block followed by `N` leaf
phases.

The forward APIs are:

```python
vector_from_theta(theta, case="real")
vector_from_theta(theta, case="complex")
```

Each tree level splits the incoming subtree amplitude with cosine and sine
factors. The complex map multiplies the final leaf magnitudes by independent
phases.

### Evidence

When `hopf_utils.py` is run directly, it constructs deterministic real and
complex coordinates, maps them to states, applies the inverse map, and checks
the reconstructed states. It also compares the analytic states with Qibo
circuits when Qibo is installed.

### Boundary

The map is a sphere parameterization, not a globally one-to-one Euclidean
chart. Boundary coordinates can be nonunique.

## 2. Explicit inverse map

### Operational statement

The inverse API is:

```python
theta_from_vector(vector, case="real")
theta_from_vector(vector, case="complex")
```

For every internal node, the magnitude angle is obtained from the norms of its
left and right leaf subtrees. In the real final layer, a signed sibling pair is
handled with a quadrant-sensitive `atan2` convention. In the complex chart, one
phase is read from each leaf.

At a zero subtree or zero sibling pair, the implementation returns angle zero.
At a zero complex leaf, its phase is set to zero. These choices select one
representative without claiming uniqueness.

### Evidence

`hopf_utils.py` checks state-level round trips rather than comparing coordinate
vectors, because multiple coordinate vectors can encode the same boundary
state.

### Boundary

The paper's `O(N)` time and `O(N)` memory statement follows from one bottom-up
subtree-weight pass and one linear phase/sign pass. The repository does not
claim a wall-clock advantage from a specific Python implementation.

## 3. Jacobian and diagonal metric

### Operational statement

The reference APIs are:

```python
jacobian(theta, case)
metric_diagonal(theta, case)
```

For real coordinates, the metric entry of a node is the incoming probability
mass of that subtree. For complex coordinates, the magnitude block is the same
and the phase block contains leaf occupations.

```math
g^{\mathbb{C}}_{i,i}=
\begin{cases}
g^{\mathbb{R}}_{i,i}, & 1\leq i\leq N-1,\\
\lvert x_{i-N}\rvert^2, & N\leq i\leq 2N-1.
\end{cases}
```

### Evidence

`hopf_utils.py` compares normalized tangent assignments with Jacobian columns.
`hopf_complex.py` additionally compares its linear-time complex coordinate
gradient with the dense complex Jacobian and checks that the metric lift equals
the projected state-space gradient.

### Boundary

The complex metric is the round-sphere Gram matrix under the real pairing. The
uniform phase direction is retained because the implementation parameterizes
normalized state vectors rather than quotienting out global phase.

## 4. Normalized tangent-state synthesis

### Operational statement

The tangent-assignment API is:

```python
theta_hopf_tangent_state(theta, i, case)
qibo_tangent_state(theta, i, case, minimal=True)
```

For a magnitude node, the construction:

1. shifts the target angle by `pi/2`;
2. clamps each ancestor to `0` or `pi/2` according to the path;
3. retains the target subtree's descendant angles;
4. retains only phases in that subtree for the complex chart; and
5. sets unrelated entries to zero.

For a complex phase coordinate, the magnitude path is clamped to the selected
leaf and its phase is shifted by `pi/2`.

### Evidence

`hopf_utils.py` forms the dense Jacobian and verifies that every nonzero
normalized column agrees with the state prepared by the corresponding tangent
assignment. Optional Qibo checks exercise the same skeleton.

### Boundary

If the metric entry is zero, the raw derivative vanishes and a normalized
coordinate tangent is not defined. The code avoids dividing by such entries in
the identity checks.

## 5. Gate schedules and circuit parity

### Operational statement

```python
gates_order(n, case)
```

returns four parallel lists:

```text
Ctrl, Anti, Targ, Index
```

They encode positive-control masks, negative-control masks, a power-of-two
target mask, and either one Hopf parameter index or a three-index complex leaf
gate. Parameter indices are one based.

```python
qibo_circuit(theta, case, minimal=False)
```

converts the schedule into Qibo gates. Negative controls are implemented by
`X` conjugation. The Hopf rotation convention sends `RY(2*theta)` to Qibo.

### Evidence

The direct utility check compares circuit statevectors with the analytic map.
The second-paper repository independently reproduces the same native schedules
and validates their initialized state columns through additional exact-logical
tests.

### Boundary

Circuit parity is a prepared-state statement from the initialized register. It
does not assert a unique full-unitary completion for every implementation.

## 6. Assigned CNOT ledger

### Operational statement

`hopf_gate_count.py` generates both schedules, counts every controlled gate
under the declared no-clean-ancilla model, and compares the totals with the
closed-form binomial formulas.

For `m` controls:

```math
c_{R_y}(m)=
\begin{cases}
0, & m=0,\\
2^{m+1}-2, & 1\leq m\leq4,\\
16(m+1)-40, & m\geq5.
\end{cases}
```

The complex promoted gate uses the same small-control values and a separate
linear large-control formula. The script reports both numerical and theoretical
totals.

### Boundary

These are assigned logical CNOT charges. They exclude device topology,
transpiler choices, approximate synthesis, state-measurement cost, and the
observable layer. The asymptotic `O(n*2**n)` conclusion is mathematical rather
than a fit to the plotted data.

## 7. Fast coordinate gradients and metric lift

### Operational statement

The real stress-test code evaluates the coordinate gradient by a bottom-up
subtree-response pass followed by a top-down incoming-amplitude pass. The
complex code uses the analogous complex contraction for magnitude angles and a
leaf-local phase formula.

The optimizer then converts the coordinate gradient to a state-sphere tangent
using the analytic diagonal metric. No dense metric inversion is performed.

### Evidence

- `finite_shot_sanity_check.py` compares its independent real tree gradient
  with a dense Jacobian reference and the signed-branch identity.
- `hopf_complex.py` compares its fast gradient with `Re(J† h)` and its metric
  lift with a direct tangent projection.
- the stress-test scripts monitor state-norm errors and nonfinite gradients.

### Boundary

The stress-test implementations use a small metric floor when lifting near a
singular chart boundary. This is an optimizer safeguard, not a redefinition of
the exact metric.

## 8. Signed branch estimator

### Operational statement

For a normalized tangent `e_i`, define

```math
\lvert\varphi_i^{(s)}\rangle
=
\frac{\lvert\psi\rangle+s\lvert e_i\rangle}{\sqrt{2}},
\qquad s\in\{+1,-1\}.
```

Then the symmetric energy identity is

```math
\partial_{\theta_i}E
=
\sqrt{g_{i,i}}
\left(E_{\varphi_i^{(+)}}-E_{\varphi_i^{(-)}}\right).
```

`finite_shot_sanity_check.py` tests this identity component by component before
adding sampling noise. `VQE_qibo.py` implements the equivalent baseline,
tangent-energy, and sign-conditioned branch form.

### Boundary

The finite-shot script allocates `S` shots to each sign for each gradient
component. It does not claim that `S` is the total full-gradient budget.
Hamiltonian grouping and device noise are outside that experiment.

## 9. Layerwise circuit organization

### Operational statement

Magnitude coordinates at depth `d` form one layer of size `2**d`. The complex
phase block forms one additional `N`-label layer. The paper's setting count is:

```math
N_{\mathrm{set}}^{\mathbb{R}}=1+2n,
\qquad
N_{\mathrm{set}}^{\mathbb{C}}=1+2(n+1).
```

`VQE_qibo.py` constructs matching layer specifications and, at `n = 4`, builds
explicit label-controlled tangent and branch circuits or an equivalent ideal
statevector sampler.

### Boundary

The explicit Qibo path is intentionally local and transparent. It is not the
asymptotically optimized general compiler. The setting-count and per-setting
resource conclusions are analytic statements in the paper.

## 10. Optimization studies

### Operational statement

The geometry-native real study implements:

- `Hopf-EGT-CG`;
- `Hopf-Riemannian-LBFGS`; and
- `Hopf-Riemannian-BB`.

All consume the same cost plus Hopf-coordinate-gradient interface. The Adam
study adds coordinate `Hopf-Adam` and a real Möttönen physical-angle baseline.
The focused complex study runs four Hopf modes at `n = 6`.

### Evidence

The generator scripts record complete task, seed, mode, and step grids.
Streaming diagnostics check missing rows, vector lengths, nonfinite values,
state normalization, threshold hits, and final rankings. Plotting scripts read
the regenerated CSVs rather than hard-coding result arrays.

### Boundary

The tasks are deterministic synthetic stress tests. They demonstrate
operational behavior under the stated designs; they do not establish universal
optimizer rankings or physical-device speedups.

## 11. What the repository does not establish

The repository should not be used as evidence for any of the following:

- that finite tests prove universality or all-size diagonal geometry;
- that assigned CNOT counts equal routed hardware counts;
- that the finite-shot plot is a complete sampling-complexity analysis;
- that the `n = 4` layerwise demo realizes the general optimized compiler;
- that the real multi-size results automatically extend to every complex size;
- that one optimizer dominates on arbitrary external applications;
- that a complex Möttönen comparison was performed; or
- that noisy hardware, measurement grouping, or error mitigation has been
  validated.

For exact-logical reverse-gradient constructions beyond the first paper's
layerwise safeguard, use the companion
[Hopf-QBP repository](https://github.com/GoGoKo699/Hopf-QBP).
