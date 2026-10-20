# Project understanding map: 30 milestones

This map contains exactly **30 project-intrinsic milestones**. It is a common
route through the scientific contribution, independent of a particular
reader's prior knowledge. The single educational anchor is **Nielsen and
Chuang, *Quantum Computation and Quantum Information*, 10th Anniversary
Edition**; [the anchor bridge](ANCHOR.md) fixes the edition, notation and
selected sections.

The allocation is **12 background/preliminaries milestones** and **18
narrative/innovation milestones**. The background joins states, controlled
circuits, differential geometry and statistical interpretation into one
operational problem. The larger narrative phase follows the chart, gradient
interface, optimizer constructions and scientific evidence. Both phases
progress from individual objects to the complete picture. Together the
milestones explain the project without replacing its complete proofs.

| Phase | Milestones | Connection established |
|---|---|---|
| Background/preliminaries | M01–M04 | State vectors, phases, basis labels and executable rotations |
| Background/preliminaries | M05–M08 | Controlled branches, objective derivatives, coordinates and tree weights |
| Background/preliminaries | M09–M12 | Geometry, statistical records and the complete optimization interface |
| Narrative/innovation | M13–M16 | Invertible tree coordinates, native circuits, metric and preparable tangents |
| Narrative/innovation | M17–M22 | Transition moments, indexed gradient access, assigned costs and a complete example |
| Narrative/innovation | M23–M26 | Differential calculus, metric lift and state-sphere optimization |
| Narrative/innovation | M27–M30 | Metrological objectives, evidence and the complete scientific contribution |

## Phase I: background and preliminaries

### M01 Normalized states and real inner products

**Central question.** What object is being optimized?

**Understanding target.** Interpret a normalized vector in $\mathbb R^N$ or
$\mathbb C^N$ as the amplitude representation, with $N=2^n$. For complex
vectors, use $\operatorname{Re}\langle u\mid v\rangle$ when treating the
sphere as a real manifold.

**Prior dependencies.** None within this map.

**Conceptual connection.** A coordinate system must preserve normalization
and carry real parameter displacements into the appropriate tangent space.

**Sources.** N&C §§2.1.1–2.1.7, 2.2.1; [notation](ANCHOR.md#notation-and-operational-conventions);
[coordinate domain](../science/chart-and-geometry.md#Thmdefinition2).

### M02 Computational labels and tensor factors

**Central question.** How does a basis string identify a leaf?

**Understanding target.** Read $\lvert q_n\cdots q_1\rangle$ from most to
least significant bit, distinguish a leaf's basis label from its tree index,
and interpret each bit as a left/right branch choice.

**Prior dependencies.** M01.

**Conceptual connection.** The same addressing convention must join amplitude
products, controlled gates, tangent labels and executable code.

**Sources.** N&C §§2.1.7, 2.2.8; [tree definition](../science/chart-and-geometry.md#Thmdefinition1);
[four-qubit indexing](../science/four-qubit-example.md#A2.SSx1).

### M03 Phase and state-vector representatives

**Central question.** What is the difference between a state vector and a physical ray?

**Understanding target.** Distinguish relative leaf phases from a common
phase, and explain why the complex chart uses $2N-1$ real parameters even
though expectation values are invariant under a common phase shift.

**Prior dependencies.** M01–M02.

**Conceptual connection.** This fixes the geometric object before interpreting
a diagonal metric or counting parameters.

**Sources.** N&C §2.2.7; [sphere versus projective bridge](ANCHOR.md#4-complex-state-vectors-and-physical-rays-use-different-metrics);
[complex chart](../science/chart-and-geometry.md#Thmdefinition2).

### M04 Rotations as local amplitude splits

**Central question.** How does one physical angle divide amplitude?

**Understanding target.** Apply $R_y(\theta)=e^{-i\theta Y}$ to
$\lvert0\rangle$, obtain $(\cos\theta,\sin\theta)$, and translate to
software's half-angle rotation convention.

**Prior dependencies.** M01–M02.

**Conceptual connection.** The parameter used in the geometry is the angle at
the designated split; a factor-of-two mismatch changes both states and gradients.

**Sources.** N&C §4.2; [local split](ANCHOR.md#1-one-split-becomes-an-entire-state);
[Hopf circuit definition](../science/chart-and-geometry.md#Thmdefinition3).

### M05 Conditional operations and coherent branches

**Central question.** How can a rotation address one subtree?

**Understanding target.** Interpret positive and negative controls, preserve
their basis pattern, and distinguish ordinary prepared-state fidelity from
relative-phase correctness inside a coherent branch construction.

**Prior dependencies.** M02–M04.

**Conceptual connection.** Addressing makes the local split executable;
branch-phase conventions later determine the sign of a transition moment.

**Sources.** N&C §§4.3–4.4, 2.2.7; [gate routines](../science/algorithms.md#A1.SSx2);
[signed branches](../science/gradient-access.md#S3.SSx1).

### M06 Expectation objectives and first derivatives

**Central question.** What does a small state change do to an observable's expectation?

**Understanding target.** Differentiate
$E=\langle\psi\mid H\mid\psi\rangle$ for Hermitian $H$ to obtain
$\partial_jE=2\operatorname{Re}\langle\partial_j\psi\mid H\mid\psi\rangle$.

**Prior dependencies.** M01, M03.

**Conceptual connection.** Differentiation produces a transition moment
between the state and a direction, setting the target for a measurement circuit.

**Sources.** N&C §§2.1.6, 2.2.3; [transition-moment bridge](ANCHOR.md#5-a-derivative-becomes-a-measurable-transition-moment);
[gradient theorem](../science/gradient-access.md#Thmtheorem2).

### M07 Coordinates, expressivity and singular representatives

**Central question.** What must a universal state chart provide besides reachability?

**Understanding target.** Separate covering every normalized vector from
having one globally unique coordinate representative. Understand why zero
subtrees and zero leaves make an inverse convention necessary.

**Prior dependencies.** M01–M04.

**Conceptual connection.** Dimension matching, reconstruction and local
regularity are distinct properties of the same representation.

**Sources.** N&C §§2.2.1, 4.5.1; [chart domain](../science/chart-and-geometry.md#Thmdefinition2);
[inverse bridge](ANCHOR.md#2-inverting-the-tree-means-reading-subtree-norms).

### M08 Orthogonal branches and subtree probability

**Central question.** Why does a split tree stay normalized?

**Understanding target.** Use orthogonal branch supports and
$\cos^2\theta+\sin^2\theta=1$ recursively. Identify incoming subtree
probability as the product of ancestor branch probabilities.

**Prior dependencies.** M02, M04, M07.

**Conceptual connection.** The tree's probability accounting later becomes
the diagonal metric and the normalization factor of each derivative direction.

**Sources.** N&C §§2.1.4, 2.1.7; [two-qubit example](ANCHOR.md#1-one-split-becomes-an-entire-state);
[amplitude product formula](../science/proofs.md#A5.SSx1).

### M09 Jacobians and pullback metrics

**Central question.** How does coordinate distance become state-space distance?

**Understanding target.** Define $J_j=\partial_j\psi$ and
$g_{ij}=\operatorname{Re}\langle J_i\mid J_j\rangle$; interpret
$d\theta^{\mathsf T}g\,d\theta$ as squared first-order state displacement.

**Prior dependencies.** M01, M06, M08.

**Conceptual connection.** Coordinate gradients describe sensitivity, while
the metric describes the amount of state motion caused by an update.

**Sources.** N&C §§2.1.2, 2.1.4; local calculus in
[metric bridge](ANCHOR.md#3-the-metric-measures-how-far-an-angle-moves-the-state);
[metric definition](../science/chart-and-geometry.md#Thmdefinition4).

### M10 Sphere geometry and objective invariance

**Central question.** Can an objective-invariant direction still have nonzero metric length?

**Understanding target.** Derive the common-phase tangent $i\psi$ and
contrast the phase blocks $\operatorname{diag}(p)$ and
$\operatorname{diag}(p)-pp^{\mathsf T}$. Identify the first as the paper's convention.

**Prior dependencies.** M03, M09.

**Conceptual connection.** The geometry of state-vector representatives and
the invariances of a scalar objective answer different questions.

**Sources.** N&C §§2.2.7, 9.2.2; [complex geometry bridge](ANCHOR.md#4-complex-state-vectors-and-physical-rays-use-different-metrics);
[state-sphere conventions](../science/optimizers.md#A3.SSx2).

### M11 Circuit settings, observations and accuracy

**Central question.** Which resource is reduced by sharing a circuit family?

**Understanding target.** Distinguish a compiled template, an observed index,
one measurement repetition, an unbiased mean and a complete-vector accuracy
requirement. Read a shot count with its per-sign or per-label convention.

**Prior dependencies.** M05–M06, M08.

**Conceptual connection.** A compact circuit family can organize many
coordinates without fixing the statistical work required for their estimates.

**Sources.** N&C §2.2.3, Appendix 1; [settings bridge](ANCHOR.md#6-a-circuit-family-is-not-a-statistical-budget);
[complete estimator](../science/gradient-access.md#S3.SSx4).

### M12 The complete arbitrary-state optimization interface

**Central question.** How do preparation, inverse coordinates, geometry and gradient access fit together?

**Understanding target.** State the full task with explicit amplitude input,
physical controls, scalar objective, derivative target and resources. Explain
why these interfaces must use consistent coordinates and conventions.

**Prior dependencies.** M01–M11.

**Conceptual connection.** This closes the background phase by joining local
quantum operations, differential geometry and measurement into one variational loop.

**Sources.** N&C §§3.2.1, 4.6; [motivation](../science/motivation.md#S1);
[bridge synthesis](ANCHOR.md#9-put-the-mechanisms-together).

## Phase II: narrative and innovation

### M13 The universal tree chart and its inverse

**Central question.** How can the same tree prepare a vector and recover its angles?

**Understanding target.** Build real and complex amplitudes from path
products, recover splits from subtree norms, handle final real signs and leaf
phases, and account for the linear-time inverse.

**Prior dependencies.** M07–M08, M12.

**Conceptual connection.** The inverse converts a target or updated state
back into directly executable coordinates, making the chart bidirectional.

**Sources.** N&C §§2.2.1, 4.2; [coordinate definition](../science/chart-and-geometry.md#Thmdefinition2);
[inverse lemma](../science/algorithms.md#Thmlemma1);
[inverse proof](../science/proofs.md#A5.SSx6).

### M14 Native real and complex gate schedules

**Central question.** How is the chart implemented with one designated physical angle per coordinate?

**Understanding target.** Follow control, anticontrol, target and parameter
index lists; explain the promoted complex leaf gate and the importance of gate order.

**Prior dependencies.** M04–M05, M13.

**Conceptual connection.** The circuit realization preserves the
coordinate-to-control correspondence needed by the geometry and derivative constructions.

**Sources.** N&C §§4.2–4.3; [ansatz definition](../science/chart-and-geometry.md#Thmdefinition3);
[gate routines](../science/algorithms.md#A1.SSx2).

### M15 The diagonal metric theorem

**Central question.** Why do off-diagonal metric entries vanish?

**Understanding target.** Organize the proof by disjoint branches and
ancestor–descendant relationships; recover subtree weights for magnitudes and
leaf probabilities for complex phases using the real inner product.

**Prior dependencies.** M08–M10, M13.

**Conceptual connection.** Orthogonality turns a dense geometric calculation
into tree probability accounting and componentwise rescaling.

**Sources.** N&C §2.1.4; [metric theorem](../science/chart-and-geometry.md#Thmtheorem1);
[proof](../science/proofs.md#A5.SSx2).

### M16 Coordinate tangents as preparable states

**Central question.** How can an active derivative direction be prepared on the same circuit skeleton?

**Understanding target.** Clamp the ancestor path, shift the chosen magnitude
split by $\pi/2$, retain the appropriate descendants, and derive the leaf-local
complex phase tangent. Preserve the positive-metric condition.

**Prior dependencies.** M05, M14–M15.

**Conceptual connection.** A differential object becomes an operational input
to transition-moment measurements.

**Sources.** N&C §§4.2–4.3; [tangent definition](../science/gradient-access.md#Thmdefinition5);
[tangent theorem](../science/gradient-access.md#Thmtheorem3);
[proof](../science/proofs.md#A5.SSx4).

### M17 Weighted transition moments give coordinate gradients

**Central question.** Which metric factor converts a tangent response into the raw derivative?

**Understanding target.** Derive
$\partial_jE=2\sqrt{g_{jj}}\operatorname{Re}\langle e_j\mid H\mid\psi\rangle$
and verify it for the two-qubit $Z\otimes I$ example.

**Prior dependencies.** M06, M15–M16.

**Conceptual connection.** The same weight that measures geometric speed also
rescales the measurable response of a unit direction.

**Sources.** N&C §§2.1.6, 2.2.3; [gradient theorem](../science/gradient-access.md#Thmtheorem2);
[proof](../science/proofs.md#A5.SSx3);
[worked example](ANCHOR.md#5-a-derivative-becomes-a-measurable-transition-moment).

### M18 Signed branches expose the transition moment

**Central question.** How can expectation measurements recover the cross term?

**Understanding target.** Expand the energies of
$(\psi\pm e_j)/\sqrt2$, check their normalization, and relate the symmetric
and single-sign-centered estimators while retaining branch-phase conventions.

**Prior dependencies.** M05, M11, M17.

**Conceptual connection.** Interference converts an off-diagonal matrix
element into ordinary expectation data with a specified statistical estimator.

**Sources.** N&C §§2.2.3, 4.3–4.4; [signed branch construction](../science/gradient-access.md#S3.SSx1).

### M19 Magnitude-layer batching

**Central question.** How does one indexed construction cover every magnitude coordinate at one depth?

**Understanding target.** Map node $j=2^d+r$ to the $d$-bit layer label,
follow derivative and branch preparation, and identify which coordinate an
observed label updates.

**Prior dependencies.** M02, M11, M16–M18.

**Conceptual connection.** Shared circuit structure reduces the number of
compiled families while preserving the label-resolved estimator.

**Sources.** N&C §§4.3–4.4; [magnitude batching](../science/gradient-access.md#S3.SSx2);
[layerwise example](../science/four-qubit-example.md#A2.SSx5).

### M20 One additional complex phase family

**Central question.** Why can all phase tangents be organized in one leaf-indexed layer?

**Understanding target.** Use the leaf-local phase differential to derive its
indexing and estimator, and obtain the real and complex setting counts with
the baseline included.

**Prior dependencies.** M03, M16, M19.

**Conceptual connection.** The phase structure adds a single indexed family,
while the statistical interpretation remains label dependent.

**Sources.** N&C §§2.2.7, 4.3–4.4; [phase batching](../science/gradient-access.md#S3.SSx3);
[complete estimator and counts](../science/gradient-access.md#S3.SSx4).

### M21 Assigned circuit costs and statistical costs

**Central question.** Under which compiler model is a per-setting cost claim valid?

**Understanding target.** Apply the declared mixed-control charges, follow
native schedule counts, and separate controlled-gate work from observable
readout and repetitions for a requested accuracy.

**Prior dependencies.** M11, M14, M19–M20.

**Conceptual connection.** The combined resource statement needs both the
circuit ledger and the statistical task, rather than a template count alone.

**Sources.** N&C §§3.2.1, 4.5.1–4.5.2;
[uniform per-setting theorem](../science/gradient-access.md#Thmtheorem4);
[counting model](../science/algorithms.md#Thmlemma2);
[native CNOT counts](../science/algorithms.md#Thmtheorem5).

### M22 The four-qubit construction as a complete instance

**Central question.** Can every abstract index be matched to a concrete operation or metric factor?

**Understanding target.** Follow one four-qubit magnitude coordinate and one
phase coordinate from the tree through the gate schedule, tangent assignment,
metric and branch circuit, using consistent basis order.

**Prior dependencies.** M13–M21.

**Conceptual connection.** This joins the local mechanisms into an executable
instance before the narrative turns to optimization.

**Sources.** N&C §§4.2–4.4; [four-qubit companion](../science/four-qubit-example.md#A2).

### M23 Jacobian formulas and inverse differentials

**Central question.** How does the analytic Jacobian connect coordinates to arbitrary tangent displacements?

**Understanding target.** Interpret its magnitude and phase columns, recover
its Gram form, and handle boundary coordinates without unstable divisions by
amplitudes or trigonometric factors.

**Prior dependencies.** M09, M13, M15–M16.

**Conceptual connection.** The explicit derivative calculus connects a
measurement interface to the algebra required by an optimizer.

**Sources.** N&C §§2.1.2, 2.1.4; [Jacobian theorem](../science/optimizers.md#Thmtheorem6);
[proof](../science/proofs.md#A5.SSx7).

### M24 Lifting coordinate gradients to the state sphere

**Central question.** How does diagonal metric rescaling produce a state-space direction?

**Understanding target.** Explain $G=Jg^{-1}\nabla_\theta E$ on the regular
chart, distinguish coordinate from state gradients, and identify what changes
when numerical damping is used near a zero weight.

**Prior dependencies.** M10, M17, M23.

**Conceptual connection.** Geometric rescaling compensates for coordinate
speed rather than creating new objective information.

**Sources.** N&C §2.1.4; [metric lift](../science/optimizers.md#A3.SSx3);
[update bridge](ANCHOR.md#7-geometry-supplies-a-state-space-update).

### M25 Exact sphere updates and vector transport

**Central question.** How can the state move along a tangent while preserving normalization?

**Understanding target.** Derive the great-circle update, interpret transport
of tangent directions along that curve, and explain why the inverse map
returns the new state to executable coordinates.

**Prior dependencies.** M13, M24.

**Conceptual connection.** Geometry makes the forward state chart part of a
closed optimization loop rather than a one-way preparation routine.

**Sources.** N&C §§2.1.4, 2.2.1; [geodesic and transport construction](../science/optimizers.md#A3.SSx4).

### M26 The geometry-native optimizer family

**Central question.** What is shared, and what changes, across the three optimization schemes?

**Understanding target.** Compare conjugate-gradient memory, transported
L-BFGS curvature pairs and Barzilai–Borwein step estimates. Identify common
state-sphere geometry and the role of line-search safeguards.

**Prior dependencies.** M24–M25.

**Conceptual connection.** The Hopf interface supports several optimizer
choices; their numerical behavior is distinct from the exact chart identities.

**Sources.** N&C §3.2.1; [optimizer family](../science/optimizers.md#A3.SSx7);
[conjugate gradients](../science/optimizers.md#A3.SSx4);
[L-BFGS](../science/optimizers.md#A3.SSx5);
[Barzilai–Borwein](../science/optimizers.md#A3.SSx6).

### M27 QFI as a local chain-rule observable

**Central question.** How does a nonlinear probe objective fit the transition-moment interface?

**Understanding target.** Differentiate pure-state QFI through its expectation
values, form the local Hermitian chain-rule observable, and hold its
coefficients fixed during branch measurements at the current state.

**Prior dependencies.** M06, M17–M18.

**Conceptual connection.** The derivative interface extends beyond a fixed
energy observable through a local chain rule, with explicit estimation assumptions.

**Sources.** N&C §§2.1.6, 2.2.3; [chain-rule form](../science/metrological-objectives.md#A4.SSx1);
[pure-state QFI](../science/metrological-objectives.md#A4.SSx2);
[worked chain rule](ANCHOR.md#8-nonlinear-objectives-require-a-local-observable).

### M28 Fixed-readout CFI and its measurement conditions

**Central question.** Which quantities must be available to turn fixed-readout CFI into a gradient observable?

**Understanding target.** Follow the probability and parameter-derivative
terms in the CFI chain rule, retain its local-estimation and support
conditions, and interpret its compiled-setting scaling.

**Prior dependencies.** M11, M20–M21, M27.

**Conceptual connection.** A metrological objective adds a specified physical
readout and local coefficients to the common state-coordinate construction.

**Sources.** N&C §§2.2.3, 2.2.5, Appendix 1;
[fixed-readout CFI](../science/metrological-objectives.md#A4.SSx3);
[setting scaling](../science/metrological-objectives.md#A4.SSx4).

### M29 Interpret the numerical and analytical evidence

**Central question.** Which scientific conclusions follow from proofs, and which describe the chosen benchmarks?

**Understanding target.** Read task definitions, exact objective/gradient
access, known optima, initialization and aggregation conventions before
interpreting performance. Connect figures and code to their specific claims.

**Prior dependencies.** M21–M28.

**Conceptual connection.** Exact identities, finite implementation checks and
optimizer experiments provide complementary evidence with different domains.

**Sources.** N&C §3.2.1, Appendix 1;
[numerical design](../science/numerical-studies.md#S4.SSx1);
[aggregate results](../science/numerical-studies.md#S4.SSx2);
[proofs](../science/proofs.md#A5); [reproduction](../verification/README.md).

### M30 The complete Hopf optimization contribution

**Central question.** How does the whole construction turn arbitrary-state preparation into a navigable optimization framework?

**Understanding target.** Give a connected account of expressivity, inverse
coordinates, diagonal geometry, preparable tangents, gradient access and
state-sphere updates. State their assumptions, evidence and resource meanings
without conflating circuit families with complete-gradient measurement cost.

**Prior dependencies.** M12–M29.

**Conceptual connection.** The local split, differential and measurement
mechanisms combine into the full scientific contribution and its interpretation.

**Sources.** N&C §§3.2.1, 4.6;
[motivation](../science/motivation.md#S1);
[conclusion](../science/conclusion.md#S5);
[bridge synthesis](ANCHOR.md#9-put-the-mechanisms-together).
