# Educational anchor and bridge

The single educational anchor is **Michael A. Nielsen and Isaac L. Chuang,
*Quantum Computation and Quantum Information*, 10th Anniversary Edition,
Cambridge University Press, ISBN 9781107002173**. The
[publisher's contents](https://assets.cambridge.org/97811070/02173/toc/9781107002173_toc.pdf)
fix the section numbering used here. The companion
[Hopf-QBP learning route](https://github.com/GoGoKo699/Hopf-QBP/blob/main/learning/ANCHOR.md)
uses the same edition and notation bridge.

This route starts from linear algebra, elementary calculus, and basic
probability. Read selected sections alongside the local explanations below;
the whole textbook is not a prerequisite. The [30 milestones](MILESTONES.md)
organize understanding of the project. Exact claims, proofs, algorithms,
figures, and attribution belong to the linked scientific account of
[*A Compass on the Quantum State Sphere: The Hopf Ansatz for Arbitrary
Pure-State Optimization*](https://arxiv.org/abs/2607.14231).

## Selected textbook route

| Nielsen–Chuang sections | Foundation | Connection developed here |
|---|---|---|
| §§2.1.1–2.1.7; §§2.2.1–2.2.2 | Bases, inner products, Hermitian operators, tensor products, states and evolution | A normalized amplitude vector becomes an executable tree of coordinates. |
| §§2.2.3, 2.2.5, 2.2.7–2.2.8 | Measurement, phase and composite systems | Distinguish phase conventions, state-vector geometry and measurable expectation values. |
| §§4.2–4.4 | Single-qubit operations, controlled operations and measurement | Each tree split is a directly programmed rotation; signed branches expose transition moments. |
| §§3.2.1, 4.5.1–4.5.2 | Resource accounting and synthesis using one-qubit gates and CNOTs | Interpret the paper's assigned controlled-gate ledger and its circuit-family counts. |
| §§9.2.1–9.2.2 | Trace distance and fidelity | Relate state comparison to optimization, while keeping the chosen sphere metric explicit. |
| Appendix 1 | Probability fundamentals | Separate expectation identities, independent repetitions and finite-shot error. |

The book supplies the quantum-mechanical and circuit foundation. The binary-tree
chart, its metric and inverse, tangent-state preparations, optimization
constructions, and metrological chain rules are developed in this repository.
Their source pointers are given below rather than attributed to the textbook.

## Notation and operational conventions

| Symbol or convention | Meaning |
|---|---|
| $n$, $N=2^n$ | Number of qubits and number of computational-basis amplitudes. |
| $\lvert q_n\cdots q_1\rangle$ | Basis label, read from most-significant to least-significant bit. |
| $j=1,\ldots,N-1$ | Breadth-first internal tree node; children are $2j$ and $2j+1$. |
| $\theta_j$ | Angle splitting the amplitude entering node $j$. |
| $\phi_\ell=\theta_{N+\ell}$ | Phase of complex leaf $\ell$, for $0\leq\ell<N$. |
| $J_j=\partial_{\theta_j}\psi$ | Raw coordinate differential, a vector in the state space. |
| $g_{ij}=\mathrm{Re}\langle J_i\mid J_j\rangle$ | Pullback of the round metric on normalized state vectors. |
| $w_j=\sqrt{g_{jj}}$ | Nonnegative incoming magnitude on the canonical chart. |
| $e_j=J_j/w_j$ | Normalized coordinate tangent when $w_j>0$. |
| $R_y(\theta)=e^{-i\theta Y}$ | Paper convention; a software gate using $e^{-i\varphi Y/2}$ receives $\varphi=2\theta$. |

The real chart has $N-1$ parameters. Its non-final magnitude angles lie in
$[0,\pi/2]$; final sibling-pair angles lie in $[0,2\pi)$ so they can encode
signs. The complex chart uses $N-1$ nonnegative magnitude splits and $N$ leaf
phases, hence $2N-1$ real parameters. These are dimensions of normalized
state-vector spheres. The complex chart retains the common phase.

## 1. One split becomes an entire state

Take orthonormal branch states $\lvert u\rangle$ and $\lvert v\rangle$.
A normalized split and its complementary direction are

```math
\lvert s(\theta)\rangle
=\cos\theta\lvert u\rangle+\sin\theta\lvert v\rangle,
\qquad
\partial_\theta\lvert s\rangle
=-\sin\theta\lvert u\rangle+\cos\theta\lvert v\rangle.
```

The derivative has norm one and is orthogonal to the split state. The
orthogonality of the two branches is what makes both statements immediate.
Replacing each branch by another orthogonal split builds a normalized binary
tree without a final normalization step.

For two qubits, put $a=\theta_1$, $b=\theta_2$, $c=\theta_3$:

```math
\psi(a,b,c)=
\begin{pmatrix}
\cos a\cos b\\
\cos a\sin b\\
\sin a\cos c\\
\sin a\sin c
\end{pmatrix}.
```

At $(a,b,c)=(\pi/4,\pi/6,\pi/3)$ this gives

```math
\psi=\frac14(\sqrt6,\sqrt2,\sqrt2,\sqrt6)^{\mathsf T}.
```

The first split chooses the total probability of the first or second leaf
pair. The next two splits choose probabilities within their own pairs. In the
complex chart, multiply component $\ell$ by $e^{i\phi_\ell}$.

**Connect to the paper:** [chart definition](../science/chart-and-geometry.md#Thmdefinition2),
[circuit definition](../science/chart-and-geometry.md#Thmdefinition3), and
[amplitude product proof](../science/proofs.md#A5.SSx1).

## 2. Inverting the tree means reading subtree norms

For a normalized real vector $x=(x_0,x_1,x_2,x_3)$, define

```math
r_L=\sqrt{x_0^2+x_1^2},\qquad
r_R=\sqrt{x_2^2+x_3^2}.
```

Then $a=\mathrm{atan2}(r_R,r_L)$,
$b=\mathrm{atan2}(x_1,x_0)$ and
$c=\mathrm{atan2}(x_3,x_2)$, with final angles wrapped into
$[0,2\pi)$. These recover the preceding example. For complex vectors, use
absolute values to recover magnitude splits and use
$\phi_\ell=\arg x_\ell$ on nonzero leaves.

A zero subtree has no preferred internal coordinates: changing its angles
cannot change the state. Choosing zero angles there is an inverse convention.
The appropriate round-trip check is reconstruction of the input vector, not
equality with every possible parameter representation. Accumulating norms
bottom-up visits each tree node once, giving $O(N)$ classical work.

**Connect to the paper:** [inverse map and cost](../science/algorithms.md#Thmlemma1)
and [its proof](../science/proofs.md#A5.SSx6).

## 3. The metric measures how far an angle moves the state

For a small parameter displacement $d\theta$, the state changes to first
order by $d\psi=J\,d\theta$. Its squared round-sphere length is

```math
\lVert d\psi\rVert^2
=\sum_{i,j}g_{ij}\,d\theta_i\,d\theta_j,
\qquad
 g_{ij}=\mathrm{Re}\langle\partial_i\psi\mid\partial_j\psi\rangle.
```

Different branch supports are orthogonal. Within a branch, its normalized
state is orthogonal to its real magnitude derivative. Recursing through these
two facts makes the Hopf metric diagonal. In the two-qubit example,

```math
g=\mathrm{diag}(1,\cos^2 a,\sin^2 a).
```

At $a=\pi/4$, a unit change of $b$ or $c$ has squared state-space speed
$1/2$, while a unit change of $a$ has speed one. The diagonal weight is the
probability entering the relevant subtree. A zero weight means the raw
coordinate differential vanishes; dividing by zero does not define a tangent.

**Connect to the paper:** [metric convention and theorem](../science/chart-and-geometry.md#Thmtheorem1),
[metric proof](../science/proofs.md#A5.SSx2), and
[Jacobian and Gram forms](../science/optimizers.md#Thmtheorem6).

## 4. Complex state vectors and physical rays use different metrics

For one complex qubit,

```math
\psi=(e^{i\phi_0}\cos\theta,e^{i\phi_1}\sin\theta)^{\mathsf T}.
```

The paper's round-sphere line element is

```math
ds_{\mathrm{sph}}^2
=d\theta^2+\cos^2\theta\,d\phi_0^2
+\sin^2\theta\,d\phi_1^2.
```

A uniform phase displacement $d\phi_0=d\phi_1=d\alpha$ produces
$d\psi=i\psi\,d\alpha$ and has squared length $d\alpha^2$.
It changes the state vector even though it does not change expectation values.
For general leaf probabilities $p_\ell=|x_\ell|^2$, the adopted phase block is
$\mathrm{diag}(p)$.

If the common phase is instead quotiented out, the Fubini–Study convention
(without an extra overall factor) subtracts
$|\langle\psi\mid d\psi\rangle|^2$. Its phase block is

```math
G_{\mathrm{ph}}^{\mathrm{FS}}
=\mathrm{diag}(p)-pp^{\mathsf T}.
```

This block annihilates the all-ones direction. It is a different geometric
object from the diagonal sphere metric used here. An expectation objective
has zero derivative along the common-phase direction in either description;
objective invariance does not make that direction null in the sphere metric.

**Connect to the paper:** [state-vector chart](../science/chart-and-geometry.md#Thmdefinition2),
[metric definition](../science/chart-and-geometry.md#Thmdefinition4), and
[state-sphere conventions](../science/optimizers.md#A3.SSx2).

## 5. A derivative becomes a measurable transition moment

Let $H=H^\dagger$ and $E(\theta)=\langle\psi\mid H\mid\psi\rangle$.
For an active coordinate with $J_j=w_j e_j$,

```math
\partial_j E
=2\mathrm{Re}\langle J_j\mid H\mid\psi\rangle
=2w_j\mathrm{Re}\langle e_j\mid H\mid\psi\rangle.
```

The normalized tangent can be prepared by a shifted split and clamped
ancestor path. For complex phase coordinate $\phi_\ell$, its raw differential
is $i x_\ell\lvert\ell\rangle$, so the normalized tangent is leaf-local when
$x_\ell\ne0$.

Because $\mathrm{Re}\langle\psi\mid e_j\rangle=0$, both signed
branches $\varphi_j^\pm=(\psi\pm e_j)/\sqrt2$ are normalized. Expanding
their expectation values gives

```math
\partial_jE
=w_j\bigl(\langle H\rangle_{\varphi_j^+}
-\langle H\rangle_{\varphi_j^-}\bigr).
```

For the two-qubit example, choose $H=Z\otimes I$. Then
$E=\cos2a$ and $\nabla E=(-2\sin2a,0,0)$. At $a=\pi/4$ the root
branches have energies $-1$ and $+1$, giving derivative $-2$.
The other derivatives vanish because the objective only distinguishes the two
subtree probabilities.

The relative branch phase matters: replacing a tangent by a phase-equivalent
state can change an interference experiment. The signed branch must therefore
be prepared with the convention appearing in the identity.

**Connect to the paper:** [gradient identity](../science/gradient-access.md#Thmtheorem2),
[tangent preparation](../science/gradient-access.md#Thmtheorem3), and
[signed branches](../science/gradient-access.md#S3.SSx1).

## 6. A circuit family is not a statistical budget

Magnitude coordinates at depth $d$ have labels
$j=2^d,\ldots,2^{d+1}-1$. A common indexed construction covers that layer.
The logical setting count is $1+2n$ for the real chart and $1+2(n+1)$ for
the complex chart: a baseline, tangent and branch constructions by magnitude
depth, and one additional complex phase family.

A setting describes a circuit template. Its repeated measurements generate
random observations and labels. In this indexed construction, one labelled
observation contributes to its selected coordinate; a count of templates does
not specify the repetitions required to estimate every coordinate accurately.
The confidence target, label frequencies, objective decomposition and variance
belong to the statistical account. The companion QBP paper changes the readout
mechanism to reuse one record across magnitude coordinates.

**Connect to the paper:** [magnitude batching](../science/gradient-access.md#S3.SSx2),
[phase batching](../science/gradient-access.md#S3.SSx3),
[complete estimator](../science/gradient-access.md#S3.SSx4), and
[controlled-gate ledger](../science/algorithms.md#A1.SSx3).

## 7. Geometry supplies a state-space update

Write $c_j=\partial_jE$. On the regular chart, the sphere gradient is

```math
G=\sum_j\frac{c_j}{g_{jj}}J_j.
```

This expression converts coordinate changes into the actual motion of the
state. For a unit vector $\psi$ and a real-sphere tangent direction $v$ with
$\mathrm{Re}\langle\psi\mid v\rangle=0$, the great-circle update is

```math
\psi(t)=\cos(t\lVert v\rVert)\psi
+\sin(t\lVert v\rVert)\frac{v}{\lVert v\rVert}
\quad (v\ne0).
```

Its norm remains one. Applying the inverse Hopf map to the updated vector
returns circuit coordinates. Conjugate-gradient, L-BFGS and Barzilai–Borwein
schemes differ in their direction memory and step choices; they share this
state-space interpretation. Metric floors and line-search safeguards are
numerical choices that must be distinguished from the exact metric identity.

**Connect to the paper:** [metric lift](../science/optimizers.md#A3.SSx3),
[geodesic transport and conjugate gradients](../science/optimizers.md#A3.SSx4),
[L-BFGS](../science/optimizers.md#A3.SSx5), and
[Barzilai–Borwein](../science/optimizers.md#A3.SSx6).

## 8. Nonlinear objectives require a local observable

Suppose $F(\psi)=f(m_1,\ldots,m_K)$ where
$m_k=\langle\psi\mid A_k\mid\psi\rangle$ and each $A_k$ is Hermitian.
At the current state, define

```math
H_F=\sum_k\left.\frac{\partial f}{\partial m_k}\right|_{m(\psi)}A_k.
```

The first derivative of $F$ is obtained by the same transition-moment identity
with $H_F$. The coefficients are evaluated at the original state and then held
fixed in those measurements. Evaluating the nonlinear function $F$ on the
branch states instead would generally measure another quantity.

For example, for a pure probe and generator $A$, the QFI is
$F_Q=4(\langle A^2\rangle-\langle A\rangle^2)$.
Its local chain-rule observable is
$H_{F_Q}=4A^2-8\langle A\rangle A$. Fixed-readout CFI uses derivatives
of outcome probabilities and its own local assumptions.

**Connect to the paper:** [chain-rule form](../science/metrological-objectives.md#A4.SSx1),
[QFI](../science/metrological-objectives.md#A4.SSx2), and
[fixed-readout CFI](../science/metrological-objectives.md#A4.SSx3).

## 9. Put the mechanisms together

The scientific chain is: orthogonal splits give a universal chart; the same
tree gives an inverse and diagonal metric; active tangents become preparable
states; transition moments supply coordinate gradients; the metric connects
those derivatives to state-sphere motion. The deterministic numerical study
examines specified optimizer behavior using exact costs and gradients. Its
measured performance and the all-size mathematical statements have different
evidence routes.

Use the [four-qubit example](../science/four-qubit-example.md#A2) to join
indexing, gates, tangents and metric factors in one concrete instance. Use
[the numerical study](../science/numerical-studies.md#S4) for its benchmark
conditions, [the proofs](../science/proofs.md#A5) for general identities, and
[reproducibility instructions](../verification/README.md) for executable checks.
The final milestone joins these sources into the complete scientific argument.
