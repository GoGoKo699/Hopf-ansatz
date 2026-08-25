# Engineering guide

This guide specifies the implementable Hopf-ansatz interface without reproducing
the paper's proofs, four-qubit walkthrough, or extended motivation.

It is organized around the questions an engineer needs to answer:

1. How are coordinates ordered and constrained?
2. How are state vectors mapped to and from coordinates?
3. How is the native gate schedule generated?
4. How are Jacobians, metrics, and normalized tangents obtained?
5. How are coordinate gradients lifted to the state sphere?
6. Which substitutions and numerical safeguards are valid?

## 1. Input and output contracts

Let:

```text
n = number of system qubits
N = 2**n
```

The core implementation is `hopf_utils.py`.

### Real chart

Input:

```text
theta.shape == (N - 1,)
```

Output:

```text
psi.shape == (N,)
psi is real and normalized
```

### Complex chart

Input:

```text
theta.shape == (2*N - 1,)
theta[:N-1]   = magnitude angles
theta[N-1:]   = N leaf phases
```

Output:

```text
psi.shape == (N,)
psi is complex and normalized
```

Always pass `case="real"` or `case="complex"` explicitly. A parameter-vector
length alone is ambiguous: a complex `n`-qubit vector has the same length as a
real `(n + 1)`-qubit vector.

```python
n = infer_n_from_theta(theta, case="real")
n = infer_n_from_theta(theta, case="complex")
```

## 2. Tree, basis, and coordinate conventions

### Breadth-first tree indices

The root is internal node `1`. For internal node `j`:

```text
left child  = 2*j
right child = 2*j + 1
```

Internal nodes are `1, ..., N - 1`. Leaves have tree indices `N, ..., 2*N - 1`.
Python arrays store an internal-node parameter at `theta[j - 1]`.

### Basis order

Basis states are written

```math
\lvert q_n\cdots q_1\rangle.
```

The path is read from most-significant bit to least-significant bit. Qibo qubit
index `0` is the most-significant basis bit in this repository's circuit
convention.

The schedule's integer masks use bit position `0` for the least-significant
bit. The helper functions in `hopf_utils.py` convert these masks to Qibo qubit
indices. Do not mix raw mask bit positions with Qibo indices.

### Coordinate order and domains

Real:

```text
(theta_1, ..., theta_(N-1))
```

Canonical ranges:

```text
nodes 1, ..., N/2 - 1:  [0, pi/2]
nodes N/2, ..., N - 1: [0, 2*pi)
```

The extended final-layer range carries arbitrary real signs.

Complex:

```text
(theta_1, ..., theta_(N-1), theta_N, ..., theta_(2N-1))
```

Canonical ranges:

```text
all N - 1 magnitude angles: [0, pi/2]
all N leaf phases:           [0, 2*pi)
```

Tangent-state assignments may intentionally use angles outside these canonical
coordinate ranges. They remain valid physical gate parameters on the same
skeleton.

## 3. Forward coordinate map

Use:

```python
psi_real = vector_from_theta(theta_real, "real")
psi_complex = vector_from_theta(theta_complex, "complex")
```

At each depth, every incoming subtree amplitude is split into a left child by a
cosine factor and a right child by a sine factor.

For a leaf with path bits `q_n ... q_1`, the real amplitude is a product of one
factor per path node:

```math
x_\ell
=
\prod_{k=1}^{n}
\begin{cases}
\cos\theta_{j_k}, & q_{n+1-k}=0,\\
\sin\theta_{j_k}, & q_{n+1-k}=1.
\end{cases}
```

The complex map multiplies the nonnegative tree magnitude by its leaf phase:

```math
x_\ell=r_\ell e^{i\theta_{N+\ell}}.
```

The implementation constructs levels iteratively and therefore evaluates the
complete state in `O(N)` arithmetic and `O(N)` output memory.

## 4. Inverse coordinate map

Use:

```python
theta_real = theta_from_vector(real_state, "real")
theta_complex = theta_from_vector(complex_state, "complex")
```

The state length must be a power of two. Normalize the state before calling the
inverse if normalization is not already guaranteed.

### Magnitude angles

For an internal node `j`, let `S_L(j)` and `S_R(j)` be the Euclidean norms of
its left and right leaf subtrees. On a nonzero subtree:

```math
\theta_j=\mathrm{atan2}\!\left(S_R(j),S_L(j)\right).
```

This produces an angle in `[0, pi/2]` when both arguments are nonnegative.

### Real final-layer signs

For the real sibling pair `(x_(2k), x_(2k+1))`:

```math
\theta_{N/2+k}
=
\mathrm{atan2}\!\left(x_{2k+1},x_{2k}\right)\bmod 2\pi.
```

This is why the real final internal layer must not be clipped to `[0, pi/2]`.

### Complex phases

For a nonzero complex leaf:

```math
\theta_{N+\ell}=\mathrm{arg}(x_\ell)\bmod 2\pi.
```

### Zero conventions

The inverse is nonunique at chart boundaries. The implementation chooses:

```text
zero magnitude subtree -> angle 0
zero real sibling pair -> final angle 0
zero complex leaf      -> phase 0
```

Validate inverse maps by reconstructing the state, not by demanding equality
with an arbitrary original coordinate vector.

### Numerical clipping for optimizer traces

The core real-chart safeguard is:

```python
theta = clip_theta_hopf_real(
    theta,
    n=n,
    eps_internal=1e-6,
    eps_sing=1e-9,
)
```

It clips non-final angles into the regular interior, wraps the final real layer
to `[0, 2*pi)`, and nudges that layer away from exact sine or cosine zeros.
The complex experiment scripts use the analogous policy for their magnitude
blocks and wrap all leaf phases.

Clipping is an optimizer-stability convention. It is not part of the exact
forward or inverse definition, and it should be disabled when deliberately
testing boundary identities or singular coordinates.

## 5. Native gate schedule

Use:

```python
ctrl, anti, targ, index = gates_order(n, case="real")
ctrl, anti, targ, index = gates_order(n, case="complex")
```

All four lists have one entry per internal-node gate. For each entry:

- `Ctrl` is the positive-control bitmask;
- `Anti` is the negative-control bitmask;
- `Targ` is a power-of-two target bitmask; and
- `Index` is a one-based parameter index or a three-index complex leaf tuple.

The schedule is deterministic. It is generated from Hamming-weight layers and
the `find_pairs` pairing rule. Do not sort or reorder the four lists
independently.

### Negative controls

A negative control triggers on `0`. The Qibo builder implements it by applying
`X` before and after the controlled operation. In the assigned CNOT model,
these single-qubit `X` gates do not change the CNOT charge.

### Real gate convention

The paper's gate is

```math
R_y(\theta)=
\begin{pmatrix}
\cos\theta & -\sin\theta\\
\sin\theta & \cos\theta
\end{pmatrix}.
```

Qibo's conventional half-angle gate must therefore receive:

```python
gates.RY(target, 2.0 * theta)
```

Using `RY(theta)` produces the wrong state.

### Complex leaf gate

The promoted final-layer gate is

```math
R_{\mathbb{C}}(\theta_a,\theta_b,\theta_c)
=
\begin{pmatrix}
e^{i\theta_b}\cos\theta_a & -e^{-i\theta_c}\sin\theta_a\\
e^{i\theta_c}\sin\theta_a & e^{-i\theta_b}\cos\theta_a
\end{pmatrix}.
```

The three one-based indices identify the parent magnitude angle and the two
sibling leaf phases. The matrix helper in `hopf_utils.py` is named `SU2`; the
paper and documentation refer to the same gate as `R_C` or
`R_{\mathbb{C}}`.

## 6. Qibo circuit construction

Use:

```python
state = qibo_circuit(theta, case="real")
state = qibo_circuit(theta, case="complex")

circuit = qibo_circuit(
    theta,
    case="real",
    return_circuit=True,
)
```

Qibo is optional. If it is not installed, the function raises a runtime error
rather than silently substituting another backend.

### Minimal tangent circuits

For tangent assignments, most schedule entries may be inactive or equal to the
identity. Use:

```python
circuit = qibo_circuit(
    tangent_theta,
    case,
    minimal=True,
    return_circuit=True,
)
```

The pruning logic is conservative:

- it tracks a superset of the current computational-basis support;
- it removes a controlled operation only when its controls cannot be reached;
- it removes `+I` operations but does not treat controlled `-I` as an identity;
- it may leave redundant gates rather than risk changing the state.

This is a state-preparation optimization, not a proof of globally minimal gate
count.

## 7. Jacobian

Use:

```python
J = jacobian(theta, case="real")
J = jacobian(theta, case="complex")
```

Columns are raw coordinate derivatives:

```math
J_{:,i}=\frac{\partial\psi}{\partial\theta_i}.
```

Magnitude columns are evaluated through the shifted tangent-state construction
multiplied by the signed incoming amplitude. This avoids unstable `tan(theta)`
or `cot(theta)` rewrites at exact branch boundaries.

Complex phase columns are leaf local:

```math
\frac{\partial\psi}{\partial\theta_{N+\ell}}
=
i x_\ell\lvert b_\ell\rangle.
```

Dense Jacobians are useful as references and at small sizes. Production
optimization scripts use tree contractions instead of materializing the full
matrix at every step.

## 8. Diagonal metric

Use:

```python
g = metric_diagonal(theta, case="real")
g = metric_diagonal(theta, case="complex")
```

For a real magnitude node, the diagonal entry is the squared amplitude entering
that node. It is the product of ancestor branch probabilities.

For the complex chart:

```math
g^{\mathbb{C}}_{i,i}=
\begin{cases}
g^{\mathbb{R}}_{i,i}, & 1\leq i\leq N-1,\\
\lvert x_{i-N}\rvert^2, & N\leq i\leq 2N-1.
\end{cases}
```

No quantum metric estimation or dense inversion is required by the reference
optimization pipeline. The complex implementation keeps the uniform global
phase as a genuine tangent direction of the normalized state-vector sphere; it
does not quotient the chart to projective Hilbert space.

### Metric zeros

A zero metric entry means the raw coordinate differential vanishes at that
representative. Do not divide by it to manufacture a normalized tangent or a
natural-gradient step. The optimization scripts use a documented small floor
only as a numerical safeguard when lifting near a chart boundary.

## 9. Normalized tangent-state assignments

Use:

```python
tangent_theta = theta_hopf_tangent_state(theta, i, case)
tangent_state = vector_from_theta(tangent_theta, case)
```

Parameter index `i` is one based.

### Magnitude coordinate

For node `i`:

1. set the target to `theta_i + pi/2`;
2. clamp a left-path ancestor to `0`;
3. clamp a right-path ancestor to `pi/2`;
4. retain strict descendant magnitude angles in the selected subtree;
5. in the complex chart, retain phases only on leaves of that subtree; and
6. set all unrelated entries to zero.

The result prepares the normalized derivative direction whenever the incoming
metric factor is nonzero.

### Complex phase coordinate

For phase coordinate `N + ell`:

1. clamp all magnitude ancestors to the path selecting leaf `ell`;
2. set unrelated magnitudes to zero;
3. set the selected phase to its original value plus `pi/2`; and
4. set all other phases to zero.

This prepares `i * exp(i*phase) |ell>` with unit norm.

### Qibo convenience wrapper

```python
state_or_circuit = qibo_tangent_state(
    theta,
    i,
    case,
    minimal=True,
    return_circuit=False,
)
```

## 10. Coordinate gradients

For a scalar objective `C(psi)` with ambient gradient `h`, the coordinate
gradient is

```math
\nabla_{\theta}C
=
\mathrm{Re}\!\left(J^{\dagger}h\right)
```

for the complex chart, with the corresponding real formula in the real chart.

The stress-test scripts provide boundary-safe linear-time routines:

- `hopf_coordinate_gradient` in `hopf_data.py` and `adam_data.py` for real
  coordinates;
- `complex_hopf_coordinate_gradient` in `hopf_complex.py` for complex
  coordinates.

These functions are application helpers rather than part of a packaged public
API. When reusing them, preserve the exact ambient-gradient convention of the
calling objective.

## 11. State-sphere metric lift

Given coordinate gradient `c_i` and metric diagonal `g_i`, the state-sphere
gradient is

```math
G
=
J g^{-1}\nabla_{\theta}C
=
\sum_i\frac{c_i}{g_{i,i}}\frac{\partial\psi}{\partial\theta_i}.
```

The real implementation evaluates this without a dense Jacobian by visiting
only the `n` path nodes contributing to each leaf. The resulting classical cost
is `O(n*N)` with `O(N)` working memory.

For the full complex chart on its regular set, the scripts use the equivalent
direct projection of the ambient gradient to the real tangent space of the
complex unit sphere.

After a state-sphere update, return to executable circuit coordinates with:

```python
next_theta = theta_from_vector(next_state, case)
```

## 12. Signed branch gradient identity

For a normalized tangent state `e_i`, define

```math
\lvert\varphi_i^{(s)}\rangle
=
\frac{\lvert\psi\rangle+s\lvert e_i\rangle}{\sqrt{2}},
\qquad s\in\{+1,-1\}.
```

For a Hermitian observable `H`:

```math
\partial_{\theta_i}E
=
2\sqrt{g_{i,i}}\,s
\left[
E_{\varphi_i^{(s)}}-
\frac{E_{\psi}+E_i^{\mathrm{tan}}}{2}
\right].
```

If both signs are estimated:

```math
\partial_{\theta_i}E
=
\sqrt{g_{i,i}}
\left(E_{\varphi_i^{(+)}}-E_{\varphi_i^{(-)}}\right).
```

`finite_shot_sanity_check.py` uses the symmetric form. `VQE_qibo.py` uses the
single-sign-centered form and averages the available sign estimates.

The relative phase between the baseline and tangent preparation branches is
load bearing. A different controlled-preparation convention can rotate or flip
the inferred transition moment.

## 13. Layerwise organization

Magnitude node indices at depth `d` are:

```text
2**d, ..., 2**(d + 1) - 1
```

The reduced layer label is `r = i - 2**d` and requires `d` index qubits. The
complex phase block is one additional `n`-qubit indexed layer.

The logical compiled-setting counts are:

```math
N_{\mathrm{set}}^{\mathbb{R}}=1+2n,
\qquad
N_{\mathrm{set}}^{\mathbb{C}}=1+2(n+1).
```

These counts include one baseline setting, one tangent setting per layer, and
one branch setting per layer. They do not specify the number of measurement
repetitions needed for a target accuracy.

`VQE_qibo.py` provides a transparent `n = 4` implementation with explicit
index-controlled preparation. For a more developed exact-logical reverse
implementation, use the companion Hopf-QBP repository.

## 14. Geometry-native optimizer pipeline

The real and complex stress-test scripts implement the following common
pipeline:

```text
state
-> inverse Hopf map
-> objective and coordinate gradient
-> analytic diagonal-metric lift
-> state-sphere search direction
-> exact great-circle update
-> inverse Hopf map for the next circuit
```

### Hopf-EGT-CG

- exact sphere exponential map;
- exact geodesic vector transport;
- hybrid conjugate-gradient memory;
- strong-Wolfe geodesic line search with bounded fallback.

### Hopf Riemannian L-BFGS

- transported curvature pairs;
- tangent-space two-loop recursion;
- positive-curvature filtering;
- geodesic line search.

### Hopf Riemannian Barzilai-Borwein

- transported spectral step estimate;
- BB1, BB2, or alternating choice;
- clipped step range;
- geodesic line-search safeguard.

The quantum-information interface assumed by these three modes is only cost and
Hopf-coordinate-gradient access. Line-search counters are not shot counts.

## 15. Assigned CNOT ledger

Use:

```bash
python hopf_gate_count.py --nmin 2 --nmax 10 --out hopf_gate_count.pdf
```

For `m` controls, the assigned controlled-`R_y` cost is:

```math
c_{R_y}(m)=
\begin{cases}
0, & m=0,\\
2^{m+1}-2, & 1\leq m\leq4,\\
16(m+1)-40, & m\geq5.
\end{cases}
```

The promoted complex gate uses:

```math
c_{R_{\mathbb{C}}}(m)=
\begin{cases}
0, & m=0,\\
2^{m+1}-2, & 1\leq m\leq4,\\
20(m+1)-38, & m\geq5\ \text{and}\ m+1\ \text{is odd},\\
20(m+1)-42, & m\geq5\ \text{and}\ m+1\ \text{is even}.
\end{cases}
```

The script compares generated totals with the closed-form schedule counts. It
counts neither Hamiltonian readout nor device-specific compilation.

## 16. Minimal examples

### Real round trip, metric, and tangent

```python
import numpy as np

from hopf_utils import (
    jacobian,
    metric_diagonal,
    theta_from_vector,
    theta_hopf_tangent_state,
    vector_from_theta,
)

n = 3
N = 2**n
rng = np.random.default_rng(7)

state = rng.normal(size=N)
state /= np.linalg.norm(state)

theta = theta_from_vector(state, "real")
reconstructed = vector_from_theta(theta, "real")
print(np.linalg.norm(reconstructed - state))

J = jacobian(theta, "real")
g = metric_diagonal(theta, "real")

i = 3  # one-based Hopf coordinate index
tangent_theta = theta_hopf_tangent_state(theta, i, "real")
normalized_tangent = vector_from_theta(tangent_theta, "real")
raw_tangent = J[:, i - 1]

if g[i - 1] > 0.0:
    print(np.linalg.norm(normalized_tangent - raw_tangent / np.sqrt(g[i - 1])))
```

### Complex round trip

```python
import numpy as np

from hopf_utils import theta_from_vector, vector_from_theta

n = 3
N = 2**n
rng = np.random.default_rng(11)

state = rng.normal(size=N) + 1j * rng.normal(size=N)
state /= np.linalg.norm(state)

theta = theta_from_vector(state, "complex")
reconstructed = vector_from_theta(theta, "complex")
print(np.linalg.norm(reconstructed - state))
```

### Build a Qibo circuit

```python
from hopf_utils import qibo_circuit, theta_from_vector

# `state` is a normalized real vector of power-of-two length.
theta = theta_from_vector(state, "real")
circuit = qibo_circuit(theta, "real", return_circuit=True)
result = circuit(nshots=0)
prepared = result.state()
```

### Assigned gate schedule

```python
from hopf_utils import gates_order

ctrl, anti, targ, index = gates_order(4, "complex")
for specification in zip(ctrl, anti, targ, index):
    print(specification)
```

## 17. Safe adaptation rules

### Replacing the backend

Preserve:

- basis-bit significance;
- mask-to-qubit translation;
- negative-control semantics;
- `RY(2*theta)` normalization;
- the exact `R_C` matrix;
- one-based `Index` entries; and
- gate order.

Add statevector parity tests for both charts and for tangent assignments.

### Replacing the inverse map

Test reconstructed states, including:

- generic interior vectors;
- signed real sibling pairs;
- zero subtrees;
- zero complex leaves; and
- global sign or phase conventions.

Do not compare coordinates at nonunique boundary points without first fixing a
common convention.

### Replacing a tangent compiler

Check both:

1. equality with the normalized Jacobian column when the metric is nonzero; and
2. the prepared support and branch phase needed by the gradient estimator.

State fidelity alone can hide a sign or global phase that matters inside an
interference circuit.

### Adding a new objective

Document the ambient-gradient convention. For expectation-value chain rules,
build the Hermitian chain-rule observable from coefficients evaluated at the
current state, then hold those coefficients fixed while estimating transition
moments. Do not substitute the nonlinear cost of a branch state for the required
linear transition observable.

### Adding shot sampling

State explicitly whether `S` means:

- shots per sign;
- shots per label;
- shots per Pauli term;
- shots per compiled setting; or
- total shots for a complete gradient.

The existing scripts use different, deliberately local conventions for their
specific safeguards.

## 18. Common failure modes

| Failure | Consequence |
|---|---|
| Inferring the chart solely from parameter-vector length | Confuses complex `n` qubits with real `n + 1` qubits. |
| Clipping the real final layer to `[0, pi/2]` | Removes arbitrary real signs. |
| Using `RY(theta)` in Qibo | Introduces a factor-of-two error in every amplitude. |
| Reversing bit significance | Corrupts target masks, controls, and leaf labels. |
| Ignoring negative controls | Applies rotations to unintended basis sectors. |
| Treating controlled `-I` as an identity | Removes a conditional phase. |
| Comparing inverse coordinates at a zero subtree | Flags harmless coordinate nonuniqueness as an error. |
| Dividing by a zero metric entry | Creates an artificial singular normalized tangent. |
| Treating a tangent state's global phase as irrelevant inside a branch interferometer | Flips or rotates the inferred transition moment. |
| Calling assigned CNOT counts hardware counts | Misstates the resource model. |
| Treating the `n = 4` Qibo toy as the general optimized compiler | Overstates the implementation scope. |
| Interpreting line-search work counters as circuit shots | Mixes classical optimizer work with measurement cost. |

## 19. Validation workflow for an engineering change

1. Run `python hopf_utils.py`.
2. Test real and complex forward-after-inverse round trips on interior and
   boundary cases.
3. Compare the proposed Jacobian or gradient routine with the dense reference
   at small `n`.
4. Compare every nonzero normalized tangent with its Jacobian column.
5. If the gate schedule changes, rerun `hopf_gate_count.py` and Qibo parity.
6. Run the small `VQE_qibo.py` smoke path with the intended backend.
7. Run `finite_shot_sanity_check.py` if branch-state or measurement conventions
   changed.
8. Run `hopf_complex.py --quick` if complex geometry changed.
9. Regenerate diagnostic data before changing any reported optimization result.
10. Record whether the change affects exact identities, numerical safeguards,
    assigned resources, or only documentation.
