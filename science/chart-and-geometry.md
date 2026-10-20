[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="S2"></a>

# II A universal Hopf variational chart

<a id="S2.p1.1"></a>

We use “Hopf” in the sense of generalized Hopf coordinates from the polyspherical-coordinate literature [[21](../references/paper-bibliography.md#bib.bib4), [22](../references/paper-bibliography.md#bib.bib3)]. In that construction, one recursively splits a Euclidean radius into two orthogonal parts, producing a binary sine–cosine tree with ordinary polar phases at the terminal two-dimensional leaves. The Hopf ansatz imports this architecture into quantum state preparation.

Throughout, $`n`$ denotes the number of qubits and $`N:=2^{n}`$ the Hilbert-space dimension, equivalently the number of leaves in the Hopf tree. We write the imaginary unit as $`\mathrm{i}`$, reserving italic letters such as $`i,j,k`$ for indices.

For $`n`$ qubits, the tree leaves are computational-basis states. Each internal node carries one angle that splits incoming probability mass between its left and right subtrees, so a root-to-leaf path gives one amplitude as a product of sine and cosine factors. The complex version adds one independent phase to each leaf. This structure yields an explicit inverse map, a diagonal pullback metric, and exactly preparable normalized coordinate tangents. Thus the Hopf ansatz is a variational chart for moving through the full real or complex pure-state sphere, with state loading as one component of the same structure.

<a id="Thmdefinition1"></a>

**Definition 1 (Hopf binary tree).**

Fix $`n\in\mathbb{N}`$ and set $`N=2^{n}`$. The Hopf binary tree is the complete rooted binary tree with root at depth $`0`$ and leaves at depth $`n`$.

Each node has a breadth–first tree index

<a id="S2.Ex1"></a>

```math
j\in\{1,2,\ldots,2N-1\},
```

with root $`j=1`$, left child $`2j`$, and right child $`2j+1`$.

- Nodes $`1\leq j\leq N-1`$ are internal nodes and carry split parameters $`\theta_{j}`$.

- Nodes $`N\leq j\leq 2N-1`$ are leaves.

A root-to-leaf path is encoded by bits $`q_{n}\cdots q_{1}`$, read from MSB to LSB: left branch $`0`$, right branch $`1`$. Let

<a id="S2.Ex2"></a>

```math
i=(q_{n}\cdots q_{1})_{2}\in\{0,\ldots,N-1\}.
```

The corresponding leaf has tree index $`j_{\mathrm{leaf}}=N+i`$, labels the computational-basis state

<a id="S2.Ex3"></a>

```math
\vert b_{i}\rangle\equiv\vert q_{n}\cdots q_{1}\rangle,
```

and carries amplitude $`x_{i}`$.

<a id="Thmdefinition2"></a>

**Definition 2 (Hopf coordinate system).**

Fix $`n\in\mathbb{N}`$, set $`N=2^{n}`$, and consider the Hopf binary tree of Definition [1](chart-and-geometry.md#Thmdefinition1). For each computational basis index $`i\in\{0,\ldots,N-1\}`$, write

<a id="S2.Ex4"></a>

```math
i=(q_{n}q_{n-1}\cdots q_{1})_{2},\qquad q_{k}\in\{0,1\},
```

and follow the unique path in the Hopf binary tree from the root to the corresponding leaf.

Initialize $`x_{i}=1`$ and $`j=1`$. For $`k=n,n-1,\ldots,1`$, update

<a id="S2.Ex5"></a>

```math
\begin{aligned}x_{i} & \leftarrow\begin{cases}x_{i}\cos\theta_{j},&q_{k}=0,\\
x_{i}\sin\theta_{j},&q_{k}=1,\end{cases}\end{aligned}
```

```math
\begin{aligned}j & \leftarrow\begin{cases}2j,&q_{k}=0,\\
2j+1,&q_{k}=1.\end{cases}\end{aligned}
```

(A) Real Hopf coordinates. The real Hopf coordinate system uses the internal-node angles $`\theta_{1},\ldots,\theta_{N-1}`$ with depth-dependent ranges

<a id="S2.Ex6"></a>

```math
\begin{aligned}\theta_{j} & \in[0,\tfrac{\pi}{2}], & 1\leq j\leq N/2-1,\end{aligned}
```

```math
\begin{aligned}\theta_{j} & \in[0,2\pi), & N/2\leq j\leq N-1.\end{aligned}
```

The final-layer range $`\theta_{j}\in[0,2\pi)`$ allows the last split in each sibling leaf pair to encode arbitrary real signs. With these ranges, the construction produces

<a id="S2.Ex7"></a>

```math
\mathbf{x}=(x_{0},\ldots,x_{N-1})\in\mathbb{R}^{N},\qquad\sum_{i=0}^{N-1}x_{i}^{2}=1.
```

Thus the real Hopf coordinates give a dimension-matched parameterization of the real unit sphere $`\mathbb{S}^{N-1}`$ by $`N-1`$ continuous parameters. The explicit inverse convention for the signed final-layer angles is stated in Lemma [1](algorithms.md#Thmlemma1).

(B) Complex Hopf coordinates. The complex Hopf coordinate system uses the same internal-node tree, but now all magnitude angles lie in

<a id="S2.Ex8"></a>

```math
\theta_{j}\in[0,\tfrac{\pi}{2}],\qquad 1\leq j\leq N-1,
```

so that the tree produces nonnegative magnitudes. These magnitudes are then augmented by independent leaf phases

<a id="S2.Ex9"></a>

```math
\theta_{N+i}\in[0,2\pi),\qquad i=0,\ldots,N-1,
```

through the assignment

<a id="S2.Ex10"></a>

```math
x_{i}\;\leftarrow\;x_{i}\,e^{\mathrm{i}\theta_{N+i}}.
```

This yields

<a id="S2.Ex11"></a>

```math
\mathbf{x}\in\mathbb{C}^{N},\qquad\sum_{i=0}^{N-1}\vert x_{i}\vert ^{2}=1.
```

The target manifold is the unit sphere in $`\mathbb{C}^{N}`$, equivalently $`\mathbb{S}^{2N-1}`$ as a real manifold. The complex Hopf representation uses $`N-1`$ magnitude angles and $`N`$ leaf phases, giving $`2N-1`$ real parameters, matching the dimension of the normalized state-vector sphere. In the circuit realization (Definition [3](chart-and-geometry.md#Thmdefinition3)), the leaf phases are implemented by promoting the final-layer rotations to $`R_{\mathbb{C}}`$ gates, each of which injects the two phases associated with a sibling leaf pair.

In both the real and complex cases, the Hopf coordinates should be understood as a full sphere parameterization rather than as a single global one-to-one chart. As usual for spherical coordinates, boundary values and zero-subtree weights lead to nonunique coordinate representatives. On the regular set where the relevant subtree weights are nonzero, the inverse map from amplitudes $`\mathbf{x}`$ to Hopf parameters $`\boldsymbol{\theta}`$ is explicit and computable in linear classical time; see Lemma [1](algorithms.md#Thmlemma1) (Appendix [A](algorithms.md#A1.SSx1)).

<a id="Thmdefinition3"></a>

**Definition 3 (Hopf ansatz).**

The Hopf ansatz is the circuit realization of Definition [2](chart-and-geometry.md#Thmdefinition2). Acting on $`\lvert 0\rangle ^{\otimes n}`$, it prepares

<a id="S2.Ex12"></a>

```math
\lvert \psi(\boldsymbol{\theta})\rangle =\sum_{i=0}^{N-1}x_{i}(\boldsymbol{\theta})\lvert b_{i}\rangle ,
```

with amplitudes given exactly by the binary-tree Hopf recursion.

The circuit is generated from the tree as four parallel lists $`\mathtt{Ctrl}`$, $`\mathtt{Anti}`$, $`\mathtt{Targ}`$ and $`\mathtt{Index}`$, encoding positive controls, negative controls, target qubits, and the relevant Hopf parameter index or indices. The masks $`\mathtt{Ctrl}`$ and $`\mathtt{Anti}`$ specify multi-qubit control conditions, $`\mathtt{Targ}`$ selects the single target qubit, and $`\mathtt{Index}`$ selects the tree parameter. Each internal-node parameter corresponds to one gate location in the fixed skeleton. In the complex case, leaf phases are added without extra layers: each final-layer rotation is promoted to a three-parameter gate carrying the parent magnitude angle and the two sibling leaf phases.

The Hopf ansatz uses the following single-qubit gates.

(A) Real Hopf ansatz. In the real case, each gate is a possibly multi-controlled $`R_{y}`$ rotation,

<a id="S2.E1"></a>

```math
R_{y}(\theta)=\begin{bmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\ \cos\theta\end{bmatrix},\tag{1}
```

which corresponds to the standard quantum gate $`R_{y}(2\theta)`$ under the conventional definition $`R_{y}(\varphi)=e^{-\mathrm{i}\varphi Y/2}`$. Its action on the computational basis is

<a id="S2.E2"></a>

<a id="S2.E2X"></a>

<a id="S2.E2Xa"></a>

```math
\begin{aligned}R_{y}(\theta)\lvert 0\rangle & =\cos\theta\lvert 0\rangle +\sin\theta\lvert 1\rangle , \\
R_{y}(\theta)\lvert 1\rangle & =-\sin\theta\lvert 0\rangle +\cos\theta\lvert 1\rangle .\end{aligned}\tag{2}
```

Given angles

<a id="S2.Ex13"></a>

```math
\boldsymbol{\theta}=(\theta_{1},\ldots,\theta_{N-1})
```

in the real Hopf coordinate ranges, Algorithm [1](algorithms.md#algorithm-1) (Appendix [A](algorithms.md#A1.SSx2)) generates the corresponding gate list. Applying the resulting circuit gives

<a id="S2.E3"></a>

```math
\mathrm{HopfReal}(\boldsymbol{\theta})\lvert 0\rangle ^{\otimes n}=\sum_{i=0}^{N-1}x_{i}\lvert b_{i}\rangle ,\tag{3}
```

with real amplitudes $`x_{i}`$.

(B) Complex Hopf ansatz. In the complex case, the magnitude-splitting structure is identical to the real case, while phases are attached at the leaf level. Accordingly, the depth-$`(n-1)`$ rotations are promoted to generalized single-qubit rotations

<a id="S2.E4"></a>

```math
R_{\mathbb{C}}(\theta_{a},\theta_{b},\theta_{c})=\begin{bmatrix}e^{\mathrm{i}\theta_{b}}\cos\theta_{a}&-e^{-\mathrm{i}\theta_{c}}\sin\theta_{a}\\
e^{\mathrm{i}\theta_{c}}\sin\theta_{a}&\ \ e^{-\mathrm{i}\theta_{b}}\cos\theta_{a}\end{bmatrix},\tag{4}
```

with action

<a id="S2.E5"></a>

<a id="S2.E5X"></a>

<a id="S2.E5Xa"></a>

```math
\begin{aligned}R_{\mathbb{C}}(\theta_{a},\theta_{b},\theta_{c})\lvert 0\rangle & =e^{\mathrm{i}\theta_{b}}\cos\theta_{a}\lvert 0\rangle +e^{\mathrm{i}\theta_{c}}\sin\theta_{a}\lvert 1\rangle , \\
R_{\mathbb{C}}(\theta_{a},\theta_{b},\theta_{c})\lvert 1\rangle & =-e^{-\mathrm{i}\theta_{c}}\sin\theta_{a}\lvert 0\rangle +e^{-\mathrm{i}\theta_{b}}\cos\theta_{a}\lvert 1\rangle .\end{aligned}\tag{5}
```

This gate reduces to $`R_{y}(\theta_{a})`$ when $`\theta_{b}=\theta_{c}=0`$. For state preparation only the action on $`\lvert 0\rangle`$ is used; the action on $`\lvert 1\rangle`$ specifies a unitary extension and supports controlled use of the same gate family.

Given parameters

<a id="S2.Ex14"></a>

```math
\boldsymbol{\theta}=(\theta_{1},\ldots,\theta_{2N-1}),
```

Algorithm [2](algorithms.md#algorithm-2) (Appendix [A](algorithms.md#A1.SSx2)) generates the corresponding gate list. Applying the resulting circuit gives

<a id="S2.E6"></a>

```math
\mathrm{HopfComplex}(\boldsymbol{\theta})\lvert 0\rangle ^{\otimes n}=\sum_{i=0}^{N-1}x_{i}\lvert b_{i}\rangle ,\tag{6}
```

with complex amplitudes satisfying $`\sum_{i}\vert x_{i}\vert ^{2}=1`$.

For tangent-state synthesis we also use gate-parameter assignments on the same Hopf skeleton outside the canonical ranges of Definition [2](chart-and-geometry.md#Thmdefinition2); for instance, some internal-node angles are shifted by $`+\pi/2`$. These remain valid physical gate settings. Thus “prepared by the Hopf ansatz” means prepared by the same ordered Hopf skeleton and gate types, not necessarily by a point in the canonical coordinate domain.

<a id="S2.SSx1"></a>

## Hopf geometry on the state sphere

<a id="S2.SSx1.p1.1"></a>

The Hopf chart exposes the optimization geometry of the full-state family in closed form. Its pullback metric is analytic and diagonal in both the real and complex ansatze, so Riemannian or state-sphere descent directions can be formed without hardware estimation or inversion of a dense quantum geometric tensor. In the real case, the Fubini–Study metric reduces here to the usual sphere metric. In the complex case, we keep the global phase as a tangent direction of the normalized state-vector sphere and use the ambient round-sphere metric. This is the metric convention naturally used by optimizer updates that move directly on normalized state vectors [[24](../references/paper-bibliography.md#bib.bib2)].

Throughout this section,

<a id="S2.Ex15"></a>

```math
\lvert \psi(\boldsymbol{\theta})\rangle =\sum_{i=0}^{N-1}x_{i}(\boldsymbol{\theta})\lvert b_{i}\rangle
```

<a id="S2.SSx1.p2.2"></a>

is the Hopf state of Definition [3](chart-and-geometry.md#Thmdefinition3). We treat the real and complex ansatze in parallel; they differ only in the metric convention: pullback Fubini–Study in the real case, and pullback ambient round-sphere metric in the complex case [[24](../references/paper-bibliography.md#bib.bib2)].

<a id="Thmdefinition4"></a>

**Definition 4 (Metrics in parameter space).**

Let $`\boldsymbol{\theta}`$ denote the Hopf parameters and define the raw coordinate tangent associated with $`\theta_{i}`$ by

<a id="S2.Ex16"></a>

```math
\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle :=\partial_{\theta_{i}}\lvert \psi(\boldsymbol{\theta})\rangle .
```

<a id="Thmdefinition4.p2.1"></a>

(A) Real-state Hopf ansatz. Assume $`x_{k}(\boldsymbol{\theta})\in\mathbb{R}`$ and $`\sum_{k=0}^{N-1}x_{k}(\boldsymbol{\theta})^{2}=1`$. We define the parameter-space metric as the pullback of the Fubini–Study metric [[23](../references/paper-bibliography.md#bib.bib6)]:

<a id="S2.E7"></a>

```math
\begin{split}g^{\mathbb{R}}_{i\ell}(\boldsymbol{\theta}):=\mathrm{Re}\!\Big[&\langle \dot{\psi}_{i}(\boldsymbol{\theta})\vert \dot{\psi}_{\ell}(\boldsymbol{\theta})\rangle \\
&-\langle \dot{\psi}_{i}(\boldsymbol{\theta})\vert \psi(\boldsymbol{\theta})\rangle \langle \psi(\boldsymbol{\theta})\vert \dot{\psi}_{\ell}(\boldsymbol{\theta})\rangle \Big].\end{split}\tag{7}
```

Because $`\sum_{k=0}^{N-1}x_{k}^{2}=1`$ implies $`\langle \psi\vert \dot{\psi}_{i}\rangle =\sum_{k=0}^{N-1}x_{k}\,\partial_{\theta_{i}}x_{k}=0`$, the projection term vanishes, and hence

<a id="S2.E8"></a>

```math
g^{\mathbb{R}}_{i\ell}(\boldsymbol{\theta})=\langle \dot{\psi}_{i}(\boldsymbol{\theta})\vert \dot{\psi}_{\ell}(\boldsymbol{\theta})\rangle =\sum_{k=0}^{N-1}\frac{\partial x_{k}}{\partial\theta_{i}}\frac{\partial x_{k}}{\partial\theta_{\ell}}.\tag{8}
```

(B) Complex-state Hopf ansatz. Assume $`x_{k}(\boldsymbol{\theta})\in\mathbb{C}`$ and $`\sum_{k=0}^{N-1}\vert x_{k}(\boldsymbol{\theta})\vert ^{2}=1,`$ so that $`\psi(\boldsymbol{\theta})\in\mathbb{S}^{2N-1}\subset\mathbb{C}^{N}.`$ We define the parameter-space metric as the pullback of the ambient round-sphere metric on $`\mathbb{S}^{2N-1}`$:

<a id="S2.E9"></a>

```math
g^{\mathbb{C}}_{i\ell}(\boldsymbol{\theta}):=\mathrm{Re}\,\langle \dot{\psi}_{i}(\boldsymbol{\theta})\vert \dot{\psi}_{\ell}(\boldsymbol{\theta})\rangle .\tag{9}
```

<a id="Thmdefinition4.p3.2"></a>

This convention treats the global phase direction as a genuine tangent direction of the unit sphere. If desired, one may instead pass to the projective Fubini–Study metric by subtracting the standard projection term $`\langle \dot{\psi}_{i}\vert \psi\rangle \langle \psi\vert \dot{\psi}_{\ell}\rangle`$. Here we use the ambient-sphere convention because the optimizer-side state updates are formulated on normalized state vectors rather than on projective Hilbert space [[24](../references/paper-bibliography.md#bib.bib2)].

Our first structural result shows that the Hopf chart diagonalizes the pullback metric in closed form.

<a id="Thmtheorem1"></a>

**Theorem 1 (Diagonal metric of the Hopf ansatz).**

Let $`\lvert \psi(\boldsymbol{\theta})\rangle`$ be prepared by the Hopf ansatz in Definition [3](chart-and-geometry.md#Thmdefinition3).

(A) Real-state version. For $`\boldsymbol{\theta}=(\theta_{1},\ldots,\theta_{N-1})`$, the metric $`\boldsymbol{g}^{\mathbb{R}}`$ of ([8](chart-and-geometry.md#S2.E8)) is diagonal and $`g^{\mathbb{R}}_{1,1}=1`$. For $`i>1`$, let $`d_{i}=\lfloor\log_{2}i\rfloor`$, and let $`q_{1}\cdots q_{d_{i}}`$ be the length-$`d_{i}`$ binary expansion of $`i-2^{d_{i}}`$, padded with leading zeros and written most-significant bit first. Let $`s_{1},\ldots,s_{d_{i}}`$ be the ancestor indices from the root $`s_{1}=1`$ to the parent of node $`i`$. Then

<a id="S2.E10"></a>

```math
g^{\mathbb{R}}_{i,i}(\boldsymbol{\theta})=\prod_{m=1}^{d_{i}}\begin{cases}\cos^{2}\theta_{s_{m}},&q_{m}=0,\\
\sin^{2}\theta_{s_{m}},&q_{m}=1.\end{cases}\tag{10}
```

(B) Complex-state version. For parameters ordered as $`\theta_{1},\ldots,\theta_{N-1}`$ (magnitudes) followed by $`\theta_{N},\ldots,\theta_{2N-1}`$ (phases), the round-sphere pullback metric $`\boldsymbol{g}^{\mathbb{C}}`$ of Definition [4](chart-and-geometry.md#Thmdefinition4)(B) is diagonal with $`g^{\mathbb{C}}_{1,1}=1`$ and

<a id="S2.E11"></a>

```math
g^{\mathbb{C}}_{i,i}(\boldsymbol{\theta})=\begin{cases}g^{\mathbb{R}}_{i,i}(\boldsymbol{\theta}),&1\leq i\leq N-1,\\[4.0pt]
\vert x_{i-N}(\boldsymbol{\theta})\vert ^{2},&N\leq i\leq 2N-1.\end{cases}\tag{11}
```

Thus, in the complex Hopf ansatz, the diagonal metric has two blocks: a magnitude block $`1\leq i\leq N-1`$, identical to the real Hopf metric, and a phase block $`N\leq i\leq 2N-1`$, whose entries are the leaf occupations $`\vert x_{i-N}(\boldsymbol{\theta})\vert ^{2}`$. Hence phase directions are weighted by local squared amplitudes, while magnitude directions inherit the real Hopf tree metric.

Proofs of the main results are collected in Appendix [E](proofs.md#A5).

Theorem [1](chart-and-geometry.md#Thmtheorem1) separates the variational geometry into independent coordinates. For a magnitude parameter, $`g_{i,i}`$ is the squared norm of the amplitudes in the subtree rooted at node $`i`$, equivalently the probability mass entering that node; for a phase parameter, it is the occupation probability of the corresponding leaf. Since the metric is analytic and diagonal, Riemannian descent requires only componentwise gradient rescaling, not hardware estimation or inversion of a dense quantum geometric tensor. Appendix [C](optimizers.md#A3) gives the complementary Jacobian formulas and state-sphere maps used by the optimizer family. The same tangent-space structure also enables gradient access: the next section shows that Hopf coordinate tangents can be prepared exactly and organized coherently across the tree.
