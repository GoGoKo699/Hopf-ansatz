[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="S3"></a>

# III Hopf gradient access and scaling

We next address gradient access. In the Hopf chart, coordinate tangents are exactly synthesizable quantum states within the same circuit family. This makes differential information directly accessible and permits coherent magnitude-layer batching, together with one indexed batch for the complex phase block.

<a id="Thmdefinition5"></a>

**Definition 5 (Normalized coordinate tangent).**

Let $`g_{i\ell}(\boldsymbol{\theta})`$ be the parameter-space metric of Definition [4](chart-and-geometry.md#Thmdefinition4). Whenever $`g_{i,i}(\boldsymbol{\theta})>0`$, define

<a id="S3.E12"></a>

```math
\lvert e_{i}(\boldsymbol{\theta})\rangle :=\frac{\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle }{\sqrt{g_{i,i}(\boldsymbol{\theta})}}.\tag{12}
```

Thus $`\lvert \dot{\psi}_{i}\rangle`$ denotes the raw coordinate tangent, whereas $`\lvert e_{i}\rangle`$ denotes its normalized tangent state, with $`\|\lvert e_{i}(\boldsymbol{\theta})\rangle \|=1`$.

If $`g_{i,i}(\boldsymbol{\theta})=0`$, the corresponding coordinate lies on a boundary stratum of the Hopf parameterization, such as a zero-weight subtree or zero-occupation leaf, and the normalized coordinate tangent is not part of the regular chart. The identities below are therefore stated on the regular set; boundary cases are handled by the inverse-map convention of Appendix [A](algorithms.md#A1.SSx1) or by continuity.

With this normalization, the gradient formula takes a uniform tangent-state form.

<a id="Thmtheorem2"></a>

**Theorem 2 (Gradient estimation for VQE).**

Let $`\lvert \psi(\boldsymbol{\theta})\rangle =U_{\boldsymbol{\theta}}\lvert 0\rangle ^{\otimes n}`$ be a differentiable family of pure states, and let

<a id="S3.Ex17"></a>

```math
E(\boldsymbol{\theta})=\langle \psi(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle
```

for a Hermitian Hamiltonian $`H`$. For each parameter $`\theta_{i}`$,

<a id="S3.E13"></a>

```math
\partial_{\theta_{i}}E(\boldsymbol{\theta})=2\sqrt{g_{i,i}(\boldsymbol{\theta})}\;\mathrm{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert \,H\,\lvert \psi(\boldsymbol{\theta})\rangle \right].\tag{13}
```

Here $`\lvert e_{i}(\boldsymbol{\theta})\rangle`$ denotes the normalized coordinate tangent of Definition [5](gradient-access.md#Thmdefinition5). The formula ([13](gradient-access.md#S3.E13)) holds for both the real-state and complex-state versions of the Hopf ansatz, with the metric $`g_{i,i}`$ chosen according to Definition [4](chart-and-geometry.md#Thmdefinition4).

The ansatz-specific step is exact preparation of these normalized tangents.

<a id="Thmtheorem3"></a>

**Theorem 3 (Normalized coordinate tangents as quantum states (Hopf)).**

Let

<a id="S3.Ex18"></a>

```math
\lvert \psi(\boldsymbol{\theta})\rangle =\begin{cases}\mathrm{HopfReal}(\boldsymbol{\theta})\lvert 0\rangle ^{\otimes n},\\
\mathrm{HopfComplex}(\boldsymbol{\theta})\lvert 0\rangle ^{\otimes n},\end{cases}
```

with parameters arranged on the Hopf binary tree as in Definitions [2](chart-and-geometry.md#Thmdefinition2) and [3](chart-and-geometry.md#Thmdefinition3). For any parameter index $`i`$ with $`g_{i,i}(\boldsymbol{\theta})>0`$, there exists an explicit gate-parameter assignment $`\boldsymbol{\theta}^{(i)}`$ on the same ordered Hopf circuit skeleton such that

<a id="S3.E14"></a>

```math
\lvert e_{i}(\boldsymbol{\theta})\rangle =\begin{cases}\mathrm{HopfReal}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n},\\
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n},\end{cases}.\tag{14}
```

Here $`\mathrm{HopfReal}`$ and $`\mathrm{HopfComplex}`$ are understood according to the off-chart gate-assignment convention stated above.

Fix a magnitude-node index $`i\in\{1,\ldots,N-1\}`$ and define $`d_{i}=\lfloor\log_{2}i\rfloor`$. Let $`s_{1},\ldots,s_{d_{i}}`$ be the ancestor indices on the unique path from the root $`s_{1}=1`$ to the parent of node $`i`$, and let $`q_{1}\cdots q_{d_{i}}\in\{0,1\}^{d_{i}}`$ be the corresponding left/right bits (MSB first), where $`q_{m}=0`$ for a left branch and $`q_{m}=1`$ for a right branch. Define

<a id="S3.Ex19"></a>

```math
A_{L}^{(i)}:=\{s_{m}\mid q_{m}=0\},\qquad A_{R}^{(i)}:=\{s_{m}\mid q_{m}=1\},
```

let $`D^{(i)}`$ denote the set of strict descendants of node $`i`$ (internal nodes only), and let $`\mathcal{L}(i)`$ be the set of leaves of the subtree rooted at $`i`$.

(A) Magnitude parameters (real and complex cases). This construction applies identically to the real Hopf gate skeleton and to the magnitude block of the complex Hopf gate skeleton. The resulting assignments $`\boldsymbol{\theta}^{(i)}`$ are gate settings for that skeleton and need not belong to the canonical Hopf chart domain.

Define $`\boldsymbol{\theta}^{(i)}`$ by

<a id="S3.EGx1"></a>

<a id="S3.E15"></a>

```math
\begin{aligned}\theta^{(i)}_{i} & =\theta_{i}+\tfrac{\pi}{2},\end{aligned}\tag{15}
```

<a id="S3.E16"></a>

```math
\begin{aligned}\theta^{(i)}_{a} & =0,\qquad a\in A_{L}^{(i)},\end{aligned}\tag{16}
```

<a id="S3.E17"></a>

```math
\begin{aligned}\theta^{(i)}_{a} & =\tfrac{\pi}{2},\qquad a\in A_{R}^{(i)},\end{aligned}\tag{17}
```

<a id="S3.E18"></a>

```math
\begin{aligned}\theta^{(i)}_{d} & =\theta_{d},\qquad d\in D^{(i)},\end{aligned}\tag{18}
```

<a id="S3.E19"></a>

```math
\begin{aligned}\theta^{(i)}_{j} & =0,\qquad\text{otherwise},\qquad 1\leq j\leq N-1.\end{aligned}\tag{19}
```

In the complex case, the phase block is further set to

<a id="S3.E20"></a>

```math
\theta^{(i)}_{N+\ell}=\begin{cases}\theta_{N+\ell},&\ell\in\mathcal{L}(i),\\
0,&\text{otherwise}.\end{cases}\tag{20}
```

Then

<a id="S3.E21"></a>

<a id="S3.E21X"></a>

<a id="S3.E21Xa"></a>

```math
\begin{aligned}\mathrm{HopfReal}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n} & =\lvert e_{i}(\boldsymbol{\theta})\rangle , \\
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n} & =\lvert e_{i}(\boldsymbol{\theta})\rangle .\end{aligned}\tag{21}
```

If node $`i`$ lies at depth $`d`$ (root at depth $`0`$), the resulting state is supported only on the $`N/2^{d}`$ computational basis states corresponding to the leaves of the subtree rooted at $`i`$. Hence its support size is at most $`N/2^{d}`$, with equality whenever its amplitude on every leaf of that subtree is nonzero.

(B) Phase parameters (complex-state Hopf ansatz only). Let $`i=N+\ell`$ with $`0\leq\ell\leq N-1`$. Write

<a id="S3.Ex20"></a>

```math
x_{\ell}(\boldsymbol{\theta})=r_{\ell}(\boldsymbol{\theta})\,e^{\mathrm{i}\theta_{N+\ell}},\qquad r_{\ell}(\boldsymbol{\theta})\geq 0.
```

Let $`(t_{1},\ldots,t_{n})`$ be the internal nodes encountered along the root-to-leaf path of leaf $`\ell`$, with corresponding bits $`(b_{1},\ldots,b_{n})`$ (MSB first).

Define $`\boldsymbol{\theta}^{(i)}`$ by clamping the magnitude block so that $`\vert x_{\ell}\vert =1`$:

<a id="S3.E22"></a>

<a id="S3.E22X"></a>

<a id="S3.E22Xa"></a>

```math
\begin{aligned}\theta^{(i)}_{t_{m}} & =\begin{cases}0,&b_{m}=0,\\
\tfrac{\pi}{2},&b_{m}=1,\end{cases} \\
\theta^{(i)}_{j} & =0,\qquad j\notin\{t_{1},\ldots,t_{n}\},\quad 1\leq j\leq N-1,\end{aligned}\tag{22}
```

and in the phase block set

<a id="S3.E23"></a>

```math
\theta^{(i)}_{N+\ell}=\theta_{N+\ell}+\tfrac{\pi}{2},\qquad\theta^{(i)}_{N+m}=0\quad\text{for }m\neq\ell.\tag{23}
```

Then

<a id="S3.Ex21"></a>

```math
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n}=\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle ,
```

since $`\lvert \dot{\psi}_{N+\ell}(\boldsymbol{\theta})\rangle =\mathrm{i}\,x_{\ell}(\boldsymbol{\theta})\lvert b_{\ell}\rangle`$ and $`g_{i,i}(\boldsymbol{\theta})=\vert x_{\ell}(\boldsymbol{\theta})\vert ^{2}`$.

In each case, the ancestor assignments isolate the relevant branch and the $`+\pi/2`$ shift produces the normalized tangent state.

In what follows, we assume coherent controlled preparation of the relevant Hopf ansatz with a fixed relative phase, so that superpositions of $`\lvert \psi(\boldsymbol{\theta})\rangle`$ and $`\lvert e_{i}(\boldsymbol{\theta})\rangle`$ are well defined.

<a id="S3.SSx1"></a>

## Gradient estimation from signed branch states

With the normalized tangent states available, each gradient component reduces to the transition moment $`\mathrm{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle \right].`$ Rather than using a separate Hadamard test for each parameter, we introduce the signed superposition states

<a id="S3.E24"></a>

```math
\lvert \varphi_{i}^{(s)}(\boldsymbol{\theta})\rangle :=\frac{1}{\sqrt{2}}\left(\lvert \psi(\boldsymbol{\theta})\rangle +s\,\lvert e_{i}(\boldsymbol{\theta})\rangle \right),\;\;s\in\{+1,-1\},\tag{24}
```

together with the energies

<a id="S3.E25"></a>

<a id="S3.E25X"></a>

<a id="S3.E25Xa"></a>

<a id="S3.E25Xb"></a>

```math
\begin{aligned}E_{\psi}:=\langle \psi(\boldsymbol{\theta})\rvert & H\lvert \psi(\boldsymbol{\theta})\rangle , \\
E_{i}^{\mathrm{tan}}:=\langle e_{i}(\boldsymbol{\theta})\rvert & H\lvert e_{i}(\boldsymbol{\theta})\rangle , \\
E_{\varphi_{i}^{(s)}}:=\langle \varphi_{i}^{(s)}(\boldsymbol{\theta})\rvert & H\lvert \varphi_{i}^{(s)}(\boldsymbol{\theta})\rangle .\end{aligned}\tag{25}
```

A direct expansion yields

<a id="S3.E26"></a>

```math
E_{\varphi_{i}^{(s)}}=\frac{1}{2}\big(E_{\psi}+E_{i}^{\mathrm{tan}}\big)+s\,\mathrm{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle \right].\tag{26}
```

Rearranging gives, for either sign $`s=\pm 1`$,

<a id="S3.E27"></a>

```math
\mathrm{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle \right]=s\!\left(E_{\varphi_{i}^{(s)}}-\tfrac{1}{2}\big(E_{\psi}+E_{i}^{\mathrm{tan}}\big)\right).\tag{27}
```

The single-branch identity is useful when $`E_{\psi}`$ and $`E_{i}^{\mathrm{tan}}`$ are estimated separately: either sign-conditioned branch average then recovers the desired matrix element, with the observed sign $`s`$ entering as a known prefactor. No postselection is required.

If both branches are combined, one obtains the equivalent symmetric identity

<a id="S3.E28"></a>

```math
\mathrm{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle \right]=\frac{1}{2}\Big(E_{\varphi_{i}^{(+)}}-E_{\varphi_{i}^{(-)}}\Big).\tag{28}
```

Substituting ([27](gradient-access.md#S3.E27)) into the general gradient formula ([13](gradient-access.md#S3.E13)) gives the single-branch estimator

<a id="S3.E29"></a>

```math
\partial_{\theta_{i}}E(\boldsymbol{\theta})=2\sqrt{g_{i,i}(\boldsymbol{\theta})}\;s\!\left(E_{\varphi_{i}^{(s)}}-\tfrac{1}{2}\big(E_{\psi}+E_{i}^{\mathrm{tan}}\big)\right),\tag{29}
```

valid for either $`s=+1`$ or $`s=-1`$. Equivalently, if both branch energies are estimated, one may write

<a id="S3.E30"></a>

```math
\partial_{\theta_{i}}E(\boldsymbol{\theta})=\sqrt{g_{i,i}(\boldsymbol{\theta})}\;\Big(E_{\varphi_{i}^{(+)}}-E_{\varphi_{i}^{(-)}}\Big).\tag{30}
```

Since $`g_{i,i}(\boldsymbol{\theta})`$ is known analytically from Theorem [1](chart-and-geometry.md#Thmtheorem1), this gives a derivative-plus-branch estimator: estimate one baseline energy $`E_{\psi}`$, and for each $`i`$, estimate $`E_{i}^{\mathrm{tan}}`$ together with one sign-conditioned branch energy. Equivalently, the symmetric identity ([30](gradient-access.md#S3.E30)) trades the separate tangent energy for the two branch signs. Both forms use the same tangent-state and branch-state constructions, and the compiled-setting organization below is unchanged at the level of asymptotic scaling.

<a id="S3.SSx2"></a>

## Magnitude-layer batching

The preceding estimators are written coordinate by coordinate. We now organize their tangent and branch constructions by tree depth, so each compiled setting serves an entire magnitude layer.

For each depth $`d\in\{0,1,\ldots,n-1\}`$, define the internal-node layer

<a id="S3.E31"></a>

<a id="S3.E31X"></a>

<a id="S3.E31Xa"></a>

```math
\begin{aligned}\mathcal{L}_{d}: & =\{\,i\in\{1,\ldots,N-1\}:\lfloor\log_{2}i\rfloor=d\,\} \\
=\{2^{d},2^{d}+1,\ldots,2^{d+1}-1\}.\end{aligned}\tag{31}
```

Introduce a reduced index register $`I`$ of $`d`$ qubits storing the layer offset

<a id="S3.E32"></a>

```math
r:=i-2^{d}\in\{0,\ldots,2^{d}-1\},\;\;i=2^{d}+r\in\mathcal{L}_{d}.\tag{32}
```

Thus $`\lvert r\rangle _{I}`$ labels internal node $`i=2^{d}+r`$.

Prepare the uniform superposition over $`\mathcal{L}_{d}`$,

<a id="S3.E33"></a>

```math
\lvert \mathrm{unif}(\mathcal{L}_{d})\rangle _{I}:=\frac{1}{\sqrt{\vert \mathcal{L}_{d}\vert }}\sum_{r=0}^{2^{d}-1}\lvert r\rangle _{I},\;\;\vert \mathcal{L}_{d}\vert =2^{d}.\tag{33}
```

For $`d\geq 1`$, this is $`\lvert +\rangle ^{\otimes d}`$; for $`d=0`$, the index register is empty and the layer consists only of the root node.

Using the explicit synthesis maps $`i\mapsto\boldsymbol{\theta}^{(i)}`$, we implement a coherent index-controlled preparation on $`I\otimes\mathcal{H}`$ that maps

<a id="S3.Ex22"></a>

```math
\lvert r\rangle _{I}\lvert 0\rangle ^{\otimes n}\longmapsto\lvert r\rangle _{I}\lvert e_{2^{d}+r}(\boldsymbol{\theta})\rangle \;\;\forall\,r\in\{0,\ldots,2^{d}-1\}.
```

Concretely, the same Hopf gate skeleton is reused throughout, and the dependence on $`r`$ enters only through coherent control logic selecting the shifted or clamped values prescribed by Theorem [3](gradient-access.md#Thmtheorem3). Because the synthesis rule $`i\mapsto\boldsymbol{\theta}^{(i)}`$ uses only a constant-size set of values at each gate location, each location is implemented by a constant number of index-conditioned multi-controlled rotations, without angle multiplexing.

Acting on ([33](gradient-access.md#S3.E33)), this preparation yields the indexed tangent-state preparation

<a id="S3.E34"></a>

```math
\lvert \Phi_{d}(\boldsymbol{\theta})\rangle :=\frac{1}{\sqrt{\vert \mathcal{L}_{d}\vert }}\sum_{r=0}^{2^{d}-1}\lvert r\rangle _{I}\lvert e_{2^{d}+r}(\boldsymbol{\theta})\rangle .\tag{34}
```

Measuring $`I`$ in the computational basis and the system register using standard Pauli measurement schemes yields label-conditioned samples of

<a id="S3.Ex23"></a>

```math
E_{i}^{\mathrm{tan}}=\langle e_{i}(\boldsymbol{\theta})\rvert H\lvert e_{i}(\boldsymbol{\theta})\rangle ,\qquad i=2^{d}+r,
```

for all $`i\in\mathcal{L}_{d}`$. In the ideal-expectation limit, the corresponding conditional averages recover these quantities for every label in the layer.

To estimate the branch energies simultaneously for all $`i\in\mathcal{L}_{d}`$, introduce one ancilla qubit $`a`$ and prepare

<a id="S3.E35"></a>

<a id="S3.E35X"></a>

<a id="S3.E35Xa"></a>

```math
\begin{aligned}\lvert \Omega_{d}(\boldsymbol{\theta})\rangle \\
:= & \frac{1}{\sqrt{\vert \mathcal{L}_{d}\vert }}\sum_{r=0}^{2^{d}-1}\lvert r\rangle _{I}\left(\frac{\lvert 0\rangle _{a}\lvert \psi(\boldsymbol{\theta})\rangle +\lvert 1\rangle _{a}\lvert e_{2^{d}+r}(\boldsymbol{\theta})\rangle }{\sqrt{2}}\right).\end{aligned}\tag{35}
```

Measuring $`a`$ in the $`X`$ basis yields an outcome $`s\in\{+1,-1\}`$, and conditioned on the observed index $`r`$, prepares the corresponding state $`\lvert \varphi_{i}^{(s)}(\boldsymbol{\theta})\rangle`$ with $`i=2^{d}+r`$. Because $`\lvert e_{i}(\boldsymbol{\theta})\rangle`$ is normalized and $`\mathrm{Re}\langle \psi(\boldsymbol{\theta})\vert e_{i}(\boldsymbol{\theta})\rangle =0`$, the state $`\lvert \varphi_{i}^{(s)}(\boldsymbol{\theta})\rangle`$ is normalized for both $`s=\pm 1`$.

Operationally, measuring the ancilla in the $`X`$ basis is implemented by a final Hadamard followed by a $`Z`$-basis measurement, yielding a bit $`b\in\{0,1\}`$ and hence a known sign $`s=(-1)^{b}`$. The observed pair $`(r,s)`$, together with the Hamiltonian measurement outcome on the system register, then contributes to the corresponding single-branch estimator ([29](gradient-access.md#S3.E29)) for $`i=2^{d}+r`$. Thus, each experimental repetition contributes directly to one gradient component, without postselection.

<a id="S3.SSx3"></a>

## Single phase-layer batching for the complex ansatz

The complex Hopf ansatz contains an additional block of leaf-phase parameters

<a id="S3.Ex24"></a>

```math
\theta_{N+\ell},\qquad\ell=0,\ldots,N-1,
```

whose normalized tangent states are supported on single computational-basis leaves (Theorem [3](gradient-access.md#Thmtheorem3)). Unlike the magnitude block, these parameters do not split by tree depth; instead they admit a single leaf-indexed batch.

Introduce an $`n`$-qubit phase-index register $`P`$, with the data register’s MSB-to-LSB ordering, and prepare

<a id="S3.E36"></a>

```math
\lvert \mathrm{unif}_{\mathrm{ph}}\rangle _{P}:=\frac{1}{\sqrt{N}}\sum_{\ell=0}^{N-1}\lvert \ell\rangle _{P}.\tag{36}
```

Using the explicit synthesis map $`\ell\mapsto\boldsymbol{\theta}^{(N+\ell)}`$ from Theorem [3](gradient-access.md#Thmtheorem3)(B), we implement a coherent index-controlled preparation on $`P\otimes\mathcal{H}`$ such that

<a id="S3.Ex25"></a>

```math
\lvert \ell\rangle _{P}\lvert 0\rangle ^{\otimes n}\longmapsto\lvert \ell\rangle _{P}\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle \;\;\forall\,\ell\in\{0,\ldots,N-1\}.
```

This yields the indexed phase-derivative state

<a id="S3.E37"></a>

```math
\lvert \Phi_{\mathrm{ph}}(\boldsymbol{\theta})\rangle :=\frac{1}{\sqrt{N}}\sum_{\ell=0}^{N-1}\lvert \ell\rangle _{P}\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle .\tag{37}
```

Measuring $`P`$ in the computational basis and the system register with standard Pauli measurement schemes yields label-conditioned samples of

<a id="S3.Ex26"></a>

```math
E_{N+\ell}^{\mathrm{tan}}=\langle e_{N+\ell}(\boldsymbol{\theta})\rvert H\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle
```

for all phase parameters. In the ideal-expectation limit, the corresponding conditional averages recover these phase-block tangent energies.

To estimate the phase-block branch energies simultaneously, introduce one ancilla qubit $`a`$ and prepare

<a id="S3.E38"></a>

<a id="S3.E38X"></a>

<a id="S3.E38Xa"></a>

```math
\begin{aligned}\lvert \Omega_{\mathrm{ph}}(\boldsymbol{\theta})\rangle \\
:= & \frac{1}{\sqrt{N}}\sum_{\ell=0}^{N-1}\lvert \ell\rangle _{P}\left(\frac{\lvert 0\rangle _{a}\lvert \psi(\boldsymbol{\theta})\rangle +\lvert 1\rangle _{a}\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle }{\sqrt{2}}\right).\end{aligned}\tag{38}
```

Measuring $`a`$ in the $`X`$ basis and $`P`$ in the computational basis yields, for the observed label $`\ell`$ and sign $`s\in\{+1,-1\}`$, label- and sign-conditioned samples of

<a id="S3.Ex27"></a>

```math
E_{\varphi_{N+\ell}^{(s)}}=\langle \varphi_{N+\ell}^{(s)}(\boldsymbol{\theta})\rvert H\lvert \varphi_{N+\ell}^{(s)}(\boldsymbol{\theta})\rangle .
```

Exactly the same single-branch estimator ([29](gradient-access.md#S3.E29)) applies to the phase block, with the observed pair $`(\ell,s)`$, together with the Hamiltonian measurement outcome on the system register, selecting the corresponding gradient component. No additional postselection step is required.

<a id="Thmtheorem4"></a>

**Theorem 4 (Uniform per-setting cost of Hopf gradient configurations).**

Under the CNOT counting model of Lemma [2](algorithms.md#Thmlemma2), every compiled circuit setting used in the Hopf gradient-access protocol has compiled execution cost $`O(nN)`$. More precisely:

(A) Real Hopf ansatz. For every depth $`d\in\{0,1,\ldots,n-1\},`$ the indexed derivative setting $`\lvert \Phi_{d}(\boldsymbol{\theta})\rangle`$ defined in ([34](gradient-access.md#S3.E34)) and the indexed branch setting $`\lvert \Omega_{d}(\boldsymbol{\theta})\rangle`$ defined in ([35](gradient-access.md#S3.E35)) each have compiled execution cost $`O(nN)`$.

(B) Complex Hopf ansatz. For every magnitude depth $`d\in\{0,1,\ldots,n-1\},`$ the magnitude-block settings $`\lvert \Phi_{d}(\boldsymbol{\theta})\rangle`$ and $`\lvert \Omega_{d}(\boldsymbol{\theta})\rangle`$ again each have compiled execution cost $`O(nN)`$.

In addition, the single phase-block derivative setting $`\lvert \Phi_{\mathrm{ph}}(\boldsymbol{\theta})\rangle`$ defined in ([37](gradient-access.md#S3.E37)) and the single phase-block branch setting $`\lvert \Omega_{\mathrm{ph}}(\boldsymbol{\theta})\rangle`$ defined in ([38](gradient-access.md#S3.E38)) each have compiled execution cost $`O(nN)`$.

Consequently, for both versions of the Hopf ansatz, every gradient-access configuration, whether derivative or branch, magnitude or phase, has the same asymptotic compiled execution cost as a forward Hopf energy-evaluation setting:

<a id="S3.Ex28"></a>

```math
C_{\mathrm{cfg}}(\text{one gradient setting})=O(nN).
```

<a id="S3.SSx4"></a>

## Complete estimator and scaling analysis

Since $`N=2^{n}`$, the magnitude tree has exactly $`n=\log_{2}N`$ layers. The real ansatz therefore requires one baseline setting, one tangent setting per magnitude layer, and one branch setting per magnitude layer:

<a id="S3.E39"></a>

```math
\boxed{N_{\mathrm{set}}^{\mathbb{R}}=1+2n=O(\log N).}\tag{39}
```

Here a setting means a fixed compiled gate pattern and indexed preparation choice; the angles $`\theta_{j}`$ are streamed as physical rotation parameters within that setting.

For the complex ansatz, the magnitude block contributes the same $`2n`$ settings and the phase block adds one tangent and one branch setting:

<a id="S3.Ex29"></a>

```math
\boxed{N_{\mathrm{set}}^{\mathbb{C}}=1+2(n+1)=O(\log N).}
```

Since the number of magnitude parameters is $`N-1`$, these setting counts are also logarithmic in the magnitude-parameter count.

In the hardware-faithful model of Lemma [2](algorithms.md#Thmlemma2), let $`C_{\mathrm{cfg}}`$ denote the compiled execution cost of one circuit setting, using CNOT count and excluding repeated measurements. A baseline energy-evaluation setting costs, by Theorem [5](algorithms.md#Thmtheorem5),

<a id="S3.E40"></a>

```math
\boxed{C_{\mathrm{cfg}}(E)=O(nN)=O(N\log N).}\tag{40}
```

The indexed magnitude and phase preparations reuse the forward Hopf structure through subtree-local tangent synthesis. By Theorem [4](gradient-access.md#Thmtheorem4), each gradient-access setting has the same asymptotic compiled cost as a forward setting. Let $`C_{\mathrm{grad}}(\nabla E)`$ denote the total compiled cost of the gradient-access settings needed for one full gradient, excluding repeated measurements. Then

<a id="S3.E41"></a>

```math
\boxed{C_{\mathrm{grad}}(\nabla E)=O(n)\,C_{\mathrm{cfg}}(E)=O(\log N)\,C_{\mathrm{cfg}}(E).}\tag{41}
```

The indexed preparations organize gradient access into $`O(\log N)`$ compiled circuit families, with the measured index specifying which gradient component each repetition informs. Accordingly, the total number of measurement repetitions required for one full gradient at a prescribed accuracy depends on the observable decomposition, label–sign sampling probabilities, estimator variances, and chosen error criterion, and is not determined by the compiled-setting count alone.

<a id="S3.SSx4.p6.1"></a>

Related approaches reduce gradient overhead through algebraic circuit structure, density-based models, or expressivity–measurement tradeoffs [[28](../references/paper-bibliography.md#bib.bib7), [29](../references/paper-bibliography.md#bib.bib22), [30](../references/paper-bibliography.md#bib.bib24), [31](../references/paper-bibliography.md#bib.bib23)]. The Hopf construction organizes the required gradient settings within the same circuit skeleton through exact tangent-state synthesis and analytic binary-tree geometry, with the attainable precision of each component set by the associated measurement statistics.

Although stated for the VQE objective $`E(\boldsymbol{\theta})=\langle H\rangle_{\boldsymbol{\theta}}`$, the construction is not tied to Hamiltonian minimization. It only requires that the objective derivative reduce to transition elements of Hermitian observables,

<a id="S3.E42"></a>

```math
\operatorname{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert O\lvert \psi(\boldsymbol{\theta})\rangle \right].\tag{42}
```

<a id="S3.SSx4.p7.2"></a>

The same framework applies, after local linearization, to metrological costs built from finitely many expectation values, probabilities, or probability phase-derivatives, including local QFI and fixed-readout CFI objectives [[32](../references/paper-bibliography.md#bib.bib33), [33](../references/paper-bibliography.md#bib.bib34), [34](../references/paper-bibliography.md#bib.bib35), [35](../references/paper-bibliography.md#bib.bib36), [36](../references/paper-bibliography.md#bib.bib37)]. Appendix [D](metrological-objectives.md#A4) gives two local examples: pure-state QFI as an ultimate probe-state benchmark and fixed-readout Ramsey CFI objective [[37](../references/paper-bibliography.md#bib.bib38)] for a specified parity readout at a calibrated operating point. In both, Hopf tangent preparation and compiled-setting organization are unchanged after replacing $`H`$ by the metrological chain-rule observable. The optimal measurement, estimator, noise model, and adaptive strategy remain problem dependent.
