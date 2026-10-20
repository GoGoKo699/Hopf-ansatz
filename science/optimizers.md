[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="A3"></a>

# Appendix C Geometry-native optimizers for the Hopf ansatz

This appendix records the optimizer layer used with the Hopf ansatz, focusing on methods whose quantum interface consists of objective evaluations and Hopf-coordinate gradient evaluations. All geometry-dependent operations below are classical post-processing: diagonal metric rescaling, Jacobian lifting, exact sphere exponential maps, exact vector transport, line searches, and the inverse Hopf map. Adam-type coordinate updates are used in the main text as baselines; the appendix focuses on optimizers that explicitly use the Hopf metric and state-sphere geometry.

The starting point is the same cost-and-gradient loop used in standard variational algorithms. At a parameter value $`\boldsymbol{\theta}`$, the quantum device supplies an estimate of

<a id="A3.E61"></a>

```math
{C(\boldsymbol{\theta})=C\big(\lvert \psi(\boldsymbol{\theta})\rangle \big)}\tag{C1}
```

and of the coordinate-gradient components

<a id="A3.E62"></a>

```math
{c_{i}:=\partial_{\theta_{i}}C(\boldsymbol{\theta}).}\tag{C2}
```

The Hopf structure then converts these coordinate data into state-sphere data using analytic information from the chart: the inverse map, diagonal metric, and raw Jacobian directions are all known explicitly.

<a id="A3.SSx1"></a>

## Jacobian formulas and inverse differential map

The main text uses only the diagonal metric, but the optimizer layer also uses the differential of the Hopf chart. The Jacobian maps parameter variations to state-space tangent vectors; its columns are the raw coordinate tangents in the ambient state-vector basis.

<a id="Thmtheorem6"></a>

**Theorem 6 (Jacobian formulas and Gram forms).**

Let $`N=2^{n}`$, and let

<a id="A3.Ex60"></a>

```math
\mathbf{x}(\boldsymbol{\theta})=(x_{0},\ldots,x_{N-1})
```

be the Hopf coordinate map.

In the real case define

<a id="A3.E63"></a>

```math
{J^{\mathbb{R}}(\boldsymbol{\theta})\in\mathbb{R}^{N\times(N-1)},\qquad(J^{\mathbb{R}})_{aj}:=\frac{\partial x_{a}}{\partial\theta_{j}},}\tag{C3}
```

with $`a=0,\ldots,N-1`$ and $`j=1,\ldots,N-1`$. In the complex case define

<a id="A3.E64"></a>

```math
{J^{\mathbb{C}}(\boldsymbol{\theta})\in\mathbb{C}^{N\times(2N-1)},\qquad(J^{\mathbb{C}})_{aj}:=\frac{\partial x_{a}}{\partial\theta_{j}},}\tag{C4}
```

with $`a=0,\ldots,N-1`$ and $`j=1,\ldots,2N-1`$.

For a fixed amplitude index $`a`$, let $`P(a)`$ be the set of internal nodes on the unique root-to-leaf path that generates $`x_{a}`$. For $`j\in P(a)`$, define $`\varepsilon_{aj}\in\{0,1\}`$ by

<a id="A3.Ex61"></a>

```math
\varepsilon_{aj}=\begin{cases}0,&\text{if the path contributes a }\cos\theta_{j}\text{ factor},\\
1,&\text{if the path contributes a }\sin\theta_{j}\text{ factor}.\end{cases}
```

The following logarithmic-derivative formulas hold literally on the regular set where the displayed $`\tan`$ and $`\cot`$ factors are finite. At zeros of the corresponding sine or cosine factor, they are understood through the direct product derivative, equivalently by continuous extension.

(A) Real case. For $`1\leq j\leq N-1`$,

<a id="A3.E65"></a>

```math
{\frac{\partial x_{a}}{\partial\theta_{j}}=\begin{cases}-\tan\theta_{j}\,x_{a},&j\in P(a),\ \varepsilon_{aj}=0,\\[2.0pt]
\ \ \cot\theta_{j}\,x_{a},&j\in P(a),\ \varepsilon_{aj}=1,\\[2.0pt]
0,&j\notin P(a).\end{cases}}\tag{C5}
```

Moreover, the real pullback metric is the Gram matrix

<a id="A3.E66"></a>

```math
{g^{\mathbb{R}}(\boldsymbol{\theta})=\big(J^{\mathbb{R}}(\boldsymbol{\theta})\big)^{\top}J^{\mathbb{R}}(\boldsymbol{\theta}).}\tag{C6}
```

(B) Complex case. Write

<a id="A3.Ex62"></a>

```math
x_{a}(\boldsymbol{\theta})=r_{a}(\boldsymbol{\theta})e^{\mathrm{i}\theta_{N+a}},\qquad r_{a}\geq 0,
```

where $`r_{a}`$ is produced by the internal-node magnitude angles. For $`1\leq j\leq N-1`$, the same path-derivative rule holds, and for the phase block there is one additional nonzero derivative:

<a id="A3.E67"></a>

```math
{\frac{\partial x_{a}}{\partial\theta_{j}}=\begin{cases}-\tan\theta_{j}\,x_{a},&j\in P(a),\ \varepsilon_{aj}=0,\\[2.0pt]
\ \ \cot\theta_{j}\,x_{a},&j\in P(a),\ \varepsilon_{aj}=1,\\[2.0pt]
\mathrm{i}\,x_{a},&j=N+a,\\[2.0pt]
0,&\text{otherwise}.\end{cases}}\tag{C7}
```

Moreover, under the ambient round-sphere convention of Definition [4](chart-and-geometry.md#Thmdefinition4), the complex pullback metric is

<a id="A3.E68"></a>

```math
{g^{\mathbb{C}}(\boldsymbol{\theta})=\mathrm{Re}\!\left[\big(J^{\mathbb{C}}(\boldsymbol{\theta})\big)^{\dagger}J^{\mathbb{C}}(\boldsymbol{\theta})\right].}\tag{C8}
```

Thus the Hopf Jacobian is explicit, and its Gram matrix is exactly the pullback metric. Since the metric is diagonal, the differential map can be inverted componentwise on the regular set, rather than by solving a dense linear system.

Define the regular parameter regions

<a id="A3.Ex63"></a>

```math
\mathcal{U}_{\mathbb{R}}:=\left\{\boldsymbol{\theta}\mid g^{\mathbb{R}}_{i,i}(\boldsymbol{\theta})>0\ \forall i\right\},\qquad\mathcal{U}_{\mathbb{C}}:=\left\{\boldsymbol{\theta}\mid g^{\mathbb{C}}_{i,i}(\boldsymbol{\theta})>0\ \forall i\right\}.
```

On $`\mathcal{U}_{\mathbb{R}}`$, the Jacobian $`J^{\mathbb{R}}`$ has full column rank and its Moore–Penrose left inverse is

<a id="A3.Ex64"></a>

```math
(J^{\mathbb{R}})^{+}=(g^{\mathbb{R}})^{-1}(J^{\mathbb{R}})^{\top},\qquad(J^{\mathbb{R}})^{+}J^{\mathbb{R}}=I_{N-1}.
```

For the complex ansatz, after realification

<a id="A3.Ex65"></a>

```math
J_{\mathbb{R}}:=\begin{bmatrix}\mathrm{Re}J^{\mathbb{C}}\\
\mathrm{Im}J^{\mathbb{C}}\end{bmatrix},
```

the corresponding left inverse on $`\mathcal{U}_{\mathbb{C}}`$ is

<a id="A3.Ex66"></a>

```math
(J_{\mathbb{R}})^{+}=(g^{\mathbb{C}})^{-1}(J_{\mathbb{R}})^{\top},\qquad(J_{\mathbb{R}})^{+}J_{\mathbb{R}}=I_{2N-1}.
```

Consequently, on the regular set, the diagonal metric not only gives componentwise Riemannian rescaling, but also gives an explicit inverse differential map from state-space tangent vectors back to Hopf-coordinate variations.

<a id="A3.SSx2"></a>

## State-sphere conventions

Let $`\mathbb{S}_{\star}`$ denote the normalized state-vector sphere on which the optimizer moves. In the real ansatz, $`\mathbb{S}_{\star}=\mathbb{S}^{N-1}\subset\mathbb{R}^{N}`$. In the complex ansatz, $`\mathbb{C}^{N}`$ is regarded as a real Euclidean space and $`\mathbb{S}_{\star}=\mathbb{S}^{2N-1}`$. The real pairing is

<a id="A3.E69"></a>

```math
{\langle u,v\rangle_{\star}=\begin{cases}\langle u\vert v\rangle ,&\text{real case},\\[2.0pt]
\mathrm{Re}\langle u\vert v\rangle ,&\text{complex round-sphere case}.\end{cases}}\tag{C9}
```

The tangent projection at $`\lvert \psi\rangle \in\mathbb{S}_{\star}`$ is

<a id="A3.E70"></a>

```math
{\Pi_{\psi}(w)=w-\langle\psi,w\rangle_{\star}\psi.}\tag{C10}
```

For $`u\in T_{\psi}\mathbb{S}_{\star}`$, the exact sphere exponential map is

<a id="A3.E71"></a>

```math
{\mathrm{Exp}_{\psi}(u)=\begin{cases}\psi,&\|u\|_{\star}=0,\\[4.0pt]
\cos(\|u\|_\star)\psi+\dfrac{\sin(\|u\|_\star)}{\|u\|_{\star}}u,&\|u\|_{\star}>0.\end{cases}}\tag{C11}
```

<a id="A3.SSx2.p1.4"></a>

Equivalently, if the geodesic is parameterized as $`\gamma(\eta)=\mathrm{Exp}_{\psi}(\eta u)`$, then its angular length is $`\eta\|u\|_{\star}`$. This is the exact-geodesic update used in Ref. [[24](../references/paper-bibliography.md#bib.bib2)].

The exact parallel transport of a tangent vector $`\zeta\in T_{\psi}\mathbb{S}_{\star}`$ along the same geodesic $`\gamma(\eta)=\mathrm{Exp}_{\psi}(\eta u)`$ is

<a id="A3.E72"></a>

```math
{\mathcal{T}_{\eta u}(\zeta)=\zeta-\sin(\eta\|u\|_\star)\frac{\langle u,\zeta\rangle_{\star}}{\|u\|_{\star}}\,\psi+\big[\cos(\eta\|u\|_\star)-1\big]\frac{\langle u,\zeta\rangle_{\star}}{\|u\|_{\star}^{2}}\,u,}\tag{C12}
```

with the continuous identity limit when $`\|u\|_{\star}=0`$. This is the transport used by Hopf-EGT-CG and by the Riemannian memory methods below. It is a vector-transport rule along the accepted geodesic; it is not a second state-update geodesic.

<a id="A3.SSx3"></a>

## Lifting Hopf coordinate gradients to the state sphere

For brevity in this appendix, write

<a id="A3.E73"></a>

```math
{\lvert \dot{\psi}_{i}\rangle :=\partial_{\theta_{i}}\lvert \psi(\boldsymbol{\theta})\rangle ,\qquad\lvert e_{i}\rangle :=\frac{\lvert \dot{\psi}_{i}\rangle }{\sqrt{g_{i,i}(\boldsymbol{\theta})}}}\tag{C13}
```

where $`g_{i,i}`$ is the diagonal Hopf metric from Theorem [1](chart-and-geometry.md#Thmtheorem1). On the regular set,

<a id="A3.E74"></a>

```math
{\langle e_{i},e_{j}\rangle_{\star}=\delta_{ij}.}\tag{C14}
```

Thus the Hopf differential has the orthonormal-frame form

<a id="A3.E75"></a>

```math
{D_{\boldsymbol{\theta}}\psi[\delta\boldsymbol{\theta}]=\sum_{i}\delta\theta_{i}\sqrt{g_{i,i}(\boldsymbol{\theta})}\,e_{i}.}\tag{C15}
```

If the coordinate gradient components $`c_{i}=\partial_{\theta_{i}}C`$ are known, the corresponding state-sphere gradient is

<a id="A3.E76"></a>

```math
{G(\boldsymbol{\theta}):=\mathrm{grad}_{\mathbb{S}_{\star}}C=\sum_{i}\frac{c_{i}}{\sqrt{g_{i,i}(\boldsymbol{\theta})}}e_{i}=\sum_{i}\frac{c_{i}}{g_{i,i}(\boldsymbol{\theta})}\dot{\psi}_{i}.}\tag{C16}
```

In matrix notation this is $`G=Jg^{-1}\nabla_{\theta}C`$, where $`J`$ is the Hopf Jacobian and $`g=J^{\top}J`$ in the real case, or the corresponding realified Gram matrix in the complex round-sphere convention.

Equation ([C16](optimizers.md#A3.E76)) is the optimizer-side use of the Hopf metric. Since $`g`$ is diagonal, no dense metric inversion is performed. In the real tree, the state-gradient lift can be evaluated without materializing $`J`$: for each leaf, only the $`n`$ internal nodes on its root-to-leaf path contribute. Thus the lift costs $`O(nN)`$ classical time and $`O(N)`$ memory for an $`n`$-qubit real Hopf state, after the coordinate gradient has been supplied.

After a state-sphere update produces $`\psi_{k+1}`$, the next circuit parameters are obtained by the inverse Hopf map,

<a id="A3.E77"></a>

```math
{\boldsymbol{\theta}_{k+1}=\mathrm{Hopf}^{-1}(\psi_{k+1}),}\tag{C17}
```

using Lemma [1](algorithms.md#Thmlemma1). This is the step that lets a classical state-sphere optimizer remain an executable Hopf-circuit optimizer.

<a id="A3.SSx4"></a>

## Hopf exact geodesic transport with conjugate gradients

The basic exact-geodesic-transport descent direction is

<a id="A3.E78"></a>

```math
{v_{k}:=-G_{k},}\tag{C18}
```

where $`G_{k}=G(\boldsymbol{\theta}_{k})`$. The first-order Hopf-EGT update is

<a id="A3.E79"></a>

```math
{\psi_{k+1}=\mathrm{Exp}_{\psi_{k}}(\eta_{k}v_{k}),\qquad\boldsymbol{\theta}_{k+1}=\mathrm{Hopf}^{-1}(\psi_{k+1}).}\tag{C19}
```

Infinitesimally, this reduces to the diagonal natural-gradient parameter step

<a id="A3.E80"></a>

```math
{\delta\theta_{i}=-\eta_{k}\frac{\partial_{\theta_{i}}C}{g_{i,i}}+O(\eta_{k}^{2}),}\tag{C20}
```

but Eq. ([C19](optimizers.md#A3.E79)) keeps the state normalized exactly and uses the full sphere geodesic rather than a flat coordinate step.

Hopf-EGT-CG adds a transported conjugate-gradient memory direction. Let $`u_{k}\in T_{\psi_{k}}\mathbb{S}_{\star}`$ denote the search direction. Initialize

<a id="A3.E81"></a>

```math
{u_{0}=v_{0}=-G_{0}.}\tag{C21}
```

Given $`u_{k}`$, choose $`\eta_{k}>0`$ and update

<a id="A3.E82"></a>

```math
{\psi_{k+1}=\mathrm{Exp}_{\psi_{k}}(\eta_{k}u_{k}).}\tag{C22}
```

Then transport the previous direction and previous gradient to $`T_{\psi_{k+1}}\mathbb{S}_{\star}`$:

<a id="A3.E83"></a>

```math
{\widetilde{u}_{k}:=\mathcal{T}_{\eta_{k}u_{k}}(u_{k}),\qquad\widetilde{G}_{k}:=\mathcal{T}_{\eta_{k}u_{k}}(G_{k}).}\tag{C23}
```

<a id="A3.SSx4.p2.4"></a>

Using the hybrid Dai–Yuan/Hestenes–Stiefel coefficient in gradient notation [[55](../references/paper-bibliography.md#bib.bib48), [56](../references/paper-bibliography.md#bib.bib52)], define

<a id="A3.E84"></a>

```math
{D_{k}:=\langle G_{k+1},\widetilde{u}_{k}\rangle_{\star}-\langle G_{k},u_{k}\rangle_{\star},}\tag{C24}
```

<a id="A3.E85"></a>

```math
{\beta_{k+1}^{\rm DY}:=\frac{\langle G_{k+1},G_{k+1}\rangle_{\star}}{D_{k}},\qquad\beta_{k+1}^{\rm HS}:=\frac{\langle G_{k+1},G_{k+1}-\widetilde{G}_{k}\rangle_{\star}}{D_{k}},}\tag{C25}
```

and

<a id="A3.E86"></a>

```math
{\beta_{k+1}:=\max\{0,\min(\beta_{k+1}^{\rm DY},\beta_{k+1}^{\rm HS})\}.}\tag{C26}
```

If $`D_{k}`$ is numerically singular or the resulting direction is not a descent direction, we restart with $`u_{k+1}=-G_{k+1}`$. Otherwise,

<a id="A3.E87"></a>

```math
{u_{k+1}=-G_{k+1}+\beta_{k+1}\widetilde{u}_{k}.}\tag{C27}
```

The line search is performed on the one-dimensional geodesic

<a id="A3.E88"></a>

```math
{\gamma_{k}(\eta)=\mathrm{Exp}_{\psi_{k}}(\eta u_{k}).}\tag{C28}
```

The strong-Wolfe conditions take the sphere form

<a id="A3.EGx3"></a>

<a id="A3.E89"></a>

```math
{\begin{aligned}C(\gamma_{k}(\eta_{k}))-C(\psi_{k}) & \leq c_{1}\eta_{k}\langle G_{k},u_{k}\rangle_{\star},\end{aligned}}\tag{C29}
```

<a id="A3.E90"></a>

```math
{\begin{aligned}\left\vert \left\langle G_{k+1},\mathcal{T}_{\eta_{k}u_{k}}(u_{k})\right\rangle_{\star}\right\vert & \leq c_{2}\left\vert \langle G_{k},u_{k}\rangle_{\star}\right\vert ,\end{aligned}}\tag{C30}
```

<a id="A3.SSx4.p3.3"></a>

with $`0<c_{1}<c_{2}<1`$. These are the cost-and-gradient tests used by the EGT-CG optimizer of Ref. [[24](../references/paper-bibliography.md#bib.bib2)]. They require additional objective and, when the curvature condition is checked, gradient evaluations; they do not require any new quantum primitive beyond the cost and Hopf-gradient queries.

<a id="A3.SSx5"></a>

## Hopf Riemannian L-BFGS

<a id="A3.SSx5.p1.1"></a>

The Riemannian L-BFGS variant [[39](../references/paper-bibliography.md#bib.bib49), [40](../references/paper-bibliography.md#bib.bib50)] keeps a limited memory of curvature pairs in state-sphere tangent spaces. At iteration $`k`$, assume the accepted update was

<a id="A3.E91"></a>

```math
{\psi_{k+1}=\mathrm{Exp}_{\psi_{k}}(\eta_{k}p_{k}),}\tag{C31}
```

where $`p_{k}\in T_{\psi_{k}}\mathbb{S}_{\star}`$ is the L-BFGS search direction. The step and gradient-difference vectors at the new point are

<a id="A3.E92"></a>

```math
{s_{k}:=\mathcal{T}_{\eta_{k}p_{k}}(\eta_{k}p_{k}),\qquad y_{k}:=G_{k+1}-\mathcal{T}_{\eta_{k}p_{k}}(G_{k}).}\tag{C32}
```

The pair is accepted only when

<a id="A3.E93"></a>

```math
{\langle s_{k},y_{k}\rangle_{\star}>0.}\tag{C33}
```

All older stored pairs are exactly transported along each newly accepted geodesic before the next two-loop recursion is applied. Thus every vector used by the recursion lies in the current tangent space.

Let the currently stored transported pairs be $`(s_{j},y_{j})`$, $`j=k-r,\ldots,k-1`$, with $`r\leq m`$. Set

<a id="A3.E94"></a>

```math
{\rho_{j}:=\frac{1}{\langle s_{j},y_{j}\rangle_{\star}}.}\tag{C34}
```

The standard two-loop recursion is then applied with the pairing $`\langle\cdot,\cdot\rangle_{\star}`$. Starting from $`q=G_{k}`$, the backward loop uses

<a id="A3.E95"></a>

```math
{\alpha_{j}=\rho_{j}\langle s_{j},q\rangle_{\star},\qquad q\leftarrow q-\alpha_{j}y_{j}.}\tag{C35}
```

The initial inverse-Hessian scale may be chosen as

<a id="A3.E96"></a>

```math
{\gamma_{k}=\frac{\langle s_{k-1},y_{k-1}\rangle_{\star}}{\langle y_{k-1},y_{k-1}\rangle_{\star}},}\tag{C36}
```

when a previous pair exists, and $`\gamma_{k}=1`$ otherwise. Setting $`r=\gamma_{k}q`$, the forward loop uses

<a id="A3.E97"></a>

```math
{\beta_{j}=\rho_{j}\langle y_{j},r\rangle_{\star},\qquad r\leftarrow r+s_{j}(\alpha_{j}-\beta_{j}).}\tag{C37}
```

The search direction is

<a id="A3.E98"></a>

```math
{p_{k}=-\Pi_{\psi_{k}}(r),}\tag{C38}
```

with restart to $`-G_{k}`$ if $`p_{k}`$ is not a descent direction. A geodesic line search, for instance the strong-Wolfe search ([C29](optimizers.md#A3.E89))–([C30](optimizers.md#A3.E90)), then chooses $`\eta_{k}`$.

This optimizer uses the transparent Hopf geometry in two places. First, $`G_{k}`$ is obtained from Hopf-coordinate gradients by the analytic diagonal metric lift ([C16](optimizers.md#A3.E76)). Second, the memory vectors are kept in the correct tangent spaces by exact sphere transport, rather than by treating the parameter space as flat.

<a id="A3.SSx6"></a>

## Hopf Riemannian Barzilai–Borwein

<a id="A3.SSx6.p1.1"></a>

The Riemannian Barzilai–Borwein method [[38](../references/paper-bibliography.md#bib.bib51)] is a spectral-gradient method on the same state sphere. It has no multi-vector memory and therefore provides a simpler third cost-and-gradient optimizer.

After a step from $`\psi_{k-1}`$ to $`\psi_{k}`$, define

<a id="A3.E99"></a>

```math
{s_{k-1}:=\mathcal{T}(\eta_{k-1}p_{k-1}),\qquad y_{k-1}:=G_{k}-\mathcal{T}(G_{k-1}),}\tag{C39}
```

where $`\mathcal{T}`$ denotes exact transport along the accepted geodesic. The two standard spectral choices are

<a id="A3.E100"></a>

```math
{\alpha_{k}^{\rm BB1}=\frac{\langle s_{k-1},s_{k-1}\rangle_{\star}}{\langle s_{k-1},y_{k-1}\rangle_{\star}},\qquad\alpha_{k}^{\rm BB2}=\frac{\langle s_{k-1},y_{k-1}\rangle_{\star}}{\langle y_{k-1},y_{k-1}\rangle_{\star}}.}\tag{C40}
```

The implementation may use BB1, BB2, or alternate between them, with clipping $`\alpha_{k}\in[\alpha_{\min},\alpha_{\max}]`$. The update is

<a id="A3.E101"></a>

```math
{\psi_{k+1}=\mathrm{Exp}_{\psi_{k}}(-\alpha_{k}G_{k}),\qquad\boldsymbol{\theta}_{k+1}=\mathrm{Hopf}^{-1}(\psi_{k+1}).}\tag{C41}
```

For stability one may add a geodesic line-search safeguard, such as Armijo backtracking or the same strong-Wolfe/bounded geodesic search used above, seeded by the Barzilai–Borwein spectral step. This safeguard uses only objective and, when the curvature test is checked, gradient evaluations, and does not change the quantum-information interface.

<a id="A3.SSx7"></a>

## Summary of the optimizer family

The optimizer family used here is therefore:

- **Hopf-EGT-CG**: exact sphere geodesic update, exact vector transport, hybrid Riemannian conjugate-gradient memory, and strong-Wolfe/bounded geodesic line search;

- **Hopf Riemannian L-BFGS**: analytic Hopf metric lift, transported limited-memory curvature pairs, tangent-space two-loop recursion, and geodesic line search;

- **Hopf Riemannian Barzilai–Borwein**: analytic Hopf metric lift, transported spectral step-size estimation, and a geodesic line-search safeguard seeded by the spectral step.

All three optimizers use the same quantum process: evaluate the cost and the Hopf-coordinate gradient at selected Hopf states. The remaining operations are classical consequences of the Hopf chart. Ritz-type and local-subspace updates are intentionally excluded from this appendix because they require stronger oracle information or inner subspace optimization beyond the cost-and-gradient interface considered here.
