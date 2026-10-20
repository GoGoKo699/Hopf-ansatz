[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="A5"></a>

# Appendix E Proofs of the main results

<a id="A5.SSx1"></a>

## Closed-form product formula for Hopf amplitudes

<a id="Thmlemma3"></a>

**Lemma 3 (Closed-form product formula for Hopf amplitudes).**

Fix $`n\in\mathbb{N}`$, set $`N=2^{n}`$, and let $`i\in\{0,\ldots,N-1\}`$ have binary expansion

<a id="A5.Ex96"></a>

```math
i=(q_{n}q_{n-1}\cdots q_{1})_{2},\qquad q_{k}\in\{0,1\},
```

using the qubit order $`(q_{n},\ldots,q_{1})`$ (MSB to LSB). Define the sequence of internal-node indices along the root-to-leaf construction by

<a id="A5.Ex97"></a>

```math
s_{t}(i):=2^{t-1}+\big(q_{n}q_{n-1}\cdots q_{n-t+2}\big)_{2},\qquad t=1,\ldots,n,
```

where for $`t=1`$ the empty bit string is interpreted as $`0`$, so that $`s_{1}(i)=1`$. Equivalently, $`s_{1}(i)=1`$ and

<a id="A5.Ex98"></a>

```math
s_{t+1}(i)=2\,s_{t}(i)+q_{n-t+1},\qquad t=1,\ldots,n-1.
```

(A) Real Hopf coordinates. The amplitude produced by Definition [2](chart-and-geometry.md#Thmdefinition2)(A) satisfies

<a id="A5.Ex99"></a>

```math
x_{i}(\boldsymbol{\theta})=\prod_{t=1}^{n}\begin{cases}\cos\theta_{s_{t}(i)},&q_{n-t+1}=0,\\[2.0pt]
\sin\theta_{s_{t}(i)},&q_{n-t+1}=1.\end{cases}
```

In particular,

<a id="A5.Ex100"></a>

```math
x_{i}(\boldsymbol{\theta})^{2}=\prod_{t=1}^{n}\begin{cases}\cos^{2}\theta_{s_{t}(i)},&q_{n-t+1}=0,\\[2.0pt]
\sin^{2}\theta_{s_{t}(i)},&q_{n-t+1}=1.\end{cases}
```

(B) Complex Hopf coordinates. Let $`r_{i}(\boldsymbol{\theta})\geq 0`$ denote the magnitude produced from the internal-node angles as in (A), using $`\theta_{1},\ldots,\theta_{N-1}`$. Then Definition [2](chart-and-geometry.md#Thmdefinition2)(B) yields

<a id="A5.Ex101"></a>

```math
x_{i}(\boldsymbol{\theta})=r_{i}(\boldsymbol{\theta})\,e^{\mathrm{i}\theta_{N+i}},\qquad\vert x_{i}(\boldsymbol{\theta})\vert ^{2}=r_{i}(\boldsymbol{\theta})^{2}.
```

<a id="A5.SSx1.1"></a>

**Proof.**

Fix $`i`$ and run the update rule in Definition [2](chart-and-geometry.md#Thmdefinition2). At the start, $`x_{i}=1`$ and $`j=1`$. At step $`t`$ (corresponding to $`k=n-t+1`$), the algorithm multiplies $`x_{i}`$ by $`\cos\theta_{j}`$ if $`q_{n-t+1}=0`$ and by $`\sin\theta_{j}`$ if $`q_{n-t+1}=1`$, and then updates $`j\leftarrow 2j`$ or $`2j+1`$ accordingly. By construction, the value of $`j`$ used at step $`t`$ is exactly $`s_{t}(i)`$. Therefore after $`n`$ steps, $`x_{i}`$ equals the stated product of $`n`$ cosine/sine factors, proving (A), and squaring gives the formula for $`x_{i}^{2}`$.

In the complex case, the same internal-node recursion produces the nonnegative magnitude $`r_{i}(\boldsymbol{\theta})`$, after which Definition [2](chart-and-geometry.md#Thmdefinition2)(B) multiplies by the leaf phase $`e^{\mathrm{i}\theta_{N+i}}`$, yielding (B). ∎

<a id="A5.SSx2"></a>

## Proof of Theorem [1](chart-and-geometry.md#Thmtheorem1)

<a id="A5.SSx2.1"></a>

**Proof.**

We write $`\lvert \psi(\boldsymbol{\theta})\rangle =\sum_{k=0}^{N-1}x_{k}(\boldsymbol{\theta})\lvert b_{k}\rangle .`$ By Lemma [3](proofs.md#Thmlemma3), each amplitude $`x_{k}`$ is a product of $`\cos\theta_{j}`$ and $`\sin\theta_{j}`$ factors along a unique root-to-leaf path, and in the complex case an additional leaf phase factor.

(A) Real case. Here $`x_{k}\in\mathbb{R}`$ and the metric is $`g^{\mathbb{R}}_{i\ell}=\langle \dot{\psi}_{i}\vert \dot{\psi}_{\ell}\rangle =\sum_{k}\partial_{\theta_{i}}x_{k}\,\partial_{\ell}x_{k}`$ (cf. Definition [4](chart-and-geometry.md#Thmdefinition4)(A)).

Step 1: Support structure of $`\lvert \dot{\psi}_{i}\rangle`$. Fix an internal node index $`i\in\{1,\ldots,N-1\}`$ and let $`\mathcal{L}(i)`$ be the set of leaf indices in the subtree rooted at $`i`$. If $`k\notin\mathcal{L}(i)`$ then the product formula for $`x_{k}`$ contains no factor $`\cos\theta_{i}`$ or $`\sin\theta_{i}`$, hence $`\partial_{\theta_{i}}x_{k}=0`$. Therefore $`\lvert \dot{\psi}_{i}\rangle`$ is supported only on $`\{\lvert b_{k}\rangle :k\in\mathcal{L}(i)\}`$.

Moreover, writing $`\mathcal{L}(i)=\mathcal{L}_{L}(i)\,\dot{\cup}\,\mathcal{L}_{R}(i)`$ as the leaves in the left and right child subtrees of $`i`$, Lemma [3](proofs.md#Thmlemma3) implies that for $`k\in\mathcal{L}_{L}(i)`$ the amplitude has the form

<a id="A5.Ex102"></a>

```math
x_{k}(\boldsymbol{\theta})=\alpha_{k}(\boldsymbol{\theta})\cos\theta_{i},
```

where $`\alpha_{k}`$ does not depend on $`\theta_{i}`$, and for $`k\in\mathcal{L}_{R}(i)`$,

<a id="A5.Ex103"></a>

```math
x_{k}(\boldsymbol{\theta})=\beta_{k}(\boldsymbol{\theta})\sin\theta_{i},
```

where $`\beta_{k}`$ does not depend on $`\theta_{i}`$. Thus

<a id="A5.E116"></a>

```math
{\partial_{\theta_{i}}x_{k}=\begin{cases}-\alpha_{k}\sin\theta_{i}=-\tan\theta_{i}\,x_{k},&k\in\mathcal{L}_{L}(i),\\[2.0pt]
\ \ \beta_{k}\cos\theta_{i}=\ \cot\theta_{i}\,x_{k},&k\in\mathcal{L}_{R}(i),\\[2.0pt]
0,&k\notin\mathcal{L}(i).\end{cases}}\tag{E1}
```

Step 2: Diagonal entries. Define the incoming mass at node $`i`$ by

<a id="A5.Ex104"></a>

```math
S_{i}(\boldsymbol{\theta}):=\sum_{k\in\mathcal{L}(i)}x_{k}(\boldsymbol{\theta})^{2}.
```

Using ([E1](proofs.md#A5.E116)),

<a id="A5.EGx5"></a>

```math
\begin{aligned}g^{\mathbb{R}}_{ii} & =\sum_{k=0}^{N-1}(\partial_{\theta_{i}}x_{k})^{2}\end{aligned}
```

```math
=\tan^{2}\theta_{i}\sum_{k\in\mathcal{L}_{L}(i)}x_{k}^{2}\;+\;\cot^{2}\theta_{i}\sum_{k\in\mathcal{L}_{R}(i)}x_{k}^{2}.
```

But from the factorization $`x_{k}=\alpha_{k}\cos\theta_{i}`$ on $`\mathcal{L}_{L}(i)`$ and $`x_{k}=\beta_{k}\sin\theta_{i}`$ on $`\mathcal{L}_{R}(i)`$ we have

<a id="A5.Ex107"></a>

```math
\sum_{k\in\mathcal{L}_{L}(i)}x_{k}^{2}=\cos^{2}\theta_{i}\sum_{k\in\mathcal{L}_{L}(i)}\alpha_{k}^{2},\qquad\sum_{k\in\mathcal{L}_{R}(i)}x_{k}^{2}=\sin^{2}\theta_{i}\sum_{k\in\mathcal{L}_{R}(i)}\beta_{k}^{2}.
```

Furthermore, by construction the descendant angles below $`i`$ parameterize unit vectors on each child subtree, so both sums of squares of $`\alpha_{k}`$ and $`\beta_{k}`$ equal the same incoming mass:

<a id="A5.Ex108"></a>

```math
\sum_{k\in\mathcal{L}_{L}(i)}\alpha_{k}^{2}=\sum_{k\in\mathcal{L}_{R}(i)}\beta_{k}^{2}=S_{i}(\boldsymbol{\theta}).
```

Substituting gives

<a id="A5.Ex109"></a>

```math
g^{\mathbb{R}}_{ii}=\tan^{2}\theta_{i}\cdot\cos^{2}\theta_{i}\,S_{i}+\cot^{2}\theta_{i}\cdot\sin^{2}\theta_{i}\,S_{i}=\sin^{2}\theta_{i}\,S_{i}+\cos^{2}\theta_{i}\,S_{i}=S_{i}(\boldsymbol{\theta}).
```

Finally, $`S_{i}`$ depends only on the ancestor split probabilities above $`i`$. Let $`d_{i}=\lfloor\log_{2}i\rfloor`$ be the depth of node $`i`$ (root at depth $`0`$), and let $`s_{1},\ldots,s_{d_{i}}`$ be the ancestor indices from the root $`s_{1}=1`$ to the parent of $`i`$, with corresponding left/right bits $`q_{1},\ldots,q_{d_{i}}`$ (MSB first) as in the theorem statement. Each ancestor contributes a factor $`\cos^{2}\theta_{s_{m}}`$ if the path goes left and $`\sin^{2}\theta_{s_{m}}`$ if it goes right, hence

<a id="A5.Ex110"></a>

```math
g^{\mathbb{R}}_{ii}(\boldsymbol{\theta})=S_{i}(\boldsymbol{\theta})=\prod_{m=1}^{d_{i}}\begin{cases}\cos^{2}\theta_{s_{m}},&q_{m}=0,\\
\sin^{2}\theta_{s_{m}},&q_{m}=1,\end{cases}
```

which is ([10](chart-and-geometry.md#S2.E10)). For $`i=1`$ the product is empty, so $`g^{\mathbb{R}}_{11}=1`$.

Step 3: Off-diagonal entries (diagonality). Let $`i\neq\ell`$ be internal nodes. If $`\mathcal{L}(i)\cap\mathcal{L}(\ell)=\varnothing`$, then $`\lvert \dot{\psi}_{i}\rangle`$ and $`\lvert \dot{\psi}_{\ell}\rangle`$ have disjoint supports, hence $`g^{\mathbb{R}}_{i\ell}=0`$.

Otherwise, one of the nodes is an ancestor of the other; assume $`\ell`$ lies in the subtree of $`i`$. Then on the entire leaf set $`\mathcal{L}(\ell)`$, the scaling factor in ([E1](proofs.md#A5.E116)) is constant (all leaves of $`\mathcal{L}(\ell)`$ lie entirely in either $`\mathcal{L}_{L}(i)`$ or $`\mathcal{L}_{R}(i)`$). Therefore there exists a real constant $`c`$ such that for all $`k\in\mathcal{L}(\ell)`$, $`\partial_{\theta_{i}}x_{k}=c\,x_{k}`$ and $`\partial_{\theta_{i}}x_{k}=0`$ for $`k\notin\mathcal{L}(i)\supseteq\mathcal{L}(\ell)`$. Hence

<a id="A5.Ex111"></a>

```math
g^{\mathbb{R}}_{i\ell}=\sum_{k}\partial_{\theta_{i}}x_{k}\,\partial_{\ell}x_{k}=c\sum_{k\in\mathcal{L}(\ell)}x_{k}\,\partial_{\ell}x_{k}.
```

But the vector of amplitudes restricted to $`\mathcal{L}(\ell)`$ has fixed squared norm $`S_{\ell}(\boldsymbol{\theta})`$, and $`\partial_{\ell}`$ varies only the internal angle at node $`\ell`$, so differentiating $`S_{\ell}=\sum_{k\in\mathcal{L}(\ell)}x_{k}^{2}`$ with respect to $`\theta_{\ell}`$ yields $`0=\partial_{\ell}S_{\ell}=2\sum_{k\in\mathcal{L}(\ell)}x_{k}\,\partial_{\ell}x_{k}.`$ Thus $`\sum_{k\in\mathcal{L}(\ell)}x_{k}\,\partial_{\ell}x_{k}=0`$ and therefore $`g^{\mathbb{R}}_{i\ell}=0`$. This proves that $`\boldsymbol{g}^{\mathbb{R}}`$ is diagonal.

(B) Complex case. Write $`x_{k}(\boldsymbol{\theta})=r_{k}(\boldsymbol{\theta})e^{\mathrm{i}\phi_{k}}`$ with $`\phi_{k}=\theta_{N+k}`$ and $`r_{k}\geq 0`$ determined by the internal-node angles. By Lemma [3](proofs.md#Thmlemma3), the magnitudes $`r_{k}`$ obey the same cosine–sine tree factorization as in the real case.

For $`1\leq i,\ell\leq N-1`$ (magnitude parameters), we have $`\partial_{\theta_{i}}x_{k}=(\partial_{\theta_{i}}r_{k})e^{\mathrm{i}\phi_{k}}`$ and therefore

<a id="A5.Ex112"></a>

```math
\langle \dot{\psi}_{i}\vert \dot{\psi}_{\ell}\rangle =\sum_{k}\overline{\partial_{\theta_{i}}x_{k}}\,\partial_{\ell}x_{k}=\sum_{k}(\partial_{\theta_{i}}r_{k})(\partial_{\ell}r_{k}),
```

which is real. The same support/cancellation argument from (A) applies to $`\{r_{k}\}`$, so the magnitude block is diagonal and $`g^{\mathbb{C}}_{ii}=g^{\mathbb{R}}_{ii}`$ for $`1\leq i\leq N-1`$.

For phase parameters $`i=N+k`$, $`\partial_{\theta_{N+k}}x_{m}=\mathrm{i}\,x_{m}\,\delta_{mk},`$ so

<a id="A5.Ex113"></a>

```math
g^{\mathbb{C}}_{N+k,\,N+k}=\mathrm{Re}\,\langle \dot{\psi}_{N+k}\vert \dot{\psi}_{N+k}\rangle =\mathrm{Re}\,(\vert \mathrm{i}x_{k}\vert ^{2})=\vert x_{k}\vert ^{2},
```

and for $`k\neq m`$ the phase derivatives are orthogonal, hence the phase block is diagonal.

Finally, for a magnitude parameter $`i\leq N-1`$ and a phase parameter $`N+k`$,

<a id="A5.Ex114"></a>

```math
\langle \dot{\psi}_{i}\vert \dot{\psi}_{N+k}\rangle =\overline{\partial_{\theta_{i}}x_{k}}\,(\mathrm{i}x_{k})=\mathrm{i}\,(\partial_{\theta_{i}}r_{k})\,r_{k},
```

which is purely imaginary; taking the real part (Definition [4](chart-and-geometry.md#Thmdefinition4)(B)) gives $`g^{\mathbb{C}}_{i,\,N+k}=0`$. Hence the full complex metric is diagonal and ([11](chart-and-geometry.md#S2.E11)) holds. ∎

<a id="A5.SSx3"></a>

## Proof of Theorem [2](gradient-access.md#Thmtheorem2)

<a id="A5.SSx3.1"></a>

**Proof.**

Let $`\lvert \psi(\boldsymbol{\theta})\rangle =U_{\boldsymbol{\theta}}\lvert 0\rangle ^{\otimes n}`$ be differentiable and define $`E(\boldsymbol{\theta})=\langle \psi(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle`$ with $`H=H^{\dagger}`$. Using the product rule,

<a id="A5.E117"></a>

```math
{\partial_{\theta_{i}}E(\boldsymbol{\theta})=\partial_{\theta_{i}}\big(\langle \psi(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle \big)\\
=\langle \dot{\psi}_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle +\langle \psi(\boldsymbol{\theta})\rvert H\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle .}\tag{E2}
```

Since $`H`$ is Hermitian, $`\langle \psi(\boldsymbol{\theta})\rvert H\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle =\overline{\langle \dot{\psi}_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle },`$ hence

<a id="A5.Ex115"></a>

```math
\partial_{\theta_{i}}E(\boldsymbol{\theta})=2\,\mathrm{Re}\!\left[\langle \dot{\psi}_{i}(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle \right].
```

Whenever $`g_{i,i}(\boldsymbol{\theta})\gt 0`$, Definition [5](gradient-access.md#Thmdefinition5) gives

<a id="A5.Ex116"></a>

```math
\lvert e_{i}(\boldsymbol{\theta})\rangle =\frac{\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle }{\sqrt{g_{i,i}(\boldsymbol{\theta})}},\qquad\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle =\sqrt{g_{i,i}(\boldsymbol{\theta})}\,\lvert e_{i}(\boldsymbol{\theta})\rangle ,
```

so substituting yields

<a id="A5.Ex117"></a>

```math
\partial_{\theta_{i}}E(\boldsymbol{\theta})=2\sqrt{g_{i,i}(\boldsymbol{\theta})}\;\mathrm{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert \,H\,\lvert \psi(\boldsymbol{\theta})\rangle \right],
```

which is ([13](gradient-access.md#S3.E13)). ∎

<a id="A5.SSx4"></a>

## Proof of Theorem [3](gradient-access.md#Thmtheorem3)

<a id="A5.SSx4.1"></a>

**Proof.**

We treat magnitude parameters and phase parameters separately. Throughout, write $`\lvert \psi(\boldsymbol{\theta})\rangle =\sum_{k=0}^{N-1}x_{k}(\boldsymbol{\theta})\lvert b_{k}\rangle .`$

In this proof, the vectors $`\boldsymbol{\theta}^{(i)}`$ are interpreted as gate-parameter assignments on the fixed Hopf circuit skeleton. They are not required to satisfy the canonical coordinate-range restrictions of Definition [2](chart-and-geometry.md#Thmdefinition2); only the resulting circuit action matters.

(A) Magnitude parameters ($`1\leq i\leq N-1`$). Fix an internal node index $`i`$ and let $`\mathcal{L}(i)`$ be the set of leaves in the subtree rooted at $`i`$. Let $`\mathcal{L}_{L}(i)`$ and $`\mathcal{L}_{R}(i)`$ denote the leaves in the left and right child subtrees of $`i`$. By Lemma [3](proofs.md#Thmlemma3), every amplitude $`x_{k}`$ factorizes as a product of trigonometric terms along the unique root-to-leaf path.

Step 1: Isolating the subtree rooted at $`i`$. Let $`A_{L}^{(i)}`$ and $`A_{R}^{(i)}`$ be the ancestor sets defined in the theorem statement. Setting $`\theta^{(i)}_{a}=0`$ for $`a\in A_{L}^{(i)}`$ forces the corresponding ancestor factor to be $`\cos\theta^{(i)}_{a}=\cos 0=1`$ and the sibling-branch factor to be $`0`$ (since $`\sin 0=0`$), thereby eliminating amplitudes supported on the “wrong” sibling subtrees. Similarly, setting $`\theta^{(i)}_{a}=\pi/2`$ for $`a\in A_{R}^{(i)}`$ forces $`\sin\theta^{(i)}_{a}=\sin(\pi/2)=1`$ and $`\cos(\pi/2)=0`$, again eliminating all sibling subtrees not on the path to $`i`$. Consequently,

<a id="A5.Ex118"></a>

```math
x_{k}\!\big(\boldsymbol{\theta}^{(i)}\big)=0\qquad\text{for all }k\notin\mathcal{L}(i),
```

so the prepared state has support only on $`\{\lvert b_{k}\rangle :k\in\mathcal{L}(i)\}`$.

Moreover, for $`k\in\mathcal{L}(i)`$ the product of all ancestor factors above $`i`$ equals $`1`$ under $`\boldsymbol{\theta}^{(i)}`$, so the remaining dependence of $`x_{k}(\boldsymbol{\theta}^{(i)})`$ comes only from the angle at $`i`$ and from strict descendants of $`i`$. Since the construction keeps all descendant angles unchanged (([18](gradient-access.md#S3.E18))), the relative amplitude pattern inside the subtree is preserved.

Step 2: A two-branch decomposition at node $`i`$. There exist vectors supported on the two child subtrees,

<a id="A5.Ex119"></a>

```math
\lvert \psi_{L}\rangle :=\sum_{k\in\mathcal{L}_{L}(i)}\alpha_{k}(\boldsymbol{\theta})\lvert b_{k}\rangle ,\qquad\lvert \psi_{R}\rangle :=\sum_{k\in\mathcal{L}_{R}(i)}\beta_{k}(\boldsymbol{\theta})\lvert b_{k}\rangle ,
```

such that $`\|\psi_{L}\|=\|\psi_{R}\|=1`$ and, on the subtree $`\mathcal{L}(i)`$,

<a id="A5.E118"></a>

```math
{\sum_{k\in\mathcal{L}(i)}x_{k}(\boldsymbol{\theta})\lvert b_{k}\rangle =\gamma_{i}(\boldsymbol{\theta})\Big(\cos\theta_{i}\,\lvert \psi_{L}\rangle +\sin\theta_{i}\,\lvert \psi_{R}\rangle \Big),}\tag{E3}
```

where the nonnegative scalar $`\gamma_{i}(\boldsymbol{\theta})`$ is the common product of ancestor trigonometric factors above node $`i`$ (hence independent of $`\theta_{i}`$). Equivalently, $`\gamma_{i}(\boldsymbol{\theta})^{2}=\sum_{k\in\mathcal{L}(i)}\vert x_{k}(\boldsymbol{\theta})\vert ^{2}.`$

Differentiating ([E3](proofs.md#A5.E118)) with respect to $`\theta_{i}`$ yields

<a id="A5.E119"></a>

```math
{\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle =\gamma_{i}(\boldsymbol{\theta})\Big(-\sin\theta_{i}\,\lvert \psi_{L}\rangle +\cos\theta_{i}\,\lvert \psi_{R}\rangle \Big),}\tag{E4}
```

and therefore

<a id="A5.Ex120"></a>

```math
\|\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle \|^{2}=\gamma_{i}(\boldsymbol{\theta})^{2}.
```

In the real case this equals $`g^{\mathbb{R}}_{i,i}(\boldsymbol{\theta})`$, and in the complex case (magnitude block) it equals $`g^{\mathbb{C}}_{i,i}(\boldsymbol{\theta})`$ (consistent with Theorem [1](chart-and-geometry.md#Thmtheorem1)). Hence, whenever $`g_{i,i}(\boldsymbol{\theta})\gt 0`$, we have

<a id="A5.E120"></a>

```math
{\lvert e_{i}(\boldsymbol{\theta})\rangle =\frac{\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle }{\sqrt{g_{i,i}(\boldsymbol{\theta})}}=-\sin\theta_{i}\,\lvert \psi_{L}\rangle +\cos\theta_{i}\,\lvert \psi_{R}\rangle .}\tag{E5}
```

Step 3: Real ansatz synthesis. By Step 1, under $`\boldsymbol{\theta}^{(i)}`$ all amplitude weight is clamped into the subtree $`\mathcal{L}(i)`$ and the ancestor prefactor becomes $`1`$. By ([18](gradient-access.md#S3.E18)), the descendant angles below $`i`$ are unchanged, so the same $`\lvert \psi_{L}\rangle`$ and $`\lvert \psi_{R}\rangle`$ appear. Finally, by ([15](gradient-access.md#S3.E15)) we shift $`\theta_{i}\mapsto\theta_{i}+\pi/2`$, and using $`\cos(\theta_i+\pi/2)=-\sin\theta_{i}`$ and $`\sin(\theta_i+\pi/2)=\cos\theta_{i},`$ the resulting subtree superposition is exactly $`-\sin\theta_{i}\,\lvert \psi_{L}\rangle +\cos\theta_{i}\,\lvert \psi_{R}\rangle .`$ Therefore,

<a id="A5.Ex121"></a>

```math
\mathrm{HopfReal}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n}=\lvert e_{i}(\boldsymbol{\theta})\rangle .
```

Step 4: Complex ansatz synthesis (magnitude block). For $`\mathrm{HopfComplex}`$, the magnitude recursion is identical to the real case, and the phase assignment ([20](gradient-access.md#S3.E20)) keeps the original leaf phases on $`\mathcal{L}(i)`$ while setting all other phases to $`0`$. Since Step 1 makes $`x_{k}(\boldsymbol{\theta}^{(i)})=0`$ for $`k\notin\mathcal{L}(i)`$, those outside phases are irrelevant. Inside the subtree, the same trigonometric shift at node $`i`$ produces the same linear combination ([E5](proofs.md#A5.E120)), now with the correct inherited leaf phases. Hence,

<a id="A5.Ex122"></a>

```math
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n}=\lvert e_{i}(\boldsymbol{\theta})\rangle .
```

Finally, if node $`i`$ has depth $`d`$ (root at depth $`0`$), its subtree contains exactly $`N/2^{d}`$ leaves. The resulting state is supported only on the corresponding computational basis states, so its support size is at most $`N/2^{d}`$. Equality holds whenever its amplitude on every leaf in $`\mathcal{L}(i)`$ is nonzero.

(B) Phase parameters ($`i=N+\ell`$). Fix $`\ell\in\{0,\ldots,N-1\}`$ and write $`x_{\ell}(\boldsymbol{\theta})=r_{\ell}(\boldsymbol{\theta})e^{\mathrm{i}\theta_{N+\ell}}`$ with $`r_{\ell}\geq 0`$. Then

<a id="A5.Ex123"></a>

```math
\lvert \dot{\psi}_{N+\ell}(\boldsymbol{\theta})\rangle =\mathrm{i}\,x_{\ell}(\boldsymbol{\theta})\,\lvert b_{\ell}\rangle ,\qquad\left\|\lvert \dot{\psi}_{N+\ell}(\boldsymbol{\theta})\rangle \right\|^{2}=\vert x_{\ell}(\boldsymbol{\theta})\vert ^{2},
```

so

<a id="A5.E121"></a>

```math
{\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle =\frac{\mathrm{i}\,x_{\ell}(\boldsymbol{\theta})}{\vert x_{\ell}(\boldsymbol{\theta})\vert }\lvert b_{\ell}\rangle .}\tag{E6}
```

Now construct $`\boldsymbol{\theta}^{(i)}`$ as in ([22](gradient-access.md#S3.E22))–([23](gradient-access.md#S3.E23)). The ancestor clamping forces the magnitude to be supported on the single leaf $`\ell`$ with $`r_{\ell}(\boldsymbol{\theta}^{(i)})=1`$ and $`r_{m}(\boldsymbol{\theta}^{(i)})=0`$ for $`m\neq\ell`$. The phase shift sets $`\theta^{(i)}_{N+\ell}=\theta_{N+\ell}+\pi/2,`$ so the prepared state is

<a id="A5.Ex124"></a>

```math
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n}=\mathrm{i}\,e^{\mathrm{i}\theta_{N+\ell}}\lvert b_{\ell}\rangle .
```

Since

<a id="A5.Ex125"></a>

```math
x_{\ell}(\boldsymbol{\theta})=r_{\ell}(\boldsymbol{\theta})e^{\mathrm{i}\theta_{N+\ell}},\qquad\vert x_{\ell}(\boldsymbol{\theta})\vert =r_{\ell}(\boldsymbol{\theta}),
```

we have

<a id="A5.Ex126"></a>

```math
\frac{\mathrm{i}\,x_{\ell}(\boldsymbol{\theta})}{\vert x_{\ell}(\boldsymbol{\theta})\vert }=\mathrm{i}\,e^{\mathrm{i}\theta_{N+\ell}}.
```

Hence

<a id="A5.Ex127"></a>

```math
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(i)}\big)\lvert 0\rangle ^{\otimes n}=\mathrm{i}\,e^{\mathrm{i}\theta_{N+\ell}}\lvert b_{\ell}\rangle =\frac{\mathrm{i}\,x_{\ell}(\boldsymbol{\theta})}{\vert x_{\ell}(\boldsymbol{\theta})\vert }\lvert b_{\ell}\rangle =\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle .
```

This completes the proof. ∎

<a id="A5.SSx5"></a>

## Proof of Theorem [4](gradient-access.md#Thmtheorem4)

<a id="A5.SSx5.1"></a>

**Proof.**

We prove the per-setting compiled execution bound separately for magnitude-layer settings and phase-layer settings.

Throughout, the forward Hopf ansatz uses a fixed binary-tree gate skeleton with $`N-1`$ single-qubit rotation locations, one for each internal node of the Hopf tree (Definition [3](chart-and-geometry.md#Thmdefinition3)). At each location, the original gate already carries some data-register controls; their number is at most $`n-1`$.

(A) Magnitude-layer settings: $`\lvert \Phi_{d}\rangle`$ and $`\lvert \Omega_{d}\rangle`$. Fix a depth $`d\in\{0,1,\ldots,n-1\}`$. The reduced index register $`I`$ contains $`d`$ qubits and labels the nodes in the layer $`\mathcal{L}_{d}=\{2^{d},\ldots,2^{d+1}-1\}`$ through the offset $`r=i-2^{d}`$.

For a fixed magnitude parameter $`i\in\mathcal{L}_{d}`$, Theorem [3](gradient-access.md#Thmtheorem3)(A) shows that the tangent-state preparation is obtained from the same Hopf skeleton by replacing each magnitude angle by one of a constant-size set of values:

<a id="A5.Ex128"></a>

```math
0,\qquad\tfrac{\pi}{2},\qquad\theta_{j},\qquad\theta_{j}+\tfrac{\pi}{2}.
```

When the preparation is made coherent over all $`i\in\mathcal{L}_{d}`$, the nonzero assignment at a gate of depth $`h`$ is selected by a single mixed-polarity prefix condition on the layer index: selection of its right subtree for $`h\lt d`$, equality with the gate for $`h=d`$, or equality with its depth-$`d`$ ancestor for $`h\gt d`$. Hence each gate location requires only $`O(1)`$ index-conditioned cases, with no angle multiplexing.

By Lemma [2](algorithms.md#Thmlemma2), the CNOT cost of a controlled single-qubit gate is linear in the total number of controls once the small-control regime is absorbed into constants. Hence each index-conditioned case costs $`O(n)`$ CNOTs. Since there are $`N-1=O(N)`$ gate locations and only $`O(1)`$ cases per location, the total compiled execution cost of the indexed derivative setting $`\lvert \Phi_{d}(\boldsymbol{\theta})\rangle`$ is

<a id="A5.Ex129"></a>

```math
O(N)\times O(n)=O(nN).
```

For the branch setting $`\lvert \Omega_{d}(\boldsymbol{\theta})\rangle`$, one ancilla qubit is added to coherently combine the baseline state and the indexed derivative preparation. This adds at most one extra control to the same family of gates and therefore does not change the asymptotic bound:

<a id="A5.Ex130"></a>

```math
\mathrm{Cost}\big(\lvert \Omega_{d}(\boldsymbol{\theta})\rangle \big)=O(nN).
```

This proves the real-case claim and the magnitude-block part of the complex-case claim.

(B) Phase-layer settings: $`\lvert \Phi_{\mathrm{ph}}\rangle`$ and $`\lvert \Omega_{\mathrm{ph}}\rangle`$. Now consider the complex Hopf ansatz and the phase parameters $`\theta_{N+\ell}`$, $`\ell=0,\ldots,N-1`$. Let $`P`$ be the $`n`$-qubit phase-index register.

For each fixed leaf label $`\ell=(b_{n}b_{n-1}\cdots b_{1})_{2},`$ Theorem [3](gradient-access.md#Thmtheorem3)(B) gives an explicit parameter vector $`\boldsymbol{\theta}^{(N+\ell)}`$ such that

<a id="A5.Ex131"></a>

```math
\mathrm{HopfComplex}\!\big(\boldsymbol{\theta}^{(N+\ell)}\big)\lvert 0\rangle ^{\otimes n}=\lvert e_{N+\ell}(\boldsymbol{\theta})\rangle .
```

Thus the indexed phase derivative preparation is obtained by implementing the map $`\ell\mapsto\boldsymbol{\theta}^{(N+\ell)}`$ coherently under control of $`P`$.

For an internal-node gate location $`j`$, the phase-case construction ([22](gradient-access.md#S3.E22)) sets its magnitude angle either to $`0`$ or to $`\pi/2`$, depending only on whether $`j`$ lies on the selected root-to-leaf path and on the next branch bit. Hence each internal-node gate location requires at most one nontrivial index-conditioned case.

For a promoted final-layer $`R_{\mathbb{C}}`$ gate corresponding to a sibling leaf pair, the phase-case construction ([23](gradient-access.md#S3.E23)) requires at most two nontrivial cases: one selecting the left leaf and one selecting the right leaf. Each such case uses only the two local leaf phases already attached to that gate location, together with the constant shift $`+\pi/2`$. Again, no angle multiplexing is required.

Therefore the indexed phase-layer implementation uses the same Hopf skeleton as the forward complex ansatz, with only $`O(1)`$ case splitting at each original gate location. Each case inherits the original data-register controls and acquires at most $`n`$ additional controls from the phase-index register $`P`$. By Lemma [2](algorithms.md#Thmlemma2), each such case costs $`O(n)`$ CNOTs. Since there are $`O(N)`$ gate locations and $`O(1)`$ cases per location, the indexed phase-derivative setting satisfies

<a id="A5.Ex132"></a>

```math
\mathrm{Cost}\big(\lvert \Phi_{\mathrm{ph}}(\boldsymbol{\theta})\rangle \big)=O(nN).
```

For the phase-branch setting $`\lvert \Omega_{\mathrm{ph}}(\boldsymbol{\theta})\rangle`$, the additional ancilla used to form the signed superposition adds at most one extra control and therefore preserves the same asymptotic bound:

<a id="A5.Ex133"></a>

```math
\mathrm{Cost}\big(\lvert \Omega_{\mathrm{ph}}(\boldsymbol{\theta})\rangle \big)=O(nN).
```

Combining parts (A) and (B), every gradient-access configuration used in the Hopf protocol—whether derivative or branch, magnitude or phase—has compiled execution cost $`O(nN)`$. ∎

<a id="A5.SSx6"></a>

## Proof of Lemma [1](algorithms.md#Thmlemma1)

<a id="A5.SSx6.1"></a>

**Proof.**

Let $`N=2^{n}`$.

**(A) Real case. Let $`\mathbf{x}\in\mathbb{R}^{N}`$ satisfy $`\sum_{\ell=0}^{N-1}x_{\ell}^{2}=1`$. Use a tree-indexed array of weights, with leaf weights**

<a id="A5.Ex134"></a>

```math
w_{N+\ell}:=x_{\ell}^{2},\qquad\ell=0,\ldots,N-1.
```

For internal nodes $`j=N-1,N-2,\ldots,1`$, define

<a id="A5.Ex135"></a>

```math
w_{j}:=w_{2j}+w_{2j+1}.
```

Then $`w_{j}`$ is the total squared mass in the subtree rooted at $`j`$. In particular,

<a id="A5.Ex136"></a>

```math
w_{j}=\sum_{\ell\in\mathcal{L}(j)}x_{\ell}^{2},\qquad w_{2j}=\sum_{\ell\in\mathcal{L}_{L}(j)}x_{\ell}^{2},\qquad w_{2j+1}=\sum_{\ell\in\mathcal{L}_{R}(j)}x_{\ell}^{2}.
```

Thus

<a id="A5.Ex137"></a>

```math
S_{L}(j)=\sqrt{w_{2j}},\qquad S_{R}(j)=\sqrt{w_{2j+1}},\qquad S(j)=\sqrt{w_{j}}.
```

For every non-final internal node $`1\leq j\leq N/2-1`$, set

<a id="A5.Ex138"></a>

```math
\theta_{j}=\begin{cases}\mathrm{atan2}\!\big(\sqrt{w_{2j+1}},\sqrt{w_{2j}}\big),&w_{j}\gt 0,\\[2.0pt]
0,&w_{j}=0.\end{cases}
```

Equivalently, when $`w_{j}\gt 0`$,

<a id="A5.Ex139"></a>

```math
\cos\theta_{j}=\frac{\sqrt{w_{2j}}}{\sqrt{w_{j}}},\qquad\sin\theta_{j}=\frac{\sqrt{w_{2j+1}}}{\sqrt{w_{j}}}.
```

Therefore the forward Hopf recursion sends the incoming magnitude $`\sqrt{w_{j}}`$ at node $`j`$ to child magnitudes

<a id="A5.Ex140"></a>

```math
\sqrt{w_{j}}\cos\theta_{j}=\sqrt{w_{2j}},\qquad\sqrt{w_{j}}\sin\theta_{j}=\sqrt{w_{2j+1}}.
```

By induction from the root down to depth $`n-1`$, the recursion reproduces the correct incoming magnitude for every final-layer sibling pair.

It remains only to recover the signs inside each final sibling pair. For a final-layer node $`j=N/2+k`$, the two descendant leaves are $`2k`$ and $`2k+1`$. Let

<a id="A5.Ex141"></a>

```math
S_{k}:=\sqrt{x_{2k}^{2}+x_{2k+1}^{2}}.
```

If $`S_{k}\gt 0`$, set

<a id="A5.Ex142"></a>

```math
\theta_{N/2+k}=\mathrm{atan2}\!\big(x_{2k+1},x_{2k}\big)\pmod{2\pi}.
```

Then

<a id="A5.Ex143"></a>

```math
S_{k}\cos\theta_{N/2+k}=x_{2k},\qquad S_{k}\sin\theta_{N/2+k}=x_{2k+1}.
```

If $`S_{k}=0`$, both leaf amplitudes vanish and the angle is arbitrary; the convention in the lemma sets it to zero. Hence the real forward Hopf map reproduces the signed vector $`\mathbf{x}`$. Nonuniqueness occurs exactly at zero-subtree or zero-leaf-pair points.

**(B) Complex case. Let $`\mathbf{x}\in\mathbb{C}^{N}`$ satisfy $`\sum_{\ell=0}^{N-1}\vert x_{\ell}\vert ^{2}=1`$, and write**

<a id="A5.Ex144"></a>

```math
x_{\ell}=r_{\ell}e^{\mathrm{i}\phi_{\ell}},\qquad r_{\ell}=\vert x_{\ell}\vert \geq 0,
```

with the convention $`\phi_{\ell}=0`$ when $`x_{\ell}=0`$.

Apply the same bottom-up construction to the nonnegative leaf weights

<a id="A5.Ex145"></a>

```math
w_{N+\ell}:=r_{\ell}^{2}=\vert x_{\ell}\vert ^{2}.
```

For every internal node $`j=1,\ldots,N-1`$, set

<a id="A5.Ex146"></a>

```math
\theta_{j}=\begin{cases}\mathrm{atan2}\!\big(\sqrt{w_{2j+1}},\sqrt{w_{2j}}\big),&w_{j}\gt 0,\\[2.0pt]
0,&w_{j}=0.\end{cases}
```

Since all magnitudes are nonnegative, these angles lie in $`[0,\pi/2]`$, and the forward magnitude recursion reproduces exactly the magnitudes $`r_{0},\ldots,r_{N-1}`$.

Finally set the leaf phases to

<a id="A5.Ex147"></a>

```math
\theta_{N+\ell}:=\phi_{\ell},\qquad\ell=0,\ldots,N-1.
```

Then the complex forward Hopf map gives

<a id="A5.Ex148"></a>

```math
r_{\ell}e^{\mathrm{i}\theta_{N+\ell}}=r_{\ell}e^{\mathrm{i}\phi_{\ell}}=x_{\ell}
```

for every leaf. If $`r_{\ell}=0`$, the corresponding phase is irrelevant, so the inverse is not unique at zero-amplitude leaves.

**Complexity. The bottom-up computation of all weights uses one traversal of a complete binary tree with $`O(N)`$ nodes. Once the weights are known, each magnitude angle is obtained in $`O(1)`$ time. In the real case, the final-layer signed angles are obtained directly from the pairs $`(x_{2k},x_{2k+1})`$. In the complex case, the leaf phases require one additional linear pass. Therefore**

<a id="A5.Ex149"></a>

```math
\text{Time}_{\mathrm{Hopf}^{-1}}(N)=O(N),\qquad\text{Memory}_{\mathrm{Hopf}^{-1}}(N)=O(N).
```

This proves the lemma. ∎

<a id="A5.SSx7"></a>

## Proof of Theorem [6](optimizers.md#Thmtheorem6)

<a id="A5.SSx7.1"></a>

**Proof.**

Fix a leaf (basis) index $`i\in\{0,\ldots,N-1\}`$. By Lemma [3](proofs.md#Thmlemma3), the Hopf amplitude $`x_{i}(\boldsymbol{\theta})`$ is a product of trigonometric factors contributed by the internal nodes on the unique root-to-leaf path of $`i`$, and (in the complex case) an additional leaf phase factor.

Let $`P(i)`$ denote the set of internal nodes on the root-to-parent path that generates $`x_{i}`$ (as in the theorem statement), and for $`j\in P(i)`$ let $`\varepsilon_{ij}\in\{0,1\}`$ indicate whether $`x_{i}`$ contains the factor $`\cos\theta_{j}`$ ($`\varepsilon_{ij}=0`$) or $`\sin\theta_{j}`$ ($`\varepsilon_{ij}=1`$).

(A) Real case. Assume $`x_{i}(\boldsymbol{\theta})\in\mathbb{R}`$. If $`j\notin P(i)`$, then the product formula for $`x_{i}`$ contains no factor depending on $`\theta_{j}`$, hence $`\partial x_{i}/\partial\theta_{j}=0`$.

Now assume $`j\in P(i)`$. Write the product formula as

<a id="A5.Ex150"></a>

```math
x_{i}(\boldsymbol{\theta})=f_{ij}(\theta_{j})\,u_{ij}(\boldsymbol{\theta}),
```

where $`u_{ij}`$ is the product of all factors in $`x_{i}`$ except the one at node $`j`$, so $`u_{ij}`$ is independent of $`\theta_{j}`$. Here

<a id="A5.Ex151"></a>

```math
f_{ij}(\theta_{j})=\begin{cases}\cos\theta_{j},&\varepsilon_{ij}=0,\\
\sin\theta_{j},&\varepsilon_{ij}=1.\end{cases}
```

The direct product derivative is always

<a id="A5.Ex152"></a>

```math
\frac{\partial x_{i}}{\partial\theta_{j}}=f^{\prime}_{ij}(\theta_{j})\,u_{ij}(\boldsymbol{\theta}).
```

On the regular set where $`f_{ij}(\theta_{j})\neq 0`$, this may be written as

<a id="A5.Ex153"></a>

```math
\frac{\partial x_{i}}{\partial\theta_{j}}=\frac{f^{\prime}_{ij}(\theta_{j})}{f_{ij}(\theta_{j})}x_{i}(\boldsymbol{\theta}).
```

Thus, if $`\varepsilon_{ij}=0`$, then

<a id="A5.Ex154"></a>

```math
\frac{\partial x_{i}}{\partial\theta_{j}}=-\tan\theta_{j}\,x_{i},
```

while if $`\varepsilon_{ij}=1`$, then

<a id="A5.Ex155"></a>

```math
\frac{\partial x_{i}}{\partial\theta_{j}}=\cot\theta_{j}\,x_{i}.
```

At points where $`f_{ij}(\theta_{j})=0`$, the displayed $`\tan`$ or $`\cot`$ formula is interpreted by the direct product derivative above, equivalently by continuous extension. This proves ([C5](optimizers.md#A3.E65)) with the stated regular-set convention.

Finally, since in the real case $`g^{\mathbb{R}}_{j\ell}(\boldsymbol{\theta})=\sum_{i=0}^{N-1}\frac{\partial x_{i}}{\partial\theta_{j}}\frac{\partial x_{i}}{\partial\theta_{\ell}}`$ (cf. ([8](chart-and-geometry.md#S2.E8))), the matrix identity $`\boldsymbol{g}^{\mathbb{R}}=(\boldsymbol{J}^{\mathbb{R}})^{\top}\boldsymbol{J}^{\mathbb{R}}`$ is exactly the statement that $`g^{\mathbb{R}}_{j\ell}`$ is the Euclidean inner product of columns $`j`$ and $`\ell`$ of $`\boldsymbol{J}^{\mathbb{R}}`$, proving ([C6](optimizers.md#A3.E66)).

(B) Complex case. Write

<a id="A5.Ex156"></a>

```math
x_{i}(\boldsymbol{\theta})=r_{i}(\boldsymbol{\theta})\,e^{\mathrm{i}\theta_{N+i}},\qquad r_{i}(\boldsymbol{\theta})\geq 0,
```

where $`r_{i}`$ is built from internal-node angles as in Lemma [3](proofs.md#Thmlemma3)(B).

For $`1\leq j\leq N-1`$ (magnitude parameters), the phase $`e^{\mathrm{i}\theta_{N+i}}`$ is independent of $`\theta_{j}`$, so the same product/log-derivative argument as in (A) applies verbatim and yields the two path-derivative cases in ([C7](optimizers.md#A3.E67)), together with the zero case when $`j\notin P(i)`$, now with $`x_{i}`$ complex-valued.

For the phase parameter $`j=N+i`$,

<a id="A5.Ex157"></a>

```math
\frac{\partial x_{i}}{\partial\theta_{N+i}}=r_{i}(\boldsymbol{\theta})\cdot\mathrm{i}e^{\mathrm{i}\theta_{N+i}}=\mathrm{i}\,x_{i}(\boldsymbol{\theta}),
```

and for $`j=N+m`$ with $`m\neq i`$ the derivative is $`0`$ since $`x_{i}`$ does not depend on that phase. This completes ([C7](optimizers.md#A3.E67)).

Finally, by Definition [4](chart-and-geometry.md#Thmdefinition4)(B),

<a id="A5.Ex158"></a>

```math
g^{\mathbb{C}}_{j\ell}(\boldsymbol{\theta})=\mathrm{Re}\,\langle\partial_{\theta_{j}}\psi,\partial_{\ell}\psi\rangle=\mathrm{Re}\sum_{i=0}^{N-1}\overline{\frac{\partial x_{i}}{\partial\theta_{j}}}\frac{\partial x_{i}}{\partial\theta_{\ell}},
```

which is precisely the $`(j,\ell)`$ entry of $`\mathrm{Re}\big((\boldsymbol{J}^{\mathbb{C}})^{\dagger}\boldsymbol{J}^{\mathbb{C}}\big)`$, proving ([C8](optimizers.md#A3.E68)). ∎
