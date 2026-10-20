[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="A4"></a>

# Appendix D Metrological objectives as chain-rule observables

This appendix states a limited transfer principle: the Hopf gradient-access framework used for VQE also applies to metrological costs whose local differentials are expectation-value differentials of Hermitian observables. The main text gives the estimator for the energy objective

<a id="A4.Ex67"></a>

```math
E(\boldsymbol{\theta})=\langle \psi(\boldsymbol{\theta})\rvert H\lvert \psi(\boldsymbol{\theta})\rangle ,
```

but the same calculus applies to smooth scalar objectives built from measured probabilities, expectation values, or local phase derivatives.

We keep QFI and fixed-readout CFI separate. A QFI objective is an ultimate local probe-state benchmark assuming an optimal measurement; a fixed-readout CFI objective is the Fisher information of a specified measurement protocol. The fixed-readout Ramsey construction below optimizes a probe for a fixed parity readout at a calibrated operating point; noisy optimality, global phase estimation, and measurement co-optimization are outside this example.

The probe is

<a id="A4.Ex68"></a>

```math
\lvert \psi(\boldsymbol{\theta})\rangle =U_{\boldsymbol{\theta}}\lvert 0\rangle ^{\otimes n},
```

and the normalized coordinate tangent is the same as in Definition [5](gradient-access.md#Thmdefinition5),

<a id="A4.Ex69"></a>

```math
\lvert e_{i}(\boldsymbol{\theta})\rangle =\frac{\lvert \dot{\psi}_{i}(\boldsymbol{\theta})\rangle }{\sqrt{g_{i,i}(\boldsymbol{\theta})}}.
```

No new tangent-state synthesis is needed; only the Hermitian observable in the transition element changes.

<a id="A4.SSx1"></a>

## A general chain-rule form

Let a scalar cost

<a id="A4.Ex70"></a>

```math
\mathcal{C}(\boldsymbol{\theta})=f\!\left(m_{1}(\boldsymbol{\theta}),\ldots,m_{K}(\boldsymbol{\theta})\right)
```

depend on finitely many expectation values

<a id="A4.Ex71"></a>

```math
m_{\alpha}(\boldsymbol{\theta})=\langle \psi(\boldsymbol{\theta})\rvert O_{\alpha}\lvert \psi(\boldsymbol{\theta})\rangle ,\qquad O_{\alpha}=O_{\alpha}^{\dagger}.
```

Define the classical chain-rule coefficients

<a id="A4.Ex72"></a>

```math
a_{\alpha}(\boldsymbol{\theta}):=\frac{\partial f}{\partial m_{\alpha}}\bigg\vert _{m=m(\boldsymbol{\theta})}
```

and the associated chain-rule observable

<a id="A4.E102"></a>

```math
O_{\mathcal{C}}(\boldsymbol{\theta}):=\sum_{\alpha=1}^{K}a_{\alpha}(\boldsymbol{\theta})\,O_{\alpha}.\tag{D1}
```

Then

<a id="A4.E103"></a>

```math
\partial_{\theta_{i}}\mathcal{C}(\boldsymbol{\theta})=2\sqrt{g_{i,i}(\boldsymbol{\theta})}\,\operatorname{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert O_{\mathcal{C}}(\boldsymbol{\theta})\lvert \psi(\boldsymbol{\theta})\rangle \right].\tag{D2}
```

Thus the VQE formula ([13](gradient-access.md#S3.E13)) is recovered by taking $`K=1`$, $`O_{1}=H`$, and $`f(m_{1})=m_{1}`$.

For each $`i`$, define the same signed branch state as in the main text,

<a id="A4.Ex73"></a>

```math
\lvert \varphi_{i}^{(s)}(\boldsymbol{\theta})\rangle =\frac{\lvert \psi(\boldsymbol{\theta})\rangle +s\lvert e_{i}(\boldsymbol{\theta})\rangle }{\sqrt{2}},\qquad s\in\{+1,-1\}.
```

For any Hermitian observable $`O`$, write

<a id="A4.Ex74"></a>

```math
B_{\chi}[O]:=\langle \chi\rvert O\lvert \chi\rangle .
```

Then

<a id="A4.E104"></a>

```math
\operatorname{Re}\!\left[\langle e_{i}\rvert O\lvert \psi\rangle \right]=s\left(B_{\varphi_{i}^{(s)}}[O]-\frac{1}{2}\left(B_{\psi}[O]+B_{e_{i}}[O]\right)\right),\tag{D3}
```

where

<a id="A4.Ex75"></a>

```math
B_{\psi}[O]=\langle \psi\rvert O\lvert \psi\rangle ,\qquad B_{e_{i}}[O]=\langle e_{i}\rvert O\lvert e_{i}\rangle .
```

Applying ([D3](metrological-objectives.md#A4.E104)) to $`O=O_{\mathcal{C}}(\boldsymbol{\theta})`$ gives

<a id="A4.E105"></a>

```math
\partial_{\theta_{i}}\mathcal{C}(\boldsymbol{\theta})=2\sqrt{g_{i,i}(\boldsymbol{\theta})}\,s\left(B_{\varphi_{i}^{(s)}}[O_{\mathcal{C}}]-\frac{1}{2}\left(B_{\psi}[O_{\mathcal{C}}]+B_{e_{i}}[O_{\mathcal{C}}]\right)\right).\tag{D4}
```

In using this identity, $`O_{\mathcal{C}}(\boldsymbol{\theta})`$ is built from chain-rule coefficients evaluated at the current probe state and then held fixed while branch-state transition elements are estimated. Thus $`B_{\varphi_{i}^{(s)}}[O_{\mathcal{C}}]`$ is a transition-moment quantity, not the nonlinear cost $`\mathcal{C}(\varphi_{i}^{(s)})`$. If the coefficients are estimated from empirical data, the resulting gradient is generally a stochastic approximation to $`\nabla\mathcal{C}`$; nonlinear functions of noisy estimates can introduce bias without debiasing. This matters especially for inverse-Fisher objectives, whose prefactors can be sensitive to small probabilities or slopes.

<a id="A4.SSx2"></a>

## Pure-state QFI objective

Consider a single unknown parameter $`\phi`$ encoded by a unitary

<a id="A4.Ex76"></a>

```math
U_{\phi}=\exp(-\mathrm{i}\phi G),
```

with Hermitian generator $`G`$. For a pure probe state, the quantum Fisher information is

<a id="A4.E106"></a>

```math
F_{Q}(\boldsymbol{\theta})=4\operatorname{Var}_{\boldsymbol{\theta}}(G)=4\left(\langle G^{2}\rangle_{\boldsymbol{\theta}}-\langle G\rangle_{\boldsymbol{\theta}}^{2}\right).\tag{D5}
```

Here and below

<a id="A4.Ex77"></a>

```math
\langle O\rangle_{\boldsymbol{\theta}}:=\langle \psi(\boldsymbol{\theta})\rvert O\lvert \psi(\boldsymbol{\theta})\rangle .
```

A natural minimized metrology cost is the regularized quantum Cramér–Rao proxy

<a id="A4.E107"></a>

```math
\mathcal{C}_{Q}(\boldsymbol{\theta})=\frac{1}{F_{Q}(\boldsymbol{\theta})+\eta},\qquad\eta>0.\tag{D6}
```

This cost defines a local probe-state benchmark associated with the optimal-measurement QFI. A fixed experimental readout realizes this benchmark only when it saturates the QFI at the operating point. By contrast, the fixed-readout Ramsey objective below is a classical Fisher-information objective for a specified parity measurement. Let

<a id="A4.Ex78"></a>

```math
\mu(\boldsymbol{\theta})=\langle G\rangle_{\boldsymbol{\theta}}.
```

From ([D5](metrological-objectives.md#A4.E106)),

<a id="A4.Ex79"></a>

```math
\partial_{\theta_{i}}F_{Q}=4\,\partial_{\theta_{i}}\langle G^{2}\rangle-8\mu\,\partial_{\theta_{i}}\langle G\rangle.
```

Therefore

<a id="A4.E108"></a>

```math
\partial_{\theta_{i}}F_{Q}=8\sqrt{g_{i,i}(\boldsymbol{\theta})}\left(T_{i}[G^{2}]-2\mu\,T_{i}[G]\right),\tag{D7}
```

where

<a id="A4.Ex80"></a>

```math
T_{i}[O]:=\operatorname{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert O\lvert \psi(\boldsymbol{\theta})\rangle \right].
```

The minimized-cost gradient is then

<a id="A4.E109"></a>

```math
\partial_{\theta_{i}}\mathcal{C}_{Q}=-\frac{\partial_{\theta_{i}}F_{Q}}{\left(F_{Q}+\eta\right)^{2}}.\tag{D8}
```

Equivalently, this is the general chain-rule formula ([D2](metrological-objectives.md#A4.E103)) with the chain-rule observable

<a id="A4.E110"></a>

```math
O_{\mathcal{C}_{Q}}=-\frac{4G^{2}-8\mu G}{\left(F_{Q}+\eta\right)^{2}}.\tag{D9}
```

Thus the QFI objective requires the same Hopf tangent-state machinery as VQE, with the Hamiltonian $`H`$ replaced by the observables $`G`$ and $`G^{2}`$, or by their linear combination ([D9](metrological-objectives.md#A4.E110)).

<a id="A4.SSx3"></a>

## Fixed-readout Ramsey CFI objective

We next record the exact chain-rule observable for a fixed-readout Ramsey CFI objective. Let

<a id="A4.Ex81"></a>

```math
U_{\phi}=\exp(-\mathrm{i}\phi G),\qquad A_{\mathrm{e}}=W^{\dagger}\Pi_{\mathrm{e}}W,
```

and define the effective even-parity POVM element

<a id="A4.Ex82"></a>

```math
M_{\mathrm{e}}(\phi)=U_{\phi}^{\dagger}A_{\mathrm{e}}U_{\phi}.
```

At a fixed operating point $`\phi_{0}`$, define

<a id="A4.Ex83"></a>

```math
p(\boldsymbol{\theta})=\langle \psi(\boldsymbol{\theta})\rvert M_{\mathrm{e}}(\phi_{0})\lvert \psi(\boldsymbol{\theta})\rangle ,
```

and

<a id="A4.Ex84"></a>

```math
q(\boldsymbol{\theta})=\left.\frac{\partial}{\partial\phi}\langle \psi(\boldsymbol{\theta})\rvert M_{\mathrm{e}}(\phi)\lvert \psi(\boldsymbol{\theta})\rangle \right\vert _{\phi=\phi_{0}}.
```

Equivalently,

<a id="A4.Ex85"></a>

```math
q(\boldsymbol{\theta})=\langle \psi(\boldsymbol{\theta})\rvert D_{\mathrm{e}}(\phi_{0})\lvert \psi(\boldsymbol{\theta})\rangle ,
```

where the derivative observable is

<a id="A4.E111"></a>

```math
D_{\mathrm{e}}(\phi_{0}):=\left.\frac{\partial M_{\mathrm{e}}(\phi)}{\partial\phi}\right\vert _{\phi=\phi_{0}}=\mathrm{i}U_{\phi_{0}}^{\dagger}\left[G,A_{\mathrm{e}}\right]U_{\phi_{0}}.\tag{D10}
```

At an interior operating point satisfying

<a id="A4.Ex86"></a>

```math
0<p(\boldsymbol{\theta})<1,
```

the fixed-readout classical Fisher information is

<a id="A4.Ex87"></a>

```math
F_{C}(\boldsymbol{\theta})=\frac{q(\boldsymbol{\theta})^{2}}{p(\boldsymbol{\theta})[1-p(\boldsymbol{\theta})]},
```

and the minimized local inverse-CFI cost is

<a id="A4.Ex88"></a>

```math
\mathcal{C}_{\mathrm{Ram}}(\boldsymbol{\theta})=\frac{1}{F_{C}(\boldsymbol{\theta})+\eta},\qquad\eta>0.
```

The derivatives below are taken on this interior domain. At the boundary points $`p(\boldsymbol{\theta})\in\{0,1\}`$, the quotient defining $`F_{C}`$ is singular and must instead be treated by a separate limiting prescription or by explicitly regularizing the Bernoulli denominator $`p(1-p)`$. The additive constant $`\eta`$ in the inverse cost does not by itself regularize this denominator.

For $`0<p(\boldsymbol{\theta})<1`$, the chain-rule coefficients with respect to $`p`$ and $`q`$ are

<a id="A4.EGx4"></a>

```math
\begin{aligned}a_{p}(\boldsymbol{\theta}) & :=\frac{\partial\mathcal{C}_{\mathrm{Ram}}}{\partial p}=\frac{q(\boldsymbol{\theta})^{2}[1-2p(\boldsymbol{\theta})]}{[F_{C}(\boldsymbol{\theta})+\eta]^{2}p(\boldsymbol{\theta})^{2}[1-p(\boldsymbol{\theta})]^{2}},\end{aligned}
```

<a id="A4.E112"></a>

```math
\begin{aligned}a_{q}(\boldsymbol{\theta}) & :=\frac{\partial\mathcal{C}_{\mathrm{Ram}}}{\partial q}=-\frac{2q(\boldsymbol{\theta})}{[F_{C}(\boldsymbol{\theta})+\eta]^{2}p(\boldsymbol{\theta})[1-p(\boldsymbol{\theta})]}.\end{aligned}\tag{D11}
```

Thus the exact fixed-readout Ramsey chain-rule observable is

<a id="A4.E113"></a>

```math
O_{\mathcal{C}_{\mathrm{Ram}}}(\boldsymbol{\theta})=a_{p}(\boldsymbol{\theta})M_{\mathrm{e}}(\phi_{0})+a_{q}(\boldsymbol{\theta})D_{\mathrm{e}}(\phi_{0}),\tag{D12}
```

and the Hopf gradient is

<a id="A4.Ex90"></a>

```math
\partial_{\theta_{i}}\mathcal{C}_{\mathrm{Ram}}=2\sqrt{g_{i,i}(\boldsymbol{\theta})}\,\operatorname{Re}\!\left[\langle e_{i}(\boldsymbol{\theta})\rvert O_{\mathcal{C}_{\mathrm{Ram}}}(\boldsymbol{\theta})\lvert \psi(\boldsymbol{\theta})\rangle \right].
```

In a centered finite-difference implementation, the derivative observable $`D_{\mathrm{e}}(\phi_{0})`$ is replaced by

<a id="A4.E114"></a>

```math
D_{\mathrm{e},\delta}=\frac{M_{\mathrm{e}}(\phi_{0}+\delta)-M_{\mathrm{e}}(\phi_{0}-\delta)}{2\delta}.\tag{D13}
```

The corresponding finite-difference chain-rule observable is

<a id="A4.Ex91"></a>

```math
O_{\mathcal{C}_{\mathrm{Ram}},\delta}(\boldsymbol{\theta})=a_{p}(\boldsymbol{\theta})M_{\mathrm{e}}(\phi_{0})+a_{q}(\boldsymbol{\theta})D_{\mathrm{e},\delta}.
```

<a id="A4.SSx4"></a>

## Compiled-setting scaling

The metrology transfer does not alter the Hopf state preparation or gradient access circuits. The indexed tangent-state preparations $`\lvert \Phi_{d}(\boldsymbol{\theta})\rangle`$, the indexed branch states $`\lvert \Omega_{d}(\boldsymbol{\theta})\rangle`$, and, for the complex ansatz, the corresponding phase-block states $`\lvert \Phi_{\mathrm{ph}}(\boldsymbol{\theta})\rangle`$ and $`\lvert \Omega_{\mathrm{ph}}(\boldsymbol{\theta})\rangle`$, are exactly the same objects as in Section [III](gradient-access.md#S3). Only the observable or readout configuration used to evaluate the transition moments changes:

<a id="A4.Ex92"></a>

```math
H\quad\longrightarrow\quad O_{\mathcal{C}}(\boldsymbol{\theta}).
```

Let $`K_{\mathrm{read}}`$ denote the number of logical readout configurations or moment types required by the chosen metrological cost. For the pure-QFI formula above, $`K_{\mathrm{read}}`$ includes the moments needed for $`G`$ and $`G^{2}`$, together with any additional grouping cost from their observable decompositions. For a centered finite-difference Ramsey implementation, the logical phase settings are

<a id="A4.Ex93"></a>

```math
\phi_{0},\qquad\phi_{0}+\delta,\qquad\phi_{0}-\delta,
```

so $`K_{\mathrm{read}}=3`$ before any further hardware-specific repetition or sampling allocation.

Thus the Hopf compiled-setting count appears per logical readout configuration. For the real ansatz,

<a id="A4.Ex94"></a>

```math
N_{\mathrm{set,met}}^{\mathbb{R}}=K_{\mathrm{read}}(1+2n),
```

and for the complex ansatz,

<a id="A4.Ex95"></a>

```math
N_{\mathrm{set,met}}^{\mathbb{C}}=K_{\mathrm{read}}\bigl[1+2(n+1)\bigr].
```

Equivalently, if $`C_{\mathrm{cfg}}^{\mathrm{read}}`$ denotes the compiled execution cost of one Hopf gradient-access setting for one logical readout configuration, then

<a id="A4.E115"></a>

```math
\boxed{C_{\mathrm{cfg}}(\nabla\mathcal{C}_{\mathrm{met}})=O\!\left(K_{\mathrm{read}}n\right)C_{\mathrm{cfg}}^{\mathrm{read}}=O\!\left(K_{\mathrm{read}}\log N\right)C_{\mathrm{cfg}}^{\mathrm{read}}.}\tag{D14}
```

When $`K_{\mathrm{read}}`$ is independent of $`N`$, this preserves the $`O(\log N)`$ compiled-setting scaling of the VQE case. If the chosen metrological observable or readout scheme requires a number of logical settings that grows with problem size, that cost enters multiplicatively.

Equation ([D14](metrological-objectives.md#A4.E115)) gives the compiled-setting and hardware-reconfiguration scaling for metrological objectives. The corresponding sampling budget is componentwise and is determined by the chosen metrological readout: by the measurements needed for $`G`$ and $`G^{2}`$ in the pure-QFI formulation, or by the Ramsey probability and local phase-slope measurements in the fixed-readout formulation above.
