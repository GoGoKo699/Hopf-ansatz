[Scientific text](README.md) · [Reviewer route](../review/README.md) · [Repository](../README.md)

<a id="A1"></a>

# Appendix A Algorithms and subroutines

<a id="A1.SSx1"></a>

## Inverse Hopf map

<a id="Thmlemma1"></a>

**Lemma 1 (Inverse Hopf map and classical cost).**

Let $`n\in\mathbb{N}`$ and $`N=2^{n}`$.

(A) Real case. Let $`\mathbf{x}=(x_{0},\ldots,x_{N-1})\in\mathbb{R}^{N}`$ satisfy $`\sum_{\ell=0}^{N-1}x_{\ell}^{2}=1`$. For each internal node $`j\in\{1,\ldots,N-1\}`$, define its left- and right-subtree leaf sets by

<a id="A1.Ex30"></a>

```math
\mathcal{L}_{L}(j):=\{\ell:\lvert b_{\ell}\rangle \text{ is a leaf in the left subtree of }j\},\qquad\mathcal{L}_{R}(j):=\{\ell:\lvert b_{\ell}\rangle \text{ is a leaf in the right subtree of }j\},
```

and the associated subtree norms

<a id="A1.E43"></a>

```math
{{S_{L}(j):=\sqrt{\sum_{\ell\in\mathcal{L}_{L}(j)}x_{\ell}^{2}},\qquad S_{R}(j):=\sqrt{\sum_{\ell\in\mathcal{L}_{R}(j)}x_{\ell}^{2}},\qquad S(j):=\sqrt{S_{L}(j)^{2}+S_{R}(j)^{2}}.}}\tag{A1}
```

For the non-final internal nodes, namely $`1\leq j\leq N/2-1`$, set

<a id="A1.E44"></a>

```math
{\theta_{j}:=\begin{cases}\mathrm{atan2}\!\big(S_{R}(j),S_{L}(j)\big),&S(j)\gt 0,\\[2.0pt]
0,&S(j)=0,\end{cases}\qquad\theta_{j}\in[0,\pi/2].}\tag{A2}
```

For the final-layer internal nodes, write $`j=N/2+k`$ with $`k=0,\ldots,N/2-1`$. These nodes split the sibling leaf pair $`(2k,2k+1)`$, and their angles are chosen with the signed quadrant convention

<a id="A1.E45"></a>

```math
{\theta_{N/2+k}:=\begin{cases}\mathrm{atan2}\!\big(x_{2k+1},x_{2k}\big)\pmod{2\pi},&x_{2k}^{2}+x_{2k+1}^{2}\gt 0,\\[2.0pt]
0,&x_{2k}=x_{2k+1}=0.\end{cases}}\tag{A3}
```

With the ranges of Definition [2](chart-and-geometry.md#Thmdefinition2)(A), the forward Hopf recursion then reproduces the signed vector $`\mathbf{x}`$ exactly on the regular set. At zero-subtree and zero-leaf-pair points, the inverse is not unique; the choices above fix one convention.

(B) Complex case. Let $`\mathbf{x}=(x_{0},\ldots,x_{N-1})\in\mathbb{C}^{N}`$ satisfy $`\sum_{\ell}\vert x_{\ell}\vert ^{2}=1`$. Write

<a id="A1.Ex31"></a>

```math
x_{\ell}=r_{\ell}e^{\mathrm{i}\phi_{\ell}},\qquad r_{\ell}=\vert x_{\ell}\vert \geq 0,
```

where

<a id="A1.Ex32"></a>

```math
\phi_{\ell}:=\arg(x_{\ell})\pmod{2\pi}\in[0,2\pi),
```

with the convention $`\phi_{\ell}:=0`$ if $`x_{\ell}=0`$. For each internal node $`j\in\{1,\ldots,N-1\}`$, define

<a id="A1.Ex33"></a>

```math
S_{L,r}(j):=\sqrt{\sum_{\ell\in\mathcal{L}_{L}(j)}r_{\ell}^{2}},\qquad S_{R,r}(j):=\sqrt{\sum_{\ell\in\mathcal{L}_{R}(j)}r_{\ell}^{2}},\qquad S_{r}(j):=\sqrt{S_{L,r}(j)^{2}+S_{R,r}(j)^{2}}.
```

Set

<a id="A1.E46"></a>

```math
{\theta_{j}:=\begin{cases}\mathrm{atan2}\!\big(S_{R,r}(j),S_{L,r}(j)\big),&S_{r}(j)\gt 0,\\[2.0pt]
0,&S_{r}(j)=0,\end{cases}\qquad j=1,\ldots,N-1.}\tag{A4}
```

The leaf phases are

<a id="A1.E47"></a>

```math
{\theta_{N+\ell}:=\phi_{\ell}\in[0,2\pi),\qquad\ell=0,\ldots,N-1.}\tag{A5}
```

Then the complex Hopf recursion of Definition [2](chart-and-geometry.md#Thmdefinition2)(B) reproduces $`\mathbf{x}`$ up to the chosen phase convention at zero entries. The inverse is not unique whenever some magnitude subtree or leaf amplitude vanishes.

We denote by $`\mathrm{Hopf}^{-1}`$ any choice of the inverse Hopf map constructed above on the regular set. There is an $`O(N)`$-time, $`O(N)`$-memory procedure to compute it. Use a tree-indexed array of weights with leaves $`w_{N+\ell}=x_{\ell}^{2}`$ in the real case and $`w_{N+\ell}=\vert x_{\ell}\vert ^{2}`$ in the complex case. A single bottom-up pass gives

<a id="A1.Ex34"></a>

```math
w_{j}=w_{2j}+w_{2j+1},\qquad j=N-1,N-2,\ldots,1.
```

The subtree norms are $`S_{L}(j)=\sqrt{w_{2j}}`$ and $`S_{R}(j)=\sqrt{w_{2j+1}}`$. The real final-layer angles are then obtained directly from the signed pairs by Eq. ([A3](algorithms.md#A1.E45)); the complex leaf phases require one additional linear pass. Thus,

<a id="A1.Ex35"></a>

```math
\boxed{\text{Time}_{\mathrm{Hopf}^{-1}}(N)=O(N),\qquad\text{Memory}_{\mathrm{Hopf}^{-1}}(N)=O(N).}
```

<a id="A1.SSx2"></a>

## Gate-schedule routines

The gate-ordering routines output four parallel lists

<a id="A1.Ex36"></a>

```math
\mathtt{Ctrl},\qquad\mathtt{Anti},\qquad\mathtt{Targ},\qquad\mathtt{Index}.
```

All entries are read position-wise. The masks $`\mathtt{Ctrl}`$ and $`\mathtt{Anti}`$ specify, respectively, controls on $`\lvert 1\rangle`$ and controls on $`\lvert 0\rangle`$. The mask $`\mathtt{Targ}`$ is a power of two specifying the single target qubit. We use the bitmask convention

<a id="A1.Ex37"></a>

```math
(q_{n},\ldots,q_{1})
```

from left to right, so the integer bit $`2^{t-1}`$ corresponds to qubit $`q_{t}`$. Thus, for example, for $`n=4`$, the mask $`1000`$ targets $`q_{4}`$, while $`0001`$ targets $`q_{1}`$.

In the real ansatz, every $`\mathtt{Index}`$ entry is a single integer $`j`$ and the corresponding gate is $`R_{y}(\theta_{j})`$. In the complex ansatz, $`\mathtt{Index}`$ entries have the following convention:

<a id="A1.Ex38"></a>

```math
\mathtt{Index}_{r}\in\{1,\ldots,N-1\}\quad\text{or}\quad\mathtt{Index}_{r}=[j,N+\ell_{L},N+\ell_{R}].
```

An integer entry denotes a magnitude-only $`R_{y}(\theta_{j})`$ gate. A length-three entry denotes a promoted final-layer gate

<a id="A1.Ex39"></a>

```math
R_{\mathbb{C}}(\theta_{j},\theta_{N+\ell_{L}},\theta_{N+\ell_{R}}).
```

We do not use length-one lists for magnitude-only complex gates.

<a id="algorithm-1"></a>

### Algorithm 1. Generate the order of quantum gates (real)

<a id="A1.SSx2.p3.3"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg0.l1"></a>Integer $`n`$ (number of qubits). |
| Output | <a id="alg0.l2"></a>Lists $`\mathtt{Ctrl}`$, $`\mathtt{Anti}`$, $`\mathtt{Targ}`$, $`\mathtt{Index}`$. |
| Complexity | <a id="alg0.l3"></a>$`O\!\left(nN\right)`$. |
| 1 | <a id="alg0.l4"></a>**procedure** HopfReal($`n`$) |
| 2 | <a id="alg0.l5"></a>&emsp;&emsp;$`N\leftarrow 2^{n}`$ |
| 3 | <a id="alg0.l6"></a>&emsp;&emsp;$`HW\leftarrow\text{{HammingWeight}}(n)`$ |
| 4 | <a id="alg0.l7"></a>&emsp;&emsp;$`\mathtt{Ctrl}\leftarrow[0,\dots,0]`$ $`\triangleright`$ $`n`$ zeros |
| 5 | <a id="alg0.l8"></a>&emsp;&emsp;$`\mathtt{Anti}\leftarrow[0,N-2^{n-1},\dots,N-2]`$ $`\triangleright`$ $`0`$; then $`N-2^{n-i}`$ for $`i=1,\ldots,n-1`$ |
| 6 | <a id="alg0.l9"></a>&emsp;&emsp;$`\mathtt{Targ}\leftarrow[2^{n-1},2^{n-2},\dots,2^{0}]`$ $`\triangleright`$ $`2^{n-i}`$ for $`i=1,\ldots,n`$ |
| 7 | <a id="alg0.l10"></a>&emsp;&emsp;$`\mathtt{Index}\leftarrow[1,2,4,\dots,2^{n-1}]`$ $`\triangleright`$ $`n`$ entries: $`\mathtt{Index}[i]=2^{i}`$ for $`i=0,\dots,n-1`$ (one per level on the initial root-to-leaf chain) |
| 8 | <a id="alg0.l11"></a>&emsp;&emsp;**for** $`k=1`$ **to** $`n-1`$ **do** |
| 9 | <a id="alg0.l12"></a>&emsp;&emsp;&emsp;$`A,B\leftarrow\text{{FindPairs}}(HW[k],HW[k+1])`$ |
| 10 | <a id="alg0.l13"></a>&emsp;&emsp;&emsp;$`\mathtt{Ctrl}.\text{extend}(A)`$ |
| 11 | <a id="alg0.l14"></a>&emsp;&emsp;&emsp;$`\mathtt{Anti}.\text{extend}(\text{{Anti}}(A))`$ |
| 12 | <a id="alg0.l15"></a>&emsp;&emsp;&emsp;$`\mathtt{Targ}.\text{extend}(B-A)`$ $`\triangleright`$ elementwise $`B_{i}-A_{i}`$ |
| 13 | <a id="alg0.l16"></a>&emsp;&emsp;&emsp;$`L\leftarrow\vert A\vert`$ |
| 14 | <a id="alg0.l17"></a>&emsp;&emsp;&emsp;**for** $`i=0`$ **to** $`L-1`$ **do** |
| 15 | <a id="alg0.l18"></a>&emsp;&emsp;&emsp;&emsp;&emsp;$`\mathtt{Index}.\text{append}(\text{{ThetaReal}}(n,A[i],B[i]))`$ |
| 16 | <a id="alg0.l19"></a>&emsp;&emsp;&emsp;**end** **for** |
| 17 | <a id="alg0.l20"></a>&emsp;&emsp;**end** **for** |
| 18 | <a id="alg0.l21"></a>&emsp;&emsp;**return** $`\mathtt{Ctrl}`$, $`\mathtt{Anti}`$, $`\mathtt{Targ}`$, $`\mathtt{Index}`$ |
| 19 | <a id="alg0.l22"></a>**end** **procedure** |

<a id="subroutine-1"></a>

### Subroutine 1. Group numbers by Hamming weight

<a id="A1.SSx2.p5.1"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg0.l1a"></a>Integer $`n`$ for the bit length. |
| Output | <a id="alg0.l2a"></a>A list $`HW`$ of length $`n+1`$, where $`HW[k]`$ contains all integers with Hamming weight $`k`$. |
| Complexity | <a id="alg0.l3a"></a>$`O\!\left(nN\right)`$. |
| 1 | <a id="alg0.l4a"></a>**procedure** HammingWeight($`n`$) |
| 2 | <a id="alg0.l5a"></a>&emsp;&emsp;$`N\leftarrow 2^{n}`$ |
| 3 | <a id="alg0.l6a"></a>&emsp;&emsp;Initialize $`HW`$ as $`[\,[\,],\,\ldots,\,[\,]\,]`$ $`\triangleright`$ $`n+1`$ empty lists |
| 4 | <a id="alg0.l7a"></a>&emsp;&emsp;**for** $`i=0`$ **to** $`N-1`$ **do** |
| 5 | <a id="alg0.l8a"></a>&emsp;&emsp;&emsp;$`k\leftarrow\text{count\_ones}(i)`$ |
| 6 | <a id="alg0.l9a"></a>&emsp;&emsp;&emsp;$`HW[k].\text{append}(i)`$ |
| 7 | <a id="alg0.l10a"></a>&emsp;&emsp;**end** **for** |
| 8 | <a id="alg0.l11a"></a>&emsp;&emsp;**return** $`HW`$ |
| 9 | <a id="alg0.l12a"></a>**end** **procedure** |

<a id="subroutine-2"></a>

### Subroutine 2. Find pairs between elements of two lists

<a id="A1.SSx2.p6.3"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg0.l1b"></a>Sorted lists $`A`$ and $`B`$ (ascending). |
| Output | <a id="alg0.l2b"></a>Two lists $`A^{\prime}`$ and $`B^{\prime}`$ (of equal length) obtained by scanning $`B`$ from largest to smallest and pairing each $`b`$ with the current largest $`a\in A`$ such that $`a\lt b`$ (with $`a`$ reusable), and emitting the resulting pairs ordered by increasing $`a`$, with $`b`$ decreasing within each group. |
| Complexity | <a id="alg0.l3b"></a>$`O(\vert A\vert +\vert B\vert )`$. |
| 1 | <a id="alg0.l4b"></a>**procedure** FindPairs($`A`$, $`B`$) |
| 2 | <a id="alg0.l5b"></a>&emsp;&emsp;$`\widetilde{A}\leftarrow[\,]`$, $`\widetilde{B}\leftarrow[\,]`$ |
| 3 | <a id="alg0.l6b"></a>&emsp;&emsp;$`i\leftarrow\vert A\vert -1`$ |
| 4 | <a id="alg0.l7b"></a>&emsp;&emsp;**for** $`b`$ **in** $`B`$ **traversed in reverse order** **do** |
| 5 | <a id="alg0.l8b"></a>&emsp;&emsp;&emsp;**while** $`i\geq 0`$ **and** $`A[i]\geq b`$ **do** |
| 6 | <a id="alg0.l9b"></a>&emsp;&emsp;&emsp;&emsp;&emsp;$`i\leftarrow i-1`$ |
| 7 | <a id="alg0.l10b"></a>&emsp;&emsp;&emsp;**end** **while** |
| 8 | <a id="alg0.l11b"></a>&emsp;&emsp;&emsp;**if** $`i\lt 0`$ **then** |
| 9 | <a id="alg0.l12b"></a>&emsp;&emsp;&emsp;&emsp;&emsp;**break** |
| 10 | <a id="alg0.l13a"></a>&emsp;&emsp;&emsp;**end** **if** |
| 11 | <a id="alg0.l14a"></a>&emsp;&emsp;&emsp;append $`A[i]`$ to $`\widetilde{A}`$ |
| 12 | <a id="alg0.l15a"></a>&emsp;&emsp;&emsp;append $`b`$ to $`\widetilde{B}`$ |
| 13 | <a id="alg0.l16a"></a>&emsp;&emsp;**end** **for** |
| 14 | <a id="alg0.l17a"></a>&emsp;&emsp;$`A^{\prime}\leftarrow[\,]`$, $`B^{\prime}\leftarrow[\,]`$ |
| 15 | <a id="alg0.l18a"></a>&emsp;&emsp;$`r\leftarrow\vert \widetilde{A}\vert -1`$ |
| 16 | <a id="alg0.l19a"></a>&emsp;&emsp;**while** $`r\geq 0`$ **do** |
| 17 | <a id="alg0.l20a"></a>&emsp;&emsp;&emsp;$`a\leftarrow\widetilde{A}[r]`$ |
| 18 | <a id="alg0.l21a"></a>&emsp;&emsp;&emsp;$`\ell\leftarrow r`$ |
| 19 | <a id="alg0.l22a"></a>&emsp;&emsp;&emsp;**while** $`\ell\gt 0`$ **and** $`\widetilde{A}[\ell-1]=a`$ **do** |
| 20 | <a id="alg0.l23"></a>&emsp;&emsp;&emsp;&emsp;&emsp;$`\ell\leftarrow\ell-1`$ |
| 21 | <a id="alg0.l24"></a>&emsp;&emsp;&emsp;**end** **while** |
| 22 | <a id="alg0.l25"></a>&emsp;&emsp;&emsp;**for** $`t=\ell`$ **to** $`r`$ **do** |
| 23 | <a id="alg0.l26"></a>&emsp;&emsp;&emsp;&emsp;&emsp;append $`\widetilde{A}[t]`$ to $`A^{\prime}`$ |
| 24 | <a id="alg0.l27"></a>&emsp;&emsp;&emsp;&emsp;&emsp;append $`\widetilde{B}[t]`$ to $`B^{\prime}`$ |
| 25 | <a id="alg0.l28"></a>&emsp;&emsp;&emsp;**end** **for** |
| 26 | <a id="alg0.l29"></a>&emsp;&emsp;&emsp;$`r\leftarrow\ell-1`$ |
| 27 | <a id="alg0.l30"></a>&emsp;&emsp;**end** **while** |
| 28 | <a id="alg0.l31"></a>&emsp;&emsp;**return** $`A^{\prime},B^{\prime}`$ |
| 29 | <a id="alg0.l32"></a>**end** **procedure** |

<a id="subroutine-3"></a>

### Subroutine 3. Find anti-controls based on control multiplicities

<a id="A1.SSx2.p7.3"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg0.l1c"></a>A list $`A`$ of integers representing control bitmasks. |
| Output | <a id="alg0.l2c"></a>A list $`\mathtt{Anti}`$ of integers, constructed as follows: for each distinct value $`a\in A`$ with multiplicity $`c`$, the list contains the values $`2^{c}-2^{i}`$ for $`i=c,c-1,\dots,1`$. The output depends only on multiplicities in $`A`$ and is ordered deterministically by ascending distinct values of $`a`$. |
| Complexity | <a id="alg0.l3c"></a>$`O\!\left(\vert A\vert \log\vert A\vert \right)`$. |
| 1 | <a id="alg0.l4c"></a>**procedure** Anti($`A`$) |
| 2 | <a id="alg0.l5c"></a>&emsp;&emsp;$`\mathtt{Anti}\leftarrow[\,]`$ |
| 3 | <a id="alg0.l6c"></a>&emsp;&emsp;$`\mathtt{Dict}\leftarrow`$ frequency map of elements in $`A`$ |
| 4 | <a id="alg0.l7c"></a>&emsp;&emsp;$`\mathtt{Keys}\leftarrow`$ sorted list of distinct elements of $`A`$ (ascending) |
| 5 | <a id="alg0.l8c"></a>&emsp;&emsp;**for** each $`a`$ **in** $`\mathtt{Keys}`$ **do** |
| 6 | <a id="alg0.l9c"></a>&emsp;&emsp;&emsp;$`c\leftarrow\mathtt{Dict}[a]`$ |
| 7 | <a id="alg0.l10c"></a>&emsp;&emsp;&emsp;**for** $`i=c`$ **down to** $`1`$ **do** |
| 8 | <a id="alg0.l11c"></a>&emsp;&emsp;&emsp;&emsp;&emsp;append $`(2^{c}-2^{i})`$ to $`\mathtt{Anti}`$ |
| 9 | <a id="alg0.l12c"></a>&emsp;&emsp;&emsp;**end** **for** |
| 10 | <a id="alg0.l13b"></a>&emsp;&emsp;**end** **for** |
| 11 | <a id="alg0.l14b"></a>&emsp;&emsp;**return** $`\mathtt{Anti}`$ |
| 12 | <a id="alg0.l15b"></a>**end** **procedure** |

<a id="subroutine-4"></a>

### Subroutine 4. Find θ index for the gate (real)

<a id="A1.SSx2.p8.3"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg0.l1d"></a>Integer $`n`$ for the bit length, integer $`a`$ encoding initial basis state, integer $`b`$ encoding target basis state. |
| Output | <a id="alg0.l2d"></a>An integer giving the $`\theta`$ index. |
| Complexity | <a id="alg0.l3d"></a>$`O(1)`$. |
| 1 | <a id="alg0.l4d"></a>**procedure** ThetaReal($`n`$, $`a`$, $`b`$) |
| 2 | <a id="alg0.l5d"></a>&emsp;&emsp;$`N\leftarrow 2^{n}`$ |
| 3 | <a id="alg0.l6d"></a>&emsp;&emsp;$`d\leftarrow b-a`$ |
| 4 | <a id="alg0.l7d"></a>&emsp;&emsp;**return** $`\dfrac{N+a}{2d}`$ |
| 5 | <a id="alg0.l8d"></a>**end** **procedure** |
| 6 | <a id="alg0.l9d"></a> |
| 7 | <a id="alg0.l10d"></a>**Remark (integrality).** For every $`(a,b)`$ produced by $`\text{FindPairs}(HW[k],HW[k+1])`$, the difference $`d=b-a`$ is a power of two, hence $`(N+a)/(2d)\in\mathbb{Z}`$. |

<a id="subroutine-5"></a>

### Subroutine 5. Find θ index for the gate (complex)

<a id="A1.SSx2.p9.3"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg0.l1e"></a>Integer $`n`$, integer $`a`$ (initial state), integer $`b`$ (target state). |
| Output | <a id="alg0.l2e"></a>Either one integer index, or a length-three list of indices. |
| Complexity | <a id="alg0.l3e"></a>$`O(1)`$. |
| 1 | <a id="alg0.l4e"></a>**procedure** ThetaComplex($`n`$, $`a`$, $`b`$) |
| 2 | <a id="alg0.l5e"></a>&emsp;&emsp;$`N\leftarrow 2^{n}`$ |
| 3 | <a id="alg0.l6e"></a>&emsp;&emsp;$`d\leftarrow b-a`$ |
| 4 | <a id="alg0.l7e"></a>&emsp;&emsp;$`j\leftarrow\dfrac{N+a}{2d}`$ |
| 5 | <a id="alg0.l8e"></a>&emsp;&emsp;**if** $`d=1`$ **then** |
| 6 | <a id="alg0.l9e"></a>&emsp;&emsp;&emsp;**return** $`[j,\,N+a,\,N+b]`$ $`\triangleright`$ leaf-level transition, requiring phase angles |
| 7 | <a id="alg0.l10e"></a>&emsp;&emsp;**else** |
| 8 | <a id="alg0.l11d"></a>&emsp;&emsp;&emsp;**return** $`j`$ |
| 9 | <a id="alg0.l12d"></a>&emsp;&emsp;**end** **if** |
| 10 | <a id="alg0.l13c"></a>**end** **procedure** |
| 11 | <a id="alg0.l14c"></a> |
| 12 | <a id="alg0.l15c"></a>**Remark (integrality).** For every $`(a,b)`$ produced by $`\text{FindPairs}(HW[k],HW[k+1])`$, the difference $`d=b-a`$ is a power of two, hence $`(N+a)/(2d)\in\mathbb{Z}`$. |

<a id="algorithm-2"></a>

### Algorithm 2. Generate the order of quantum gates

<a id="A1.SSx2.p10.3"></a>

| Line | Algorithm |
| ---: | --- |
| Input | <a id="alg1.l1"></a>Integer $`n`$. |
| Output | <a id="alg1.l2"></a>Lists $`\mathtt{Ctrl}`$, $`\mathtt{Anti}`$, $`\mathtt{Targ}`$, $`\mathtt{Index}`$, where each $`\mathtt{Index}`$ entry is either an integer or a length-$`3`$ list. |
| Complexity | <a id="alg1.l3"></a>$`O\!\left(nN\right)`$. |
| 1 | <a id="alg1.l4"></a>**procedure** HopfComplex($`n`$) |
| 2 | <a id="alg1.l5"></a>&emsp;&emsp;**if** $`n=1`$ **then** |
| 3 | <a id="alg1.l6"></a>&emsp;&emsp;&emsp;**return** $`[0]`$, $`[0]`$, $`[1]`$, $`[[1,2,3]]`$ |
| 4 | <a id="alg1.l7"></a>&emsp;&emsp;**end** **if** |
| 5 | <a id="alg1.l8"></a>&emsp;&emsp;$`N\leftarrow 2^{n}`$ |
| 6 | <a id="alg1.l9"></a>&emsp;&emsp;$`HW\leftarrow\text{{HammingWeight}}(n)`$ |
| 7 | <a id="alg1.l10"></a>&emsp;&emsp;$`\mathtt{Ctrl}\leftarrow[0,\dots,0]`$ |
| 8 | <a id="alg1.l11"></a>&emsp;&emsp;$`\mathtt{Anti}\leftarrow[0,N-2^{n-1},\dots,N-2]`$ |
| 9 | <a id="alg1.l12"></a>&emsp;&emsp;$`\mathtt{Targ}\leftarrow[2^{n-1},2^{n-2},\dots,2^{0}]`$ |
| 10 | <a id="alg1.l13"></a>&emsp;&emsp;$`\mathtt{Index}\leftarrow[1,2,4,\dots,2^{n-2},\,[2^{n-1},\,N,\,N+1]]`$ $`\triangleright`$ $`n`$ initial gates: the first $`n-1`$ are magnitude-only; the last is promoted to $`R_{\mathbb{C}}`$ to attach the first sibling-pair phases |
| 11 | <a id="alg1.l14"></a>&emsp;&emsp;**for** $`k=1`$ **to** $`n-1`$ **do** |
| 12 | <a id="alg1.l15"></a>&emsp;&emsp;&emsp;$`A,B\leftarrow\text{{FindPairs}}(HW[k],HW[k+1])`$ |
| 13 | <a id="alg1.l16"></a>&emsp;&emsp;&emsp;$`\mathtt{Ctrl}.\text{extend}(A)`$ |
| 14 | <a id="alg1.l17"></a>&emsp;&emsp;&emsp;$`\mathtt{Anti}.\text{extend}(\text{{Anti}}(A))`$ |
| 15 | <a id="alg1.l18"></a>&emsp;&emsp;&emsp;$`\mathtt{Targ}.\text{extend}(B-A)`$ |
| 16 | <a id="alg1.l19"></a>&emsp;&emsp;&emsp;$`L\leftarrow\vert A\vert`$ |
| 17 | <a id="alg1.l20"></a>&emsp;&emsp;&emsp;**for** $`i=0`$ **to** $`L-1`$ **do** |
| 18 | <a id="alg1.l21"></a>&emsp;&emsp;&emsp;&emsp;&emsp;$`\mathtt{Index}.\text{append}(\text{{ThetaComplex}}(n,A[i],B[i]))`$ |
| 19 | <a id="alg1.l22"></a>&emsp;&emsp;&emsp;**end** **for** |
| 20 | <a id="alg1.l23"></a>&emsp;&emsp;**end** **for** |
| 21 | <a id="alg1.l24"></a>&emsp;&emsp;**return** $`\mathtt{Ctrl}`$, $`\mathtt{Anti}`$, $`\mathtt{Targ}`$, $`\mathtt{Index}`$ |
| 22 | <a id="alg1.l25"></a>**end** **procedure** |

<a id="A1.SSx3"></a>

## CNOT counting model

<a id="Thmlemma2"></a>

**Lemma 2 (CNOT counting model for mixed-sign multi-controls).**

We count only CNOT gates and assume no clean ancilla.

Consider a single-qubit gate acting on one target and conditioned on a set of control qubits, where each control is specified either as a positive control (trigger on $`\lvert 1\rangle`$) or a negative control (trigger on $`\lvert 0\rangle`$). Let $`m`$ denote the total number of controls (positive plus negative). In the implementation, controls are specified by two disjoint bitmasks $`\mathtt{Ctrl}`$ and $`\mathtt{Anti}`$, and

<a id="A1.Ex40"></a>

```math
m=\mathrm{popcount}(\mathtt{Ctrl}\lor\mathtt{Anti}).
```

- **Sign-independence.** Negative controls are implemented by conjugating the corresponding control qubits by Pauli-$`X`$ (before and after the controlled operation). Since $`X`$ can be taken as a CNOT-free primitive in our accounting, the CNOT cost depends only on $`m`$ and not on how many controls are negative.

- **Controlled-$`R_{y}`$.** Let $`C_{R}(m)`$ be the CNOT cost of an $`m`$-controlled $`R_{y}`$ gate (no ancilla). We use the piecewise model



  <a id="A1.Ex41"></a>

  ```math
  C_{R}(m)=\begin{cases}0,&m=0,\\
  2^{m+1}-2,&1\leq m\leq 4,\\
  16(m+1)-40,&m\geq 5,\end{cases}
  ```



  <a id="A1.I1.i2.p1.2"></a>

  where the small-$`m`$ values reflect best-known no-ancilla decompositions [[51](../references/paper-bibliography.md#bib.bib10), [52](../references/paper-bibliography.md#bib.bib9)], and the linear regime for $`m\geq 5`$ follows the explicit constructions of [[53](../references/paper-bibliography.md#bib.bib8)] (with $`n_{g}=m+1`$ total qubits involved).

- **Controlled-$`R_{\mathbb{C}}`$ (special unitary single-qubit gate).** Let $`C_{U}(m)`$ be the CNOT cost of an $`m`$-controlled special unitary single-qubit gate (no ancilla), which in our application is the promoted $`R_{\mathbb{C}}`$ gate. We use



  <a id="A1.Ex42"></a>

  ```math
  C_{U}(m)=\begin{cases}0,&m=0,\\
  2^{m+1}-2,&1\leq m\leq 4,\\
  \begin{cases}20(m+1)-38,&m+1\ \mathrm{odd},\\
  20(m+1)-42,&m+1\ \mathrm{even},\end{cases}&m\geq 5,\end{cases}
  ```



  <a id="A1.I1.i3.p1.2"></a>

  again using standard small-control decompositions [[51](../references/paper-bibliography.md#bib.bib10), [52](../references/paper-bibliography.md#bib.bib9)] and the linear no-ancilla constructions of [[53](../references/paper-bibliography.md#bib.bib8)].

All rotation parameters are implemented directly as physical gate parameters; no multiplexing or angle recombination is used.

<a id="Thmtheorem5"></a>

**Theorem 5 (CNOT counts for Hopf ansatze under the piecewise no-ancilla model).**

Let $`n`$ be the number of data qubits, and let $`m`$ denote the number of controls on a multi-controlled single-qubit gate. Under the explicit no-ancilla decompositions used in this work, we count CNOTs using the piecewise model

<a id="A1.Ex43"></a>

```math
c_{\mathrm{small}}(m):=\begin{cases}0,&m=0,\\
2^{m+1}-2,&1\leq m\leq 4,\end{cases}
```

and for $`m\geq 5`$ the linear-cost formulas of Lemma [2](algorithms.md#Thmlemma2),

<a id="A1.Ex44"></a>

```math
c_{R_{y}}(m):=16(m+1)-40,\qquad c_{R_{\mathbb{C}}}(m):=\begin{cases}20(m+1)-38,&m+1\ \mathrm{odd},\\
20(m+1)-42,&m+1\ \mathrm{even}.\end{cases}
```

(Equivalently, $`c_{R_{y}}(m)=c_{\mathrm{small}}(m)`$ for $`m\leq 4`$ and $`c_{R_{\mathbb{C}}}(m)=c_{\mathrm{small}}(m)`$ for $`m\leq 4`$.)

For $`n\in\{2,3,4,5\}`$, the real-state and complex-state versions of the Hopf ansatz generally differ at the gate level, because in the complex ansatz the depth-$`(n-1)`$ rotations are promoted to $`R_{\mathbb{C}}`$ gates that inject leaf phases (Definition [3](chart-and-geometry.md#Thmdefinition3) and Algorithm [2](algorithms.md#algorithm-2)). Nevertheless, under the counting model above, their total CNOT counts coincide for these sizes:

<a id="A1.Ex45"></a>

```math
\begin{array}{c|c|c}n&G_{\mathbb{R}}(n)&G_{\mathbb{C}}(n)\\
\hline 2&4&4\\
3&24&24\\
4&100&100\\
5&360&360\end{array}
```

These values are obtained by applying the CNOT-counting model above to the full generated gate lists, including the promoted $`R_{\mathbb{C}}`$ gates in the complex case.

For general $`n\geq 2`$, the control-number distribution induced by the Hopf ordering is as follows.

For the real ansatz, for each $`m\in\{0,1,\dots,n-1\}`$, there are exactly $`\binom{n}{m}`$ $`m`$-controlled $`R_{y}`$ gates. Hence

<a id="A1.Ex46"></a>

```math
G_{\mathbb{R}}(n)=\sum_{m=0}^{n-1}\binom{n}{m}\,c_{R_{y}}(m).
```

For the complex ansatz, the gate types split as follows:

- for each $`m\in\{0,1,\dots,n-2\}`$, there are exactly $`\binom{n-1}{m}`$ $`m`$-controlled $`R_{y}`$ gates;

- for each $`m\in\{1,2,\dots,n-2\}`$, there are exactly $`\binom{n-1}{m-1}`$ $`m`$-controlled $`R_{\mathbb{C}}`$ gates;

- at the top control number $`m=n-1`$, there are exactly $`n`$ $`(n-1)`$-controlled $`R_{\mathbb{C}}`$ gates and no $`(n-1)`$-controlled $`R_{y}`$ gate.

Therefore the exact complex CNOT count is

<a id="A1.Ex47"></a>

```math
G_{\mathbb{C}}(n)=\sum_{m=0}^{n-1}\left[\binom{n-1}{m}c_{R_{y}}(m)+\binom{n-1}{m-1}c_{R_{\mathbb{C}}}(m)\right]+c_{R_{\mathbb{C}}}(n-1)-c_{R_{y}}(n-1),
```

where binomial coefficients with an invalid lower index are interpreted as zero.

In particular, $`\mathrm{HopfComplex}`$ promotes exactly

<a id="A1.Ex48"></a>

```math
\sum_{m=1}^{n-2}\binom{n-1}{m-1}+n=\sum_{r=0}^{n-3}\binom{n-1}{r}+n=2^{n-1}=N/2
```

gates to $`R_{\mathbb{C}}`$ gates, in agreement with the circuit construction.

Consequently, both ansatze satisfy the asymptotic bounds

<a id="A1.Ex49"></a>

```math
G_{\mathbb{R}}(n)=O(nN),\qquad G_{\mathbb{C}}(n)=O(nN).
```

Moreover, under the no-ancilla counting model above,

<a id="A1.Ex50"></a>

```math
G_{\mathbb{R}}(n)\lt 8nN,\qquad G_{\mathbb{C}}(n)\lt 9nN
```

for all sufficiently large $`n`$.
