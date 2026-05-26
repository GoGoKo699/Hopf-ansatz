# Hopf-ansatz
Reference implementation for the Hopf ansatz: a binary-tree variational chart for arbitrary real and complex pure quantum states, with inverse coordinates, diagonal metric, tangent-state synthesis, CNOT counting, and six-qubit VQE/Ramsey examples.

This repository contains the reference code for the **Hopf ansatz**, a
binary-tree variational chart for arbitrary normalized real and complex pure
quantum states.

The code accompanies the paper draft on the Hopf ansatz.  It is intended to
make the definitions, circuit conventions, gate-count estimates, tangent-state
construction, and numerical examples reproducible.

## Contents

| File | Role |
|---|---|
| `hopf_utils.py` | Core implementation of the Hopf ansatz. |
| `hopf_gate_count.py` | CNOT-count checks for the Hopf circuit skeleton. |
| `VQE_qibo.py` | Qibo safeguard for the Hopf state-preparation circuit. |
| `VQE_Layerwise_ADAM_EGTCG.py` | Six-qubit VQE example with layerwise Hopf-gradient estimates. |
| `MET_Layerwise_ADAM_EGTCG.py` | Six-qubit fixed-readout Ramsey metrology example. |

The central file is `hopf_utils.py`.  The other scripts are numerical or
circuit-level checks built around the same conventions.

## Main idea

For \(n\) qubits, the Hopf ansatz uses a complete binary tree with \(2^n\)
leaves.  Each internal node carries one angle that splits probability mass
between the left and right subtrees.  A root-to-leaf path gives one
computational-basis amplitude as a product of sine and cosine factors.

For complex states, the same magnitude tree is supplemented by one phase at
each leaf.

The implementation provides:

- real and complex Hopf coordinate maps;
- inverse maps from amplitudes to Hopf parameters;
- diagonal metric entries;
- Jacobian-related routines;
- normalized tangent-state construction;
- layerwise gradient-access utilities;
- CNOT-count checks;
- VQE and Ramsey numerical demonstrations.

## Usage

Run the CNOT-count check:

```bash
python hopf_gate_count.py
```

Run the Qibo circuit safeguard:

```bash
python VQE_qibo.py
```

Run the VQE demonstration:

```bash
python VQE_Layerwise_ADAM_EGTCG.py
```

Run the fixed-readout Ramsey demonstration:

```bash
python MET_Layerwise_ADAM_EGTCG.py
```

## Requirements

The core code uses standard scientific Python packages:

```bash
pip install numpy scipy matplotlib
```

The Qibo safeguard additionally requires:

```bash
pip install qibo
```

If Qibo is not installed, `VQE_qibo.py` can be skipped.

## What the examples show

The VQE script tests the Hopf gradient-access construction on a six-qubit
transverse-field Ising Hamiltonian.

The Ramsey script tests the same tangent-state machinery on a six-qubit
fixed-readout metrology objective, where the cost depends on measurement
probabilities and local phase slopes rather than on an energy expectation.

Both examples compare exact state-vector gradients with finite-shot layerwise
Hopf-gradient estimates, using Adam and EGT-CG optimization.

## Scope

The Hopf ansatz is a full-state ansatz.  It does not remove the exponential
dimension of arbitrary pure-state optimization, and the layerwise gradient
protocol does not claim logarithmic total shot complexity.

Its purpose is structural: the coordinates are explicit, the inverse map is
available, the metric is diagonal, tangent states are exactly preparable, and
the number of distinct compiled gradient configurations grows only
logarithmically with the magnitude-tree size.
