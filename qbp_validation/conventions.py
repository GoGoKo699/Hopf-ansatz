"""Shared indexing and tree conventions for the exact-logical validation suite.

This module deliberately contains no Qibo imports.  It fixes the manuscript's
big-endian basis labels and the combinatorial data used independently by the
reference and circuit implementations.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterator, Mapping, TypeAlias


@dataclass(frozen=True)
class FrameGateSpec:
    node: int
    depth: int
    position: int
    anchor: int
    marker: int
    target: int  # Qibo/system-list index, with 0 the most-significant qubit.


@dataclass(frozen=True)
class DepthGateSpec:
    node: int
    depth: int
    position: int
    target: int  # Qibo/system-list index.


@dataclass(frozen=True)
class PolyLeaf:
    """A leaf encoded by one computational-basis label."""

    label: int


@dataclass(frozen=True)
class PolyBranch:
    """An ordered full-binary-tree internal node.

    ``key`` identifies the angle in an external mapping.  Keeping the topology
    separate from the numerical angles makes it possible for the circuit and
    reference paths to consume the same immutable test case without sharing
    state-vector code.
    """

    key: str
    left: "PolyTree"
    right: "PolyTree"


PolyTree: TypeAlias = PolyLeaf | PolyBranch


def infer_n_from_magnitudes(theta: object) -> int:
    length = len(theta)  # type: ignore[arg-type]
    n_float = math.log2(length + 1)
    n = int(round(n_float))
    if n < 1 or (1 << n) - 1 != length:
        raise ValueError("Magnitude vector must have length 2**n - 1.")
    return n


def node_depth_position(node: int) -> tuple[int, int]:
    if node < 1:
        raise ValueError("Internal-node indices start at 1.")
    depth = node.bit_length() - 1
    return depth, node - (1 << depth)


def anchor_label(node: int, n: int) -> int:
    depth, position = node_depth_position(node)
    if depth >= n:
        raise ValueError("Node is not internal for the requested qubit count.")
    return position << (n - depth)


def marker_label(node: int, n: int) -> int:
    depth, position = node_depth_position(node)
    if depth >= n:
        raise ValueError("Node is not internal for the requested qubit count.")
    return (2 * position + 1) << (n - depth - 1)


def frame_gate_specs(n: int) -> tuple[FrameGateSpec, ...]:
    specs: list[FrameGateSpec] = []
    for depth in range(n):
        for position in range(1 << depth):
            node = (1 << depth) + position
            specs.append(
                FrameGateSpec(
                    node=node,
                    depth=depth,
                    position=position,
                    anchor=anchor_label(node, n),
                    marker=marker_label(node, n),
                    target=depth,
                )
            )
    return tuple(specs)


def depth_gate_specs(n: int, depth: int | None = None) -> tuple[DepthGateSpec, ...]:
    depths = range(n) if depth is None else (depth,)
    specs: list[DepthGateSpec] = []
    for d in depths:
        if not 0 <= d < n:
            raise ValueError("Depth must lie in 0, ..., n-1.")
        for position in range(1 << d):
            specs.append(
                DepthGateSpec(
                    node=(1 << d) + position,
                    depth=d,
                    position=position,
                    target=d,
                )
            )
    return tuple(specs)


def bit_at(label: int, qubit: int, n: int) -> int:
    """Return the bit on big-endian qubit ``qubit`` (0 is most significant)."""

    if not 0 <= qubit < n:
        raise ValueError("Qubit index out of range.")
    return (label >> (n - 1 - qubit)) & 1


def parity(label: int, outcome: int) -> int:
    """Return the mod-two inner product of two basis labels."""

    return (label & outcome).bit_count() & 1


def marker_map(n: int) -> dict[int, int]:
    return {spec.node: spec.marker for spec in frame_gate_specs(n)}


def controlled_ry_cnot_charge(num_controls: int) -> int:
    if num_controls < 0:
        raise ValueError("Number of controls must be nonnegative.")
    if num_controls == 0:
        return 0
    if num_controls <= 4:
        return (1 << (num_controls + 1)) - 2
    return 16 * (num_controls + 1) - 40


def frame_cnot_charge(n: int) -> int:
    return ((1 << n) - 1) * controlled_ry_cnot_charge(n - 1)


def depth_layer_cnot_charge(depth: int) -> int:
    if depth < 0:
        raise ValueError("Depth must be nonnegative.")
    return (1 << depth) * controlled_ry_cnot_charge(depth)


def depth_preparation_cnot_charge(n: int) -> int:
    return sum(depth_layer_cnot_charge(d) for d in range(n))


def checkpoint_cnot_charge_without_observable(n: int, depth: int) -> int:
    if not 0 <= depth < n:
        raise ValueError("Depth must lie in 0, ..., n-1.")
    suffix = sum(depth_layer_cnot_charge(d) for d in range(depth + 1, n))
    return depth_preparation_cnot_charge(n) + suffix


def poly_preorder(tree: PolyTree) -> Iterator[PolyBranch]:
    if isinstance(tree, PolyLeaf):
        return
    yield tree
    yield from poly_preorder(tree.left)
    yield from poly_preorder(tree.right)


def poly_leaves(tree: PolyTree) -> tuple[int, ...]:
    if isinstance(tree, PolyLeaf):
        return (tree.label,)
    return poly_leaves(tree.left) + poly_leaves(tree.right)


def poly_anchor(tree: PolyTree) -> int:
    if isinstance(tree, PolyLeaf):
        return tree.label
    return poly_anchor(tree.left)


def poly_marker(node: PolyBranch) -> int:
    return poly_anchor(node.right)


def poly_relative_markers(tree: PolyTree) -> dict[str, int]:
    root_anchor = poly_anchor(tree)
    return {node.key: root_anchor ^ poly_marker(node) for node in poly_preorder(tree)}


def validate_poly_tree(tree: PolyTree, n: int, angles: Mapping[str, float]) -> None:
    labels = poly_leaves(tree)
    if len(labels) != len(set(labels)):
        raise ValueError("Polyspherical leaf labels must be distinct.")
    if any(label < 0 or label >= (1 << n) for label in labels):
        raise ValueError("Polyspherical leaf label outside the encoded space.")
    keys = [node.key for node in poly_preorder(tree)]
    if len(keys) != len(set(keys)):
        raise ValueError("Polyspherical internal-node keys must be unique.")
    missing = set(keys).difference(angles)
    if missing:
        raise ValueError(f"Missing angles for internal nodes: {sorted(missing)}")
