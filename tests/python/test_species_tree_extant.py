#!/usr/bin/env python3
"""Test that simulate_species_tree with extinction produces a valid tree.

With mu > 0, the tree includes both extant and extinct species as leaves.
The parameter n controls the number of extant species at the end of the simulation,
but extinct lineages are also retained in the tree.
"""

import rustree


def test_species_tree_structure_with_extinction():
    """Tree with extinction has at least n leaves and valid binary structure."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.5, seed=42)

    # With extinction, total leaves >= n (extinct lineages are also leaves)
    assert tree.num_leaves() >= n, (
        f"Expected at least {n} leaves, got {tree.num_leaves()}"
    )

    # Binary tree: internal nodes = leaves - 1, total = 2*leaves - 1
    assert tree.num_nodes() == 2 * tree.num_leaves() - 1, (
        f"Binary tree should have 2n-1 nodes, got {tree.num_nodes()}"
    )


def test_species_tree_no_extinction():
    """Without extinction, the tree should have exactly n leaves."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.0, seed=42)

    assert tree.num_leaves() == n, (
        f"Expected exactly {n} leaves with no extinction, got {tree.num_leaves()}"
    )
    assert tree.num_nodes() == 2 * n - 1


def test_leaves_have_no_children():
    """Every structural leaf (node with no children) has both child indices as None."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    num_leaves_found = 0
    for node in tree:
        if node.left_child is None and node.right_child is None:
            num_leaves_found += 1
        else:
            # A non-leaf must have exactly two children (binary tree)
            assert node.left_child is not None and node.right_child is not None, (
                f"Node {node.index} has only one child: "
                f"left={node.left_child}, right={node.right_child}"
            )

    assert num_leaves_found == tree.num_leaves(), (
        f"Counted {num_leaves_found} leaves by iteration, "
        f"but num_leaves() returns {tree.num_leaves()}"
    )


def test_internal_nodes_have_exactly_two_children():
    """Every non-leaf node has exactly two children (strictly binary tree)."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    for node in tree:
        is_leaf = node.left_child is None and node.right_child is None
        if not is_leaf:
            assert node.left_child is not None, (
                f"Internal node {node.index} has no left child"
            )
            assert node.right_child is not None, (
                f"Internal node {node.index} has no right child"
            )


def test_root_has_no_parent():
    """The root node has no parent."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    root = tree.get_node(tree.root_index())
    assert root.parent is None, (
        f"Root node {root.index} should have no parent, got {root.parent}"
    )


def test_all_non_root_nodes_have_parent():
    """Every node except the root has a parent."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    root_idx = tree.root_index()
    for node in tree:
        if node.index == root_idx:
            assert node.parent is None, f"Root should have no parent"
        else:
            assert node.parent is not None, (
                f"Non-root node {node.index} has no parent"
            )


def test_parent_child_consistency():
    """For every non-root node, its parent's child list contains this node's index."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    root_idx = tree.root_index()
    for node in tree:
        if node.index == root_idx:
            continue
        parent = tree.get_node(node.parent)
        children = {parent.left_child, parent.right_child}
        assert node.index in children, (
            f"Node {node.index} claims parent {node.parent}, "
            f"but parent's children are {children}"
        )


def test_tree_height_positive():
    """The tree height is a positive float."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    height = tree.tree_height()
    assert isinstance(height, float), f"tree_height() returned {type(height)}, expected float"
    assert height > 0.0, f"Expected positive tree height, got {height}"


def test_leaf_names_count_matches_num_leaves():
    """leaf_names() returns exactly as many names as num_leaves()."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    leaf_names = tree.leaf_names()
    assert len(leaf_names) == tree.num_leaves(), (
        f"leaf_names() has {len(leaf_names)} entries but num_leaves() = {tree.num_leaves()}"
    )
    # Names must be unique
    assert len(set(leaf_names)) == len(leaf_names), "Duplicate leaf names found"
