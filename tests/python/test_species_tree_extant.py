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
