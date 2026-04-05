#!/usr/bin/env python3
"""Compatibility check for species tree simulation API."""

import rustree


def test_simulate_species_tree_returns_requested_extant_count():
    """simulate_species_tree(n, ...) should expose exactly n extant leaves."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)
    assert tree.num_leaves() == 20
