#!/usr/bin/env python3
"""Test to understand extant vs extinct species in simulate_species_tree."""

import rustree


def test_extant_count_equals_n():
    """The number of Leaf events (extant species) exactly equals the requested n."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()
    event_types = events["event_type"]

    n_leaves = event_types.count("Leaf")
    print(f"Requested: n={n} extant species")
    print(f"Leaf events (extant): {n_leaves}")

    assert n_leaves == n, (
        f"Expected exactly {n} extant species (Leaf events), got {n_leaves}"
    )


def test_tree_includes_extinct_lineages():
    """With mu > 0, total leaves include both extant and extinct species."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()
    event_types = events["event_type"]

    n_leaves = event_types.count("Leaf")
    n_extinctions = event_types.count("Extinction")

    print(f"Tree has: {tree.num_leaves()} total leaves")
    print(f"Total nodes: {tree.num_nodes()} nodes\n")
    print(f"Birth-death events:")
    print(f"  Leaves (extant): {n_leaves} events")
    print(f"  Extinctions: {n_extinctions} events")
    print(f"\nThe tree includes BOTH extant AND extinct species as leaves!")
    print(f"Extant leaves: {n_leaves}")
    print(f"Extinct leaves: {n_extinctions}")
    print(f"Total leaves: {tree.num_leaves()}")
    print(f"Verification: {n_leaves + n_extinctions} = {tree.num_leaves()}")

    # With mu=0.8 and seed=42 there are extinctions
    assert n_extinctions > 0, (
        f"Expected extinctions with mu=0.8, got {n_extinctions}"
    )
    # Total structural leaves = extant + extinct
    assert n_leaves + n_extinctions == tree.num_leaves(), (
        f"Leaf ({n_leaves}) + Extinction ({n_extinctions}) events "
        f"should equal num_leaves ({tree.num_leaves()})"
    )
    # Total leaves must be at least n
    assert tree.num_leaves() >= n, (
        f"Total leaves {tree.num_leaves()} should be >= n={n}"
    )


def test_no_extinction_means_only_extant_leaves():
    """With mu=0, every leaf is extant: Leaf events == num_leaves and zero Extinctions."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.0, seed=42)
    events = tree.get_bd_events()
    event_types = events["event_type"]

    n_leaves_events = event_types.count("Leaf")
    n_extinctions = event_types.count("Extinction")

    assert n_extinctions == 0, f"With mu=0, expected 0 extinctions, got {n_extinctions}"
    assert n_leaves_events == n, (
        f"With mu=0, expected {n} Leaf events, got {n_leaves_events}"
    )
    assert tree.num_leaves() == n, (
        f"With mu=0, expected exactly {n} structural leaves, got {tree.num_leaves()}"
    )


def test_speciations_equal_leaves_minus_one():
    """Number of Speciation events equals (total leaves - 1) for a binary tree."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()
    event_types = events["event_type"]

    n_speciations = event_types.count("Speciation")
    total_leaves = tree.num_leaves()

    assert n_speciations == total_leaves - 1, (
        f"Expected {total_leaves - 1} Speciation events, got {n_speciations}"
    )
