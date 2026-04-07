#!/usr/bin/env python3
"""Check event types in the BD tree."""

import rustree
from collections import Counter


def test_event_types_keys():
    """get_bd_events() returns a dict with exactly the expected keys."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    assert set(events.keys()) == {"time", "node_name", "event_type", "child1_name", "child2_name"}, (
        f"Unexpected keys: {set(events.keys())}"
    )


def test_event_types_all_lists_same_length():
    """All lists in the events dict have the same length (one entry per node)."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    lengths = {k: len(v) for k, v in events.items()}
    assert len(set(lengths.values())) == 1, f"Inconsistent lengths: {lengths}"
    assert list(lengths.values())[0] == tree.num_nodes(), (
        f"Event list length {list(lengths.values())[0]} != num_nodes {tree.num_nodes()}"
    )


def test_event_types_only_known_types():
    """event_type list contains only the three known event types."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    known = {"Leaf", "Extinction", "Speciation"}
    actual = set(events["event_type"])
    assert actual.issubset(known), f"Unknown event types: {actual - known}"


def test_event_types_counts():
    """Event counts satisfy the structural constraints of a binary tree."""
    n = 20
    tree = rustree.simulate_species_tree(n=n, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    event_counts = Counter(events["event_type"])
    n_leaves = event_counts["Leaf"]
    n_extinctions = event_counts["Extinction"]
    n_speciations = event_counts["Speciation"]

    print(f"Requested: n={n} extant species")
    print(f"Tree has: {tree.num_leaves()} total leaves\n")
    print("Event type counts:")
    for event_type, count in sorted(event_counts.items()):
        print(f"  {event_type}: {count}")
    print(f"\nEvent breakdown:")
    print(f"  Leaf (extant): {n_leaves}")
    print(f"  Extinction: {n_extinctions}")
    print(f"  Speciation: {n_speciations}")
    print(f"  Total: {n_leaves + n_extinctions + n_speciations}")

    # n controls the number of extant species
    assert n_leaves == n, f"Expected {n} Leaf events (extant species), got {n_leaves}"

    # Total leaves (extant + extinct) must equal structural leaf count
    assert n_leaves + n_extinctions == tree.num_leaves(), (
        f"Leaf + Extinction events ({n_leaves + n_extinctions}) "
        f"!= structural leaves ({tree.num_leaves()})"
    )

    # Binary tree: number of speciations == number of leaves - 1
    total_leaves = tree.num_leaves()
    assert n_speciations == total_leaves - 1, (
        f"Expected {total_leaves - 1} Speciation events, got {n_speciations}"
    )

    # Total events == total nodes
    assert n_leaves + n_extinctions + n_speciations == tree.num_nodes(), (
        f"Total events {n_leaves + n_extinctions + n_speciations} != num_nodes {tree.num_nodes()}"
    )


def test_event_types_leaf_names():
    """leaf_names() returns all structural leaves (extant + extinct)."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    leaf_names = tree.leaf_names()

    print(f"\nTotal leaf names: {len(leaf_names)}")
    print(f"First 10 leaf names: {leaf_names[:10]}")

    assert len(leaf_names) == tree.num_leaves(), (
        f"leaf_names length {len(leaf_names)} != num_leaves {tree.num_leaves()}"
    )
    # All names should be non-empty strings
    assert all(isinstance(name, str) and name != "" for name in leaf_names), (
        "Some leaf names are empty or not strings"
    )
    # All names should be unique
    assert len(set(leaf_names)) == len(leaf_names), "Leaf names are not all unique"


def test_speciation_nodes_have_children():
    """Speciation events must reference two non-empty child names."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    for et, c1, c2 in zip(events["event_type"], events["child1_name"], events["child2_name"]):
        if et == "Speciation":
            assert c1 != "" and c2 != "", (
                f"Speciation node has empty child names: child1='{c1}', child2='{c2}'"
            )
        else:
            assert c1 == "" and c2 == "", (
                f"Non-speciation event '{et}' has non-empty children: '{c1}', '{c2}'"
            )
