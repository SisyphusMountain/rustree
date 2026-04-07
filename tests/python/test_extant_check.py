#!/usr/bin/env python3
"""Check the structure of bd_events."""

import rustree


def test_bd_events_keys():
    """bd_events dict has exactly the expected keys."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    print("Keys in events dict:")
    print(events.keys())
    print("\nFull events dict:")
    for key, value in events.items():
        print(f"{key}: {len(value) if isinstance(value, list) else value}")

    expected_keys = {"time", "node_name", "event_type", "child1_name", "child2_name"}
    assert set(events.keys()) == expected_keys, (
        f"Expected keys {expected_keys}, got {set(events.keys())}"
    )


def test_bd_events_all_lists():
    """All values in the events dict are lists."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    for key, value in events.items():
        assert isinstance(value, list), f"events['{key}'] is not a list, got {type(value)}"


def test_bd_events_lengths_match_num_nodes():
    """All event lists have length equal to the number of nodes in the tree."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    num_nodes = tree.num_nodes()
    for key, value in events.items():
        assert len(value) == num_nodes, (
            f"events['{key}'] has length {len(value)}, expected {num_nodes}"
        )


def test_bd_events_times_positive():
    """All event times are non-negative floats."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    for t in events["time"]:
        assert isinstance(t, float), f"Event time {t!r} is not a float"
        assert t >= 0.0, f"Event time {t} is negative"


def test_bd_events_node_names_match_tree():
    """Every node_name in bd_events corresponds to a real node in the tree."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    # Collect all node names from the tree via iteration
    tree_node_names = {node.name for node in tree}
    event_node_names = set(events["node_name"])

    assert event_node_names == tree_node_names, (
        f"Mismatch between event node names and tree node names.\n"
        f"Only in events: {event_node_names - tree_node_names}\n"
        f"Only in tree: {tree_node_names - event_node_names}"
    )


def test_bd_events_event_types_are_strings():
    """Each event_type entry is a non-empty string."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    for et in events["event_type"]:
        assert isinstance(et, str) and et != "", f"Invalid event_type: {et!r}"
