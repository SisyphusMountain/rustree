#!/usr/bin/env python3
"""Compatibility check for birth-death events dictionary format."""

import rustree


def test_get_bd_events_has_expected_columns():
    """get_bd_events() should expose tabular columns used by the main test suite."""
    tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)
    events = tree.get_bd_events()

    expected_keys = {"time", "node_name", "event_type", "child1_name", "child2_name"}
    assert expected_keys.issubset(events.keys())
    assert len(events["time"]) == len(events["event_type"]) > 0
