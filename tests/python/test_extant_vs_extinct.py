#!/usr/bin/env python3
"""Test to understand extant vs extinct species in simulate_species_tree."""

import rustree

# Simulate with HIGH extinction rate to see the effect
tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.8, seed=42)

print(f"Requested: n=20 extant species")
print(f"Tree has: {tree.num_leaves()} total leaves")
print(f"Total nodes: {tree.num_nodes()} nodes\n")

# Get BD events to distinguish extant from extinct
events = tree.get_bd_events()
event_types = events['event_type']
n_leaves = event_types.count('Leaf')
n_extinctions = event_types.count('Extinction')
n_speciations = event_types.count('Speciation')

print(f"Birth-death events:")
print(f"  Leaves (extant): {n_leaves} events")
print(f"  Extinctions: {n_extinctions} events")
print(f"  Speciations: {n_speciations} events")

print(f"\nThe tree includes BOTH extant AND extinct species as leaves!")
print(f"Extant leaves: {n_leaves}")
print(f"Extinct leaves: {n_extinctions}")
print(f"Total leaves: {tree.num_leaves()}")
print(f"Verification: {n_leaves + n_extinctions} = {tree.num_leaves()}")
