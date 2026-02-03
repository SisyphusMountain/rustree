# rustree Python Tutorial

This tutorial explains how to use the rustree Python bindings for simulating species trees, gene trees with DTL events, and visualizing reconciled trees.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Building the Library](#building-the-library)
- [Quick Start](#quick-start)
- [Complete Example](#complete-example)
- [Full Workflow Example](#full-workflow-example)
- [Advanced Features](#advanced-features)
  - [Birth-Death Event Analysis](#birth-death-event-analysis-)
  - [Pairwise Distance Analysis](#pairwise-distance-analysis-)
  - [Subtree Extraction by Names](#subtree-extraction-by-names)
  - [Gene Tree Distance Analysis](#gene-tree-distance-analysis-)
- [Complete Workflows](#complete-workflows)
  - [Workflow 1: Simulate and Analyze BD Events](#workflow-1-simulate-tree-and-analyze-birth-death-events)
  - [Workflow 2: Distance-Based Phylogenetic Analysis](#workflow-2-compute-distances-and-phylogenetic-analysis)
  - [Workflow 3: Subtree Sampling](#workflow-3-sample-tree-for-focused-analysis)
  - [Workflow 4: Combined Pipeline](#workflow-4-combined-analysis-pipeline)
- [API Reference](#api-reference)
- [Parameters](#parameters)
- [Notes](#notes)
- [Implementation Guide](#implementation-guide-for-missing-features)

## Prerequisites

- **Python** (>= 3.8)
- **Rust** (>= 1.70)
- **maturin** (`pip install maturin`)

For SVG visualization, install thirdkind:
```bash
cargo install thirdkind
```

## Building the Library

```bash
cd /path/to/rustree

# Build and install with maturin
maturin develop --release

# Or build wheel for distribution
maturin build --release
```

## Quick Start

```python
import rustree

# Simulate a species tree with 50 extant species
sp_tree = rustree.simulate_species_tree(50, 1.0, 0.3, seed=42)

# Check the tree
print(f"Number of leaves: {sp_tree.num_leaves()}")
print(f"Newick: {sp_tree.to_newick()[:100]}...")

# Simulate gene trees with DTL events
gene_trees = sp_tree.simulate_dtl_batch(10, 0.5, 0.2, 0.3, seed=123)

# Work with the first gene tree
gt = gene_trees[0]
print(f"Number of extant genes: {gt.num_extant()}")
```

## Complete Example

### 1. Simulate a Species Tree

```python
import rustree

# Birth-death species tree
# n: number of extant species
# lambda: speciation (birth) rate
# mu: extinction (death) rate (must be < lambda)
# seed: random seed for reproducibility (optional)

sp_tree = rustree.simulate_species_tree(
    n=100,
    lambda_=1.0,
    mu=0.5,
    seed=42
)

# Inspect the tree
print(sp_tree.num_leaves())       # 100
print(sp_tree.leaf_names())       # List of leaf names
print(sp_tree.to_newick())        # Newick string
print(sp_tree.tree_height())      # Total height
```

### 2. Parse an Existing Newick Tree

```python
# Parse from string
newick_str = "((A:1,B:1):1,(C:1,D:1):1):0;"
sp_tree = rustree.parse_species_tree(newick_str)

# Or read from file
with open("species_tree.nwk") as f:
    sp_tree = rustree.parse_species_tree(f.read())
```

### 3. Simulate Gene Trees with DTL Events

```python
# Single gene tree (uniform random transfers)
gene_tree = sp_tree.simulate_dtl(
    lambda_d=0.5,    # Duplication rate
    lambda_t=0.2,    # Transfer rate
    lambda_l=0.3,    # Loss rate
    seed=123
)

# Batch of gene trees (more efficient)
gene_trees = sp_tree.simulate_dtl_batch(
    n=100,           # Number of gene trees
    lambda_d=0.5,
    lambda_t=0.2,
    lambda_l=0.3,
    seed=123
)

# Access individual trees
gt1 = gene_trees[0]
gt2 = gene_trees[1]
```

### 4. Assortative (Distance-Dependent) Transfers

By default, transfer recipients are chosen uniformly at random from contemporary species. You can enable **assortative transfers** where closer species are more likely to receive transfers:

```python
# Assortative transfers: P(recipient) ∝ exp(-alpha * distance)
# Higher alpha = more local transfers (closer species preferred)
# alpha = 0 is equivalent to uniform random

# Single gene tree with assortative transfers
gene_tree = sp_tree.simulate_dtl(
    lambda_d=0.5,
    lambda_t=0.5,
    lambda_l=0.3,
    transfer_alpha=1.0,  # Distance decay parameter
    seed=123
)

# Batch with assortative transfers
gene_trees = sp_tree.simulate_dtl_batch(
    n=100,
    lambda_d=0.5,
    lambda_t=0.5,
    lambda_l=0.3,
    transfer_alpha=2.0,  # Stronger distance preference
    seed=123
)
```

The distance between species A and B at time t is computed as:
```
d(A, B, t) = 2 × (t - depth_of_LCA(A, B))
```

### 5. Inspect Gene Trees

```python
gt = gene_trees[0]

# Number of extant genes (surviving genes on extant species)
# Only counts gene leaves mapped to extant species (not extinctions)
print(gt.num_extant())

# Convert to Newick
print(gt.to_newick())

# Count events by type: (speciations, duplications, transfers, losses, leaves)
s, d, t, l, leaves = gt.count_events()
print(f"S={s}, D={d}, T={t}, L={l}, Leaves={leaves}")

# Get extant gene names
print(gt.extant_gene_names())

# Export to pandas DataFrame
df = gt.to_csv()  # Returns DataFrame
print(df.head())
```

### 6. Sample Extant Genes

```python
# Extract induced subtree containing only extant genes
# (removes internal nodes leading only to losses)
sampled_gt = gt.sample_extant()

print(sampled_gt.to_newick())
print(f"Sampled tree has {sampled_gt.num_nodes()} nodes")

# Sample by specific gene names
sampled_gt = gt.sample_by_names(["G1_A", "G2_B", "G3_C"])
```

### 7. Export Results

```python
# Save species tree to Newick
sp_tree.save_newick("species_tree.nwk")

# Save gene tree to Newick
gt.save_newick("gene_tree.nwk")

# Save gene tree to RecPhyloXML (for visualization)
gt.save_xml("gene_tree.recphyloxml")

# Save gene tree to CSV
df = gt.to_csv("gene_tree.csv")  # Also saves to file
```

### 8. Visualize with thirdkind

```python
# Generate SVG (requires thirdkind to be installed)
# Returns SVG as string
svg_content = gt.to_svg()

# Save to file
gt.to_svg("gene_tree.svg")

# Save and open in browser
gt.to_svg("gene_tree.svg", open_browser=True)

# Display in Jupyter notebook
gt.display()  # Shows SVG inline
```

## Full Workflow Example

```python
import rustree

# 1. Simulate species tree
print("Simulating species tree...")
sp_tree = rustree.simulate_species_tree(50, 1.0, 0.3, seed=42)
print(f"Species tree has {sp_tree.num_leaves()} leaves")

# 2. Simulate 100 gene families with assortative transfers
print("Simulating 100 gene families...")
gene_trees = sp_tree.simulate_dtl_batch(
    100, 0.5, 0.2, 0.3,
    transfer_alpha=1.0,  # Distance-dependent transfers
    seed=123
)

# 3. Analyze results
extant_counts = [gt.num_extant() for gt in gene_trees]
print(f"Extant genes per family:")
print(f"  Min: {min(extant_counts)}")
print(f"  Max: {max(extant_counts)}")
print(f"  Mean: {sum(extant_counts)/len(extant_counts):.1f}")

# 4. Save results
sp_tree.save_newick("species.nwk")

for i, gt in enumerate(gene_trees):
    gt.save_newick(f"gene_family_{i:03d}.nwk")

# 5. Visualize first gene family
gene_trees[0].save_xml("gene_family_000.xml")
gene_trees[0].to_svg("gene_family_000.svg")

print("Done! Results saved.")
```

## Advanced Features

> **Implementation Status:** The features documented below are part of the rustree API. Features marked with ⚠️ are currently available in the R bindings but need to be added to the Python bindings. If a feature is not yet available, you'll see an `AttributeError`. To implement these features, add the corresponding methods to `PySpeciesTree` and `PyGeneTree` in `/home/enzo/Documents/git/WP2/rustree/src/python.rs` following the patterns in `/home/enzo/Documents/git/WP2/rustree/src/r.rs`.

### Birth-Death Event Analysis ⚠️

Birth-death simulations track all speciation and extinction events. You can analyze and export these events:

```python
import rustree

# Simulate a species tree with birth-death process
sp_tree = rustree.simulate_species_tree(20, 1.0, 0.3, seed=42)

# Get birth-death events as a list
# Returns list of TreeEvent objects with: time, node_id, event_type, child1, child2
events = sp_tree.get_bd_events()

# Print events chronologically
for event in events:
    print(f"Time {event.time:.2f}: {event.event_type} at node {event.node_id}")

# Save events to CSV file
sp_tree.save_bd_events_csv("bd_events.csv")
# CSV columns: time, node_name, event_type, child1_name, child2_name
```

The CSV file contains:
- `time`: When the event occurred (going backwards from present at 0)
- `node_name`: Name of the node where event occurred
- `event_type`: "Speciation", "Extinction", or "Leaf"
- `child1_name`, `child2_name`: Names of child nodes (for speciation events)

### Pairwise Distance Analysis ⚠️

Calculate distances between nodes in a tree for phylogenetic analysis:

```python
import rustree

# Parse or simulate a tree
sp_tree = rustree.simulate_species_tree(10, 1.0, 0.3, seed=42)

# Get pairwise distances as a list of PairwiseDistance objects
# distance_type can be "topological" or "metric"
# leaves_only=True restricts to leaf-to-leaf distances

# Topological distance (number of edges between nodes)
topo_distances = sp_tree.pairwise_distances("topological", leaves_only=True)

for dist in topo_distances[:5]:  # Show first 5
    print(f"{dist.node1} <-> {dist.node2}: {dist.distance} edges")

# Metric distance (sum of branch lengths)
metric_distances = sp_tree.pairwise_distances("metric", leaves_only=True)

for dist in metric_distances[:5]:
    print(f"{dist.node1} <-> {dist.node2}: {dist.distance:.4f} units")

# Save to CSV file
sp_tree.save_pairwise_distances_csv(
    "leaf_distances.csv",
    distance_type="metric",
    leaves_only=True
)
# CSV columns: node1, node2, distance

# Include all nodes (internal + leaves)
sp_tree.save_pairwise_distances_csv(
    "all_distances.csv",
    distance_type="topological",
    leaves_only=False
)
```

**Distance Types:**
- `"topological"`: Counts edges between nodes (integer values)
- `"metric"`: Sums branch lengths along the path (real values)

### Subtree Extraction by Names

Extract induced subtrees containing only specified taxa:

```python
import rustree

# Simulate or load a large tree
sp_tree = rustree.simulate_species_tree(50, 1.0, 0.3, seed=42)
print(f"Original tree: {sp_tree.num_leaves()} leaves")

# Get all leaf names
all_leaves = sp_tree.leaf_names()

# Select a subset of species to keep
selected_species = all_leaves[:10]  # First 10 species

# Extract induced subtree
subtree = sp_tree.extract_induced_subtree_by_names(selected_species)
print(f"Subtree: {subtree.num_leaves()} leaves")
print(f"Subtree height: {subtree.tree_height():.2f}")

# Save the subtree
subtree.save_newick("subtree.nwk")

# This also works with gene trees
gt = sp_tree.simulate_dtl(0.5, 0.2, 0.3, seed=123)
gene_names = gt.extant_gene_names()
sampled_genes = gene_names[:20]  # Keep first 20 genes
sampled_gt = gt.sample_by_names(sampled_genes)  # Different method for gene trees

# You can also sample species and automatically filter gene tree
species_subset = sp_tree.leaf_names()[:5]
sampled_gt2 = gt.sample_species_leaves(species_subset)
print(f"Gene tree after species sampling: {sampled_gt2.num_extant()} extant genes")
```

### Gene Tree Distance Analysis ⚠️

Gene trees also support distance computations, useful for analyzing gene family evolution:

```python
import rustree

# Simulate species tree and gene tree
sp_tree = rustree.simulate_species_tree(10, 1.0, 0.3, seed=42)
gt = sp_tree.simulate_dtl(0.5, 0.2, 0.3, seed=123)

print(f"Gene tree has {gt.num_extant()} extant genes")

# Compute pairwise distances between genes
gene_distances = gt.pairwise_distances("metric", leaves_only=True)
print(f"Computed {len(gene_distances)} pairwise gene distances")

# Save to CSV for analysis
gt.save_pairwise_distances_csv(
    "gene_distances.csv",
    distance_type="metric",
    leaves_only=True
)

# Compare gene tree distances to species tree distances
# This can reveal signatures of duplication, transfer, and loss events
import pandas as pd

gene_df = pd.read_csv("gene_distances.csv")
print(f"Mean gene distance: {gene_df['distance'].mean():.4f}")
print(f"Max gene distance: {gene_df['distance'].max():.4f}")
```

## Complete Workflows

### Workflow 1: Simulate Tree and Analyze Birth-Death Events

```python
import rustree
import pandas as pd

# 1. Simulate species tree
print("Simulating species tree...")
sp_tree = rustree.simulate_species_tree(30, lambda_=1.2, mu=0.4, seed=42)

# 2. Analyze the tree
print(f"Extant species: {sp_tree.num_leaves()}")
print(f"Total nodes: {sp_tree.num_nodes()}")
print(f"Tree height: {sp_tree.tree_height():.2f}")

# 3. Extract and analyze birth-death events
events = sp_tree.get_bd_events()
speciations = [e for e in events if e.event_type == "Speciation"]
extinctions = [e for e in events if e.event_type == "Extinction"]

print(f"\nEvent summary:")
print(f"  Speciations: {len(speciations)}")
print(f"  Extinctions: {len(extinctions)}")
print(f"  Extant leaves: {sp_tree.num_leaves()}")

# 4. Save events for further analysis
sp_tree.save_bd_events_csv("species_tree_events.csv")
df = pd.read_csv("species_tree_events.csv")
print(f"\nEvents saved: {len(df)} records")

# 5. Visualize event times
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.hist([e.time for e in events], bins=20, alpha=0.7)
plt.xlabel("Time (backwards from present)")
plt.ylabel("Number of events")
plt.title("Distribution of Birth-Death Events")
plt.savefig("event_times.png")
```

### Workflow 2: Compute Distances and Phylogenetic Analysis

```python
import rustree
import pandas as pd
import numpy as np

# 1. Parse or simulate a tree
sp_tree = rustree.simulate_species_tree(15, 1.0, 0.3, seed=42)

# 2. Compute pairwise distances
print("Computing pairwise distances...")

# Metric distances (branch lengths)
metric_dists = sp_tree.pairwise_distances("metric", leaves_only=True)
sp_tree.save_pairwise_distances_csv(
    "metric_distances.csv",
    "metric",
    leaves_only=True
)

# Topological distances (edge counts)
topo_dists = sp_tree.pairwise_distances("topological", leaves_only=True)
sp_tree.save_pairwise_distances_csv(
    "topo_distances.csv",
    "topological",
    leaves_only=True
)

# 3. Load and analyze distance matrix
df_metric = pd.read_csv("metric_distances.csv")
df_topo = pd.read_csv("topo_distances.csv")

print(f"Number of pairwise comparisons: {len(df_metric)}")
print(f"Mean metric distance: {df_metric['distance'].mean():.4f}")
print(f"Max metric distance: {df_metric['distance'].max():.4f}")
print(f"Mean topological distance: {df_topo['distance'].mean():.2f}")

# 4. Convert to distance matrix for clustering/visualization
leaf_names = sp_tree.leaf_names()
n = len(leaf_names)
dist_matrix = np.zeros((n, n))

name_to_idx = {name: i for i, name in enumerate(leaf_names)}
for _, row in df_metric.iterrows():
    i = name_to_idx[row['node1']]
    j = name_to_idx[row['node2']]
    dist_matrix[i, j] = row['distance']

# 5. Perform hierarchical clustering
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

# Create linkage for dendrogram
condensed_dist = dist_matrix[np.triu_indices(n, k=1)]
linkage_matrix = linkage(condensed_dist, method='average')

plt.figure(figsize=(12, 8))
dendrogram(linkage_matrix, labels=leaf_names, leaf_rotation=90)
plt.title("Hierarchical Clustering from Pairwise Distances")
plt.xlabel("Species")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig("clustering.png")
```

### Workflow 3: Sample Tree for Focused Analysis

```python
import rustree

# 1. Simulate large species tree
print("Simulating large species tree...")
sp_tree = rustree.simulate_species_tree(100, 1.0, 0.3, seed=42)
print(f"Full tree: {sp_tree.num_leaves()} species")

# 2. Get all species names
all_species = sp_tree.leaf_names()

# 3. Select species of interest (e.g., every 5th species)
selected_species = [all_species[i] for i in range(0, len(all_species), 5)]
print(f"Selected {len(selected_species)} species for analysis")

# 4. Extract subtree
print("Extracting subtree...")
subtree = sp_tree.extract_induced_subtree_by_names(selected_species)
print(f"Subtree: {subtree.num_leaves()} leaves, height: {subtree.tree_height():.2f}")

# 5. Simulate gene families on the subtree (faster than full tree)
print("Simulating gene families on subtree...")
gene_trees = subtree.simulate_dtl_batch(
    n=50,
    lambda_d=0.5,
    lambda_t=0.2,
    lambda_l=0.3,
    transfer_alpha=1.0,
    seed=123
)

# 6. Analyze gene families
extant_counts = [gt.num_extant() for gt in gene_trees]
print(f"\nGene family statistics:")
print(f"  Mean extant genes: {sum(extant_counts)/len(extant_counts):.1f}")
print(f"  Min: {min(extant_counts)}, Max: {max(extant_counts)}")

# 7. Count events across all families
total_events = {'S': 0, 'D': 0, 'T': 0, 'L': 0}
for gt in gene_trees:
    s, d, t, l, leaves = gt.count_events()
    total_events['S'] += s
    total_events['D'] += d
    total_events['T'] += t
    total_events['L'] += l

print(f"\nTotal events across {len(gene_trees)} families:")
for event_type, count in total_events.items():
    print(f"  {event_type}: {count}")

# 8. Export selected gene family
print("\nExporting first gene family...")
gene_trees[0].save_newick("gene_family_0.nwk")
gene_trees[0].save_xml("gene_family_0.xml")
gene_trees[0].to_csv("gene_family_0.csv")

# 9. Compute distances on gene tree
gene_dists = gene_trees[0].pairwise_distances("metric", leaves_only=True)
print(f"Computed {len(gene_dists)} pairwise distances for gene tree")
```

### Workflow 4: Combined Analysis Pipeline

```python
import rustree
import pandas as pd

# 1. Simulate species tree and save events
sp_tree = rustree.simulate_species_tree(50, 1.0, 0.5, seed=42)
sp_tree.save_newick("species.nwk")
sp_tree.save_bd_events_csv("species_events.csv")

# 2. Analyze birth-death process
events_df = pd.read_csv("species_events.csv")
print("Birth-death process summary:")
print(events_df['event_type'].value_counts())

# 3. Compute and save species distances
sp_tree.save_pairwise_distances_csv(
    "species_distances.csv",
    distance_type="metric",
    leaves_only=True
)

# 4. Simulate gene families
print("\nSimulating 100 gene families...")
gene_trees = sp_tree.simulate_dtl_batch(
    100, 0.5, 0.2, 0.3,
    transfer_alpha=1.0,
    seed=123
)

# 5. Export gene tree data
for i, gt in enumerate(gene_trees[:10]):  # Export first 10
    gt.save_newick(f"output/gene_family_{i:03d}.nwk")
    gt.to_csv(f"output/gene_family_{i:03d}.csv")

# 6. Sample a subset of species
species_subset = sp_tree.leaf_names()[:20]
subtree = sp_tree.extract_induced_subtree_by_names(species_subset)

# 7. Recompute distances on subtree
subtree.save_pairwise_distances_csv(
    "subtree_distances.csv",
    distance_type="metric",
    leaves_only=True
)

print("\nPipeline complete! Results saved to output directory.")
```

## API Reference

### Module Functions

| Function | Description |
|----------|-------------|
| `simulate_species_tree(n, lambda_, mu, seed=None)` | Simulate birth-death tree |
| `parse_species_tree(newick_str)` | Parse Newick string |
| `parse_recphyloxml(filepath)` | Parse RecPhyloXML file |

### PySpeciesTree Methods

| Method | Status | Description |
|--------|--------|-------------|
| `to_newick()` | ✅ | Convert to Newick string |
| `save_newick(filepath)` | ✅ | Save to Newick file |
| `num_nodes()` | ✅ | Count all nodes |
| `num_leaves()` | ✅ | Count leaves |
| `tree_height()` | ✅ | Get total height |
| `leaf_names()` | ✅ | Get list of leaf names |
| `root_index()` | ✅ | Get root node index |
| `extract_induced_subtree_by_names(names)` | ✅ | Extract subtree with specified leaves |
| `simulate_dtl(...)` | ✅ | Simulate single gene tree |
| `simulate_dtl_batch(...)` | ✅ | Simulate batch of gene trees |
| `get_bd_events()` | ⚠️ | Get birth-death events as list of TreeEvent |
| `save_bd_events_csv(filepath)` | ⚠️ | Save BD events to CSV file |
| `pairwise_distances(distance_type, leaves_only)` | ⚠️ | Compute pairwise distances |
| `save_pairwise_distances_csv(filepath, distance_type, leaves_only)` | ⚠️ | Save distances to CSV |

**Status Legend:** ✅ Available | ⚠️ Needs implementation (see R bindings for reference)

### PyGeneTree Methods

| Method | Status | Description |
|--------|--------|-------------|
| `to_newick()` | ✅ | Convert to Newick string |
| `save_newick(filepath)` | ✅ | Save to Newick file |
| `num_nodes()` | ✅ | Count all nodes |
| `num_extant()` | ✅ | Count extant genes (on extant species only) |
| `count_events()` | ✅ | Return (S, D, T, L, Leaves) counts |
| `extant_gene_names()` | ✅ | Get names of extant genes |
| `sample_extant()` | ✅ | Extract extant-only subtree |
| `sample_by_names(names)` | ✅ | Extract subtree by gene names |
| `sample_species_leaves(species_leaf_names)` | ✅ | Sample species and filter gene tree |
| `to_xml()` | ✅ | Convert to RecPhyloXML string |
| `save_xml(filepath)` | ✅ | Save to RecPhyloXML file |
| `to_csv(filepath=None)` | ✅ | Export to DataFrame (optionally save) |
| `to_svg(filepath=None, open_browser=False)` | ✅ | Generate SVG visualization |
| `display()` | ✅ | Display in Jupyter notebook |
| `pairwise_distances(distance_type, leaves_only)` | ⚠️ | Compute pairwise distances |
| `save_pairwise_distances_csv(filepath, distance_type, leaves_only)` | ⚠️ | Save distances to CSV |

## Parameters

### DTL Rates
- `lambda_d`: Duplication rate per unit time along branches
- `lambda_t`: Transfer rate per unit time along branches
- `lambda_l`: Loss rate per unit time along branches

### Assortative Transfers
- `transfer_alpha`: Distance decay parameter for transfer recipient selection
  - `None` (default): Uniform random selection among contemporary species
  - `0`: Equivalent to uniform random
  - `> 0`: Closer species are more likely to receive transfers
  - Higher values = stronger preference for nearby species

### Require Extant
- `require_extant`: Ensure at least one gene survives on an extant species
  - `False` (default): Return tree even if all genes went extinct
  - `True`: Retry simulation until at least one gene survives on an extant species
  - Useful when loss rates are high and you need functional gene families

### Distance Types
- `distance_type`: Type of distance to compute between nodes
  - `"topological"` or `"topo"`: Count of edges between nodes (integer)
  - `"metric"`, `"branch"`, or `"patristic"`: Sum of branch lengths (real)

### Distance Scope
- `leaves_only`: Whether to restrict distance computation
  - `True`: Only compute distances between leaf nodes (extant species/genes)
  - `False`: Compute distances between all nodes (including internal nodes)

## Notes

- Use `seed` parameter for reproducible simulations
- DTL rates are per unit time along branches
- thirdkind must be installed for SVG visualization
- For Jupyter notebook display, IPython must be available

### Extant Genes Definition

The `num_extant()` method counts gene leaves that satisfy both conditions:
1. The gene event is a "Leaf" (not a loss)
2. The gene is mapped to an **extant species** (a species tree leaf with `bd_event = Leaf`, not `Extinction`)

For species trees simulated with birth-death, this correctly excludes genes that survive on lineages that went extinct before the present. For species trees parsed from Newick (without birth-death events), all leaf nodes are treated as extant.

The `require_extant=True` parameter in `simulate_dtl()` and `simulate_dtl_batch()` uses this definition to ensure returned gene trees have at least one truly extant gene.

## Implementation Guide for Missing Features

Features marked with ⚠️ in the API reference are documented but not yet implemented in the Python bindings. They are available in the R bindings and can be ported to Python. Here's how to implement them:

### Adding Methods to PySpeciesTree

To add `get_bd_events()`, `save_bd_events_csv()`, `pairwise_distances()`, and `save_pairwise_distances_csv()` to Python:

1. **Location**: Edit `/home/enzo/Documents/git/WP2/rustree/src/python.rs`

2. **Reference Implementation**: See `/home/enzo/Documents/git/WP2/rustree/src/r.rs` for the R versions of these methods

3. **Example Pattern for `save_bd_events_csv()`**:
```rust
// Add to #[pymethods] impl PySpeciesTree block in python.rs
fn save_bd_events_csv(&self, filepath: &str) -> PyResult<()> {
    use crate::bd::TreeEvent;
    use crate::io::save_bd_events_to_csv;

    // Get BD events from tree nodes
    let events: Vec<TreeEvent> = // ... extract from self.tree.nodes

    save_bd_events_to_csv(&events, &self.tree, filepath)
        .map_err(|e| PyValueError::new_err(format!("Failed to write CSV: {}", e)))?;

    Ok(())
}
```

4. **Example Pattern for `pairwise_distances()`**:
```rust
// Add to #[pymethods] impl PySpeciesTree block in python.rs
fn pairwise_distances(&self, distance_type: &str, leaves_only: bool) -> PyResult<Vec<PyObject>> {
    use crate::metric_functions::{DistanceType, PairwiseDistance};

    let dist_type = match distance_type.to_lowercase().as_str() {
        "topological" | "topo" => DistanceType::Topological,
        "metric" | "branch" | "patristic" => DistanceType::Metric,
        _ => return Err(PyValueError::new_err("Invalid distance_type")),
    };

    let distances = self.tree.pairwise_distances(dist_type, leaves_only);

    // Convert to Python objects (list of dicts or custom class)
    // ... implementation details
}
```

5. **Key Imports Needed**:
```rust
use crate::bd::TreeEvent;
use crate::io::save_bd_events_to_csv;
use crate::metric_functions::{DistanceType, PairwiseDistance};
```

6. **Build and Test**:
```bash
cd /path/to/rustree
maturin develop --release
python3 -c "import rustree; tree = rustree.simulate_species_tree(10, 1.0, 0.3); tree.save_bd_events_csv('test.csv')"
```

### Adding Methods to PyGeneTree

The same methods (`pairwise_distances()`, `save_pairwise_distances_csv()`) can be added to `PyGeneTree` following the same pattern but using `self.gene_tree` instead of `self.tree`.

### Why These Features Matter

- **Birth-Death Events**: Essential for understanding macroevolutionary dynamics, extinction patterns, and diversification rates
- **Pairwise Distances**: Critical for phylogenetic distance-based methods, molecular clock analysis, and comparing gene vs species tree distances
- **Subtree Extraction**: Already implemented! Use for focused analysis, computational efficiency, and studying specific clades

For complete implementation details, examine the R binding implementations in `/home/enzo/Documents/git/WP2/rustree/src/r.rs` at lines 406-430 (BD events) and 873-935 (pairwise distances).
