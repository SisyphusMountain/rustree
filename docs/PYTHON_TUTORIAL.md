# rustree Python Tutorial

This tutorial explains how to use the rustree Python bindings for simulating species trees, gene trees with DTL events, and visualizing reconciled trees.

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

## API Reference

### Module Functions

| Function | Description |
|----------|-------------|
| `simulate_species_tree(n, lambda_, mu, seed=None)` | Simulate birth-death tree |
| `parse_species_tree(newick_str)` | Parse Newick string |

### PySpeciesTree Methods

| Method | Description |
|--------|-------------|
| `to_newick()` | Convert to Newick string |
| `save_newick(filepath)` | Save to Newick file |
| `num_nodes()` | Count all nodes |
| `num_leaves()` | Count leaves |
| `tree_height()` | Get total height |
| `leaf_names()` | Get list of leaf names |
| `root_index()` | Get root node index |
| `simulate_dtl(lambda_d, lambda_t, lambda_l, transfer_alpha=None, require_extant=False, seed=None)` | Simulate single gene tree |
| `simulate_dtl_batch(n, lambda_d, lambda_t, lambda_l, transfer_alpha=None, require_extant=False, seed=None)` | Simulate batch |

### PyGeneTree Methods

| Method | Description |
|--------|-------------|
| `to_newick()` | Convert to Newick string |
| `save_newick(filepath)` | Save to Newick file |
| `num_nodes()` | Count all nodes |
| `num_extant()` | Count extant genes (on extant species only) |
| `count_events()` | Return (S, D, T, L, Leaves) counts |
| `extant_gene_names()` | Get names of extant genes |
| `sample_extant()` | Extract extant-only subtree |
| `sample_by_names(names)` | Extract subtree by gene names |
| `to_xml()` | Convert to RecPhyloXML string |
| `save_xml(filepath)` | Save to RecPhyloXML file |
| `to_csv(filepath=None)` | Export to DataFrame (optionally save) |
| `to_svg(filepath=None, open_browser=False)` | Generate SVG visualization |
| `display()` | Display in Jupyter notebook |

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
