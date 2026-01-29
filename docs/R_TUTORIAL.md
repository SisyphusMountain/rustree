# rustree R Package Tutorial

This tutorial explains how to build and use the rustree R bindings for simulating species trees, gene trees with DTL events, and visualizing reconciled trees.

## Prerequisites

- **R** (>= 4.0)
- **Rust** (>= 1.70)
- **Cargo**

For SVG visualization, install thirdkind:
```bash
cargo install thirdkind
```

## Building the Library

```bash
cd /path/to/rustree

# Build with R feature only (important: exclude Python to avoid symbol conflicts)
cargo build --release --no-default-features --features r

# The library is at: target/release/librustree.so (Linux)
#                    target/release/librustree.dylib (macOS)
#                    target/release/rustree.dll (Windows)
```

## Quick Start

```r
# Load the library (adjust path as needed)
dyn.load("target/release/librustree.so")

# Source the wrapper functions
source("R/rustree.R")

# Simulate a species tree with 50 extant species
sp_tree <- simulate_species_tree(50L, 1.0, 0.3, 42L)

# Check the tree
cat("Number of leaves:", tree_num_leaves(sp_tree), "\n")
cat("Newick:", tree_to_newick(sp_tree), "\n")

# Simulate gene trees with DTL events
gene_trees <- simulate_dtl_batch(sp_tree, 10L, 0.5, 0.2, 0.3, 123L)

# Work with the first gene tree
gt <- gene_trees[[1]]
cat("Number of extant genes:", gene_tree_num_extant(gt), "\n")
```

## Complete Example

### 1. Simulate a Species Tree

```r
# Birth-death species tree
# n: number of extant species
# lambda: speciation (birth) rate
# mu: extinction (death) rate (must be < lambda)
# seed: random seed for reproducibility (use NULL for random)

sp_tree <- simulate_species_tree(
  n = 100L,
  lambda = 1.0,
  mu = 0.5,
  seed = 42L
)

# Inspect the tree
tree_num_leaves(sp_tree)      # 100
tree_leaf_names(sp_tree)      # Vector of leaf names
tree_to_newick(sp_tree)       # Newick string
```

### 2. Parse an Existing Newick Tree

```r
# Parse from string
newick_str <- "((A:1,B:1):1,(C:1,D:1):1):0;"
sp_tree <- parse_newick(newick_str)

# Or read from file
newick_str <- readLines("species_tree.nwk")
sp_tree <- parse_newick(newick_str)
```

### 3. Simulate Gene Trees with DTL Events

```r
# Single gene tree
gene_tree <- simulate_dtl(
  species_tree = sp_tree,
  lambda_d = 0.5,    # Duplication rate
  lambda_t = 0.2,    # Transfer rate
  lambda_l = 0.3,    # Loss rate
  seed = 123L
)

# Batch of gene trees (more efficient)
gene_trees <- simulate_dtl_batch(
  species_tree = sp_tree,
  n = 100L,          # Number of gene trees
  lambda_d = 0.5,
  lambda_t = 0.2,
  lambda_l = 0.3,
  seed = 123L
)

# Access individual trees
gt1 <- gene_trees[[1]]
gt2 <- gene_trees[[2]]
```

### 4. Inspect Gene Trees

```r
gt <- gene_trees[[1]]

# Number of extant genes (surviving, non-lost)
gene_tree_num_extant(gt)

# Convert to Newick
gene_tree_to_newick(gt)

# The gene tree is a list with columns:
# - name: node names
# - parent: parent indices (NA for root)
# - left_child, right_child: child indices (NA for leaves)
# - length: branch lengths
# - depth: node depths
# - species_node: mapped species name
# - event: "Speciation", "Duplication", "Transfer", "Loss", "Leaf"

# Access as data frame
df <- as.data.frame(gt[c("name", "parent", "species_node", "event")])
head(df)
```

### 5. Sample Extant Genes

```r
# Extract induced subtree containing only extant genes
# (removes internal nodes leading only to losses)
sampled_gt <- sample_extant(gt)

gene_tree_to_newick(sampled_gt)
```

### 6. Export Results

```r
# Save species tree to Newick
save_newick(sp_tree, "species_tree.nwk")

# Save gene tree to Newick
save_newick(gt, "gene_tree.nwk")

# Save gene tree to RecPhyloXML (for visualization)
save_xml(gt, "gene_tree.recphyloxml")

# Save gene tree to CSV
save_csv(gt, "gene_tree.csv")
```

### 7. Visualize with thirdkind

```r
# Generate SVG (requires thirdkind to be installed)
# Returns SVG as string
svg_content <- gene_tree_to_svg(gt)

# Save to file
gene_tree_to_svg(gt, "gene_tree.svg")

# Save and open in browser
gene_tree_to_svg(gt, "gene_tree.svg", open_browser = TRUE)
```

## Full Workflow Example

```r
# Complete simulation workflow

# Load library
dyn.load("target/release/librustree.so")
source("R/rustree.R")

# 1. Simulate species tree
cat("Simulating species tree...\n")
sp_tree <- simulate_species_tree(50L, 1.0, 0.3, 42L)
cat("Species tree has", tree_num_leaves(sp_tree), "leaves\n")

# 2. Simulate 100 gene families
cat("Simulating 100 gene families...\n")
gene_trees <- simulate_dtl_batch(sp_tree, 100L, 0.5, 0.2, 0.3, 123L)

# 3. Analyze results
extant_counts <- sapply(gene_trees, gene_tree_num_extant)
cat("Extant genes per family:\n")
cat("  Min:", min(extant_counts), "\n")
cat("  Max:", max(extant_counts), "\n")
cat("  Mean:", mean(extant_counts), "\n")

# 4. Save results
save_newick(sp_tree, "species.nwk")

for (i in seq_along(gene_trees)) {
  save_newick(gene_trees[[i]], sprintf("gene_family_%03d.nwk", i))
}

# 5. Visualize first gene family
save_xml(gene_trees[[1]], "gene_family_001.xml")
gene_tree_to_svg(gene_trees[[1]], "gene_family_001.svg")

cat("Done! Results saved.\n")
```

## Function Reference

### Species Tree Functions

| Function | Description |
|----------|-------------|
| `simulate_species_tree(n, lambda, mu, seed)` | Simulate birth-death tree |
| `parse_newick(newick_str)` | Parse Newick string |
| `tree_to_newick(tree)` | Convert to Newick string |
| `tree_num_leaves(tree)` | Count leaves |
| `tree_leaf_names(tree)` | Get leaf names |

### Gene Tree Functions

| Function | Description |
|----------|-------------|
| `simulate_dtl(sp_tree, lambda_d, lambda_t, lambda_l, seed)` | Simulate single gene tree |
| `simulate_dtl_batch(sp_tree, n, lambda_d, lambda_t, lambda_l, seed)` | Simulate batch |
| `gene_tree_num_extant(gt)` | Count extant genes |
| `gene_tree_to_newick(gt)` | Convert to Newick |
| `gene_tree_to_xml(gt)` | Convert to RecPhyloXML |
| `sample_extant(gt)` | Extract extant-only subtree |

### I/O Functions

| Function | Description |
|----------|-------------|
| `save_newick(tree, filepath)` | Save to Newick file |
| `save_xml(gt, filepath)` | Save to RecPhyloXML file |
| `save_csv(gt, filepath)` | Save to CSV file |
| `gene_tree_to_svg(gt, filepath, open_browser)` | Generate SVG |

## Notes

- All integer arguments must use `L` suffix in R (e.g., `50L` not `50`)
- Use `NULL` for optional seed to get random behavior
- DTL rates are per unit time along branches
- The `seed` parameter ensures reproducibility
