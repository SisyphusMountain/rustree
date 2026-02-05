# rustree R Package Tutorial

This tutorial explains how to build and use the rustree R bindings for simulating species trees, gene trees with DTL (Duplication-Transfer-Loss) events, and visualizing reconciled trees.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Building the Library](#building-the-library)
3. [Quick Start](#quick-start)
4. [Complete Examples](#complete-example)
   - Simulate Species Trees
   - Parse Newick Trees
   - Simulate Gene Trees with DTL
   - Assortative Transfers
   - Per-Species DTL Model
   - Inspect Gene Trees
   - Sample Extant Genes
   - Export Birth-Death Events
   - Pairwise Distances
   - Export Gene Tree Results
   - Visualize with thirdkind
5. [Analyzing Birth-Death Events](#analyzing-birth-death-events)
6. [Pairwise Distances](#pairwise-distances)
7. [Full Workflow Example](#full-workflow-example)
8. [Advanced Examples](#advanced-comparative-analysis-example)
9. [Function Reference](#function-reference)
10. [Understanding Data Structures](#understanding-the-data-structures)
11. [Parameters and Configuration](#parameters-and-configuration)
12. [Important Notes](#important-notes)

## Quick Reference Card

```r
# === LIBRARY SETUP ===
dyn.load("target/release/librustree.so")  # Linux
# dyn.load("target/release/librustree.dylib")  # macOS
source("R/rustree.R")

# === SPECIES TREE ===
sp_tree <- simulate_species_tree(50L, lambda=1.0, mu=0.3, seed=42L)
sp_tree <- parse_newick("((A:1,B:1):1,C:1):0;")
newick <- tree_to_newick(sp_tree)
n_leaves <- tree_num_leaves(sp_tree)
leaves <- tree_leaf_names(sp_tree)

# === BIRTH-DEATH EVENTS (NEW!) ===
save_bd_events_csv(sp_tree, "bd_events.csv")
events <- sp_tree$events  # time, node_id, event_type, child1, child2

# === GENE TREES ===
gt <- simulate_dtl(sp_tree, lambda_d=0.5, lambda_t=0.2, lambda_l=0.3, seed=123L)
gts <- simulate_dtl_batch(sp_tree, n=100L, lambda_d=0.5, lambda_t=0.2, lambda_l=0.3)
# Ensure trees have at least one extant gene (retry until success):
gt <- simulate_dtl(sp_tree, 0.5, 0.2, 0.3, require_extant=TRUE, seed=123L)
n_extant <- gene_tree_num_extant(gt)
sampled <- sample_extant(gt)

# === TREE SAMPLING ===
subset <- extract_induced_subtree_by_names(sp_tree, c("A", "B", "C"))  # Keep specific leaves

# === ASSORTATIVE TRANSFERS ===
gt <- simulate_dtl(sp_tree, 0.5, 0.5, 0.3, transfer_alpha=1.0)  # Local transfers

# === PER-SPECIES DTL MODEL (Zombi-style) ===
gt <- simulate_dtl_per_species(sp_tree, 0.5, 0.2, 0.3, seed=123L)
gts <- simulate_dtl_per_species_batch(sp_tree, n=100L, 0.5, 0.2, 0.3, seed=123L)

# === PAIRWISE DISTANCES ===
dists <- pairwise_distances(sp_tree, "metric", leaves_only = TRUE)
save_pairwise_distances_csv(sp_tree, "distances.csv", "metric")

# === EXPORT ===
save_newick(sp_tree, "species.nwk")  # Species or gene tree
save_bd_events_csv(sp_tree, "bd_events.csv")  # Birth-death events
save_csv(gt, "gene_tree.csv")  # Gene tree with events
save_xml(gt, "gene_tree.xml")  # RecPhyloXML
gene_tree_to_svg(gt, "gene_tree.svg")  # SVG visualization

# === ANALYSIS ===
table(gt$event)  # Count DTL event types
table(sp_tree$events$event_type)  # Count BD event types
```

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

# NEW: Access birth-death events
cat("BD Events:", nrow(as.data.frame(sp_tree$events)), "\n")
table(sp_tree$events$event_type)

# Export birth-death events
save_bd_events_csv(sp_tree, "bd_events.csv")

# Simulate gene trees with DTL events
gene_trees <- simulate_dtl_batch(sp_tree, 10L, 0.5, 0.2, 0.3, seed = 123L)

# Work with the first gene tree
gt <- gene_trees[[1]]
cat("Number of extant genes:", gene_tree_num_extant(gt), "\n")

# Export gene tree events
save_csv(gt, "gene_tree.csv")
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
# Single gene tree (uniform random transfers)
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

# Ensure all gene trees have at least one extant gene
# (useful when losses are high and some trees might go completely extinct)
gene_tree <- simulate_dtl(
  species_tree = sp_tree,
  lambda_d = 0.5,
  lambda_t = 0.2,
  lambda_l = 0.8,    # High loss rate
  require_extant = TRUE,  # Retry until we get an extant gene
  seed = 123L
)

# Works with batch too - only keeps trees with extant genes
gene_trees <- simulate_dtl_batch(
  species_tree = sp_tree,
  n = 100L,
  lambda_d = 0.5,
  lambda_t = 0.2,
  lambda_l = 0.8,
  require_extant = TRUE,  # All returned trees have >= 1 extant gene
  seed = 123L
)
```

### 4. Assortative (Distance-Dependent) Transfers

By default, transfer recipients are chosen uniformly at random from contemporary species. You can enable **assortative transfers** where closer species are more likely to receive transfers:

```r
# Assortative transfers: P(recipient) ∝ exp(-alpha * distance)
# Higher alpha = more local transfers (closer species preferred)
# alpha = 0 is equivalent to uniform random

# Single gene tree with assortative transfers
gene_tree <- simulate_dtl(
  species_tree = sp_tree,
  lambda_d = 0.5,
  lambda_t = 0.5,
  lambda_l = 0.3,
  transfer_alpha = 1.0,  # Distance decay parameter
  seed = 123L
)

# Batch with assortative transfers
gene_trees <- simulate_dtl_batch(
  species_tree = sp_tree,
  n = 100L,
  lambda_d = 0.5,
  lambda_t = 0.5,
  lambda_l = 0.3,
  transfer_alpha = 2.0,  # Stronger distance preference
  seed = 123L
)
```

The distance between species A and B at time t is computed as:
```
d(A, B, t) = 2 × (t - depth_of_LCA(A, B))
```

### 4a. Per-Species DTL Model (Zombi-Style)

The standard DTL model (`simulate_dtl`) draws event rates proportional to the number of **gene copies** -- more gene copies means more events per unit time. The per-species model (`simulate_dtl_per_species`) instead draws event rates proportional to the number of **alive species**. When an event is triggered, a random alive species is chosen uniformly. If that species has no gene copies, the event "fails": time advances but nothing happens. This means duplications do not increase the overall event rate, producing dynamics that more closely match the Zombi simulator.

```r
# Single gene tree with per-species model
gene_tree <- simulate_dtl_per_species(
  species_tree = sp_tree,
  lambda_d = 0.5,    # Duplication rate per species per unit time
  lambda_t = 0.2,    # Transfer rate per species per unit time
  lambda_l = 0.3,    # Loss rate per species per unit time
  seed = 123L
)

# Batch of gene trees (more efficient, precomputes shared data once)
gene_trees <- simulate_dtl_per_species_batch(
  species_tree = sp_tree,
  n = 100L,
  lambda_d = 0.5,
  lambda_t = 0.2,
  lambda_l = 0.3,
  seed = 123L
)

# All the same options are supported: assortative transfers and require_extant
gene_tree <- simulate_dtl_per_species(
  species_tree = sp_tree,
  lambda_d = 0.5,
  lambda_t = 0.3,
  lambda_l = 0.5,
  transfer_alpha = 1.0,   # Distance-dependent transfers
  require_extant = TRUE,   # Retry until at least one extant gene
  seed = 123L
)
```

**Key difference from `simulate_dtl`:** In the standard model, if a gene duplicates 10 times in one species, the total event rate increases because there are now 10 gene copies each independently experiencing events. In the per-species model, the event rate stays the same because it depends only on the number of alive species. A species with 10 gene copies is no more likely to be chosen for the next event than a species with 1 copy.

```r
# Compare the two models on the same species tree
sp_tree <- simulate_species_tree(30L, 1.0, 0.4, seed = 42L)

# Standard per-gene-copy model
gt_standard <- simulate_dtl_batch(
  sp_tree, 50L,
  lambda_d = 0.8, lambda_t = 0.2, lambda_l = 0.3,
  seed = 1L
)

# Per-species model (Zombi-style)
gt_per_species <- simulate_dtl_per_species_batch(
  sp_tree, 50L,
  lambda_d = 0.8, lambda_t = 0.2, lambda_l = 0.3,
  seed = 1L
)

# Per-species model typically produces smaller gene families
# because duplications don't cause explosive growth
cat("Standard model - mean extant:", mean(sapply(gt_standard, gene_tree_num_extant)), "\n")
cat("Per-species model - mean extant:", mean(sapply(gt_per_species, gene_tree_num_extant)), "\n")
```

### 5. Inspect Gene Trees

```r
gt <- gene_trees[[1]]

# Number of extant genes (surviving genes on extant species)
# Only counts gene leaves mapped to extant species (not extinctions)
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

### 6. Sample Extant Genes

```r
# Extract induced subtree containing only extant genes
# (removes internal nodes leading only to losses)
sampled_gt <- sample_extant(gt)

gene_tree_to_newick(sampled_gt)
```

### 6a. Extract Induced Subtree by Leaf Names

You can extract an induced subtree containing only specified leaves. This is useful for focusing on particular species or genes of interest:

```r
# Parse a species tree
tree <- parse_newick("((A:1,B:1):1,(C:1,D:1):1):0;")
cat("Original tree leaves:", paste(tree_leaf_names(tree), collapse=", "), "\n")
# Output: A, B, C, D

# Extract subset with only A and C
# The induced subtree keeps only these leaves and their MRCA (root)
# Intermediate nodes with only one kept descendant are collapsed
subset_AC <- extract_induced_subtree_by_names(tree, c("A", "C"))
cat("Subset tree:", tree_to_newick(subset_AC), "\n")
# Output: (A:2.0,C:2.0):0;
# Note: Branch lengths are accumulated (1.0 to parent + 1.0 to root = 2.0)

# Extract siblings A and B
# Their parent node is preserved since both subtrees have descendants
subset_AB <- extract_induced_subtree_by_names(tree, c("A", "B"))
cat("Siblings tree:", tree_to_newick(subset_AB), "\n")
# Output: (A:1.0,B:1.0):1.0;

# Extract single leaf
# All branch lengths to root are accumulated
single <- extract_induced_subtree_by_names(tree, c("B"))
cat("Single leaf:", tree_to_newick(single), "\n")
# Output: B:2.0;

# Works with gene trees too (reconciliation info not preserved in output)
gt <- simulate_dtl(sp_tree, 0.5, 0.2, 0.3, seed = 42L)
all_leaves <- tree_leaf_names(gt)
subset_genes <- extract_induced_subtree_by_names(gt, all_leaves[1:5])
cat("Gene subset leaves:", tree_num_leaves(subset_genes), "\n")
```

### 7. Export Birth-Death Events

When you simulate a species tree, the simulation events (speciations, extinctions) are automatically captured and included in the tree object:

```r
# Simulate species tree with events
sp_tree <- simulate_species_tree(
  n = 20L,
  lambda = 1.0,
  mu = 0.5,
  seed = 42L
)

# The tree object now includes an 'events' field
# Inspect the events structure
str(sp_tree$events)

# Events include:
# - time: when the event occurred (going backward from present at 0)
# - node_id: the node where the event occurred
# - event_type: "Speciation", "Extinction", or "Leaf"
# - child1, child2: child node IDs (for speciation events)

# Save birth-death events to CSV
save_bd_events_csv(sp_tree, "bd_events.csv")

# Read and analyze the events
events_df <- read.csv("bd_events.csv")
head(events_df)

# Count event types
table(events_df$event_type)

# Plot event times
hist(events_df$time[events_df$event_type == "Speciation"],
     main = "Speciation Events Over Time",
     xlab = "Time (backward from present)",
     col = "skyblue")
```

### 8. Export Gene Tree Results

```r
# Save species tree to Newick
save_newick(sp_tree, "species_tree.nwk")

# Save gene tree to Newick
save_newick(gt, "gene_tree.nwk")

# Save gene tree to RecPhyloXML (for visualization)
save_xml(gt, "gene_tree.recphyloxml")

# Save gene tree events to CSV
save_csv(gt, "gene_tree.csv")
```

### 9. Visualize with thirdkind

```r
# Generate SVG (requires thirdkind to be installed)
# Returns SVG as string
svg_content <- gene_tree_to_svg(gt)

# Save to file
gene_tree_to_svg(gt, "gene_tree.svg")

# Save and open in browser
gene_tree_to_svg(gt, "gene_tree.svg", open_browser = TRUE)
```

## Analyzing Birth-Death Events

The birth-death simulation tracks all events that occur during the tree-building process. This information is useful for understanding the dynamics of speciation and extinction:

```r
# Simulate a species tree
sp_tree <- simulate_species_tree(50L, 1.0, 0.5, seed = 42L)

# Access the events
events <- sp_tree$events

# Convert to data frame for analysis
events_df <- data.frame(
  time = events$time,
  node_id = events$node_id,
  event_type = events$event_type,
  child1 = ifelse(is.na(events$child1), NA, events$child1),
  child2 = ifelse(is.na(events$child2), NA, events$child2)
)

# Count events by type
table(events_df$event_type)
# Speciation Extinction      Leaf
#        49         50        50

# Analyze temporal patterns
# (Note: time goes backward from present at 0)
speciation_times <- events_df$time[events_df$event_type == "Speciation"]
extinction_times <- events_df$time[events_df$event_type == "Extinction"]

cat("Tree height (time to root):", max(events_df$time), "\n")
cat("Mean speciation time:", mean(speciation_times), "\n")
cat("Mean extinction time:", mean(extinction_times), "\n")

# Plot event distribution over time
library(ggplot2)  # Optional, for nicer plots

# Simple histogram
par(mfrow = c(1, 2))
hist(speciation_times, breaks = 20,
     main = "Speciation Events", xlab = "Time (backward)",
     col = "lightblue")
hist(extinction_times, breaks = 20,
     main = "Extinction Events", xlab = "Time (backward)",
     col = "salmon")
par(mfrow = c(1, 1))

# Save events for further analysis
save_bd_events_csv(sp_tree, "bd_events.csv")

# Read back and analyze
bd_events <- read.csv("bd_events.csv")

# Calculate diversification metrics
total_time <- max(bd_events$time)
n_speciation <- sum(bd_events$event_type == "Speciation")
n_extinction <- sum(bd_events$event_type == "Extinction")

observed_lambda <- n_speciation / total_time
observed_mu <- n_extinction / total_time

cat("Observed speciation rate:", observed_lambda, "\n")
cat("Observed extinction rate:", observed_mu, "\n")
cat("Net diversification rate:", observed_lambda - observed_mu, "\n")
```

## Pairwise Distances

You can compute pairwise distances between nodes in a species tree. Two types of distances are supported:

- **Metric distance** (also called patristic distance): Sum of branch lengths along the path between two nodes
- **Topological distance**: Number of edges along the path between two nodes

```r
# Parse or simulate a species tree
sp_tree <- parse_newick("((A:1,B:1):1,C:2):0;")

# Compute metric (branch length) distances between leaves
metric_dists <- pairwise_distances(sp_tree, "metric", leaves_only = TRUE)
print(metric_dists)
#   node1 node2 distance
# 1     A     A        0
# 2     A     B        2   # A to AB (1) + AB to B (1)
# 3     A     C        4   # A to AB (1) + AB to root (1) + root to C (2)
# 4     B     A        2
# 5     B     B        0
# 6     B     C        4
# 7     C     A        4
# 8     C     B        4
# 9     C     C        0

# Compute topological (edge count) distances
topo_dists <- pairwise_distances(sp_tree, "topological", leaves_only = TRUE)
print(topo_dists)
#   node1 node2 distance
# 1     A     A        0
# 2     A     B        2   # 2 edges: A-AB, AB-B
# 3     A     C        3   # 3 edges: A-AB, AB-root, root-C
# ...

# Include internal nodes as well
all_dists <- pairwise_distances(sp_tree, "metric", leaves_only = FALSE)

# Save to CSV for external analysis
save_pairwise_distances_csv(sp_tree, "distances.csv", "metric", leaves_only = TRUE)
```

### Distance Matrix Analysis

```r
# Convert to a distance matrix for clustering or other analyses
sp_tree <- simulate_species_tree(20L, 1.0, 0.5, seed = 42L)
dists <- pairwise_distances(sp_tree, "metric", leaves_only = TRUE)

# Create a square distance matrix
leaves <- unique(dists$node1)
n <- length(leaves)
dist_matrix <- matrix(0, n, n, dimnames = list(leaves, leaves))
for (i in 1:nrow(dists)) {
  dist_matrix[dists$node1[i], dists$node2[i]] <- dists$distance[i]
}

# Use for hierarchical clustering
hc <- hclust(as.dist(dist_matrix))
plot(hc, main = "Hierarchical Clustering from Pairwise Distances")

# Or compute summary statistics
cat("Mean pairwise distance:", mean(dists$distance[dists$node1 != dists$node2]), "\n")
cat("Max pairwise distance:", max(dists$distance), "\n")
```

## Full Workflow Example

```r
# Complete simulation workflow with comprehensive analysis

# Load library
dyn.load("target/release/librustree.so")
source("R/rustree.R")

# ============================================================
# 1. SIMULATE SPECIES TREE
# ============================================================
cat("Simulating species tree...\n")
sp_tree <- simulate_species_tree(
  n = 50L,
  lambda = 1.0,  # Speciation rate
  mu = 0.3,      # Extinction rate
  seed = 42L
)
cat("Species tree has", tree_num_leaves(sp_tree), "leaves\n")

# ============================================================
# 2. ANALYZE SPECIES TREE BIRTH-DEATH EVENTS
# ============================================================
cat("\nAnalyzing birth-death events...\n")

# Save events to CSV
save_bd_events_csv(sp_tree, "bd_events.csv")

# Read and summarize
bd_events <- read.csv("bd_events.csv")
cat("Event counts:\n")
print(table(bd_events$event_type))

cat("Tree age (time to root):", max(bd_events$time), "\n")

# ============================================================
# 3. SIMULATE GENE FAMILIES
# ============================================================
cat("\nSimulating 100 gene families with DTL events...\n")
gene_trees <- simulate_dtl_batch(
  species_tree = sp_tree,
  n = 100L,
  lambda_d = 0.5,      # Duplication rate
  lambda_t = 0.2,      # Transfer rate
  lambda_l = 0.3,      # Loss rate
  transfer_alpha = 1.0,  # Distance-dependent transfers
  seed = 123L
)

# ============================================================
# 4. ANALYZE GENE TREE RESULTS
# ============================================================
cat("\nAnalyzing gene trees...\n")

extant_counts <- sapply(gene_trees, gene_tree_num_extant)
cat("Extant genes per family:\n")
cat("  Min:", min(extant_counts), "\n")
cat("  Max:", max(extant_counts), "\n")
cat("  Mean:", mean(extant_counts), "\n")
cat("  Median:", median(extant_counts), "\n")

# Count DTL events across all gene families
count_events <- function(gt) {
  table(gt$event)
}

all_events <- lapply(gene_trees, count_events)
event_summary <- Reduce("+", all_events)
cat("\nTotal DTL events across all gene families:\n")
print(event_summary)

# ============================================================
# 5. SAVE RESULTS
# ============================================================
cat("\nSaving results...\n")

# Save species tree and events
save_newick(sp_tree, "species.nwk")
save_bd_events_csv(sp_tree, "species_bd_events.csv")

# Save all gene families
dir.create("gene_families", showWarnings = FALSE)
for (i in seq_along(gene_trees)) {
  # Newick format
  save_newick(gene_trees[[i]],
              sprintf("gene_families/family_%03d.nwk", i))

  # CSV with DTL events
  save_csv(gene_trees[[i]],
           sprintf("gene_families/family_%03d.csv", i))
}

# ============================================================
# 6. VISUALIZE REPRESENTATIVE GENE FAMILY
# ============================================================
cat("\nGenerating visualization for family 1...\n")

# Save as RecPhyloXML
save_xml(gene_trees[[1]], "gene_families/family_001.xml")

# Generate SVG (requires thirdkind)
tryCatch({
  gene_tree_to_svg(gene_trees[[1]], "gene_families/family_001.svg")
  cat("  SVG saved to gene_families/family_001.svg\n")
}, error = function(e) {
  cat("  Note: thirdkind not installed, skipping SVG generation\n")
  cat("  Install with: cargo install thirdkind\n")
})

# ============================================================
# 7. SUMMARY STATISTICS
# ============================================================
cat("\n" , "=" , rep("=", 50), "\n", sep = "")
cat("SIMULATION SUMMARY\n")
cat(rep("=", 51), "\n", sep = "")
cat("Species tree:\n")
cat("  Extant species:", tree_num_leaves(sp_tree), "\n")
cat("  Total nodes:", length(sp_tree$name), "\n")
cat("  Tree age:", max(bd_events$time), "\n")
cat("\nGene families:\n")
cat("  Number of families:", length(gene_trees), "\n")
cat("  Avg extant genes:", mean(extant_counts), "\n")
cat("  Families with 0 extant:", sum(extant_counts == 0), "\n")
cat("  Families with >10 extant:", sum(extant_counts > 10), "\n")
cat("\nOutput files:\n")
cat("  species.nwk - Species tree\n")
cat("  species_bd_events.csv - Birth-death events\n")
cat("  gene_families/*.nwk - Gene family trees\n")
cat("  gene_families/*.csv - Gene family events\n")
cat("  gene_families/family_001.xml - RecPhyloXML\n")
cat("  gene_families/family_001.svg - Visualization\n")
cat(rep("=", 51), "\n", sep = "")
cat("\nDone! All results saved.\n")
```

## Advanced: Comparative Analysis Example

```r
# Compare different transfer models

# Load library
dyn.load("target/release/librustree.so")
source("R/rustree.R")

# Simulate species tree
sp_tree <- simulate_species_tree(30L, 1.0, 0.4, seed = 100L)

# Simulate with different transfer models
n_families <- 50L

# Model 1: No transfers
cat("Simulating without transfers...\n")
no_transfer <- simulate_dtl_batch(
  sp_tree, n_families,
  lambda_d = 0.5, lambda_t = 0.0, lambda_l = 0.3,
  seed = 1L
)

# Model 2: Uniform random transfers
cat("Simulating with uniform transfers...\n")
uniform_transfer <- simulate_dtl_batch(
  sp_tree, n_families,
  lambda_d = 0.5, lambda_t = 0.3, lambda_l = 0.3,
  transfer_alpha = NULL,  # or 0
  seed = 1L
)

# Model 3: Assortative transfers (local)
cat("Simulating with local transfers...\n")
local_transfer <- simulate_dtl_batch(
  sp_tree, n_families,
  lambda_d = 0.5, lambda_t = 0.3, lambda_l = 0.3,
  transfer_alpha = 2.0,  # Strong local preference
  seed = 1L
)

# Compare results
compare_models <- function(trees, model_name) {
  extant <- sapply(trees, gene_tree_num_extant)
  transfers <- sapply(trees, function(gt) sum(gt$event == "Transfer"))

  cat("\n", model_name, ":\n", sep = "")
  cat("  Mean extant genes:", mean(extant), "\n")
  cat("  Mean transfers:", mean(transfers), "\n")
  cat("  Families with transfers:", sum(transfers > 0), "/", length(trees), "\n")
}

compare_models(no_transfer, "No transfer")
compare_models(uniform_transfer, "Uniform transfer")
compare_models(local_transfer, "Local transfer")

# Visualize comparison
extant_no <- sapply(no_transfer, gene_tree_num_extant)
extant_uni <- sapply(uniform_transfer, gene_tree_num_extant)
extant_loc <- sapply(local_transfer, gene_tree_num_extant)

boxplot(extant_no, extant_uni, extant_loc,
        names = c("No transfer", "Uniform", "Local"),
        ylab = "Extant genes per family",
        main = "Effect of Transfer Model on Gene Family Size",
        col = c("lightblue", "lightgreen", "lightyellow"))
```

## Function Reference

### Species Tree Functions

| Function | Description | Returns |
|----------|-------------|---------|
| `simulate_species_tree(n, lambda, mu, seed)` | Simulate birth-death tree with n extant species | List with tree structure and events |
| `parse_newick(newick_str)` | Parse Newick-formatted string into tree | List with tree structure |
| `tree_to_newick(tree)` | Convert tree to Newick format | String |
| `tree_num_leaves(tree)` | Count number of leaves in tree | Integer |
| `tree_leaf_names(tree)` | Get names of all leaves | Character vector |

**Parameters for `simulate_species_tree`:**
- `n`: Number of extant species (Integer with `L` suffix, e.g., `50L`)
- `lambda`: Speciation/birth rate (Numeric, must be > 0)
- `mu`: Extinction/death rate (Numeric, must be >= 0 and < lambda)
- `seed`: Random seed for reproducibility (Integer with `L` or `NULL` for random)

**Returns:** A list containing:
- `name`: Node names (character vector)
- `parent`: Parent node indices (integer vector, NA for root)
- `left_child`, `right_child`: Child indices (integer vectors, NA for leaves)
- `length`: Branch lengths (numeric vector)
- `depth`: Node depths from root (numeric vector)
- `root`: Root node index (integer)
- `events`: Birth-death events (nested list with time, node_id, event_type, child1, child2)

### Gene Tree Functions

| Function | Description | Returns |
|----------|-------------|---------|
| `simulate_dtl(sp_tree, lambda_d, lambda_t, lambda_l, transfer_alpha, require_extant, seed)` | Simulate single gene tree along species tree with DTL events | List with gene tree and reconciliation |
| `simulate_dtl_batch(sp_tree, n, lambda_d, lambda_t, lambda_l, transfer_alpha, require_extant, seed)` | Simulate batch of n gene trees (more efficient) | List of gene tree lists |
| `simulate_dtl_per_species(sp_tree, lambda_d, lambda_t, lambda_l, transfer_alpha, require_extant, seed)` | Simulate single gene tree using per-species DTL model (Zombi-style) | List with gene tree and reconciliation |
| `simulate_dtl_per_species_batch(sp_tree, n, lambda_d, lambda_t, lambda_l, transfer_alpha, require_extant, seed)` | Simulate batch of n gene trees using per-species model | List of gene tree lists |
| `gene_tree_num_extant(gt)` | Count extant genes (on extant species only) | Integer |
| `gene_tree_to_newick(gt)` | Convert gene tree to Newick format | String |
| `gene_tree_to_xml(gt)` | Convert to RecPhyloXML format for visualization | String (XML) |

**Parameters for `simulate_dtl` / `simulate_dtl_batch`:**
- `species_tree`: Species tree object from `simulate_species_tree()` or `parse_newick()`
- `lambda_d`: Duplication rate per gene copy per unit time (Numeric, >= 0)
- `lambda_t`: Horizontal transfer rate per gene copy per unit time (Numeric, >= 0)
- `lambda_l`: Loss rate per gene copy per unit time (Numeric, >= 0)
- `transfer_alpha`: Distance decay for assortative transfers (Numeric or `NULL` for uniform)
  - `NULL`: Uniform random transfer recipients
  - `0`: Equivalent to uniform
  - `> 0`: Preference for closer species (higher = stronger preference)
- `require_extant`: Whether to ensure trees have extant genes (Logical, default `FALSE`)
  - `FALSE`: Return tree even if all genes are lost
  - `TRUE`: Retry simulation until at least one gene survives
- `n`: Number of gene trees to simulate (for batch, Integer with `L`)
- `seed`: Random seed (Integer with `L` or `NULL`)

**Returns (gene tree):** A list containing:
- `name`: Gene node names
- `parent`, `left_child`, `right_child`: Tree topology
- `length`: Branch lengths
- `depth`: Node depths
- `species_node`: Mapped species node name for each gene node
- `event`: Event type ("Speciation", "Duplication", "Transfer", "Loss", "Leaf")
- `species_tree`: The original species tree (nested list)

**Parameters for `simulate_dtl_per_species` / `simulate_dtl_per_species_batch`:**

These functions accept the same parameters as the standard DTL functions above, with identical types and defaults. The difference is in the underlying model, not the interface:

- `species_tree`: Species tree object from `simulate_species_tree()` or `parse_newick()`
- `lambda_d`: Duplication rate per alive species per unit time (Numeric, >= 0)
- `lambda_t`: Transfer rate per alive species per unit time (Numeric, >= 0)
- `lambda_l`: Loss rate per alive species per unit time (Numeric, >= 0)
- `transfer_alpha`: Distance decay for assortative transfers (Numeric or `NULL` for uniform)
- `require_extant`: Whether to ensure trees have extant genes (Logical, default `FALSE`)
- `n`: Number of gene trees to simulate (for batch, Integer with `L`)
- `seed`: Random seed (Integer with `L` or `NULL`)

**Model difference:** In the standard model, the total event rate at time t is `(lambda_d + lambda_t + lambda_l) * n_gene_copies(t)`. In the per-species model, the total event rate is `(lambda_d + lambda_t + lambda_l) * n_alive_species(t)`. When an event fires, a random alive species is chosen; if it has no gene copies, the event is discarded. This prevents duplications from causing runaway gene family expansion.

### Tree Sampling and Analysis Functions

| Function | Description | Returns |
|----------|-------------|---------|
| `sample_extant(gt)` | Extract induced subtree with only extant genes | List with sampled gene tree |
| `extract_induced_subtree_by_names(tree, leaf_names)` | Extract induced subtree keeping only specified leaves | List with pruned tree |
| `pairwise_distances(tree, distance_type, leaves_only)` | Compute all pairwise distances between nodes | Data frame (node1, node2, distance) |
| `save_pairwise_distances_csv(tree, filepath, distance_type, leaves_only)` | Save pairwise distances to CSV file | NULL (writes file) |

**`extract_induced_subtree_by_names` parameters:**
- `tree`: A tree list (species tree or gene tree)
- `leaf_names`: Character vector of leaf names to keep

**How it works:**
- Extracts an induced subtree containing only the specified leaves and their most recent common ancestors (MRCAs)
- Internal nodes with only one kept descendant lineage are collapsed (their branch lengths are added to the descendant)
- Works on both species trees and gene trees
- For gene trees, reconciliation information is not preserved in the output

**Example usage:**
```r
# Parse a tree
tree <- parse_newick("((A:1,B:1):1,(C:1,D:1):1):0;")

# Extract subset with A and C
subset <- extract_induced_subtree_by_names(tree, c("A", "C"))
# Result: (A:2.0,C:2.0):0;  (branch lengths collapsed)

# Extract siblings A and B
siblings <- extract_induced_subtree_by_names(tree, c("A", "B"))
# Result: (A:1.0,B:1.0):1.0;  (parent preserved)
```

**`pairwise_distances` parameters:**
- `tree`: A tree list (species tree or gene tree)
- `distance_type`: Type of distance to compute:
  - `"metric"` (or `"patristic"`, `"branch"`): Sum of branch lengths along path
  - `"topological"` (or `"topo"`): Number of edges along path
- `leaves_only`: If TRUE (default), only compute distances between leaf nodes; if FALSE, include all nodes

**Returns:** A data frame with columns:
- `node1`: Name of the first node
- `node2`: Name of the second node
- `distance`: Distance between the two nodes

The output includes all pairs (A,B), (B,A), and self-distances (A,A) for completeness.

### I/O and Export Functions

| Function | Description | Returns |
|----------|-------------|---------|
| `save_newick(tree, filepath)` | Save species or gene tree to Newick file | NULL (writes file) |
| `save_xml(gt, filepath)` | Save gene tree to RecPhyloXML file | NULL (writes file) |
| `save_csv(gt, filepath)` | Save gene tree with events to CSV | NULL (writes file) |
| `save_bd_events_csv(sp_tree, filepath)` | Save birth-death simulation events to CSV | NULL (writes file) |
| `gene_tree_to_svg(gt, filepath, open_browser)` | Generate SVG visualization using thirdkind | String (SVG content) |

**CSV File Formats:**

Birth-Death Events CSV (from `save_bd_events_csv`):
```
time,node_id,event_type,child1,child2
0.0,0,Leaf,,
0.0,1,Leaf,,
0.523,12,Speciation,0,1
1.234,15,Extinction,,
...
```

Gene Tree CSV (from `save_csv`):
```
node_id,name,parent,left_child,right_child,length,depth,species_node,event
0,g0,,,,1.234,0.0,sp5,Leaf
1,g1,5,,,0.876,1.2,sp3,Leaf
5,g5,,0,1,0.543,2.1,sp1,Duplication
...
```

## Understanding the Data Structures

### Species Tree Structure

When you simulate or parse a species tree, you get a list with the following fields:

```r
sp_tree <- simulate_species_tree(10L, 1.0, 0.5, seed = 42L)

# Tree structure (flat representation)
sp_tree$name           # Character vector: node names
sp_tree$parent         # Integer vector: parent node indices (NA for root)
sp_tree$left_child     # Integer vector: left child indices (NA for leaves)
sp_tree$right_child    # Integer vector: right child indices (NA for leaves)
sp_tree$length         # Numeric vector: branch lengths
sp_tree$depth          # Numeric vector: node depths from root
sp_tree$root           # Integer: index of root node

# Birth-death events (NEW!)
sp_tree$events$time         # Numeric vector: event times (backward from 0)
sp_tree$events$node_id      # Integer vector: node where event occurred
sp_tree$events$event_type   # Character vector: "Speciation", "Extinction", "Leaf"
sp_tree$events$child1       # Integer vector: first child (NA for non-speciation)
sp_tree$events$child2       # Integer vector: second child (NA for non-speciation)

# Example: Access leaf nodes
leaf_indices <- which(is.na(sp_tree$left_child) & is.na(sp_tree$right_child))
leaf_names <- sp_tree$name[leaf_indices]

# Example: Find the root
root_idx <- sp_tree$root
root_name <- sp_tree$name[root_idx]
```

### Gene Tree Structure

Gene trees include additional reconciliation information:

```r
gt <- simulate_dtl(sp_tree, 0.5, 0.2, 0.3, seed = 123L)

# Tree structure (same as species tree)
gt$name, gt$parent, gt$left_child, gt$right_child
gt$length, gt$depth, gt$root

# Reconciliation information
gt$species_node    # Character vector: mapped species node name for each gene node
gt$event           # Character vector: "Speciation", "Duplication", "Transfer", "Loss", "Leaf"

# The original species tree is embedded
gt$species_tree    # Nested list with species tree structure

# Example: Find all duplications
dup_indices <- which(gt$event == "Duplication")
dup_names <- gt$name[dup_indices]
dup_species <- gt$species_node[dup_indices]

cat("Found", length(dup_indices), "duplication events\n")
for (i in dup_indices) {
  cat("  Gene", gt$name[i], "duplicated in species", gt$species_node[i], "\n")
}

# Example: Count events
event_counts <- table(gt$event)
print(event_counts)

# Example: Find extant genes
extant_indices <- which(gt$event == "Leaf" & is.na(gt$left_child))
extant_names <- gt$name[extant_indices]
cat("Extant genes:", paste(extant_names, collapse = ", "), "\n")
```

### Working with Tree Indices

The trees use 0-based indexing internally but R uses 1-based indexing. The wrapper functions handle this conversion:

```r
# Access node by index (R uses 1-based indexing)
node_idx <- 5
node_name <- sp_tree$name[node_idx]
parent_idx <- sp_tree$parent[node_idx]

# Navigate the tree
if (!is.na(parent_idx)) {
  parent_name <- sp_tree$name[parent_idx]
  cat("Node", node_name, "has parent", parent_name, "\n")
}

# Get children
left_idx <- sp_tree$left_child[node_idx]
right_idx <- sp_tree$right_child[node_idx]

if (!is.na(left_idx)) {
  left_name <- sp_tree$name[left_idx]
  right_name <- sp_tree$name[right_idx]
  cat("Node", node_name, "has children", left_name, "and", right_name, "\n")
}
```

## Parameters and Configuration

### Birth-Death Parameters

- `n`: Number of extant species to simulate (Integer)
- `lambda`: Speciation/birth rate per lineage per unit time (Numeric, must be > 0)
- `mu`: Extinction/death rate per lineage per unit time (Numeric, must be >= 0 and < lambda)
- `seed`: Random seed for reproducibility (Integer or `NULL`)

The birth-death process simulates a tree with exactly `n` extant species using the backward algorithm from Stadler (2011). The tree age and number of extinction events are stochastic outcomes.

### DTL Parameters

- `lambda_d`: Duplication rate per unit time along branches (Numeric, >= 0)
- `lambda_t`: Transfer rate per unit time along branches (Numeric, >= 0)
- `lambda_l`: Loss rate per unit time along branches (Numeric, >= 0)

**Note:** For the standard model (`simulate_dtl`), these rates are per gene copy. For the per-species model (`simulate_dtl_per_species`), these rates are per alive species.
- `transfer_alpha`: Distance decay parameter for assortative transfers
  - `NULL` (default): Uniform random selection among contemporary species
  - `0`: Equivalent to uniform random
  - `> 0`: Closer species are more likely to receive transfers
  - Higher values = stronger preference for nearby species (e.g., 1.0-3.0)
- `require_extant`: Require at least one extant gene in the result (Logical)
  - `FALSE` (default): Return the simulated tree even if all genes went extinct
  - `TRUE`: Retry simulation until at least one gene survives on an extant species
  - Useful when loss rates are high and you need functional gene families

### Extant Genes Definition

The `gene_tree_num_extant()` function counts gene leaves that satisfy both conditions:
1. The gene event is a "Leaf" (not a loss)
2. The gene is mapped to an **extant species** (a species tree leaf with birth-death event type "Leaf", not "Extinction")

For species trees simulated with `simulate_species_tree()`, this correctly excludes genes that survive on lineages that went extinct before the present. For species trees parsed from Newick (without birth-death events), all leaf nodes are treated as extant.

The `require_extant = TRUE` parameter uses this definition to ensure returned gene trees have at least one truly extant gene.

DTL rates are per unit time along branches. Higher rates lead to more events. Both models support:
- **Duplications**: Gene copies within the same species
- **Transfers**: Horizontal gene transfer between contemporary species
- **Losses**: Gene lineage goes extinct
- **Speciations**: Gene follows the species tree split

Two DTL models are available:
- **Per-gene-copy model** (`simulate_dtl`): Total event rate = `(lambda_d + lambda_t + lambda_l) * number_of_gene_copies`. Duplications increase the rate because each new copy independently experiences events.
- **Per-species model** (`simulate_dtl_per_species`): Total event rate = `(lambda_d + lambda_t + lambda_l) * number_of_alive_species`. A random alive species is chosen; if it has no genes, the event fails silently. Duplications do not increase the rate.

### Assortative Transfers

Distance-dependent transfers model the observation that gene transfer is more likely between closely related or geographically proximate species:

```
P(transfer to species B | donor in species A at time t) ∝ exp(-alpha × distance(A, B, t))

where distance(A, B, t) = 2 × (t - depth_of_LCA(A, B))
```

Examples:
- `transfer_alpha = NULL`: Equal probability for all contemporary species
- `transfer_alpha = 0.5`: Slight preference for closer species
- `transfer_alpha = 1.0`: Moderate local preference
- `transfer_alpha = 2.0`: Strong local preference
- `transfer_alpha = 5.0`: Very strong local preference (rare long-distance transfers)

## Important Notes

### R-Specific Considerations

- **Integer literals**: All integer arguments must use `L` suffix in R (e.g., `50L` not `50`)
- **NULL values**: Use `NULL` for optional parameters (`seed`, `transfer_alpha`) to get default behavior
- **Index base**: R uses 1-based indexing; the wrapper handles conversion automatically
- **NA values**: Missing/optional values are represented as `NA` in integer/numeric vectors

### Performance Tips

- Use `simulate_dtl_batch()` / `simulate_dtl_per_species_batch()` instead of looping the single-tree versions - they are much faster
- For large species trees (>100 species), DTL simulation can be slow
- Setting a seed ensures reproducibility but doesn't affect performance
- Birth-death simulation time scales with the number of extinction events

### File Formats

- **Newick**: Standard phylogenetic tree format, supported by most tools
- **RecPhyloXML**: Reconciliation format for gene-species tree relationships
- **CSV**: Tabular format for detailed analysis in R, Python, or spreadsheets
- **SVG**: Vector graphics for publication-quality figures

### Troubleshooting

Common issues:

1. **"library not loaded" error**: Check that you've built with `--features r` and the path to the .so/.dylib is correct
2. **Integer/numeric type errors**: Make sure to use `L` suffix for integers (e.g., `50L`)
3. **thirdkind not found**: Install with `cargo install thirdkind` for SVG generation
4. **Simulation hangs**: Check that `lambda > mu` for birth-death; reduce rates or tree size for DTL

### References

- Stadler, T. (2011). Simulating trees with a fixed number of extant species. *Systematic Biology*, 60(5), 676-684.
- For DTL models: See the rustree paper or documentation for mathematical details

## Common Workflows

### Workflow 1: Basic Simulation for Testing

```r
# Quick simulation for testing or demonstration
library_path <- "target/release/librustree.so"
dyn.load(library_path)
source("R/rustree.R")

# Small species tree
sp <- simulate_species_tree(10L, 1.0, 0.5, seed = 1L)
cat("Species:", tree_num_leaves(sp), "\n")

# Few gene families
gts <- simulate_dtl_batch(sp, 5L, 0.5, 0.2, 0.3, seed = 1L)
sapply(gts, gene_tree_num_extant)

# Cleanup
dyn.unload(library_path)
```

### Workflow 2: Large-Scale Simulation Study

```r
# Simulate many gene families for statistical analysis
dyn.load("target/release/librustree.so")
source("R/rustree.R")

# Parameters
n_species <- 100L
n_families <- 1000L

# Simulate species tree
cat("Simulating species tree with", n_species, "species...\n")
sp_tree <- simulate_species_tree(n_species, 1.0, 0.4, seed = 42L)

# Save species tree and events
save_newick(sp_tree, "species_tree.nwk")
save_bd_events_csv(sp_tree, "species_events.csv")

# Simulate gene families in batches (for memory efficiency)
batch_size <- 100L
n_batches <- as.integer(n_families / batch_size)

all_extant_counts <- numeric(n_families)

for (batch in 1:n_batches) {
  cat("Batch", batch, "of", n_batches, "...\n")

  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- batch * batch_size

  gts <- simulate_dtl_batch(
    sp_tree, batch_size,
    lambda_d = 0.5, lambda_t = 0.2, lambda_l = 0.3,
    seed = as.integer(batch)
  )

  # Save gene trees
  for (i in seq_along(gts)) {
    idx <- start_idx + i - 1
    save_csv(gts[[i]], sprintf("gene_family_%04d.csv", idx))
    all_extant_counts[idx] <- gene_tree_num_extant(gts[[i]])
  }
}

# Analyze results
cat("\nSimulation complete!\n")
cat("Mean extant genes:", mean(all_extant_counts), "\n")
cat("Extinct families:", sum(all_extant_counts == 0), "\n")

# Save summary
write.csv(
  data.frame(family_id = 1:n_families, extant_genes = all_extant_counts),
  "summary.csv",
  row.names = FALSE
)
```

### Workflow 3: Parameter Exploration

```r
# Test different parameter combinations
dyn.load("target/release/librustree.so")
source("R/rustree.R")

# Fixed species tree
sp_tree <- simulate_species_tree(50L, 1.0, 0.4, seed = 1L)

# Parameter grid
dup_rates <- c(0.2, 0.5, 0.8)
transfer_rates <- c(0.0, 0.2, 0.4)
loss_rates <- c(0.2, 0.5, 0.8)

results <- expand.grid(
  dup = dup_rates,
  transfer = transfer_rates,
  loss = loss_rates
)
results$mean_extant <- NA
results$mean_transfers <- NA

# Run simulations
for (i in 1:nrow(results)) {
  cat("Parameter set", i, "of", nrow(results), "\n")

  gts <- simulate_dtl_batch(
    sp_tree, 100L,
    lambda_d = results$dup[i],
    lambda_t = results$transfer[i],
    lambda_l = results$loss[i],
    seed = as.integer(i)
  )

  extant <- sapply(gts, gene_tree_num_extant)
  n_transfers <- sapply(gts, function(gt) sum(gt$event == "Transfer"))

  results$mean_extant[i] <- mean(extant)
  results$mean_transfers[i] <- mean(n_transfers)
}

# Save results
write.csv(results, "parameter_exploration.csv", row.names = FALSE)

# Visualize
library(ggplot2)
ggplot(results, aes(x = dup, y = mean_extant, color = factor(transfer))) +
  geom_point(size = 3) +
  facet_wrap(~loss, labeller = label_both) +
  labs(title = "Gene Family Size vs Parameters",
       x = "Duplication Rate",
       y = "Mean Extant Genes",
       color = "Transfer Rate")
```

### Workflow 4: Reproducible Research

```r
# Fully reproducible simulation with documentation
dyn.load("target/release/librustree.so")
source("R/rustree.R")

# Document parameters
params <- list(
  n_species = 50L,
  lambda = 1.0,
  mu = 0.3,
  n_families = 100L,
  lambda_d = 0.5,
  lambda_t = 0.2,
  lambda_l = 0.3,
  transfer_alpha = 1.0,
  species_seed = 42L,
  gene_seed = 123L
)

# Save parameters
saveRDS(params, "simulation_parameters.rds")
write.csv(as.data.frame(params), "simulation_parameters.csv", row.names = FALSE)

# Simulate with documented seeds
sp_tree <- simulate_species_tree(
  params$n_species,
  params$lambda,
  params$mu,
  seed = params$species_seed
)

gene_trees <- simulate_dtl_batch(
  sp_tree,
  params$n_families,
  params$lambda_d,
  params$lambda_t,
  params$lambda_l,
  transfer_alpha = params$transfer_alpha,
  seed = params$gene_seed
)

# Save all outputs
save_newick(sp_tree, "output/species_tree.nwk")
save_bd_events_csv(sp_tree, "output/species_events.csv")

for (i in seq_along(gene_trees)) {
  save_csv(gene_trees[[i]], sprintf("output/gene_family_%03d.csv", i))
}

# Create analysis report
extant <- sapply(gene_trees, gene_tree_num_extant)

report <- list(
  parameters = params,
  species_tree_leaves = tree_num_leaves(sp_tree),
  total_gene_families = length(gene_trees),
  mean_extant_genes = mean(extant),
  median_extant_genes = median(extant),
  extinct_families = sum(extant == 0),
  timestamp = Sys.time()
)

saveRDS(report, "output/analysis_report.rds")

cat("Simulation complete and fully documented!\n")
cat("All outputs saved to output/ directory\n")
cat("Parameters saved in simulation_parameters.rds\n")
```

## Summary of All 21 Available Functions

| # | Function | Category | Description |
|---|----------|----------|-------------|
| 1 | `simulate_species_tree_r` | Species Tree | Simulate birth-death tree |
| 2 | `parse_newick_r` | Species Tree | Parse Newick string |
| 3 | `tree_to_newick_r` | Species Tree | Convert tree to Newick |
| 4 | `tree_num_leaves_r` | Species Tree | Count leaves |
| 5 | `tree_leaf_names_r` | Species Tree | Get leaf names |
| 6 | `simulate_dtl_r` | Gene Tree | Simulate single gene tree with DTL |
| 7 | `simulate_dtl_batch_r` | Gene Tree | Simulate batch of gene trees |
| 8 | `simulate_dtl_per_species_r` | Gene Tree | Simulate single gene tree with per-species DTL |
| 9 | `simulate_dtl_per_species_batch_r` | Gene Tree | Simulate batch of gene trees with per-species DTL |
| 10 | `gene_tree_num_extant_r` | Gene Tree | Count extant genes |
| 11 | `gene_tree_to_newick_r` | Gene Tree | Convert gene tree to Newick |
| 12 | `gene_tree_to_xml_r` | Gene Tree | Convert to RecPhyloXML |
| 13 | `save_newick_r` | Export | Save tree to Newick file |
| 14 | `save_xml_r` | Export | Save gene tree to RecPhyloXML |
| 15 | `save_csv_r` | Export | Save gene tree to CSV |
| 16 | `save_bd_events_csv_r` | Export | Save birth-death events to CSV |
| 17 | `save_pairwise_distances_csv_r` | Export | Save pairwise distances to CSV |
| 18 | `gene_tree_to_svg_r` | Visualization | Generate SVG with thirdkind |
| 19 | `sample_extant_r` | Analysis | Extract extant-only subtree |
| 20 | `extract_induced_subtree_by_names_r` | Analysis | Extract induced subtree by leaf names |
| 21 | `pairwise_distances_r` | Analysis | Compute pairwise node distances |

**Note:** The R wrapper functions in `R/rustree.R` remove the `_r` suffix for cleaner syntax (e.g., `simulate_species_tree` instead of `simulate_species_tree_r`).

---

**For questions, issues, or contributions, please visit the rustree GitHub repository.**
