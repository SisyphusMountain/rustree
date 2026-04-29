library(ape)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
repo_root <- normalizePath(script_dir)

lib_candidates <- file.path(
  repo_root, "target", "release",
  c("librustree.dylib", "librustree.so", "rustree.dll")
)
lib_path <- lib_candidates[file.exists(lib_candidates)][[1L]]
dyn.load(lib_path)
source(file.path(repo_root, "R", "rustree.R"))

sp_tree <- simulate_species_tree(20L, 1.0, 0.3, seed = 42L)
gene_tree <- simulate_dtl(sp_tree, 0.5, 0.2, 0.3, seed = 123L)
gene_trees <- simulate_dtl_batch(sp_tree, 10L, 0.5, 0.2, 0.3, seed = 123L)

sp_phylo <- as_ape_phylo(sp_tree)
gt_phylo <- gene_tree_to_ape(gene_tree)
gt_many <- gene_trees_to_ape(gene_trees)

attr(sp_phylo, "order")
head(sp_phylo$edge)

# Fastest path when ape output is all you need:
sp_phylo2 <- simulate_species_tree_ape(50L, 1.0, 0.3, seed = 42L)
gt_phylo2 <- simulate_dtl_ape(sp_tree, 0.5, 0.2, 0.3, seed = 123L)
gt_many2 <- simulate_dtl_batch_ape(sp_tree, 100L, 0.5, 0.2, 0.3, seed = 123L)

pdf(file.path("./", "rustree_ape_example.pdf"))
plot(sp_phylo)
plot(gt_many[[1]])
