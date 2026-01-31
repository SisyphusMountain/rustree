# rustree R wrapper functions
# Source this file after dyn.load() to get easy-to-use functions

# Species Tree Functions

simulate_species_tree <- function(n, lambda, mu, seed = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  .Call("wrap__simulate_species_tree_r", as.integer(n), as.double(lambda), as.double(mu), as.integer(seed))
}

parse_newick <- function(newick_str) {
  .Call("wrap__parse_newick_r", as.character(newick_str))
}

tree_to_newick <- function(tree) {
  .Call("wrap__tree_to_newick_r", tree)
}

tree_num_leaves <- function(tree) {
  .Call("wrap__tree_num_leaves_r", tree)
}

tree_leaf_names <- function(tree) {
  .Call("wrap__tree_leaf_names_r", tree)
}

# Gene Tree Functions

simulate_dtl <- function(species_tree, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL, seed = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  .Call("wrap__simulate_dtl_r", species_tree, as.double(lambda_d), as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha), as.integer(seed))
}

simulate_dtl_batch <- function(species_tree, n, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL, seed = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  .Call("wrap__simulate_dtl_batch_r", species_tree, as.integer(n), as.double(lambda_d), as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha), as.integer(seed))
}

gene_tree_num_extant <- function(gene_tree) {
  .Call("wrap__gene_tree_num_extant_r", gene_tree)
}

gene_tree_to_newick <- function(gene_tree) {
  .Call("wrap__gene_tree_to_newick_r", gene_tree)
}

gene_tree_to_xml <- function(gene_tree) {
  .Call("wrap__gene_tree_to_xml_r", gene_tree)
}

sample_extant <- function(gene_tree) {
  .Call("wrap__sample_extant_r", gene_tree)
}

extract_induced_subtree_by_names <- function(tree, leaf_names) {
  .Call("wrap__extract_induced_subtree_by_names_r", tree, as.character(leaf_names))
}

# I/O Functions

save_newick <- function(tree_or_gene_tree, filepath) {
  .Call("wrap__save_newick_r", tree_or_gene_tree, as.character(filepath))
}

save_xml <- function(gene_tree, filepath) {
  .Call("wrap__save_xml_r", gene_tree, as.character(filepath))
}

save_csv <- function(gene_tree, filepath) {
  .Call("wrap__save_csv_r", gene_tree, as.character(filepath))
}

save_bd_events_csv <- function(species_tree, filepath) {
  .Call("wrap__save_bd_events_csv_r", species_tree, as.character(filepath))
}

gene_tree_to_svg <- function(gene_tree, filepath = NULL, open_browser = FALSE) {
  if (is.null(filepath)) filepath <- NA_character_
  .Call("wrap__gene_tree_to_svg_r", gene_tree, as.character(filepath), as.logical(open_browser))
}
