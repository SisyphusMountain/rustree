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

name_internal_nodes <- function(tree) {
  .Call("wrap__name_internal_nodes_r", tree)
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

# Ape conversion helpers

.is_rustree_tree <- function(tree) {
  is.list(tree) &&
    all(c("name", "parent", "left_child", "right_child", "length", "root") %in% names(tree))
}

#' Convert a rustree tree to an ape phylo object.
#'
#' Both species trees and reconciled gene-tree lists use rustree's flat tree
#' columns, so this converter works for either object type. The ape layout is
#' constructed in Rust without serializing through Newick.
#'
#' @param tree A rustree species tree or gene tree list.
#' @param order Edge row order for the returned phylo object.
#' @param use_node_labels Include non-empty internal node names as node.label.
#' @param include_root_edge Store a non-zero rustree root branch length as root.edge.
#' @return An object of class `"phylo"` compatible with ape.
#' @export
as_ape_phylo <- function(tree,
                         order = c("cladewise", "postorder", "pruningwise"),
                         use_node_labels = TRUE,
                         include_root_edge = TRUE) {
  order <- match.arg(order)

  if (inherits(tree, "phylo")) {
    return(tree)
  }
  if (!.is_rustree_tree(tree)) {
    stop("tree must be a rustree tree list with name, parent, child, length, and root fields",
         call. = FALSE)
  }

  .Call("wrap__tree_to_ape_r", tree, order, as.logical(use_node_labels),
        as.logical(include_root_edge))
}

#' Convert a list of rustree trees to an ape multiPhylo object.
#'
#' @param trees A list of rustree tree lists or phylo objects.
#' @param order Edge row order for the returned phylo objects.
#' @param use_node_labels Include non-empty internal node names as node.label.
#' @param include_root_edge Store non-zero rustree root branch lengths as root.edge.
#' @return An object of class `"multiPhylo"` compatible with ape.
#' @export
as_ape_multiPhylo <- function(trees,
                              order = c("cladewise", "postorder", "pruningwise"),
                              use_node_labels = TRUE,
                              include_root_edge = TRUE) {
  order <- match.arg(order)
  if (.is_rustree_tree(trees)) {
    stop("trees is a single rustree tree; use as_ape_phylo() instead", call. = FALSE)
  }
  if (!is.list(trees)) {
    stop("trees must be a list of rustree tree lists or phylo objects", call. = FALSE)
  }

  is_rustree <- vapply(trees, .is_rustree_tree, logical(1))
  if (all(is_rustree)) {
    return(.Call("wrap__trees_to_ape_multi_r", trees, order,
                 as.logical(use_node_labels), as.logical(include_root_edge)))
  }

  converted <- lapply(trees, function(tree) {
    if (inherits(tree, "phylo")) {
      tree
    } else {
      as_ape_phylo(tree, order = order, use_node_labels = use_node_labels,
                   include_root_edge = include_root_edge)
    }
  })
  names(converted) <- names(trees)
  class(converted) <- "multiPhylo"
  converted
}

tree_to_ape <- as_ape_phylo
gene_tree_to_ape <- as_ape_phylo
gene_trees_to_ape <- as_ape_multiPhylo

simulate_species_tree_ape <- function(n, lambda, mu, seed = NULL,
                                      order = c("cladewise", "postorder", "pruningwise"),
                                      use_node_labels = TRUE,
                                      include_root_edge = TRUE) {
  if (is.null(seed)) seed <- NA_integer_
  order <- match.arg(order)
  .Call("wrap__simulate_species_tree_ape_r", as.integer(n), as.double(lambda),
        as.double(mu), as.integer(seed), order, as.logical(use_node_labels),
        as.logical(include_root_edge))
}

# Gene Tree Functions

simulate_dtl <- function(species_tree, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL,
                         require_extant = FALSE, seed = NULL, replacement_transfer = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_r", species_tree, as.double(lambda_d), as.double(lambda_t),
        as.double(lambda_l), as.double(transfer_alpha), as.double(replacement_transfer),
        as.logical(require_extant), as.integer(seed))
}

simulate_dtl_batch <- function(species_tree, n, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL,
                               require_extant = FALSE, seed = NULL, replacement_transfer = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_batch_r", species_tree, as.integer(n), as.double(lambda_d),
        as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha),
        as.double(replacement_transfer), as.logical(require_extant), as.integer(seed))
}

simulate_dtl_per_species <- function(species_tree, lambda_d, lambda_t, lambda_l,
                                     transfer_alpha = NULL, require_extant = FALSE,
                                     seed = NULL, replacement_transfer = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_per_species_r", species_tree, as.double(lambda_d),
        as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha),
        as.double(replacement_transfer), as.logical(require_extant), as.integer(seed))
}

simulate_dtl_per_species_batch <- function(species_tree, n, lambda_d, lambda_t, lambda_l,
                                           transfer_alpha = NULL, require_extant = FALSE,
                                           seed = NULL, replacement_transfer = NULL) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_per_species_batch_r", species_tree, as.integer(n),
        as.double(lambda_d), as.double(lambda_t), as.double(lambda_l),
        as.double(transfer_alpha), as.double(replacement_transfer),
        as.logical(require_extant), as.integer(seed))
}

simulate_dtl_ape <- function(species_tree, lambda_d, lambda_t, lambda_l,
                             transfer_alpha = NULL, require_extant = FALSE,
                             seed = NULL, replacement_transfer = NULL,
                             order = c("cladewise", "postorder", "pruningwise"),
                             use_node_labels = TRUE,
                             include_root_edge = TRUE) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  order <- match.arg(order)
  .Call("wrap__simulate_dtl_ape_r", species_tree, as.double(lambda_d),
        as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha),
        as.double(replacement_transfer), as.logical(require_extant),
        as.integer(seed), order, as.logical(use_node_labels),
        as.logical(include_root_edge))
}

simulate_dtl_batch_ape <- function(species_tree, n, lambda_d, lambda_t, lambda_l,
                                   transfer_alpha = NULL, require_extant = FALSE,
                                   seed = NULL, replacement_transfer = NULL,
                                   order = c("cladewise", "postorder", "pruningwise"),
                                   use_node_labels = TRUE,
                                   include_root_edge = TRUE) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  order <- match.arg(order)
  .Call("wrap__simulate_dtl_batch_ape_r", species_tree, as.integer(n),
        as.double(lambda_d), as.double(lambda_t), as.double(lambda_l),
        as.double(transfer_alpha), as.double(replacement_transfer),
        as.logical(require_extant), as.integer(seed), order,
        as.logical(use_node_labels), as.logical(include_root_edge))
}

simulate_dtl_per_species_ape <- function(species_tree, lambda_d, lambda_t, lambda_l,
                                         transfer_alpha = NULL,
                                         require_extant = FALSE, seed = NULL,
                                         replacement_transfer = NULL,
                                         order = c("cladewise", "postorder", "pruningwise"),
                                         use_node_labels = TRUE,
                                         include_root_edge = TRUE) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  order <- match.arg(order)
  .Call("wrap__simulate_dtl_per_species_ape_r", species_tree,
        as.double(lambda_d), as.double(lambda_t), as.double(lambda_l),
        as.double(transfer_alpha), as.double(replacement_transfer),
        as.logical(require_extant), as.integer(seed), order,
        as.logical(use_node_labels), as.logical(include_root_edge))
}

simulate_dtl_per_species_batch_ape <- function(species_tree, n, lambda_d, lambda_t, lambda_l,
                                               transfer_alpha = NULL,
                                               require_extant = FALSE, seed = NULL,
                                               replacement_transfer = NULL,
                                               order = c("cladewise", "postorder", "pruningwise"),
                                               use_node_labels = TRUE,
                                               include_root_edge = TRUE) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  order <- match.arg(order)
  .Call("wrap__simulate_dtl_per_species_batch_ape_r", species_tree,
        as.integer(n), as.double(lambda_d), as.double(lambda_t),
        as.double(lambda_l), as.double(transfer_alpha),
        as.double(replacement_transfer), as.logical(require_extant),
        as.integer(seed), order, as.logical(use_node_labels),
        as.logical(include_root_edge))
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

sample_leaves <- function(gene_tree, species_leaf_names) {
  .Call("wrap__sample_leaves_r", gene_tree, as.character(species_leaf_names))
}

# Induced Transfers

#' Get DTL events attached to a gene tree.
#'
#' @param gene_tree A gene tree from simulate_dtl or simulate_dtl_batch
#' @return A list with columns: event_type, time, gene_id, species,
#'         from_species, to_species. Returns NULL if no events attached.
get_dtl_events <- function(gene_tree) {
  attr(gene_tree, "dtl_events")
}

#' Compute induced transfers by projecting onto a sampled subtree.
#'
#' @param species_tree The complete species tree
#' @param sampled_leaf_names Character vector of species to keep
#' @param dtl_events DTL events list (from get_dtl_events())
#' @param mode Algorithm mode: "projection" (default) or "damien"
#' @param remove_undetectable Logical; only used in `mode="damien"`.
#' @return A data.frame with columns: time, gene_id, from_complete,
#'         to_complete, from_sampled, to_sampled
induced_transfers <- function(species_tree, sampled_leaf_names,
                              dtl_events,
                              mode = "projection",
                              remove_undetectable = FALSE) {
  result <- .Call(
    "wrap__induced_transfers_r",
    species_tree,
    as.character(sampled_leaf_names),
    dtl_events,
    as.character(mode),
    as.logical(remove_undetectable)
  )
  as.data.frame(result, stringsAsFactors = FALSE)
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

parse_recphyloxml <- function(filepath) {
  .Call("wrap__parse_recphyloxml_r", as.character(filepath))
}

gene_tree_to_svg <- function(gene_tree, filepath = NULL, open_browser = FALSE) {
  if (is.null(filepath)) filepath <- NA_character_
  .Call("wrap__gene_tree_to_svg_r", gene_tree, as.character(filepath), as.logical(open_browser))
}

# Distance Functions

#' Compute pairwise distances between nodes in a tree.
#'
#' @param tree A tree list from simulate_species_tree or parse_newick
#' @param distance_type Type of distance: "topological" (number of edges) or "metric" (sum of branch lengths)
#' @param leaves_only If TRUE, only compute distances between leaf nodes (default TRUE)
#' @return A data.frame with columns: node1, node2, distance
#' @examples
#' sp_tree <- parse_newick("((A:1,B:1):1,C:2):0;")
#' dists <- pairwise_distances(sp_tree, "metric", leaves_only = TRUE)
#' head(dists)
pairwise_distances <- function(tree, distance_type = "metric", leaves_only = TRUE) {
  result <- .Call("wrap__pairwise_distances_r", tree, as.character(distance_type), as.logical(leaves_only))
  as.data.frame(result, stringsAsFactors = FALSE)
}

#' Save pairwise distances to a CSV file.
#'
#' @param tree A tree list from simulate_species_tree or parse_newick
#' @param filepath Path to save the CSV file
#' @param distance_type Type of distance: "topological" or "metric"
#' @param leaves_only If TRUE, only compute distances between leaf nodes (default TRUE)
#' @examples
#' sp_tree <- parse_newick("((A:1,B:1):1,C:2):0;")
#' save_pairwise_distances_csv(sp_tree, "distances.csv", "metric")
save_pairwise_distances_csv <- function(tree, filepath, distance_type = "metric", leaves_only = TRUE) {
  .Call("wrap__save_pairwise_distances_csv_r", tree, as.character(filepath), as.character(distance_type), as.logical(leaves_only))
}

# Streaming Functions (memory-efficient batch output)

simulate_dtl_stream_xml <- function(newick, n, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL, replacement_transfer = NULL, require_extant = FALSE, seed = NULL, output_dir) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_stream_xml_r", as.character(newick), as.integer(n), as.double(lambda_d), as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha), as.double(replacement_transfer), as.logical(require_extant), as.integer(seed), as.character(output_dir))
}

simulate_dtl_stream_newick <- function(newick, n, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL, replacement_transfer = NULL, require_extant = FALSE, seed = NULL, output_dir) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_stream_newick_r", as.character(newick), as.integer(n), as.double(lambda_d), as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha), as.double(replacement_transfer), as.logical(require_extant), as.integer(seed), as.character(output_dir))
}

simulate_dtl_per_species_stream_xml <- function(newick, n, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL, replacement_transfer = NULL, require_extant = FALSE, seed = NULL, output_dir) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_per_species_stream_xml_r", as.character(newick), as.integer(n), as.double(lambda_d), as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha), as.double(replacement_transfer), as.logical(require_extant), as.integer(seed), as.character(output_dir))
}

simulate_dtl_per_species_stream_newick <- function(newick, n, lambda_d, lambda_t, lambda_l, transfer_alpha = NULL, replacement_transfer = NULL, require_extant = FALSE, seed = NULL, output_dir) {
  if (is.null(seed)) seed <- NA_integer_
  if (is.null(transfer_alpha)) transfer_alpha <- NA_real_
  if (is.null(replacement_transfer)) replacement_transfer <- NA_real_
  .Call("wrap__simulate_dtl_per_species_stream_newick_r", as.character(newick), as.integer(n), as.double(lambda_d), as.double(lambda_t), as.double(lambda_l), as.double(transfer_alpha), as.double(replacement_transfer), as.logical(require_extant), as.integer(seed), as.character(output_dir))
}
