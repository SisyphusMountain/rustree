# Tests for direct conversion from rustree R lists to ape-compatible objects.
# Run with: Rscript test_ape_conversion.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]))
} else {
  normalizePath("tests/r/test_ape_conversion.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
lib_candidates <- unique(file.path(
  repo_root, "target", "release",
  paste0("librustree", c(.Platform$dynlib.ext, ".dylib", ".so", ".dll"))
))
lib_path <- lib_candidates[file.exists(lib_candidates)][[1L]]
dyn.load(lib_path)
source(file.path(repo_root, "R", "rustree.R"))

test_passed <- 0
test_failed <- 0

run_test <- function(name, expr) {
  tryCatch({
    eval(expr)
    cat("PASS:", name, "\n")
    test_passed <<- test_passed + 1
  }, error = function(e) {
    cat("FAIL:", name, "-", conditionMessage(e), "\n")
    test_failed <<- test_failed + 1
  })
}

assert_true <- function(cond, msg = "") {
  if (!isTRUE(cond)) {
    stop(paste("Assertion failed:", msg))
  }
}

assert_identical <- function(actual, expected, msg = "") {
  if (!identical(actual, expected)) {
    stop(paste("Expected", paste(expected, collapse = ","), "but got",
               paste(actual, collapse = ","), msg))
  }
}

example_tree <- function() {
  list(
    name = c("root", "left", "A", "B", "C"),
    parent = c(NA_integer_, 0L, 1L, 1L, 0L),
    left_child = c(1L, 2L, NA_integer_, NA_integer_, NA_integer_),
    right_child = c(4L, 3L, NA_integer_, NA_integer_, NA_integer_),
    length = c(0.0, 0.5, 1.0, 1.2, 2.0),
    depth = c(0.0, 0.5, 1.5, 1.7, 2.0),
    root = 0L
  )
}

reference_as_ape_phylo <- function(tree,
                                   order = c("cladewise", "postorder", "pruningwise"),
                                   use_node_labels = TRUE,
                                   include_root_edge = TRUE) {
  order <- match.arg(order)
  attr_order <- if (identical(order, "pruningwise")) "postorder" else order
  n <- length(tree$name)
  root <- as.integer(tree$root[[1L]])

  stack <- root
  preorder <- integer(0)
  visited <- rep(FALSE, n)
  while (length(stack) > 0L) {
    idx <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (visited[[idx + 1L]]) {
      stop("reference converter saw a repeated node")
    }
    visited[[idx + 1L]] <- TRUE
    preorder <- c(preorder, idx)

    right <- tree$right_child[[idx + 1L]]
    left <- tree$left_child[[idx + 1L]]
    if (!is.na(right)) stack <- c(stack, as.integer(right))
    if (!is.na(left)) stack <- c(stack, as.integer(left))
  }

  left_children <- tree$left_child[preorder + 1L]
  right_children <- tree$right_child[preorder + 1L]
  is_tip <- is.na(left_children) & is.na(right_children)
  tip_nodes <- preorder[is_tip]
  internal_nodes <- preorder[!is_tip]

  ape_id <- integer(n)
  ape_id[tip_nodes + 1L] <- seq_along(tip_nodes)
  ape_id[internal_nodes + 1L] <- length(tip_nodes) + seq_along(internal_nodes)

  if (identical(attr_order, "postorder")) {
    edge_order <- integer(0)
    for (parent in rev(internal_nodes)) {
      for (child in c(tree$left_child[[parent + 1L]],
                      tree$right_child[[parent + 1L]])) {
        if (!is.na(child)) {
          edge_order <- c(edge_order, as.integer(child))
        }
      }
    }
  } else {
    edge_order <- preorder[preorder != root]
  }
  edge <- matrix(integer(length(edge_order) * 2L), ncol = 2L)
  if (length(edge_order) > 0L) {
    parents <- tree$parent[edge_order + 1L]
    edge[, 1L] <- ape_id[parents + 1L]
    edge[, 2L] <- ape_id[edge_order + 1L]
  }

  phy <- list(
    edge = edge,
    edge.length = as.numeric(tree$length[edge_order + 1L]),
    tip.label = as.character(tree$name[tip_nodes + 1L]),
    Nnode = as.integer(length(internal_nodes))
  )

  node_labels <- as.character(tree$name[internal_nodes + 1L])
  if (isTRUE(use_node_labels) && length(node_labels) > 0L && any(nzchar(node_labels))) {
    phy$node.label <- node_labels
  }

  root_length <- as.numeric(tree$length[[root + 1L]])
  if (isTRUE(include_root_edge) && is.finite(root_length)) {
    phy$root.edge <- root_length
  }

  class(phy) <- "phylo"
  attr(phy, "order") <- attr_order
  phy
}

assert_phylo_equal <- function(actual, expected, msg = "") {
  assert_identical(class(actual), class(expected), paste(msg, "class"))
  assert_true(identical(unname(actual$edge), unname(expected$edge)),
              paste(msg, "edge"))
  assert_identical(actual$edge.length, expected$edge.length,
                   paste(msg, "edge.length"))
  assert_identical(actual$tip.label, expected$tip.label,
                   paste(msg, "tip.label"))
  assert_identical(actual$Nnode, expected$Nnode, paste(msg, "Nnode"))
  assert_identical(actual$node.label, expected$node.label,
                   paste(msg, "node.label"))
  assert_identical(actual$root.edge, expected$root.edge,
                   paste(msg, "root.edge"))
  assert_identical(attr(actual, "order"), attr(expected, "order"),
                   paste(msg, "order"))
}

assert_phylo_equal_tolerant <- function(actual, expected, msg = "",
                                        tolerance = 1e-6) {
  assert_numeric_equal_abs <- function(x, y, field) {
    if (is.null(x) || is.null(y)) {
      assert_identical(x, y, paste(msg, field))
    } else {
      assert_identical(length(x), length(y), paste(msg, field, "length"))
      if (length(x) > 0L) {
        assert_true(max(abs(x - y)) <= tolerance, paste(msg, field))
      }
    }
  }

  assert_identical(class(actual), class(expected), paste(msg, "class"))
  assert_true(identical(unname(actual$edge), unname(expected$edge)),
              paste(msg, "edge"))
  assert_numeric_equal_abs(actual$edge.length, expected$edge.length,
                           "edge.length")
  assert_identical(actual$tip.label, expected$tip.label,
                   paste(msg, "tip.label"))
  assert_identical(actual$Nnode, expected$Nnode, paste(msg, "Nnode"))
  assert_identical(actual$node.label, expected$node.label,
                   paste(msg, "node.label"))
  assert_numeric_equal_abs(actual$root.edge, expected$root.edge,
                           "root.edge")
  assert_identical(attr(actual, "order"), attr(expected, "order"),
                   paste(msg, "order"))
}

cat("=== ape Conversion Tests ===\n\n")

run_test("as_ape_phylo: class and required fields", {
  phy <- as_ape_phylo(example_tree())
  assert_identical(class(phy), "phylo")
  assert_true(is.matrix(phy$edge), "edge should be a matrix")
  assert_true(is.numeric(phy$edge.length), "edge.length should be numeric")
  assert_identical(phy$tip.label, c("A", "B", "C"))
  assert_identical(phy$Nnode, 2L)
})

run_test("as_ape_phylo: ape numbering and cladewise edge order", {
  phy <- as_ape_phylo(example_tree())
  expected_edge <- matrix(
    c(
      4L, 5L,
      5L, 1L,
      5L, 2L,
      4L, 3L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  assert_true(identical(unname(phy$edge), expected_edge), "edge matrix mismatch")
  assert_identical(phy$edge.length, c(0.5, 1.0, 1.2, 2.0))
  assert_identical(phy$node.label, c("root", "left"))
  assert_identical(attr(phy, "order"), "cladewise")
})

run_test("as_ape_phylo: postorder edge rows put children before parents", {
  phy <- as_ape_phylo(example_tree(), order = "postorder")
  expected_edge <- matrix(
    c(
      5L, 1L,
      5L, 2L,
      4L, 5L,
      4L, 3L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  assert_true(identical(unname(phy$edge), expected_edge), "postorder edge matrix mismatch")
  assert_identical(phy$edge.length, c(1.0, 1.2, 0.5, 2.0))
  assert_identical(attr(phy, "order"), "postorder")
})

run_test("as_ape_phylo: single-tip tree", {
  tree <- list(
    name = "A",
    parent = NA_integer_,
    left_child = NA_integer_,
    right_child = NA_integer_,
    length = 0.0,
    root = 0L
  )
  phy <- as_ape_phylo(tree)
  assert_identical(dim(phy$edge), c(0L, 2L))
  assert_identical(phy$tip.label, "A")
  assert_identical(phy$Nnode, 0L)
})

run_test("as_ape_multiPhylo: converts a list of rustree trees", {
  trees <- list(first = example_tree(), second = example_tree())
  multi <- as_ape_multiPhylo(trees)
  assert_identical(class(multi), "multiPhylo")
  assert_identical(names(multi), c("first", "second"))
  assert_true(inherits(multi[[1]], "phylo"), "first tree should be phylo")
  assert_true(inherits(multi[[2]], "phylo"), "second tree should be phylo")
})

run_test("aliases: tree_to_ape and gene_tree_to_ape", {
  tree <- example_tree()
  assert_true(inherits(tree_to_ape(tree), "phylo"), "tree_to_ape should return phylo")
  assert_true(inherits(gene_tree_to_ape(tree), "phylo"), "gene_tree_to_ape should return phylo")
  assert_true(inherits(gene_trees_to_ape(list(tree)), "multiPhylo"),
              "gene_trees_to_ape should return multiPhylo")
})

run_test("Rust ape conversion matches pure-R reference conversion", {
  tree <- example_tree()
  for (order in c("cladewise", "postorder", "pruningwise")) {
    actual <- as_ape_phylo(tree, order = order)
    expected <- reference_as_ape_phylo(tree, order = order)
    assert_phylo_equal(actual, expected, paste("order", order))
  }
})

run_test("Rust ape conversion matches R reference for simulated trees", {
  sp_tree <- simulate_species_tree(12L, 1.0, 0.3, seed = 99L)
  assert_phylo_equal(
    as_ape_phylo(sp_tree),
    reference_as_ape_phylo(sp_tree),
    "species tree"
  )

  gene_tree <- simulate_dtl(sp_tree, 0.5, 0.2, 0.3, seed = 100L)
  assert_phylo_equal(
    gene_tree_to_ape(gene_tree, order = "postorder"),
    reference_as_ape_phylo(gene_tree, order = "postorder"),
    "gene tree"
  )
})

run_test("Rust ape conversion matches ape Newick roundtrip", {
  if (!requireNamespace("ape", quietly = TRUE)) {
    cat("SKIP: ape package is not installed\n")
  } else {
    check_ape_roundtrip <- function(tree, label) {
      via_newick <- ape::read.tree(text = tree_to_newick(tree))
      assert_phylo_equal_tolerant(
        as_ape_phylo(tree),
        via_newick,
        paste(label, "cladewise")
      )
      assert_phylo_equal_tolerant(
        as_ape_phylo(tree, order = "postorder"),
        ape::reorder.phylo(via_newick, order = "postorder"),
        paste(label, "postorder")
      )
    }

    check_ape_roundtrip(example_tree(), "example tree")

    for (seed in 1:5) {
      sp_tree <- simulate_species_tree(20L, 1.0, 0.3, seed = seed)
      check_ape_roundtrip(sp_tree, paste("species tree", seed))

      gene_tree <- simulate_dtl(sp_tree, 0.5, 0.2, 0.3,
                                require_extant = TRUE, seed = seed + 100L)
      check_ape_roundtrip(gene_tree, paste("gene tree", seed))
    }
  }
})

cat("\n=== Test Summary ===\n")
cat("Passed:", test_passed, "\n")
cat("Failed:", test_failed, "\n")

if (test_failed > 0) {
  quit(status = 1)
}
