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
      4L, 3L,
      5L, 2L,
      5L, 1L,
      4L, 5L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  assert_true(identical(unname(phy$edge), expected_edge), "postorder edge matrix mismatch")
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

cat("\n=== Test Summary ===\n")
cat("Passed:", test_passed, "\n")
cat("Failed:", test_failed, "\n")

if (test_failed > 0) {
  quit(status = 1)
}
