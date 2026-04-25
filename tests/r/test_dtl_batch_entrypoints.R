# Regression tests for R DTL batch entry points.
# Run with: Rscript tests/r/test_dtl_batch_entrypoints.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]))
} else {
  normalizePath("tests/r/test_dtl_batch_entrypoints.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
lib_candidates <- unique(file.path(
  repo_root, "target", "release",
  paste0("librustree", c(.Platform$dynlib.ext, ".dylib", ".so", ".dll"))
))
lib_candidates <- lib_candidates[file.exists(lib_candidates)]
if (length(lib_candidates) == 0L) {
  stop("No release rustree shared library found; run `cargo build --release --features r`")
}
dyn.load(lib_candidates[[1L]])
source(file.path(repo_root, "R", "rustree.R"))

test_passed <- 0L
test_failed <- 0L

run_test <- function(name, expr) {
  tryCatch({
    force(expr)
    cat("PASS:", name, "\n")
    test_passed <<- test_passed + 1L
  }, error = function(e) {
    cat("FAIL:", name, "-", conditionMessage(e), "\n")
    test_failed <<- test_failed + 1L
  })
}

assert_true <- function(cond, msg = "") {
  if (!isTRUE(cond)) {
    stop(paste("Assertion failed:", msg))
  }
}

assert_equal <- function(actual, expected, msg = "") {
  if (!identical(actual, expected)) {
    stop(paste("Expected", expected, "but got", actual, msg))
  }
}

cat("=== DTL Batch R Entrypoint Tests ===\n\n")

run_test("native symbols are loaded", {
  assert_true(is.loaded("wrap__simulate_dtl_batch_r"),
              "wrap__simulate_dtl_batch_r should be loaded")
  assert_true(is.loaded("wrap__simulate_dtl_batch_ape_r"),
              "wrap__simulate_dtl_batch_ape_r should be loaded")
})

sptr <- simulate_species_tree(8L, lambda = 1.0, mu = 0.3, seed = 42L)
n_gene_trees <- 4L
d_val <- 0.3
t_val <- 0.2
l_val <- 0.3
seed <- 123L

run_test("simulate_dtl_batch accepts observed-style named arguments", {
  gene_trees <- simulate_dtl_batch(
    species_tree = sptr,
    n = n_gene_trees,
    lambda_d = d_val,
    lambda_t = t_val,
    lambda_l = l_val,
    transfer_alpha = 0,
    require_extant = TRUE,
    seed = seed
  )

  assert_true(is.list(gene_trees), "simulate_dtl_batch should return a list")
  assert_equal(length(gene_trees), n_gene_trees, "batch length")
  extant_counts <- vapply(gene_trees, gene_tree_num_extant, integer(1L))
  assert_true(all(extant_counts >= 1L), "all require_extant trees should have extant genes")
})

run_test("simulate_dtl_batch_ape accepts matching named arguments", {
  gene_trees <- simulate_dtl_batch_ape(
    species_tree = sptr,
    n = n_gene_trees,
    lambda_d = d_val,
    lambda_t = t_val,
    lambda_l = l_val,
    transfer_alpha = 0,
    require_extant = TRUE,
    seed = seed
  )

  assert_true(inherits(gene_trees, "multiPhylo"), "simulate_dtl_batch_ape should return multiPhylo")
  assert_equal(length(gene_trees), n_gene_trees, "ape batch length")
  assert_true(all(vapply(gene_trees, inherits, logical(1L), what = "phylo")),
              "all ape batch elements should be phylo")
})

run_test("simulate_dtl_batch works for provided snippet tree", {
  snippet_newick <- paste0(
    "((((e0:0.364409,e11:0.364409):1.13153,e6:1.495939):0.556356,",
    "(((e7:0.401854,e3:0.401854):0.770008,(e16:0.976252,e",
    "\n",
    "18:0.976252):0.195611):0.270193,e12:1.442055):0.610239):1.001781,",
    "((((((e21:0.030921,(e22:0.025842,(e23:0.009523,e1:0.009523):0.016318):0.00508):0.097624,e9",
    "\n",
    ":0.128545):0.067287,e19:0.195832):0.194472,(e5:0.112536,e14:0.112536):0.277768):1.244906,",
    "((e10:0.741156,(e4:0.536745,e13:0.536745):0.204412):0.728835,e20:1",
    "\n",
    ".469991):0.165219):1.01824,((e2:0.005197,e15:0.005197):1.40283,",
    "(e24:1.014774,(e17:0.963473,e8:0.963473):0.051301):0.393253):1.245423):0.400626);"
  )
  augmented_sp_tree_ape <- ape::read.tree(text = snippet_newick)
  augmented_sp_tree_rustree <- parse_newick(ape::write.tree(augmented_sp_tree_ape))
  dtl <- c(d = 0, t = 0, l = 0.1)
  n_gene_trees <- 4L

  gene_trees <- simulate_dtl_batch(
    species_tree = augmented_sp_tree_rustree,
    n = n_gene_trees,
    lambda_d = dtl[1],
    lambda_t = dtl[2],
    lambda_l = dtl[3],
    transfer_alpha = 0,
    require_extant = TRUE,
    seed = 123L
  )

  assert_true(is.list(gene_trees), "snippet simulate_dtl_batch should return a list")
  assert_equal(length(gene_trees), n_gene_trees, "snippet batch length")
  extant_counts <- vapply(gene_trees, gene_tree_num_extant, integer(1L))
  assert_true(all(extant_counts >= 1L), "snippet trees should all have extant genes")
})

cat("\nResults:", test_passed, "passed,", test_failed, "failed\n")
if (test_failed > 0L) {
  quit(status = 1L)
}
