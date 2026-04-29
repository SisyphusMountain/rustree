library(ape)
dyn.load("target/release/librustree.dylib")
source("R/rustree.R")
n_gene_trees <- 4L
augmented_sp_tree_ape <- read.tree(text="((((e0:0.364409,e11:0.364409):1.13153,e6:1.495939):0.556356,(((e7:0.401854,e3:0.401854):0.770008,(e16:0.976252,e
18:0.976252):0.195611):0.270193,e12:1.442055):0.610239):1.001781,((((((e21:0.030921,(e22:0.025842,(e23:0.009523,e1:0.009523):0.016318):0.00508):0.097624,e9
:0.128545):0.067287,e19:0.195832):0.194472,(e5:0.112536,e14:0.112536):0.277768):1.244906,((e10:0.741156,(e4:0.536745,e13:0.536745):0.204412):0.728835,e20:1
.469991):0.165219):1.01824,((e2:0.005197,e15:0.005197):1.40283,(e24:1.014774,(e17:0.963473,e8:0.963473):0.051301):0.393253):1.245423):0.400626);")
augmented_sp_tree_rustree <- parse_newick(write.tree(augmented_sp_tree_ape))
dtl <- c(d=0, t=0, l=0.1)
gene_trees <- simulate_dtl_batch ( #_obs is for "observed"
    species_tree = augmented_sp_tree_rustree,
    n = n_gene_trees,          # Number of gene trees
    lambda_d = dtl[1],
    lambda_t = dtl[2],
    lambda_l = dtl[3],
    transfer_alpha = 0,
    require_extant = TRUE,
    seed = 123L
  )
