// rustree/src/lib.rs
#[macro_use]
extern crate pest_derive;

// Core modules
pub mod node;   // Tree data structures and traversal
pub mod newick; // Newick parsing
pub mod io;     // I/O utilities (CSV export)

// Simulation modules
pub mod bd;     // Birth-death tree simulation
pub mod dtl;    // DTL (Duplication-Transfer-Loss) simulation

// Tree operations
pub mod sampling;
pub mod surgery;
pub mod comparison;
pub mod metric_functions;
pub mod debug;

// Analysis
pub mod induced_transfers;

// External tool integration
pub mod alerax;

// Language bindings (feature-gated)
#[cfg(feature = "python")]
pub mod python;

#[cfg(feature = "r")]
pub mod r;

// Re-export key types for easier access
pub use newick::newick::parse_newick;
pub use node::{Node, FlatNode, FlatTree, TraversalOrder, node_to_flat, flat_to_node, obtain_mapping, rename_reconciled_tree};
pub use node::{RecTree, RecTreeOwned, Event, parse_recphyloxml, parse_recphyloxml_file};
pub use metric_functions::{DistanceType, PairwiseDistance, LcaTable};
pub use sampling::{
    extract_induced_subtree, extract_induced_subtree_by_names,
    compute_lca, build_leaf_pair_lca_map, build_sampled_to_original_mapping,
    get_descendant_leaf_names, find_all_leaf_indices, find_leaf_indices_by_names,
};
pub use alerax::{run_alerax, AleRaxConfig, AleRaxFamilyResult, GeneFamily, ModelType, validate_inputs};