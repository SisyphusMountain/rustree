// rustree/src/lib.rs
//! Phylogenetic tree simulation library.
//!
//! Error handling: the library is migrating from `Result<_, String>` to
//! `Result<_, RustreeError>`. Both forms coexist during the transition.

#[macro_use]
extern crate pest_derive;

// Unified error type
pub mod error;
pub use error::RustreeError;

// Core modules
pub mod io;
pub mod newick; // Newick parsing
pub mod node; // Tree data structures and traversal // I/O utilities (CSV export)

// Simulation modules
pub mod simulation;
pub use simulation::bd;
pub use simulation::dtl;

// Tree operations
pub mod comparison;
pub mod debug;
pub mod metric_functions;
pub mod robinson_foulds;
pub mod sampling;
pub mod surgery;
pub mod topology;

// Analysis
pub mod induced_transfers;

// Shared bindings utilities
pub mod bindings_common;

// External tool integration
pub mod external;
pub use external::alerax;

// Language bindings (feature-gated)
#[cfg(feature = "python")]
pub mod python;

#[cfg(feature = "r")]
pub mod r;

// Re-export core types at crate root.
// For other items, use qualified paths (e.g., `rustree::sampling::compute_lca`).
pub use newick::parse_newick;
pub use node::{parse_recphyloxml, parse_recphyloxml_file, Event, GeneForest, RecTree};
pub use node::{FlatNode, FlatTree, Node, TraversalOrder};
pub use topology::{UnlabeledShape, UnrootedShapeKey, UnrootedTopologyKey};
