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

// Python bindings
pub mod python;

// Re-export key types for easier access
pub use newick::newick::parse_newick;
pub use node::{Node, FlatNode, FlatTree, TraversalOrder, node_to_flat, flat_to_node};