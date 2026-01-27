// rustree/src/lib.rs
#[macro_use]
extern crate pest_derive;


// Expose the modules
pub mod newick; // Contains newick parsing logic (newick.rs, newick.pest, etc.)
pub mod node;   // Contains the Node, FlatNode, and traversal logic
pub mod metric_functions;  // added metric functions
pub mod sampling; // added sampling functions
pub mod surgery; // added surgery functions
pub mod comparison; // added comparison functions
pub mod debug; // added debug functions
pub mod bd; // birth-death tree simulation
pub mod bd_optimized; // optimized birth-death tree simulation
pub mod dtl; // DTL (Duplication-Transfer-Loss) simulation

// Re-export key functions for easier access:
pub use newick::newick::parse_newick;
pub use node::{Node, FlatNode, FlatTree, TraversalOrder, node_to_flat, flat_to_node};