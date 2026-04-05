//! Birth-death process simulation for phylogenetic trees.
//!
//! Provides forward- and backward-time birth-death tree generation,
//! event extraction (speciation, extinction), and CSV export of
//! the resulting event sequences.

mod events;
mod simulation;
mod types;

// Re-export public API
pub use events::{generate_events_from_tree, generate_events_with_extinction};
pub use simulation::simulate_bd_tree_bwd;
pub use types::{BDEvent, TreeEvent};

// Re-export from io module for backward compatibility
pub use crate::io::save_bd_events_to_csv as save_events_to_csv;
