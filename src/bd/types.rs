// birth-death process types

use crate::node::FlatTree;

/// Event types in a birth-death process
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BDEvent {
    /// Speciation event - lineage splits into two
    Speciation,
    /// Extinction event - lineage goes extinct
    Extinction,
    /// Leaf node - extant species at present time
    Leaf,
}

impl BDEvent {
    /// Convert to string representation for display/CSV
    pub fn as_str(&self) -> &'static str {
        match self {
            BDEvent::Speciation => "Speciation",
            BDEvent::Extinction => "Extinction",
            BDEvent::Leaf => "Leaf",
        }
    }

    /// Parse from string representation
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "Speciation" => Some(BDEvent::Speciation),
            "Extinction" => Some(BDEvent::Extinction),
            "Leaf" => Some(BDEvent::Leaf),
            _ => None,
        }
    }
}

/// Represents an event in the birth-death process
#[derive(Clone, Debug)]
pub struct TreeEvent {
    /// Time when the event occurred (going backwards from present at 0)
    pub time: f64,
    /// Node ID where the event occurred
    pub node_id: usize,
    /// Type of event
    pub event_type: BDEvent,
    /// First child node ID (for speciation events)
    pub child1: Option<usize>,
    /// Second child node ID (for speciation events)
    pub child2: Option<usize>,
}

impl TreeEvent {
    /// Convert event to CSV row format, resolving node names from the tree
    pub fn to_csv_row(&self, tree: &FlatTree) -> String {
        format!(
            "{},{},{},{},{}",
            self.time,
            tree.nodes[self.node_id].name,
            self.event_type.as_str(),
            self.child1.map_or(String::new(), |c| tree.nodes[c].name.clone()),
            self.child2.map_or(String::new(), |c| tree.nodes[c].name.clone())
        )
    }

    /// CSV header for event data
    pub fn csv_header() -> &'static str {
        "time,node_name,event_type,child1_name,child2_name"
    }
}
