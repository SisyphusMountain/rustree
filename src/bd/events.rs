// Event generation functions for birth-death trees

use crate::node::FlatTree;

use super::types::{BDEvent, TreeEvent};

/// Generate birth-death events from an existing tree.
///
/// This function respects the `bd_event` field when present (from BD simulation),
/// allowing it to distinguish extinctions from extant leaves. For trees without
/// annotations (e.g., from Newick), it infers events from tree structure.
///
/// # Arguments
/// * `tree` - A FlatTree with depths assigned
///
/// # Returns
/// A vector of TreeEvent representing the tree's evolutionary history, sorted by time
///
/// # Panics
/// Panics if any node doesn't have a depth assigned
///
/// # Example
/// ```ignore
/// let tree = parse_newick("((A:1,B:1):1,C:2):0;").unwrap();
/// let mut flat_tree = tree.to_flat_tree();
/// flat_tree.assign_depths();
/// let events = generate_events_from_tree(&flat_tree);
/// ```
pub fn generate_events_from_tree(tree: &FlatTree) -> Vec<TreeEvent> {
    let mut events = Vec::with_capacity(tree.nodes.len());

    for (idx, node) in tree.nodes.iter().enumerate() {
        let depth = node.depth.expect("Tree must have depths assigned. Call assign_depths() first.");

        // Use bd_event annotation if available, otherwise infer from structure
        let event_type = if let Some(bd_event) = node.bd_event {
            // Use the annotation (can distinguish Extinction from Leaf)
            bd_event
        } else {
            // Fall back to structure-based inference
            if node.left_child.is_none() && node.right_child.is_none() {
                BDEvent::Leaf  // Assume extant if no annotation
            } else {
                BDEvent::Speciation
            }
        };

        events.push(TreeEvent {
            time: depth,
            node_id: idx,
            event_type,
            child1: node.left_child,
            child2: node.right_child,
        });
    }

    // Sort by time (depth) for consistency
    events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());

    events
}
