//! Reconciled tree structure for DTL (Duplication-Transfer-Loss) model.

use super::{FlatNode, FlatTree};
use crate::dtl::DTLEvent;
use crate::error::RustreeError;
use std::sync::Arc;

/// Events that can occur during DTL reconciliation.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Event {
    /// Speciation event - gene tree lineage follows species tree split
    Speciation,
    /// Duplication event - gene duplicates within a species
    Duplication,
    /// Horizontal gene transfer event
    Transfer,
    /// Gene loss event
    Loss,
    /// Leaf node (extant gene)
    Leaf,
}

/// Reconciled tree structure for DTL model.
///
/// Represents a gene tree reconciled with a species tree, storing the mapping
/// between gene tree nodes and species tree nodes, as well as the evolutionary
/// events (Duplication, Transfer, Loss) at each node.
///
/// The species tree is shared via `Arc<FlatTree>`, allowing multiple gene trees
/// to reference the same species tree without cloning.
#[derive(Clone, Debug)]
pub struct RecTree {
    /// The species tree (shared via Arc)
    pub species_tree: Arc<FlatTree>,
    /// The gene tree
    pub gene_tree: FlatTree,
    /// Maps each gene tree node index to its corresponding species tree node index.
    /// `None` means the mapping is unknown (e.g., after pruning/sampling).
    pub node_mapping: Vec<Option<usize>>,
    /// Maps each gene tree node index to its evolutionary event
    pub event_mapping: Vec<Event>,
    /// Detailed DTL events from simulation (if available).
    /// `None` for trees parsed from files or after sampling operations.
    pub dtl_events: Option<Vec<DTLEvent>>,
}

fn validate_mappings(
    species_tree: &FlatTree,
    gene_tree: &FlatTree,
    node_mapping: &[Option<usize>],
    event_mapping: &[Event],
) -> Result<(), RustreeError> {
    if gene_tree.nodes.len() != node_mapping.len() {
        return Err(RustreeError::Validation(format!(
            "node_mapping length {} must match gene_tree node count {}",
            node_mapping.len(),
            gene_tree.nodes.len()
        )));
    }
    if gene_tree.nodes.len() != event_mapping.len() {
        return Err(RustreeError::Validation(format!(
            "event_mapping length {} must match gene_tree node count {}",
            event_mapping.len(),
            gene_tree.nodes.len()
        )));
    }

    let sp_len = species_tree.nodes.len();
    for (i, mapping) in node_mapping.iter().enumerate() {
        if let Some(sp_idx) = mapping {
            if *sp_idx >= sp_len {
                return Err(RustreeError::Index(format!(
                    "node_mapping[{i}] = {sp_idx} is out of bounds for species tree with {sp_len} nodes"
                )));
            }
        }
    }

    Ok(())
}

impl RecTree {
    /// Creates a new RecTree with a shared species tree.
    ///
    /// # Panics
    /// Panics if the mappings don't match the gene tree size or if
    /// node_mapping contains indices out of bounds for the species tree.
    pub fn new(
        species_tree: Arc<FlatTree>,
        gene_tree: FlatTree,
        node_mapping: Vec<Option<usize>>,
        event_mapping: Vec<Event>,
    ) -> Self {
        Self::try_new(species_tree, gene_tree, node_mapping, event_mapping)
            .expect("invalid RecTree mappings")
    }

    /// Creates a new RecTree with checked reconciliation mappings.
    pub fn try_new(
        species_tree: Arc<FlatTree>,
        gene_tree: FlatTree,
        node_mapping: Vec<Option<usize>>,
        event_mapping: Vec<Event>,
    ) -> Result<Self, RustreeError> {
        validate_mappings(&species_tree, &gene_tree, &node_mapping, &event_mapping)?;

        Ok(RecTree {
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
            dtl_events: None,
        })
    }

    /// Creates a new RecTree with DTL events.
    ///
    /// # Panics
    /// Panics if the mappings don't match the gene tree size or if
    /// node_mapping contains indices out of bounds for the species tree.
    pub fn with_dtl_events(
        species_tree: Arc<FlatTree>,
        gene_tree: FlatTree,
        node_mapping: Vec<Option<usize>>,
        event_mapping: Vec<Event>,
        dtl_events: Vec<DTLEvent>,
    ) -> Self {
        Self::try_with_dtl_events(
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
            dtl_events,
        )
        .expect("invalid RecTree mappings")
    }

    /// Creates a new RecTree with checked reconciliation mappings and DTL events.
    pub fn try_with_dtl_events(
        species_tree: Arc<FlatTree>,
        gene_tree: FlatTree,
        node_mapping: Vec<Option<usize>>,
        event_mapping: Vec<Event>,
        dtl_events: Vec<DTLEvent>,
    ) -> Result<Self, RustreeError> {
        validate_mappings(&species_tree, &gene_tree, &node_mapping, &event_mapping)?;

        Ok(RecTree {
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
            dtl_events: Some(dtl_events),
        })
    }

    /// Creates a new RecTree, wrapping the species tree in an Arc.
    ///
    /// Convenience constructor for when you have an owned FlatTree.
    pub fn new_owned(
        species_tree: FlatTree,
        gene_tree: FlatTree,
        node_mapping: Vec<Option<usize>>,
        event_mapping: Vec<Event>,
    ) -> Self {
        Self::new(
            Arc::new(species_tree),
            gene_tree,
            node_mapping,
            event_mapping,
        )
    }

    /// Creates a new RecTree with DTL events, wrapping the species tree in an Arc.
    pub fn new_owned_with_events(
        species_tree: FlatTree,
        gene_tree: FlatTree,
        node_mapping: Vec<Option<usize>>,
        event_mapping: Vec<Event>,
        dtl_events: Vec<DTLEvent>,
    ) -> Self {
        Self::with_dtl_events(
            Arc::new(species_tree),
            gene_tree,
            node_mapping,
            event_mapping,
            dtl_events,
        )
    }

    /// Gets the species tree node index for a given gene tree node.
    /// Returns `None` inside `Ok` if the mapping is unknown (e.g., after pruning).
    ///
    /// # Errors
    /// Returns an error if `gene_node_idx` is out of bounds.
    pub fn species_node_for(&self, gene_node_idx: usize) -> Result<Option<usize>, RustreeError> {
        self.node_mapping
            .get(gene_node_idx)
            .copied()
            .ok_or_else(|| {
                RustreeError::Index(format!(
                    "gene_node_idx {} is out of bounds (gene tree has {} nodes)",
                    gene_node_idx,
                    self.node_mapping.len()
                ))
            })
    }

    /// Gets the event type for a given gene tree node.
    ///
    /// # Errors
    /// Returns an error if `gene_node_idx` is out of bounds.
    pub fn event_for(&self, gene_node_idx: usize) -> Result<&Event, RustreeError> {
        self.event_mapping.get(gene_node_idx).ok_or_else(|| {
            RustreeError::Index(format!(
                "gene_node_idx {} is out of bounds (gene tree has {} nodes)",
                gene_node_idx,
                self.event_mapping.len()
            ))
        })
    }

    /// Gets the gene tree node and its corresponding species node and event.
    ///
    /// # Errors
    /// Returns an error if `gene_node_idx` is out of bounds.
    pub fn get_full_info(
        &self,
        gene_node_idx: usize,
    ) -> Result<(&FlatNode, Option<usize>, &Event), RustreeError> {
        let gene_node = self.gene_tree.nodes.get(gene_node_idx).ok_or_else(|| {
            RustreeError::Index(format!(
                "gene_node_idx {} is out of bounds (gene tree has {} nodes)",
                gene_node_idx,
                self.gene_tree.nodes.len()
            ))
        })?;
        let species_idx = self
            .node_mapping
            .get(gene_node_idx)
            .copied()
            .ok_or_else(|| {
                RustreeError::Index(format!(
                    "gene_node_idx {} is out of bounds for node_mapping (len {})",
                    gene_node_idx,
                    self.node_mapping.len()
                ))
            })?;
        let event = self.event_mapping.get(gene_node_idx).ok_or_else(|| {
            RustreeError::Index(format!(
                "gene_node_idx {} is out of bounds for event_mapping (len {})",
                gene_node_idx,
                self.event_mapping.len()
            ))
        })?;
        Ok((gene_node, species_idx, event))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_node_tree(name: &str) -> FlatTree {
        FlatTree {
            nodes: vec![FlatNode {
                name: name.to_string(),
                left_child: None,
                right_child: None,
                parent: None,
                depth: Some(0.0),
                length: 0.0,
                bd_event: None,
            }],
            root: 0,
        }
    }

    #[test]
    fn try_new_rejects_mapping_length_mismatch() {
        let species_tree = Arc::new(one_node_tree("sp"));
        let gene_tree = one_node_tree("gene");
        let err =
            RecTree::try_new(species_tree, gene_tree, Vec::new(), vec![Event::Leaf]).unwrap_err();

        assert!(matches!(err, RustreeError::Validation(_)));
    }

    #[test]
    fn try_new_rejects_species_index_out_of_bounds() {
        let species_tree = Arc::new(one_node_tree("sp"));
        let gene_tree = one_node_tree("gene");
        let err = RecTree::try_new(species_tree, gene_tree, vec![Some(1)], vec![Event::Leaf])
            .unwrap_err();

        assert!(matches!(err, RustreeError::Index(_)));
    }
}
