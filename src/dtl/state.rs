// Unified simulation state for both per-gene and per-species DTL models

use std::collections::HashMap;
use crate::node::{FlatNode, Event};
use super::event::DTLEvent;

/// Represents an active gene copy being simulated (per-gene model only)
#[derive(Clone, Debug)]
pub(crate) struct GeneCopy {
    pub gene_node_idx: usize,
    pub species_node_idx: usize,
    pub current_time: f64,
}

/// Holds the mutable state during DTL simulation.
///
/// Used by both per-gene and per-species simulation models.
/// The `genes_per_species` field is only populated in per-species mode.
pub(crate) struct SimulationState {
    pub gene_nodes: Vec<FlatNode>,
    pub node_mapping: Vec<usize>,
    pub event_mapping: Vec<Event>,
    pub events: Vec<DTLEvent>,
    /// Maps species_idx -> Vec of gene node indices in that species.
    /// Only used in per-species mode (None in per-gene mode).
    pub genes_per_species: Option<HashMap<usize, Vec<usize>>>,
}

impl SimulationState {
    /// Create state for per-gene mode (no gene-per-species tracking)
    pub fn new(capacity: usize) -> Self {
        Self {
            gene_nodes: Vec::with_capacity(capacity),
            node_mapping: Vec::with_capacity(capacity),
            event_mapping: Vec::with_capacity(capacity),
            events: Vec::with_capacity(capacity),
            genes_per_species: None,
        }
    }

    /// Create state for per-species mode (with gene-per-species tracking)
    pub fn new_per_species(capacity: usize) -> Self {
        Self {
            gene_nodes: Vec::with_capacity(capacity),
            node_mapping: Vec::with_capacity(capacity),
            event_mapping: Vec::with_capacity(capacity),
            events: Vec::with_capacity(capacity),
            genes_per_species: Some(HashMap::new()),
        }
    }

    /// Creates a new gene node and returns its index.
    pub fn create_gene_node(&mut self, parent: Option<usize>, species_idx: usize, event: Event, time: f64) -> usize {
        let idx = self.gene_nodes.len();
        self.gene_nodes.push(FlatNode {
            name: format!("{}_{}", species_idx, idx),
            left_child: None,
            right_child: None,
            parent,
            depth: Some(time),
            length: 0.0,
            bd_event: None,
        });
        self.node_mapping.push(species_idx);
        self.event_mapping.push(event);
        idx
    }

    /// Update gene's branch length and depth up to the given time
    fn update_gene_to_time(&mut self, gene_idx: usize, current_time: f64) {
        let gene_start = self.gene_nodes[gene_idx].depth.unwrap_or(0.0);
        self.gene_nodes[gene_idx].length = current_time - gene_start;
        self.gene_nodes[gene_idx].depth = Some(current_time);
    }

    /// Add gene to a species (per-species mode only, no-op otherwise)
    pub fn add_gene_to_species(&mut self, species_idx: usize, gene_idx: usize) {
        if let Some(ref mut map) = self.genes_per_species {
            map.entry(species_idx).or_default().push(gene_idx);
        }
    }

    /// Remove gene from a species (per-species mode only, no-op otherwise)
    fn remove_gene_from_species(&mut self, species_idx: usize, gene_idx: usize) {
        if let Some(ref mut map) = self.genes_per_species {
            if let Some(genes) = map.get_mut(&species_idx) {
                if let Some(pos) = genes.iter().position(|&g| g == gene_idx) {
                    genes.swap_remove(pos);
                }
            }
        }
    }

    /// Handles a duplication event: creates two children and records the event.
    pub fn handle_duplication(
        &mut self,
        parent_idx: usize,
        species_idx: usize,
        event_time: f64,
    ) -> (usize, usize) {
        self.update_gene_to_time(parent_idx, event_time);
        self.remove_gene_from_species(species_idx, parent_idx);
        self.event_mapping[parent_idx] = Event::Duplication;

        let child1_idx = self.create_gene_node(Some(parent_idx), species_idx, Event::Duplication, event_time);
        let child2_idx = self.create_gene_node(Some(parent_idx), species_idx, Event::Duplication, event_time);

        self.gene_nodes[parent_idx].left_child = Some(child1_idx);
        self.gene_nodes[parent_idx].right_child = Some(child2_idx);

        self.add_gene_to_species(species_idx, child1_idx);
        self.add_gene_to_species(species_idx, child2_idx);

        self.events.push(DTLEvent::Duplication {
            time: event_time,
            gene_id: parent_idx,
            species_id: species_idx,
            child1: child1_idx,
            child2: child2_idx,
        });

        (child1_idx, child2_idx)
    }

    /// Handles a transfer event: creates donor and recipient children.
    pub fn handle_transfer(
        &mut self,
        parent_idx: usize,
        donor_species: usize,
        recipient_species: usize,
        event_time: f64,
    ) -> (usize, usize) {
        self.update_gene_to_time(parent_idx, event_time);
        self.remove_gene_from_species(donor_species, parent_idx);
        self.event_mapping[parent_idx] = Event::Transfer;

        let donor_child_idx = self.create_gene_node(Some(parent_idx), donor_species, Event::Transfer, event_time);
        let recipient_child_idx = self.create_gene_node(Some(parent_idx), recipient_species, Event::Transfer, event_time);

        self.gene_nodes[parent_idx].left_child = Some(donor_child_idx);
        self.gene_nodes[parent_idx].right_child = Some(recipient_child_idx);

        self.add_gene_to_species(donor_species, donor_child_idx);
        self.add_gene_to_species(recipient_species, recipient_child_idx);

        self.events.push(DTLEvent::Transfer {
            time: event_time,
            gene_id: parent_idx,
            from_species: donor_species,
            to_species: recipient_species,
            donor_child: donor_child_idx,
            recipient_child: recipient_child_idx,
        });

        (donor_child_idx, recipient_child_idx)
    }

    /// Handles a loss event: marks the gene as lost.
    pub fn handle_loss(&mut self, gene_idx: usize, species_idx: usize, event_time: f64) {
        self.update_gene_to_time(gene_idx, event_time);
        self.remove_gene_from_species(species_idx, gene_idx);
        self.event_mapping[gene_idx] = Event::Loss;

        self.events.push(DTLEvent::Loss {
            time: event_time,
            gene_id: gene_idx,
            species_id: species_idx,
        });
    }

    /// Handles reaching a leaf species: gene survives as extant.
    pub fn handle_leaf(&mut self, gene_idx: usize, species_idx: usize, event_time: f64) {
        self.update_gene_to_time(gene_idx, event_time);
        self.event_mapping[gene_idx] = Event::Leaf;

        self.events.push(DTLEvent::Leaf {
            time: event_time,
            gene_id: gene_idx,
            species_id: species_idx,
        });
    }

    /// Handles a speciation event: gene follows both children species.
    pub fn handle_speciation(
        &mut self,
        parent_idx: usize,
        parent_species: usize,
        left_species: usize,
        right_species: usize,
        event_time: f64,
    ) -> (usize, usize) {
        self.update_gene_to_time(parent_idx, event_time);
        self.event_mapping[parent_idx] = Event::Speciation;

        let left_gene_idx = self.create_gene_node(Some(parent_idx), left_species, Event::Speciation, event_time);
        let right_gene_idx = self.create_gene_node(Some(parent_idx), right_species, Event::Speciation, event_time);

        self.gene_nodes[parent_idx].left_child = Some(left_gene_idx);
        self.gene_nodes[parent_idx].right_child = Some(right_gene_idx);

        self.add_gene_to_species(left_species, left_gene_idx);
        self.add_gene_to_species(right_species, right_gene_idx);

        self.events.push(DTLEvent::Speciation {
            time: event_time,
            gene_id: parent_idx,
            species_id: parent_species,
            left_child: left_gene_idx,
            right_child: right_gene_idx,
        });

        (left_gene_idx, right_gene_idx)
    }
}
