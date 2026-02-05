// Utility functions for DTL simulation

use crate::bd::BDEvent;
use crate::node::{FlatTree, FlatNode, Event, RecTree};
use rand::Rng;

/// Validate that all DTL rates are non-negative
#[inline]
pub(crate) fn validate_rates(lambda_d: f64, lambda_t: f64, lambda_l: f64) {
    assert!(lambda_d >= 0.0, "Duplication rate must be non-negative");
    assert!(lambda_t >= 0.0, "Transfer rate must be non-negative");
    assert!(lambda_l >= 0.0, "Loss rate must be non-negative");
}

/// Precompute LCA depths if transfer_alpha is provided
#[inline]
pub(crate) fn precompute_lca(
    species_tree: &FlatTree,
    transfer_alpha: Option<f64>,
) -> Option<Vec<Vec<f64>>> {
    transfer_alpha.map(|_| species_tree.precompute_lca_depths())
}

// ============================================================================
// Unified Recipient Selection (eliminates duplication between per-gene/per-species)
// ============================================================================

/// Get the list of potential transfer recipients (excluding donor) at a given time
pub(crate) fn get_contemporaneous_recipients(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    event_time: f64,
    donor: usize,
) -> Vec<usize> {
    let time_idx = find_time_index(depths, event_time);
    let contemporaries = &contemporaneity[time_idx];
    contemporaries.iter()
        .filter(|&&sp| sp != donor)
        .copied()
        .collect()
}

/// Transfer recipient selection (uniform random)
pub(crate) fn select_transfer_recipient<R: Rng>(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    event_time: f64,
    donor: usize,
    rng: &mut R,
) -> Option<usize> {
    let recipients = get_contemporaneous_recipients(depths, contemporaneity, event_time, donor);
    if recipients.is_empty() {
        None
    } else {
        Some(recipients[rng.gen_range(0..recipients.len())])
    }
}

/// Transfer recipient selection with assortative preference (distance-weighted)
pub(crate) fn select_transfer_recipient_assortative<R: Rng>(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    lca_depths: &[Vec<f64>],
    event_time: f64,
    donor: usize,
    alpha: f64,
    rng: &mut R,
) -> Option<usize> {
    let recipients = get_contemporaneous_recipients(depths, contemporaneity, event_time, donor);
    if recipients.is_empty() {
        return None;
    }

    // Compute weights for each potential recipient
    let mut weights: Vec<f64> = Vec::with_capacity(recipients.len());
    let mut total_weight = 0.0;

    for &species in &recipients {
        let lca_depth = lca_depths[donor][species];
        let distance = 2.0 * (event_time - lca_depth);
        let weight = (-alpha * distance).exp();
        weights.push(weight);
        total_weight += weight;
    }

    if total_weight <= 0.0 {
        return None;
    }

    // Sample from weighted distribution
    let threshold: f64 = rng.gen::<f64>() * total_weight;
    let mut cumulative = 0.0;

    for (i, &weight) in weights.iter().enumerate() {
        cumulative += weight;
        if cumulative >= threshold {
            return Some(recipients[i]);
        }
    }

    // Fallback to last recipient
    Some(recipients[recipients.len() - 1])
}

/// Find the index in the depths array corresponding to a given time
pub(crate) fn find_time_index(depths: &[f64], time: f64) -> usize {
    match depths.binary_search_by(|probe| probe.partial_cmp(&time).unwrap()) {
        Ok(idx) => idx,
        Err(idx) => {
            if idx == 0 {
                0
            } else if idx >= depths.len() {
                depths.len() - 1
            } else {
                // Return the index of the interval containing this time
                idx
            }
        }
    }
}

/// Draws an exponential waiting time for the next event.
#[inline]
pub(crate) fn draw_waiting_time<R: Rng>(total_rate: f64, rng: &mut R) -> f64 {
    if total_rate > 0.0 {
        let u: f64 = rng.gen();
        -u.ln() / total_rate
    } else {
        f64::INFINITY
    }
}

/// Finalizes simulation by finding the root and handling edge cases.
pub(crate) fn finalize_simulation<'a>(
    species_tree: &'a FlatTree,
    gene_nodes: Vec<FlatNode>,
    node_mapping: Vec<usize>,
    mut event_mapping: Vec<Event>,
    origin_species: usize,
    origin_depth: f64,
) -> RecTree<'a> {
    let mut final_gene_nodes = gene_nodes;
    let mut final_node_mapping = node_mapping;

    // Handle edge case: gene lost before any events
    // Node name format: {species_idx}_{gene_idx}
    if final_gene_nodes.is_empty() {
        final_gene_nodes.push(FlatNode {
            name: format!("{}_0", origin_species),
            left_child: None,
            right_child: None,
            parent: None,
            depth: Some(origin_depth),
            length: 0.0,
            bd_event: None,
        });
        final_node_mapping.push(origin_species);
        event_mapping.push(Event::Loss);
    }

    let root_idx = final_gene_nodes
        .iter()
        .position(|n| n.parent.is_none())
        .unwrap_or(0);

    let gene_tree = FlatTree {
        nodes: final_gene_nodes,
        root: root_idx,
    };

    RecTree::new(species_tree, gene_tree, final_node_mapping, event_mapping)
}

/// Counts the number of extant genes (gene leaves mapped to extant species leaves).
///
/// A gene is considered extant if:
/// 1. It is a gene leaf (event_mapping is Event::Leaf)
/// 2. It is mapped to an extant species leaf:
///    - If bd_event is set: must be BDEvent::Leaf (not Extinction)
///    - If bd_event is None (e.g., tree from Newick): falls back to checking if species node is a leaf (no children)
///
/// This excludes gene leaves that are mapped to extinct species (extinction events).
pub fn count_extant_genes(rec_tree: &RecTree) -> usize {
    rec_tree
        .gene_tree
        .nodes
        .iter()
        .enumerate()
        .filter(|(i, _node)| {
            // Check if this gene node is a Leaf event
            if rec_tree.event_mapping[*i] != Event::Leaf {
                return false;
            }
            // Check if the corresponding species node is an extant leaf (not an extinction)
            let species_idx = rec_tree.node_mapping[*i];
            let species_node = &rec_tree.species_tree.nodes[species_idx];

            // If bd_event is set, check for Leaf (not Extinction)
            // If bd_event is None (tree from Newick), fall back to checking no children
            match species_node.bd_event {
                Some(BDEvent::Leaf) => true,
                Some(BDEvent::Extinction) => false,
                Some(BDEvent::Speciation) => false, // Internal nodes aren't leaves
                None => {
                    // Fallback for trees without bd_event: check if node has no children
                    species_node.left_child.is_none() && species_node.right_child.is_none()
                }
            }
        })
        .count()
}

/// Counts events by type in a RecTree
pub fn count_events(rec_tree: &RecTree) -> (usize, usize, usize, usize, usize) {
    let mut speciations = 0;
    let mut duplications = 0;
    let mut transfers = 0;
    let mut losses = 0;
    let mut leaves = 0;

    for event in &rec_tree.event_mapping {
        match event {
            Event::Speciation => speciations += 1,
            Event::Duplication => duplications += 1,
            Event::Transfer => transfers += 1,
            Event::Loss => losses += 1,
            Event::Leaf => leaves += 1,
        }
    }

    (speciations, duplications, transfers, losses, leaves)
}
