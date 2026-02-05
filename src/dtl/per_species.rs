// Per-species DTL (Duplication-Transfer-Loss) simulation (Zombi-style)
//
// This module contains the per-species model where DTL event rates are
// proportional to the number of alive species, not the number of gene copies.

use crate::bd::{BDEvent, TreeEvent, generate_events_from_tree};
use crate::node::{FlatTree, Event, RecTree};
use rand::Rng;

use super::event::DTLEvent;
use super::state::SimulationState;
use super::utils::{
    validate_rates, precompute_lca, find_time_index,
    select_transfer_recipient, select_transfer_recipient_assortative,
    count_extant_genes,
};

/// Simulates a gene tree using Zombi-style per-species rates.
///
/// In this model, the DTL event rate is proportional to the number of ALIVE species,
/// NOT the number of gene copies. When an event is drawn, a random alive species is
/// chosen uniformly. If that species has no gene copies, the event "fails" (no action
/// is taken, but time advances). This means the effective event rate depends on
/// which species have genes, but the waiting times are drawn based on species count.
pub fn simulate_dtl_per_species<'a, R: Rng>(
    species_tree: &'a FlatTree,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    require_extant: bool,
    rng: &mut R,
) -> (RecTree<'a>, Vec<DTLEvent>) {
    validate_rates(lambda_d, lambda_t, lambda_l);

    // Precompute species events, depths, contemporaneity, and LCA depths
    let species_events = generate_events_from_tree(species_tree);
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, transfer_alpha);

    loop {
        let (rec_tree, events) = simulate_dtl_per_species_internal(
            species_tree,
            &species_events,
            &depths,
            &contemporaneity,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            lca_depths.as_ref(),
            rng,
        );

        if !require_extant || count_extant_genes(&rec_tree) > 0 {
            return (rec_tree, events);
        }
        // Retry if no extant genes
    }
}

/// Simulates multiple gene trees using the Zombi-style per-species DTL model.
pub fn simulate_dtl_per_species_batch<'a, R: Rng>(
    species_tree: &'a FlatTree,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    n_simulations: usize,
    require_extant: bool,
    rng: &mut R,
) -> (Vec<RecTree<'a>>, Vec<Vec<DTLEvent>>) {
    validate_rates(lambda_d, lambda_t, lambda_l);

    // Precompute species events, depths, contemporaneity, and LCA depths once for all simulations
    let species_events = generate_events_from_tree(species_tree);
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, transfer_alpha);

    let mut rec_trees = Vec::with_capacity(n_simulations);
    let mut all_events = Vec::with_capacity(n_simulations);

    // Simulate all gene trees using the shared precomputed data
    while rec_trees.len() < n_simulations {
        let (rec_tree, events) = simulate_dtl_per_species_internal(
            species_tree,
            &species_events,
            &depths,
            &contemporaneity,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            lca_depths.as_ref(),
            rng,
        );

        if !require_extant || count_extant_genes(&rec_tree) > 0 {
            rec_trees.push(rec_tree);
            all_events.push(events);
        }
        // Otherwise, discard and try again
    }

    (rec_trees, all_events)
}

/// Internal implementation of per-species DTL simulation
pub fn simulate_dtl_per_species_internal<'a, R: Rng>(
    species_tree: &'a FlatTree,
    species_events: &[TreeEvent],
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    lca_depths: Option<&Vec<Vec<f64>>>,
    rng: &mut R,
) -> (RecTree<'a>, Vec<DTLEvent>) {
    let total_dtl_rate = lambda_d + lambda_t + lambda_l;
    let dup_threshold = if total_dtl_rate > 0.0 { lambda_d / total_dtl_rate } else { 0.0 };
    let transfer_threshold = if total_dtl_rate > 0.0 { (lambda_d + lambda_t) / total_dtl_rate } else { 0.0 };

    // Initialize state (per-species mode with gene tracking)
    let estimated_capacity = species_tree.nodes.len() * 4;
    let mut state = SimulationState::new_per_species(estimated_capacity);

    // Track current position in species events
    let mut species_event_idx = 0;

    // Find when origin species starts
    let origin_start_time = species_tree.nodes[origin_species].depth.unwrap_or(0.0)
        - species_tree.nodes[origin_species].length;

    let mut current_time = origin_start_time;

    // Skip species events that happen before origin time
    while species_event_idx < species_events.len() && species_events[species_event_idx].time < origin_start_time {
        species_event_idx += 1;
    }

    // Create initial gene at origin
    let initial_gene_idx = state.create_gene_node(None, origin_species, Event::Speciation, origin_start_time);
    state.add_gene_to_species(origin_species, initial_gene_idx);

    // Main simulation loop
    loop {
        // Find next species event time
        let next_species_event_time = if species_event_idx < species_events.len() {
            species_events[species_event_idx].time
        } else {
            f64::INFINITY
        };

        // Check if all genes are lost
        let gps = state.genes_per_species.as_ref().unwrap();
        if gps.values().map(|g| g.len()).sum::<usize>() == 0 {
            break;
        }

        // Get alive species at current time using contemporaneity
        let time_idx = find_time_index(depths, current_time);
        let alive_species = &contemporaneity[time_idx];
        let n_alive_species = alive_species.len();

        // Draw time to next DTL event
        let dtl_total_rate = n_alive_species as f64 * total_dtl_rate;
        let dtl_waiting_time = if dtl_total_rate > 0.0 {
            -rng.gen::<f64>().ln() / dtl_total_rate
        } else {
            f64::INFINITY
        };
        let next_dtl_time = current_time + dtl_waiting_time;

        // Check if simulation is complete (no more events possible)
        if next_species_event_time == f64::INFINITY && next_dtl_time == f64::INFINITY {
            break;
        }

        if next_species_event_time <= next_dtl_time && species_event_idx < species_events.len() {
            // Process species event first
            let sp_event = &species_events[species_event_idx].clone();
            current_time = sp_event.time;

            match sp_event.event_type {
                BDEvent::Speciation => {
                    let (child1, child2) = (sp_event.child1.unwrap(), sp_event.child2.unwrap());
                    let gps = state.genes_per_species.as_mut().unwrap();
                    if let Some(genes) = gps.remove(&sp_event.node_id) {
                        for gene_idx in genes {
                            state.handle_speciation(gene_idx, sp_event.node_id, child1, child2, current_time);
                        }
                    }
                }
                BDEvent::Extinction => {
                    let gps = state.genes_per_species.as_mut().unwrap();
                    if let Some(genes) = gps.remove(&sp_event.node_id) {
                        for gene_idx in genes {
                            state.handle_loss(gene_idx, sp_event.node_id, current_time);
                        }
                    }
                }
                BDEvent::Leaf => {
                    let gps = state.genes_per_species.as_mut().unwrap();
                    if let Some(genes) = gps.remove(&sp_event.node_id) {
                        for gene_idx in genes {
                            state.handle_leaf(gene_idx, sp_event.node_id, current_time);
                        }
                    }
                }
            }

            species_event_idx += 1;
        } else {
            // Process DTL event
            current_time = next_dtl_time;

            // Choose random species from ALL alive species (using contemporaneity)
            let affected_species = alive_species[rng.gen_range(0..n_alive_species)];

            // Check if this species has any genes - if not, event fails but time advances
            let gps = state.genes_per_species.as_ref().unwrap();
            if let Some(genes) = gps.get(&affected_species) {
                if !genes.is_empty() {
                    let gene_local_idx = rng.gen_range(0..genes.len());
                    let affected_gene = genes[gene_local_idx];

                    let event_prob: f64 = rng.gen();

                    if event_prob < dup_threshold {
                        state.handle_duplication(affected_gene, affected_species, current_time);
                    } else if event_prob < transfer_threshold {
                        let recipient = match (transfer_alpha, lca_depths) {
                            (Some(alpha), Some(lca)) => select_transfer_recipient_assortative(
                                depths, contemporaneity, lca, current_time, affected_species, alpha, rng
                            ),
                            _ => select_transfer_recipient(depths, contemporaneity, current_time, affected_species, rng),
                        };

                        if let Some(recipient_species) = recipient {
                            state.handle_transfer(affected_gene, affected_species, recipient_species, current_time);
                        }
                    } else {
                        state.handle_loss(affected_gene, affected_species, current_time);
                    }
                }
            }
        }

        // Check if simulation is complete
        let gps = state.genes_per_species.as_ref().unwrap();
        if species_event_idx >= species_events.len() && gps.values().map(|g| g.len()).sum::<usize>() == 0 {
            break;
        }
    }

    // Build final tree
    let root_idx = state.gene_nodes
        .iter()
        .position(|n| n.parent.is_none())
        .unwrap_or(0);

    let gene_tree = FlatTree {
        nodes: state.gene_nodes,
        root: root_idx,
    };

    (RecTree::new(species_tree, gene_tree, state.node_mapping, state.event_mapping), state.events)
}
