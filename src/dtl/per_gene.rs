// Per-gene DTL (Duplication-Transfer-Loss) simulation
//
// This module contains the per-gene-copy model where DTL event rates scale
// with the number of gene copies.

use crate::node::{FlatTree, Event, RecTree};
use rand::Rng;

use super::event::DTLEvent;
use super::state::{SimulationState, GeneCopy};
use super::utils::{
    validate_rates, precompute_lca, select_transfer_recipient,
    select_transfer_recipient_assortative, draw_waiting_time, finalize_simulation,
    count_extant_genes,
};

/// Simulates a gene tree along a species tree using the DTL model (internal version with pre-computed data)
pub fn simulate_dtl_internal<'a, R: Rng>(
    species_tree: &'a FlatTree,
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    lca_depths: Option<&[Vec<f64>]>,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    rng: &mut R,
) -> (RecTree<'a>, Vec<DTLEvent>) {
    let origin_depth = species_tree.nodes[origin_species]
        .depth
        .expect("Species tree must have depths assigned");

    // Initialize simulation state
    let estimated_capacity = species_tree.nodes.len() * 2;
    let mut state = SimulationState::new(estimated_capacity);

    // Create the initial gene copy at the origin
    let initial_gene_idx = state.create_gene_node(None, origin_species, Event::Speciation, origin_depth);

    // Active gene copies being simulated
    let mut active_copies = vec![GeneCopy {
        gene_node_idx: initial_gene_idx,
        current_time: origin_depth,
        species_node_idx: origin_species,
    }];

    // Total event rate and rate thresholds
    let total_rate = lambda_d + lambda_t + lambda_l;
    let dup_threshold = lambda_d / total_rate;
    let transfer_threshold = (lambda_d + lambda_t) / total_rate;

    // Reusable buffers
    let mut new_copies: Vec<GeneCopy> = Vec::with_capacity(estimated_capacity / 10);
    let mut copies_to_remove: Vec<usize> = Vec::with_capacity(estimated_capacity / 10);

    // Main simulation loop
    while !active_copies.is_empty() {
        new_copies.clear();
        copies_to_remove.clear();

        for (copy_idx, copy) in active_copies.iter().enumerate() {
            let species_node = &species_tree.nodes[copy.species_node_idx];
            let species_end_time = species_node.depth.unwrap();

            let mut current_time = copy.current_time;
            let current_gene_idx = copy.gene_node_idx;
            let mut terminated = false;

            // Simulate D/T/L events along the branch
            while current_time < species_end_time && !terminated {
                let event_time = current_time + draw_waiting_time(total_rate, rng);

                if event_time >= species_end_time {
                    // No event before branch end
                    current_time = species_end_time;
                } else {
                    current_time = event_time;
                    let event_prob: f64 = rng.gen();

                    if event_prob < dup_threshold {
                        // Duplication
                        let (child1, child2) = state.handle_duplication(
                            current_gene_idx, copy.species_node_idx, event_time
                        );
                        new_copies.push(GeneCopy {
                            gene_node_idx: child1,
                            species_node_idx: copy.species_node_idx,
                            current_time: event_time,
                        });
                        new_copies.push(GeneCopy {
                            gene_node_idx: child2,
                            species_node_idx: copy.species_node_idx,
                            current_time: event_time,
                        });
                        terminated = true;

                    } else if event_prob < transfer_threshold {
                        // Transfer - find valid recipient
                        let recipient = match (transfer_alpha, lca_depths) {
                            (Some(alpha), Some(lca)) => select_transfer_recipient_assortative(
                                depths, contemporaneity, lca, event_time, copy.species_node_idx, alpha, rng
                            ),
                            _ => select_transfer_recipient(depths, contemporaneity, event_time, copy.species_node_idx, rng),
                        };
                        if let Some(recipient) = recipient {
                            let (donor_child, recipient_child) = state.handle_transfer(
                                current_gene_idx, copy.species_node_idx, recipient, event_time
                            );
                            new_copies.push(GeneCopy {
                                gene_node_idx: donor_child,
                                species_node_idx: copy.species_node_idx,
                                current_time: event_time,
                            });
                            new_copies.push(GeneCopy {
                                gene_node_idx: recipient_child,
                                species_node_idx: recipient,
                                current_time: event_time,
                            });
                            terminated = true;
                        }
                        // If no valid recipient, gene continues (failed transfer)

                    } else {
                        // Loss
                        state.handle_loss(current_gene_idx, copy.species_node_idx, event_time);
                        terminated = true;
                    }
                }
            }

            // Process species tree node if gene survived the branch
            if !terminated {
                if species_node.left_child.is_none() && species_node.right_child.is_none() {
                    // Leaf: gene survives as extant
                    state.handle_leaf(current_gene_idx, copy.species_node_idx, species_end_time);
                } else {
                    // Speciation: gene follows both children
                    let left_species = species_node.left_child.unwrap();
                    let right_species = species_node.right_child.unwrap();

                    let (left_gene, right_gene) = state.handle_speciation(
                        current_gene_idx, copy.species_node_idx,
                        left_species, right_species, species_end_time
                    );

                    new_copies.push(GeneCopy {
                        gene_node_idx: left_gene,
                        species_node_idx: left_species,
                        current_time: species_end_time,
                    });
                    new_copies.push(GeneCopy {
                        gene_node_idx: right_gene,
                        species_node_idx: right_species,
                        current_time: species_end_time,
                    });
                }
            }

            copies_to_remove.push(copy_idx);
        }

        // Update active copies
        for &idx in copies_to_remove.iter().rev() {
            active_copies.swap_remove(idx);
        }
        active_copies.extend(new_copies.drain(..));
    }

    // Finalize and return
    let rec_tree = finalize_simulation(
        species_tree,
        state.gene_nodes,
        state.node_mapping,
        state.event_mapping,
        origin_species,
        origin_depth,
    );
    (rec_tree, state.events)
}

/// Simulates a gene tree along a species tree using the DTL model
///
/// This is the main public interface. It computes depths and contemporaneity,
/// then calls the internal simulation function.
pub fn simulate_dtl<'a, R: Rng>(
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

    // Compute depths and contemporaneity
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);

    // Compute LCA depths if using assortative transfers
    let lca_depths = precompute_lca(species_tree, transfer_alpha);
    let lca_ref = lca_depths.as_ref().map(|v| v.as_slice());

    loop {
        let (rec_tree, events) = simulate_dtl_internal(
            species_tree,
            &depths,
            &contemporaneity,
            lca_ref,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            rng,
        );

        if !require_extant || count_extant_genes(&rec_tree) > 0 {
            return (rec_tree, events);
        }
        // Otherwise, loop and try again
    }
}

/// Simulates multiple gene trees efficiently with shared pre-computed data
pub fn simulate_dtl_batch<'a, R: Rng>(
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

    // Compute depths and contemporaneity once for all simulations
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);

    // Compute LCA depths once if using assortative transfers
    let lca_depths = precompute_lca(species_tree, transfer_alpha);
    let lca_ref = lca_depths.as_ref().map(|v| v.as_slice());

    let mut rec_trees = Vec::with_capacity(n_simulations);
    let mut all_events = Vec::with_capacity(n_simulations);

    // Simulate all gene trees using the shared precomputed data
    while rec_trees.len() < n_simulations {
        let (rec_tree, events) = simulate_dtl_internal(
            species_tree,
            &depths,
            &contemporaneity,
            lca_ref,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
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
