#![allow(clippy::too_many_arguments)]
// Shared Gillespie-style DTL simulation loop
// Both per-gene and per-species DTL models use the same chronological
// event loop, differing only in how the total event rate is computed
// and how the affected gene copy is selected.

use crate::bd::{BDEvent, TreeEvent};
use crate::error::RustreeError;
use crate::node::{Event, FlatTree, RecTree};
use rand::Rng;
use std::sync::Arc;

use super::event::DTLEvent;
use super::state::SimulationState;
use super::utils::{
    finalize_simulation, find_time_index, select_transfer_recipient,
    select_transfer_recipient_assortative,
};
use super::DTLConfig;
use crate::simulation::utils::draw_waiting_time;

/// Determines how DTL event rates scale with the simulation state.
#[derive(Clone, Copy)]
pub(crate) enum DTLMode {
    /// Rate proportional to number of gene copies (per-gene model)
    PerGene,
    /// Rate proportional to number of alive species (per-species/Zombi model)
    PerSpecies,
}

impl DTLMode {
    fn total_event_rate(
        self,
        state: &SimulationState<'_>,
        depths: &[f64],
        contemporaneity: &[Vec<usize>],
        current_time: f64,
        config: &DTLConfig,
    ) -> f64 {
        match self {
            DTLMode::PerGene => {
                if config.uses_branch_rates() {
                    state.total_gene_weighted_rate(|species| config.branch_total_rate(species))
                } else {
                    state.total_gene_copies() as f64 * config.branch_total_rate(0)
                }
            }
            DTLMode::PerSpecies => {
                let time_idx = find_time_index(depths, current_time);
                if config.uses_branch_rates() {
                    contemporaneity[time_idx]
                        .iter()
                        .map(|&species| config.branch_total_rate(species))
                        .sum()
                } else {
                    contemporaneity[time_idx].len() as f64 * config.branch_total_rate(0)
                }
            }
        }
    }

    fn select_affected_gene<R: Rng>(
        self,
        state: &SimulationState<'_>,
        depths: &[f64],
        contemporaneity: &[Vec<usize>],
        current_time: f64,
        config: &DTLConfig,
        rng: &mut R,
    ) -> Option<(usize, usize)> {
        match self {
            DTLMode::PerGene => {
                if config.uses_branch_rates() {
                    state
                        .random_gene_copy_weighted(|species| config.branch_total_rate(species), rng)
                } else {
                    state.random_gene_copy(rng)
                }
            }
            DTLMode::PerSpecies => {
                let time_idx = find_time_index(depths, current_time);
                let alive_species = &contemporaneity[time_idx];
                if alive_species.is_empty() {
                    return None;
                }

                let species = if config.uses_branch_rates() {
                    let total_rate: f64 = alive_species
                        .iter()
                        .map(|&species| config.branch_total_rate(species))
                        .sum();
                    if total_rate <= 0.0 {
                        return None;
                    }

                    let mut threshold = rng.gen::<f64>() * total_rate;
                    let mut species = *alive_species.last()?;
                    for &candidate in alive_species {
                        let rate = config.branch_total_rate(candidate);
                        if rate <= 0.0 {
                            continue;
                        }
                        if threshold < rate {
                            species = candidate;
                            break;
                        }
                        threshold -= rate;
                    }
                    species
                } else {
                    alive_species[rng.gen_range(0..alive_species.len())]
                };

                match state.genes_per_species.get(&species) {
                    Some(genes) if !genes.is_empty() => {
                        let gene = genes[rng.gen_range(0..genes.len())];
                        Some((species, gene))
                    }
                    _ => None,
                }
            }
        }
    }
}

/// Shared Gillespie-style DTL simulation.
///
/// Both per-gene and per-species models use this function. The `mode` parameter
/// determines:
/// - **PerGene**: total rate = n_copies × (λ_D + λ_T + λ_L), random copy is selected
/// - **PerSpecies**: total rate = n_alive_species × (λ_D + λ_T + λ_L), random species
///   is selected then random gene in that species (event fails if species has no genes)
///
/// For now, rates are the same for all species in the tree. We can change it later if we want.
/// Doing so will involve changing the rates at which events are proposed,
/// as well as the probability of a given event affecting each gene or species at time t.
/// We can have auxiliary functions to flexibly handle the different generalizations without loss of performance.
pub(crate) fn simulate_dtl_gillespie<R: Rng>(
    mode: DTLMode,
    species_tree: &Arc<FlatTree>,
    species_events: &[TreeEvent],
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    lca_depths: Option<&[Vec<f64>]>,
    origin_species: usize,
    config: &DTLConfig,
    rng: &mut R,
) -> Result<(RecTree, Vec<DTLEvent>), RustreeError> {
    let transfer_alpha = config.transfer_alpha;
    let replacement_transfer = config.replacement_transfer;
    let origin_species = config.sample_origin(origin_species, rng);

    // Preallocate as many gene nodes as species nodes (for default case where there are no transfers, no duplications, no losses)
    let estimated_capacity = species_tree.nodes.len();
    let branch_total_rates = if config.uses_branch_rates() {
        Some(
            (0..species_tree.nodes.len())
                .map(|species| config.branch_total_rate(species))
                .collect::<Vec<_>>(),
        )
    } else {
        None
    };
    let mut state = SimulationState::with_branch_total_rates(
        estimated_capacity,
        species_tree,
        branch_total_rates.as_deref(),
    );

    // Find when origin species starts (beginning of its branch)
    let origin_node = species_tree.nodes.get(origin_species).ok_or_else(|| {
        RustreeError::Index(format!(
            "origin_species {origin_species} is out of bounds (species tree has {} nodes)",
            species_tree.nodes.len()
        ))
    })?;
    let origin_start_time = origin_node.depth.ok_or_else(|| {
        RustreeError::missing_depth("simulate_dtl_gillespie", origin_species, &origin_node.name)
    })? - origin_node.length;

    let mut current_time = origin_start_time;

    // Skip species events before origin time
    // The following assumes species_events are sorted by time, which should always be the case.
    let mut species_event_idx = 0;
    while species_event_idx < species_events.len()
        && species_events[species_event_idx].time < origin_start_time
    {
        species_event_idx += 1;
    }

    // For more flexibility, maybe just create a function handle_origination in
    // SimulationState that we can call here.
    // This would avoid hardcoding the logic, and would make the code more uniform.
    // Create initial gene at origin
    // We don't have an "Origination" like in ALE_dated or ALE_undated, so we just use a speciation event with no parent,
    // on the origination species.
    let initial_gene_idx =
        state.create_gene_node(None, origin_species, Event::Speciation, origin_start_time);
    state.add_gene_to_species(origin_species, initial_gene_idx);

    // Nudge past origin so the rate computation sees the stem interval
    // (same half-open boundary issue as speciations).
    current_time = current_time.next_up();

    // Main Gillespie loop
    loop {
        // We can enter here and have species_event_idx be out of bounds. In this case, next event time will be infinity.
        // what does it mean to be out of bounds?

        let next_species_event_time = species_events
            .get(species_event_idx)
            .map_or(f64::INFINITY, |e| e.time);

        // Check if all genes are lost
        if state.total_gene_copies() == 0 {
            break;
        }

        let dtl_total_rate =
            mode.total_event_rate(&state, depths, contemporaneity, current_time, config);

        let next_dtl_time = current_time + draw_waiting_time(dtl_total_rate, rng);

        if next_species_event_time == f64::INFINITY && next_dtl_time == f64::INFINITY {
            break;
        }
        // next DTL event was supposed to occur AFTER
        // the next species event, so we actually don't do the
        // DTL event and instead process a species-level event which will affect genes.
        if next_species_event_time <= next_dtl_time && species_event_idx < species_events.len() {
            // === Process species event ===
            let sp_event = species_events[species_event_idx].clone();
            current_time = sp_event.time;

            match sp_event.event_type {
                BDEvent::Speciation => {
                    let child1 = sp_event.child1.ok_or_else(|| {
                        RustreeError::Simulation(format!(
                            "Speciation event for node {} has no child1",
                            sp_event.node_id
                        ))
                    })?;
                    let child2 = sp_event.child2.ok_or_else(|| {
                        RustreeError::Simulation(format!(
                            "Speciation event for node {} has no child2",
                            sp_event.node_id
                        ))
                    })?;
                    if let Some(genes) = state.take_genes_for_species(sp_event.node_id) {
                        for gene_idx in genes {
                            state.handle_speciation(
                                gene_idx,
                                sp_event.node_id,
                                child1,
                                child2,
                                current_time,
                            );
                        }
                    }
                }
                BDEvent::Extinction => {
                    if let Some(genes) = state.take_genes_for_species(sp_event.node_id) {
                        for gene_idx in genes {
                            state.handle_loss(gene_idx, sp_event.node_id, current_time);
                        }
                    }
                }
                BDEvent::Leaf => {
                    if let Some(genes) = state.take_genes_for_species(sp_event.node_id) {
                        for gene_idx in genes {
                            state.handle_leaf(gene_idx, sp_event.node_id, current_time);
                        }
                    }
                }
            }

            species_event_idx += 1;

            // Nudge time past the event boundary so the next rate computation
            // uses post-event contemporaneity (children instead of parent after
            // speciations, excludes extinct/leaf species after those events).
            current_time = current_time.next_up();
        } else {
            // === Process DTL event ===
            current_time = next_dtl_time;

            let selection = mode.select_affected_gene(
                &state,
                depths,
                contemporaneity,
                current_time,
                config,
                rng,
            );

            let (affected_species, affected_gene) = match selection {
                Some(s) => s,
                None => continue, // Event failed, time still advances. This can happen in PerSpecies mode if we select a species with no genes, or if there are no alive species at this time.
            };

            // Determine event type
            let (lambda_d, lambda_t, lambda_l) = config.branch_event_rates(affected_species);
            let branch_total_rate = lambda_d + lambda_t + lambda_l;
            if branch_total_rate <= 0.0 {
                continue;
            }

            let event_draw = rng.gen::<f64>() * branch_total_rate;
            if event_draw < lambda_d {
                // Duplication
                state.handle_duplication(affected_gene, affected_species, current_time);
            } else if event_draw < lambda_d + lambda_t {
                // Transfer
                let recipient = match (transfer_alpha, lca_depths) {
                    (Some(alpha), Some(lca)) => select_transfer_recipient_assortative(
                        depths,
                        contemporaneity,
                        lca,
                        current_time,
                        affected_species,
                        alpha,
                        rng,
                    ),
                    _ => select_transfer_recipient(
                        depths,
                        contemporaneity,
                        current_time,
                        affected_species,
                        rng,
                    ),
                };

                if let Some(recipient_species) = recipient {
                    // Replacement transfer: find and remove victim BEFORE adding transfer
                    let is_replacement =
                        replacement_transfer.is_some_and(|p| p > 0.0 && rng.gen::<f64>() < p);
                    if is_replacement {
                        if let Some(victim) = state.random_gene_in_species(recipient_species, rng) {
                            state.handle_loss(victim, recipient_species, current_time);
                        }
                    }
                    state.handle_transfer(
                        affected_gene,
                        affected_species,
                        recipient_species,
                        current_time,
                    );
                }
            } else if lambda_l > 0.0 {
                // Loss
                state.handle_loss(affected_gene, affected_species, current_time);
            }
        }

        // Check if simulation is complete
        if species_event_idx >= species_events.len() && state.total_gene_copies() == 0 {
            break;
        }
    }

    // Finalize tree
    let rec_tree = finalize_simulation(
        species_tree.clone(),
        state.gene_nodes,
        state.node_mapping,
        state.event_mapping,
        origin_species,
        origin_start_time,
    )?;
    Ok((rec_tree, state.events))
}
