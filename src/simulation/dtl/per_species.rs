// Per-species DTL (Duplication-Transfer-Loss) simulation (Zombi-style)
//
// This module contains the per-species model where DTL event rates are
// proportional to the number of alive species, not the number of gene copies.
// Uses the shared Gillespie loop with PerSpecies mode.

use crate::bd::{TreeEvent, generate_events_from_tree};
use crate::node::{FlatTree, RecTree};
use rand::Rng;
use std::sync::Arc;

use super::event::DTLEvent;
use super::gillespie::{DTLMode, simulate_dtl_gillespie};
use super::utils::{validate_rates, precompute_lca, count_extant_genes};

/// Simulates a gene tree using Zombi-style per-species rates.
///
/// In this model, the DTL event rate is proportional to the number of ALIVE species,
/// NOT the number of gene copies. When an event is drawn, a random alive species is
/// chosen uniformly. If that species has no gene copies, the event "fails" (no action
/// is taken, but time advances). This means the effective event rate depends on
/// which species have genes, but the waiting times are drawn based on species count.
pub fn simulate_dtl_per_species<R: Rng>(
    species_tree: &FlatTree,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    replacement_transfer: Option<f64>,
    require_extant: bool,
    rng: &mut R,
) -> Result<(RecTree, Vec<DTLEvent>), String> {
    validate_rates(lambda_d, lambda_t, lambda_l)?;

    let species_arc = Arc::new(species_tree.clone());
    let species_events = generate_events_from_tree(species_tree)?;
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, transfer_alpha);

    loop {
        let (rec_tree, events) = simulate_dtl_per_species_internal(
            &species_arc,
            &species_events,
            &depths,
            &contemporaneity,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
            lca_depths.as_ref(),
            rng,
        )?;

        if !require_extant || count_extant_genes(&rec_tree) > 0 {
            return Ok((rec_tree, events));
        }
        // Retry if no extant genes
    }
}

/// Simulates multiple gene trees using the Zombi-style per-species DTL model.
pub fn simulate_dtl_per_species_batch<R: Rng>(
    species_tree: &FlatTree,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    replacement_transfer: Option<f64>,
    n_simulations: usize,
    require_extant: bool,
    rng: &mut R,
) -> Result<(Vec<RecTree>, Vec<Vec<DTLEvent>>), String> {
    validate_rates(lambda_d, lambda_t, lambda_l)?;

    let species_arc = Arc::new(species_tree.clone());
    let species_events = generate_events_from_tree(species_tree)?;
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, transfer_alpha);

    let mut rec_trees = Vec::with_capacity(n_simulations);
    let mut all_events = Vec::with_capacity(n_simulations);

    while rec_trees.len() < n_simulations {
        let (rec_tree, events) = simulate_dtl_per_species_internal(
            &species_arc,
            &species_events,
            &depths,
            &contemporaneity,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
            lca_depths.as_ref(),
            rng,
        )?;

        if !require_extant || count_extant_genes(&rec_tree) > 0 {
            rec_trees.push(rec_tree);
            all_events.push(events);
        }
        // Otherwise, discard and try again
    }

    Ok((rec_trees, all_events))
}

/// Internal implementation of per-species DTL simulation
pub fn simulate_dtl_per_species_internal<R: Rng>(
    species_tree: &Arc<FlatTree>,
    species_events: &[TreeEvent],
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    replacement_transfer: Option<f64>,
    lca_depths: Option<&Vec<Vec<f64>>>,
    rng: &mut R,
) -> Result<(RecTree, Vec<DTLEvent>), String> {
    simulate_dtl_gillespie(
        DTLMode::PerSpecies,
        species_tree,
        species_events,
        depths,
        contemporaneity,
        lca_depths.map(|v| v.as_slice()),
        origin_species,
        lambda_d,
        lambda_t,
        lambda_l,
        transfer_alpha,
        replacement_transfer,
        rng,
    )
}
