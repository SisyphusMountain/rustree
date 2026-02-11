// Per-gene DTL (Duplication-Transfer-Loss) simulation
//
// This module contains the per-gene-copy model where DTL event rates scale
// with the number of gene copies. Uses the shared Gillespie loop with PerGene mode.

use crate::bd::{TreeEvent, generate_events_from_tree};
use crate::node::{FlatTree, RecTreeOwned};
use rand::Rng;

use super::event::DTLEvent;
use super::gillespie::{DTLMode, simulate_dtl_gillespie};
use super::utils::{validate_rates, precompute_lca, count_extant_genes};

/// Simulates a gene tree along a species tree using the DTL model (internal version with pre-computed data)
pub fn simulate_dtl_internal<R: Rng>(
    species_tree: &FlatTree,
    species_events: &[TreeEvent],
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    lca_depths: Option<&[Vec<f64>]>,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    replacement_transfer: Option<f64>,
    rng: &mut R,
) -> Result<(RecTreeOwned, Vec<DTLEvent>), String> {
    simulate_dtl_gillespie(
        DTLMode::PerGene,
        species_tree,
        species_events,
        depths,
        contemporaneity,
        lca_depths,
        origin_species,
        lambda_d,
        lambda_t,
        lambda_l,
        transfer_alpha,
        replacement_transfer,
        rng,
    )
}

/// Simulates a gene tree along a species tree using the DTL model
///
/// This is the main public interface. It computes depths, contemporaneity,
/// and species events, then calls the internal simulation function.
pub fn simulate_dtl<R: Rng>(
    species_tree: &FlatTree,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Option<f64>,
    replacement_transfer: Option<f64>,
    require_extant: bool,
    rng: &mut R,
) -> Result<(RecTreeOwned, Vec<DTLEvent>), String> {
    validate_rates(lambda_d, lambda_t, lambda_l)?;

    let species_events = generate_events_from_tree(species_tree)?;
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, transfer_alpha);
    let lca_ref = lca_depths.as_ref().map(|v| v.as_slice());

    loop {
        let (rec_tree, events) = simulate_dtl_internal(
            species_tree,
            &species_events,
            &depths,
            &contemporaneity,
            lca_ref,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
            rng,
        )?;

        if !require_extant || count_extant_genes(&rec_tree) > 0 {
            return Ok((rec_tree, events));
        }
        // Otherwise, loop and try again
    }
}

/// Simulates multiple gene trees efficiently with shared pre-computed data
pub fn simulate_dtl_batch<R: Rng>(
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
) -> Result<(Vec<RecTreeOwned>, Vec<Vec<DTLEvent>>), String> {
    validate_rates(lambda_d, lambda_t, lambda_l)?;

    let species_events = generate_events_from_tree(species_tree)?;
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, transfer_alpha);
    let lca_ref = lca_depths.as_ref().map(|v| v.as_slice());

    let mut rec_trees = Vec::with_capacity(n_simulations);
    let mut all_events = Vec::with_capacity(n_simulations);

    while rec_trees.len() < n_simulations {
        let (rec_tree, events) = simulate_dtl_internal(
            species_tree,
            &species_events,
            &depths,
            &contemporaneity,
            lca_ref,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
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
