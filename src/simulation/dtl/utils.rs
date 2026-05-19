// Utility functions for DTL simulation

use crate::bd::{generate_events_from_tree, BDEvent, TreeEvent};
use crate::error::RustreeError;
use crate::node::{Event, FlatNode, FlatTree, RecTree};
use crate::simulation::utils::draw_waiting_time;
use rand::Rng;
use std::sync::Arc;

use super::event::DTLEvent;
use super::gillespie::DTLMode;
use super::state::SimulationState;
use super::DTLConfig;

/// Shared precomputed state used by DTL iterators.
pub(crate) struct PreparedDtlSimulation {
    pub species_tree: Arc<FlatTree>,
    pub species_events: Vec<TreeEvent>,
    pub depths: Vec<f64>,
    pub contemporaneity: Vec<Vec<usize>>,
    pub lca_depths: Option<Vec<Vec<f64>>>,
    pub runtime: PreparedDtlRuntime,
}

/// Precomputed data used by repeated Gillespie simulations.
pub(crate) struct PreparedDtlRuntime {
    default_origin_start_time: f64,
    branch_start_times: Option<Vec<f64>>,
    branch_total_rates: Option<Vec<f64>>,
    per_species_interval_total_rates: Option<Vec<f64>>,
}

impl PreparedDtlRuntime {
    fn from_parts(
        default_origin_start_time: f64,
        branch_start_times: Option<Vec<f64>>,
        config: &DTLConfig,
        contemporaneity: &[Vec<usize>],
    ) -> Self {
        let branch_total_rates = config.branch_rates.as_ref().map(|rates| {
            (0..rates.lambda_d.len())
                .map(|species| rates.total_rate(species))
                .collect::<Vec<_>>()
        });

        let per_species_interval_total_rates = branch_total_rates.as_ref().map(|rates| {
            contemporaneity
                .iter()
                .map(|species| species.iter().map(|&sp| rates[sp]).sum())
                .collect()
        });

        Self {
            default_origin_start_time,
            branch_start_times,
            branch_total_rates,
            per_species_interval_total_rates,
        }
    }

    #[inline]
    fn origin_start_time(&self, origin_species: usize) -> f64 {
        self.branch_start_times
            .as_ref()
            .map_or(self.default_origin_start_time, |starts| {
                starts[origin_species]
            })
    }
}

/// Prepare shared, reusable DTL simulation state.
pub(crate) fn prepare_simulation(
    species_tree: &FlatTree,
    origin_species: usize,
    config: &DTLConfig,
) -> Result<PreparedDtlSimulation, RustreeError> {
    const OPERATION: &str = "prepare_dtl_simulation";

    config.validate_for_tree(species_tree.nodes.len())?;

    let species_events = generate_events_from_tree(species_tree)?;

    // Include possible origin start times in depths so stem periods are covered
    // by contemporaneity (make_subdivision only has node depths, not branch starts).
    let origin_start_time = branch_start_time(species_tree, origin_species, OPERATION)?;
    let branch_start_times = if config.uses_branch_rates() {
        Some(precompute_branch_start_times(species_tree, OPERATION)?)
    } else {
        None
    };

    let mut depths = species_tree.make_subdivision();
    if let Some(starts) = branch_start_times.as_ref() {
        depths.extend(starts.iter().copied());
    } else {
        depths.push(origin_start_time);
    }
    depths.sort_by(|a, b| a.total_cmp(b));
    depths.dedup();

    let contemporaneity = species_tree.find_contemporaneity(&depths);
    let lca_depths = precompute_lca(species_tree, config.transfer_alpha)?;
    let runtime = PreparedDtlRuntime::from_parts(
        origin_start_time,
        branch_start_times,
        config,
        &contemporaneity,
    );

    Ok(PreparedDtlSimulation {
        species_tree: Arc::new(species_tree.clone()),
        species_events,
        depths,
        contemporaneity,
        lca_depths,
        runtime,
    })
}

fn branch_start_time(
    species_tree: &FlatTree,
    species_idx: usize,
    operation: &'static str,
) -> Result<f64, RustreeError> {
    let node = species_tree.nodes.get(species_idx).ok_or_else(|| {
        RustreeError::Index(format!(
            "origin_species {species_idx} is out of bounds (species tree has {} nodes)",
            species_tree.nodes.len()
        ))
    })?;

    let depth = node
        .depth
        .ok_or_else(|| RustreeError::missing_depth(operation, species_idx, &node.name))?;
    if !depth.is_finite() {
        return Err(RustreeError::invalid_depth(
            operation,
            species_idx,
            &node.name,
            depth,
        ));
    }
    if !node.length.is_finite() {
        return Err(RustreeError::invalid_length(
            operation,
            species_idx,
            &node.name,
            node.length,
        ));
    }

    let start = depth - node.length;
    if !start.is_finite() {
        return Err(RustreeError::InvalidDepth {
            operation,
            node_index: species_idx,
            node_name: node.name.clone(),
            depth: start,
        });
    }

    Ok(start)
}

fn precompute_branch_start_times(
    species_tree: &FlatTree,
    operation: &'static str,
) -> Result<Vec<f64>, RustreeError> {
    (0..species_tree.nodes.len())
        .map(|idx| branch_start_time(species_tree, idx, operation))
        .collect()
}

/// Precompute LCA depths if transfer_alpha is provided
#[inline]
pub(crate) fn precompute_lca(
    species_tree: &FlatTree,
    transfer_alpha: Option<f64>,
) -> Result<Option<Vec<Vec<f64>>>, RustreeError> {
    transfer_alpha
        .map(|_| {
            species_tree
                .precompute_lca_depths()
                .map_err(RustreeError::Simulation)
        })
        .transpose()
}

fn total_event_rate(
    mode: DTLMode,
    state: &SimulationState<'_>,
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    current_time: f64,
    config: &DTLConfig,
    runtime: &PreparedDtlRuntime,
) -> f64 {
    match mode {
        DTLMode::PerGene => {
            if let Some(branch_total_rates) = runtime.branch_total_rates.as_ref() {
                state.total_gene_weighted_rate(|species| branch_total_rates[species])
            } else {
                state.total_gene_copies() as f64 * config.branch_total_rate(0)
            }
        }
        DTLMode::PerSpecies => {
            let time_idx = find_time_index(depths, current_time);
            if let Some(interval_rates) = runtime.per_species_interval_total_rates.as_ref() {
                interval_rates[time_idx]
            } else {
                contemporaneity[time_idx].len() as f64 * config.branch_total_rate(0)
            }
        }
    }
}

fn select_affected_gene<R: Rng>(
    mode: DTLMode,
    state: &SimulationState<'_>,
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    current_time: f64,
    runtime: &PreparedDtlRuntime,
    rng: &mut R,
) -> Option<(usize, usize)> {
    match mode {
        DTLMode::PerGene => {
            if let Some(branch_total_rates) = runtime.branch_total_rates.as_ref() {
                state.random_gene_copy_weighted(|species| branch_total_rates[species], rng)
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

            let species = if let Some(branch_total_rates) = runtime.branch_total_rates.as_ref() {
                let total_rate = runtime
                    .per_species_interval_total_rates
                    .as_ref()
                    .map_or_else(
                        || {
                            alive_species
                                .iter()
                                .map(|&sp| branch_total_rates[sp])
                                .sum::<f64>()
                        },
                        |rates| rates[time_idx],
                    );
                if total_rate <= 0.0 {
                    return None;
                }

                let mut threshold = rng.gen::<f64>() * total_rate;
                let mut species = *alive_species.last()?;
                for &candidate in alive_species {
                    let rate = branch_total_rates[candidate];
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

/// Shared Gillespie-style DTL simulation using prepared runtime data.
#[allow(clippy::too_many_arguments)]
pub(crate) fn simulate_dtl_gillespie_prepared<R: Rng>(
    mode: DTLMode,
    species_tree: &Arc<FlatTree>,
    species_events: &[TreeEvent],
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    lca_depths: Option<&[Vec<f64>]>,
    origin_species: usize,
    config: &DTLConfig,
    runtime: &PreparedDtlRuntime,
    rng: &mut R,
) -> Result<(RecTree, Vec<DTLEvent>), RustreeError> {
    let transfer_alpha = config.transfer_alpha;
    let replacement_transfer = config.replacement_transfer;
    let origin_species = config.sample_origin(origin_species, rng);

    let estimated_capacity = species_tree.nodes.len();
    let mut state = SimulationState::with_branch_total_rates(
        estimated_capacity,
        species_tree,
        runtime.branch_total_rates.as_deref(),
    );

    let origin_start_time = runtime.origin_start_time(origin_species);
    let mut current_time = origin_start_time;

    let mut species_event_idx = 0;
    while species_event_idx < species_events.len()
        && species_events[species_event_idx].time < origin_start_time
    {
        species_event_idx += 1;
    }

    let initial_gene_idx =
        state.create_gene_node(None, origin_species, Event::Speciation, origin_start_time);
    state.add_gene_to_species(origin_species, initial_gene_idx);

    current_time = current_time.next_up();

    loop {
        let next_species_event_time = species_events
            .get(species_event_idx)
            .map_or(f64::INFINITY, |e| e.time);

        if state.total_gene_copies() == 0 {
            break;
        }

        let dtl_total_rate = total_event_rate(
            mode,
            &state,
            depths,
            contemporaneity,
            current_time,
            config,
            runtime,
        );

        let next_dtl_time = current_time + draw_waiting_time(dtl_total_rate, rng);

        if next_species_event_time == f64::INFINITY && next_dtl_time == f64::INFINITY {
            break;
        }
        if next_species_event_time <= next_dtl_time && species_event_idx < species_events.len() {
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
            current_time = current_time.next_up();
        } else {
            current_time = next_dtl_time;

            let selection = select_affected_gene(
                mode,
                &state,
                depths,
                contemporaneity,
                current_time,
                runtime,
                rng,
            );

            let (affected_species, affected_gene) = match selection {
                Some(s) => s,
                None => continue,
            };

            let (lambda_d, lambda_t, lambda_l) = config.branch_event_rates(affected_species);
            let branch_total_rate = lambda_d + lambda_t + lambda_l;
            if branch_total_rate <= 0.0 {
                continue;
            }

            let event_draw = rng.gen::<f64>() * branch_total_rate;
            if event_draw < lambda_d {
                state.handle_duplication(affected_gene, affected_species, current_time);
            } else if event_draw < lambda_d + lambda_t {
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
                state.handle_loss(affected_gene, affected_species, current_time);
            }
        }

        if species_event_idx >= species_events.len() && state.total_gene_copies() == 0 {
            break;
        }
    }

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

// ============================================================================
// Unified Recipient Selection (eliminates duplication between per-gene/per-species)
// ============================================================================

/// Returns the list of potential horizontal transfer recipients at a given time,
/// excluding the donor species.
///
/// This is the main entry point for determining which species can receive a
/// horizontal gene transfer event. It works by:
///
/// 1. Calling [`find_time_index`] to map `event_time` to the correct time
///    subdivision index `j`.
/// 2. Looking up `contemporaneity[j]` to get all species alive during that interval.
/// 3. Filtering out the `donor` species (a lineage cannot transfer to itself).
///
/// The result is used by [`select_transfer_recipient`] (uniform random selection)
/// and [`select_transfer_recipient_assortative`] (distance-weighted selection).
///
/// # Arguments
///
/// * `depths` - Sorted time subdivision boundaries from [`FlatTree::make_subdivision`](crate::node::FlatTree::make_subdivision).
/// * `contemporaneity` - Per-interval species lists from [`FlatTree::find_contemporaneity`](crate::node::FlatTree::find_contemporaneity).
/// * `event_time` - The time at which the transfer event occurs.
/// * `donor` - The index of the donor species to exclude from recipients.
///
/// # Returns
///
/// A vector of species indices that are alive at `event_time` and are not the `donor`.
/// May be empty if the donor is the only species alive at that time.
#[cfg(test)]
fn get_contemporaneous_recipients(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    event_time: f64,
    donor: usize,
) -> Vec<usize> {
    let time_idx = find_time_index(depths, event_time);
    let contemporaries = &contemporaneity[time_idx];
    contemporaries
        .iter()
        .filter(|&&sp| sp != donor)
        .copied()
        .collect()
}

/// Transfer recipient selection (uniform random).
///
/// Selects directly from the contemporaneity slice without allocating.
pub(crate) fn select_transfer_recipient<R: Rng>(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    event_time: f64,
    donor: usize,
    rng: &mut R,
) -> Option<usize> {
    let time_idx = find_time_index(depths, event_time);
    let contemporaries = &contemporaneity[time_idx];

    // Count eligible recipients (all contemporaries except donor)
    let count = contemporaries.iter().filter(|&&sp| sp != donor).count();
    if count == 0 {
        return None;
    }

    // Pick the k-th eligible recipient (0-indexed)
    let k = rng.gen_range(0..count);
    contemporaries
        .iter()
        .filter(|&&sp| sp != donor)
        .nth(k)
        .copied()
}

/// Transfer recipient selection with assortative preference (distance-weighted).
///
/// Selects directly from the contemporaneity slice without allocating a recipient Vec.
pub(crate) fn select_transfer_recipient_assortative<R: Rng>(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    lca_depths: &[Vec<f64>],
    event_time: f64,
    donor: usize,
    alpha: f64,
    rng: &mut R,
) -> Option<usize> {
    let time_idx = find_time_index(depths, event_time);
    let contemporaries = &contemporaneity[time_idx];

    // Compute total weight over eligible recipients
    let mut total_weight = 0.0;
    for &sp in contemporaries {
        if sp != donor {
            let lca_depth = lca_depths[donor][sp];
            let distance = 2.0 * (event_time - lca_depth);
            total_weight += (-alpha * distance).exp();
        }
    }

    if total_weight <= 0.0 {
        return None;
    }

    // Sample from weighted distribution
    let threshold: f64 = rng.gen::<f64>() * total_weight;
    let mut cumulative = 0.0;
    let mut last_eligible = None;

    for &sp in contemporaries {
        if sp != donor {
            let lca_depth = lca_depths[donor][sp];
            let distance = 2.0 * (event_time - lca_depth);
            cumulative += (-alpha * distance).exp();
            last_eligible = Some(sp);
            if cumulative >= threshold {
                return Some(sp);
            }
        }
    }

    // Fallback to last eligible recipient
    last_eligible
}

/// Maps a time value to the index of the time subdivision interval containing it.
///
/// Given a sorted `depths` array representing time subdivision boundaries, this function
/// uses binary search to find the index `j` such that `depths[j-1] < time <= depths[j]`.
/// In other words, it returns the index of the interval whose upper bound is at or just
/// after `time`.
///
/// # Convention
///
/// The returned index `j` satisfies `depths[j-1] < time <= depths[j]`, which is consistent
/// with how [`FlatTree::find_contemporaneity`](crate::node::FlatTree::find_contemporaneity)
/// indexes its data: `contemporaneity[j]` contains the species alive during the interval
/// `(depths[j-1], depths[j]]`. Therefore:
///
/// ```text
/// contemporaneity[find_time_index(depths, t)]
/// ```
///
/// gives the set of species alive at time `t`.
///
/// # Boundary behavior
///
/// - **Exact match**: if `time` equals `depths[j]`, returns `j` directly.
/// - **Before first depth**: if `time < depths[0]`, returns `0`.
/// - **After last depth**: if `time > depths[depths.len() - 1]`, returns `depths.len() - 1`.
///
/// # Difference from `find_closest_index`
///
/// Unlike [`FlatTree::find_closest_index`](crate::node::FlatTree) in `metric_functions.rs`,
/// which uses nearest-neighbor matching (picking whichever of `depths[j-1]` or `depths[j]`
/// is closer to `time`), this function always returns the higher index `j` when `time` falls
/// between two depths. This "ceiling" behavior is the correct choice for looking up which
/// species are alive at a given time in the contemporaneity structure.
///
/// # Arguments
///
/// * `depths` - A sorted slice of time subdivision boundaries.
/// * `time` - The time value to locate within the subdivision.
///
/// # Returns
///
/// The index into `depths` (and `contemporaneity`) corresponding to the interval containing `time`.
pub(crate) fn find_time_index(depths: &[f64], time: f64) -> usize {
    match depths.binary_search_by(|probe| probe.total_cmp(&time)) {
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

/// Finalizes simulation by finding the root and handling edge cases.
pub(crate) fn finalize_simulation(
    species_tree: Arc<FlatTree>,
    gene_nodes: Vec<FlatNode>,
    node_mapping: Vec<Option<usize>>,
    mut event_mapping: Vec<Event>,
    origin_species: usize,
    origin_depth: f64,
) -> Result<RecTree, RustreeError> {
    let mut final_gene_nodes = gene_nodes;
    let mut final_node_mapping = node_mapping;

    // Handle edge case: gene lost before any events
    // Node name format: {species_name}_{gene_idx}
    if final_gene_nodes.is_empty() {
        let sp_name = &species_tree.nodes[origin_species].name;
        final_gene_nodes.push(FlatNode {
            name: format!("{}_0", sp_name),
            left_child: None,
            right_child: None,
            parent: None,
            depth: Some(origin_depth),
            length: 0.0,
            bd_event: None,
        });
        final_node_mapping.push(Some(origin_species));
        event_mapping.push(Event::Loss);
    }

    let root_idx = final_gene_nodes
        .iter()
        .position(|n| n.parent.is_none())
        .ok_or_else(|| {
            RustreeError::Simulation("gene tree must have at least one root node".to_string())
        })?;

    let gene_tree = FlatTree {
        nodes: final_gene_nodes,
        root: root_idx,
    };

    RecTree::try_new(species_tree, gene_tree, final_node_mapping, event_mapping)
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
            let species_idx = match rec_tree.node_mapping[*i] {
                Some(idx) => idx,
                None => return false, // unmapped node cannot be extant
            };
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: standard subdivision used across all find_time_index tests.
    fn test_depths() -> Vec<f64> {
        vec![0.0, 1.0, 2.0, 3.0, 5.0]
    }

    // ---------------------------------------------------------------
    // find_time_index: exact matches
    // ---------------------------------------------------------------

    #[test]
    fn find_time_index_exact_match_at_interior_point() {
        let depths = test_depths();
        // time == depths[2] => should return 2
        assert_eq!(find_time_index(&depths, 2.0), 2);
    }

    #[test]
    fn find_time_index_exact_match_at_each_depth() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, 0.0), 0);
        assert_eq!(find_time_index(&depths, 1.0), 1);
        assert_eq!(find_time_index(&depths, 2.0), 2);
        assert_eq!(find_time_index(&depths, 3.0), 3);
        assert_eq!(find_time_index(&depths, 5.0), 4);
    }

    // ---------------------------------------------------------------
    // find_time_index: between points (higher-index / ceiling behavior)
    // ---------------------------------------------------------------

    #[test]
    fn find_time_index_between_points_returns_higher_index() {
        let depths = test_depths();
        // 1.5 is between depths[1]=1.0 and depths[2]=2.0 => should return 2
        assert_eq!(find_time_index(&depths, 1.5), 2);
    }

    #[test]
    fn find_time_index_between_various_points() {
        let depths = test_depths();
        // Between depths[0]=0.0 and depths[1]=1.0
        assert_eq!(find_time_index(&depths, 0.5), 1);
        // Between depths[2]=2.0 and depths[3]=3.0
        assert_eq!(find_time_index(&depths, 2.5), 3);
        // Between depths[3]=3.0 and depths[4]=5.0
        assert_eq!(find_time_index(&depths, 4.0), 4);
    }

    // ---------------------------------------------------------------
    // find_time_index: boundary clamping
    // ---------------------------------------------------------------

    #[test]
    fn find_time_index_before_first_depth_returns_zero() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, -0.5), 0);
    }

    #[test]
    fn find_time_index_well_before_first_depth() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, -100.0), 0);
    }

    #[test]
    fn find_time_index_after_last_depth_returns_last_index() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, 10.0), 4);
    }

    #[test]
    fn find_time_index_well_after_last_depth() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, 1_000.0), 4);
    }

    // ---------------------------------------------------------------
    // find_time_index: first and last boundary exact matches
    // ---------------------------------------------------------------

    #[test]
    fn find_time_index_at_first_boundary() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, 0.0), 0);
    }

    #[test]
    fn find_time_index_at_last_boundary() {
        let depths = test_depths();
        assert_eq!(find_time_index(&depths, 5.0), 4);
    }

    // ---------------------------------------------------------------
    // find_time_index: edge cases with small arrays
    // ---------------------------------------------------------------

    #[test]
    fn find_time_index_single_element() {
        let depths = vec![3.0];
        assert_eq!(find_time_index(&depths, 3.0), 0);
        assert_eq!(find_time_index(&depths, 1.0), 0);
        assert_eq!(find_time_index(&depths, 5.0), 0);
    }

    #[test]
    fn find_time_index_two_elements() {
        let depths = vec![1.0, 4.0];
        assert_eq!(find_time_index(&depths, 1.0), 0);
        assert_eq!(find_time_index(&depths, 4.0), 1);
        assert_eq!(find_time_index(&depths, 2.5), 1);
        assert_eq!(find_time_index(&depths, 0.0), 0);
        assert_eq!(find_time_index(&depths, 5.0), 1);
    }

    // ---------------------------------------------------------------
    // find_time_index: values very close to boundaries
    // ---------------------------------------------------------------

    #[test]
    fn find_time_index_just_above_depth() {
        let depths = test_depths();
        // Just above depths[1]=1.0 should return index 2
        assert_eq!(find_time_index(&depths, 1.0 + 1e-15), 2);
    }

    #[test]
    fn find_time_index_just_below_depth() {
        let depths = test_depths();
        // Just below depths[2]=2.0 should return index 2 (ceiling behavior)
        assert_eq!(find_time_index(&depths, 2.0 - 1e-15), 2);
    }

    // ---------------------------------------------------------------
    // get_contemporaneous_recipients: integration with find_time_index
    // ---------------------------------------------------------------

    #[test]
    fn get_contemporaneous_recipients_excludes_donor() {
        let depths = vec![0.0, 1.0, 2.0];
        // contemporaneity[1] has species 0, 1, 2 alive in interval (0.0, 1.0]
        let contemporaneity = vec![
            vec![],        // interval ending at 0.0
            vec![0, 1, 2], // interval (0.0, 1.0]
            vec![1, 2],    // interval (1.0, 2.0]
        ];
        let recipients = get_contemporaneous_recipients(&depths, &contemporaneity, 0.5, 1);
        assert_eq!(recipients, vec![0, 2]);
    }

    #[test]
    fn get_contemporaneous_recipients_empty_when_donor_is_only_species() {
        let depths = vec![0.0, 1.0];
        let contemporaneity = vec![
            vec![],
            vec![3], // only species 3 alive
        ];
        let recipients = get_contemporaneous_recipients(&depths, &contemporaneity, 0.5, 3);
        assert!(recipients.is_empty());
    }
}
