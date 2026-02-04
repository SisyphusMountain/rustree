// DTL (Duplication-Transfer-Loss) simulation along a species tree
//
// This module simulates gene tree evolution within a species tree using the DTL model.
// Events: Speciation (S), Duplication (D), Transfer (T), Loss (L)

use crate::bd::BDEvent;
use crate::node::{FlatTree, FlatNode, Event, RecTree};
use rand::Rng;

// Static strings to avoid allocations
const EVENT_DUPLICATION: &str = "Duplication";
const EVENT_TRANSFER: &str = "Transfer";
const EVENT_LOSS: &str = "Loss";
const EVENT_LEAF: &str = "Leaf";
const EVENT_SPECIATION: &str = "Speciation";

/// Represents a DTL event during gene tree simulation
#[derive(Clone, Debug)]
pub struct DTLEvent {
    /// Time when the event occurred (absolute time from root)
    pub time: f64,
    /// Gene node ID where the event occurred
    pub gene_node_id: usize,
    /// Type of event: static string reference to avoid allocations
    pub event_type: &'static str,
    /// Species tree node index where the event occurred
    pub species_node_idx: usize,
    /// For transfers: donor species node index
    pub donor_species_idx: Option<usize>,
    /// For transfers: recipient species node index
    pub recipient_species_idx: Option<usize>,
    /// First child gene node ID (for Speciation, Duplication, Transfer)
    pub child1: Option<usize>,
    /// Second child gene node ID (for Speciation, Duplication, Transfer)
    pub child2: Option<usize>,
}

impl DTLEvent {
    /// Convert event to CSV row format, resolving names from trees
    ///
    /// # Arguments
    /// * `species_tree` - The species tree to resolve species names
    /// * `gene_tree` - The gene tree to resolve gene node names
    pub fn to_csv_row(&self, species_tree: &FlatTree, gene_tree: &FlatTree) -> String {
        format!(
            "{},{},{},{},{},{},{},{}",
            self.time,
            gene_tree.nodes[self.gene_node_id].name,
            self.event_type,
            species_tree.nodes[self.species_node_idx].name,
            self.donor_species_idx.map_or(String::new(), |idx| species_tree.nodes[idx].name.clone()),
            self.recipient_species_idx.map_or(String::new(), |idx| species_tree.nodes[idx].name.clone()),
            self.child1.map_or(String::new(), |c| gene_tree.nodes[c].name.clone()),
            self.child2.map_or(String::new(), |c| gene_tree.nodes[c].name.clone())
        )
    }

    /// CSV header for event data
    pub fn csv_header() -> &'static str {
        "time,gene_node_name,event_type,species_node,donor_species,recipient_species,child1_name,child2_name"
    }
}

/// Represents an active gene copy being simulated
#[derive(Clone, Debug)]
struct GeneCopy {
    /// Index of this gene in the gene tree nodes vector
    gene_node_idx: usize,
    /// Index of the species tree branch this gene is on
    species_node_idx: usize,
    /// Current time position on the branch (absolute time from root)
    current_time: f64,
}

/// Holds the mutable state during DTL simulation.
/// Extracted to simplify the main simulation loop and enable event-specific methods.
struct SimulationState {
    gene_nodes: Vec<FlatNode>,
    node_mapping: Vec<usize>,
    event_mapping: Vec<Event>,
    events: Vec<DTLEvent>,
}

impl SimulationState {
    fn new(capacity: usize) -> Self {
        Self {
            gene_nodes: Vec::with_capacity(capacity),
            node_mapping: Vec::with_capacity(capacity),
            event_mapping: Vec::with_capacity(capacity),
            events: Vec::with_capacity(capacity),
        }
    }

    /// Creates a new gene node and returns its index.
    /// Node name format: {species_idx}_{gene_idx}
    fn create_gene_node(&mut self, parent: Option<usize>, species_idx: usize, event: Event) -> usize {
        let idx = self.gene_nodes.len();
        self.gene_nodes.push(FlatNode {
            name: format!("{}_{}", species_idx, idx),
            left_child: None,
            right_child: None,
            parent,
            depth: None,
            length: 0.0,
            bd_event: None,
        });
        self.node_mapping.push(species_idx);
        self.event_mapping.push(event);
        idx
    }

    /// Handles a duplication event: creates two children and records the event.
    /// Returns the indices of the two child gene nodes.
    fn handle_duplication(
        &mut self,
        parent_idx: usize,
        species_idx: usize,
        event_time: f64,
    ) -> (usize, usize) {
        let child1_idx = self.create_gene_node(Some(parent_idx), species_idx, Event::Duplication);
        let child2_idx = self.create_gene_node(Some(parent_idx), species_idx, Event::Duplication);

        // Update parent node
        self.gene_nodes[parent_idx].left_child = Some(child1_idx);
        self.gene_nodes[parent_idx].right_child = Some(child2_idx);
        self.gene_nodes[parent_idx].depth = Some(event_time);
        self.event_mapping[parent_idx] = Event::Duplication;

        // Record event
        self.events.push(DTLEvent {
            time: event_time,
            gene_node_id: parent_idx,
            event_type: EVENT_DUPLICATION,
            species_node_idx: species_idx,
            donor_species_idx: None,
            recipient_species_idx: None,
            child1: Some(child1_idx),
            child2: Some(child2_idx),
        });

        (child1_idx, child2_idx)
    }

    /// Handles a transfer event: creates donor and recipient children.
    /// Returns the indices of (donor_child, recipient_child).
    fn handle_transfer(
        &mut self,
        parent_idx: usize,
        donor_species: usize,
        recipient_species: usize,
        event_time: f64,
    ) -> (usize, usize) {
        let donor_child_idx = self.create_gene_node(Some(parent_idx), donor_species, Event::Transfer);
        let recipient_child_idx = self.create_gene_node(Some(parent_idx), recipient_species, Event::Transfer);

        // Update parent node
        self.gene_nodes[parent_idx].left_child = Some(donor_child_idx);
        self.gene_nodes[parent_idx].right_child = Some(recipient_child_idx);
        self.gene_nodes[parent_idx].depth = Some(event_time);
        self.event_mapping[parent_idx] = Event::Transfer;

        // Record event
        self.events.push(DTLEvent {
            time: event_time,
            gene_node_id: parent_idx,
            event_type: EVENT_TRANSFER,
            species_node_idx: donor_species,
            donor_species_idx: Some(donor_species),
            recipient_species_idx: Some(recipient_species),
            child1: Some(donor_child_idx),
            child2: Some(recipient_child_idx),
        });

        (donor_child_idx, recipient_child_idx)
    }

    /// Handles a loss event: marks the gene as lost.
    fn handle_loss(&mut self, gene_idx: usize, species_idx: usize, event_time: f64) {
        self.gene_nodes[gene_idx].depth = Some(event_time);
        self.event_mapping[gene_idx] = Event::Loss;

        self.events.push(DTLEvent {
            time: event_time,
            gene_node_id: gene_idx,
            event_type: EVENT_LOSS,
            species_node_idx: species_idx,
            donor_species_idx: None,
            recipient_species_idx: None,
            child1: None,
            child2: None,
        });
    }

    /// Handles reaching a leaf species: gene survives as extant.
    fn handle_leaf(&mut self, gene_idx: usize, species_idx: usize, event_time: f64) {
        self.gene_nodes[gene_idx].depth = Some(event_time);
        self.event_mapping[gene_idx] = Event::Leaf;

        self.events.push(DTLEvent {
            time: event_time,
            gene_node_id: gene_idx,
            event_type: EVENT_LEAF,
            species_node_idx: species_idx,
            donor_species_idx: None,
            recipient_species_idx: None,
            child1: None,
            child2: None,
        });
    }

    /// Handles a speciation event: gene follows both children species.
    /// Returns the indices of (left_gene, right_gene).
    fn handle_speciation(
        &mut self,
        parent_idx: usize,
        parent_species: usize,
        left_species: usize,
        right_species: usize,
        event_time: f64,
    ) -> (usize, usize) {
        let left_gene_idx = self.create_gene_node(Some(parent_idx), left_species, Event::Speciation);
        let right_gene_idx = self.create_gene_node(Some(parent_idx), right_species, Event::Speciation);

        // Update parent node
        self.gene_nodes[parent_idx].left_child = Some(left_gene_idx);
        self.gene_nodes[parent_idx].right_child = Some(right_gene_idx);
        self.gene_nodes[parent_idx].depth = Some(event_time);
        self.event_mapping[parent_idx] = Event::Speciation;

        // Record event
        self.events.push(DTLEvent {
            time: event_time,
            gene_node_id: parent_idx,
            event_type: EVENT_SPECIATION,
            species_node_idx: parent_species,
            donor_species_idx: None,
            recipient_species_idx: None,
            child1: Some(left_gene_idx),
            child2: Some(right_gene_idx),
        });

        (left_gene_idx, right_gene_idx)
    }
}

/// Simulates a gene tree along a species tree using the DTL model (internal version with pre-computed data)
///
/// # Arguments
/// * `species_tree` - The species tree (must have depths assigned)
/// * `depths` - Pre-computed time subdivision
/// * `contemporaneity` - Pre-computed contemporaneity information
/// * `lca_depths` - Optional precomputed LCA depths matrix for assortative transfers
/// * `origin_species` - Species node index where the gene family originates
/// * `lambda_d` - Duplication rate per unit time
/// * `lambda_t` - Transfer rate per unit time
/// * `lambda_l` - Loss rate per unit time
/// * `transfer_alpha` - Optional distance decay parameter for assortative transfers
/// * `rng` - Random number generator
///
/// # Returns
/// A tuple containing:
/// - A `RecTree` with the simulated gene tree and mappings to the species tree
/// - A `Vec<DTLEvent>` with all events that occurred during simulation
fn simulate_dtl_internal<'a, R: Rng>(
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
    let initial_gene_idx = state.create_gene_node(None, origin_species, Event::Speciation);
    state.gene_nodes[initial_gene_idx].depth = Some(origin_depth);

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
                    state.gene_nodes[current_gene_idx].length += species_end_time - current_time;
                    current_time = species_end_time;
                } else {
                    // Event occurs
                    state.gene_nodes[current_gene_idx].length += event_time - current_time;
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
                        // Transfer - find valid recipient (uniform or assortative)
                        let recipient = match (transfer_alpha, lca_depths) {
                            (Some(alpha), Some(lca)) => find_transfer_recipient_assortative(
                                lca, depths, contemporaneity, event_time,
                                copy.species_node_idx, alpha, rng
                            ),
                            _ => find_transfer_recipient(
                                depths, contemporaneity, event_time,
                                copy.species_node_idx, rng
                            ),
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
    finalize_simulation(species_tree, state, origin_species, origin_depth)
}

/// Draws an exponential waiting time for the next event.
#[inline]
fn draw_waiting_time<R: Rng>(total_rate: f64, rng: &mut R) -> f64 {
    if total_rate > 0.0 {
        let u: f64 = rng.gen();
        -u.ln() / total_rate
    } else {
        f64::INFINITY
    }
}

/// Finds a valid transfer recipient from contemporary species (uniform random).
/// Returns None if no valid recipient exists (donor is the only contemporary species).
fn find_transfer_recipient<R: Rng>(
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    event_time: f64,
    donor_species: usize,
    rng: &mut R,
) -> Option<usize> {
    let time_idx = find_time_index(depths, event_time);
    let contemporaries = &contemporaneity[time_idx];

    // If donor is alone, no valid recipient
    if contemporaries.len() <= 1 {
        return None;
    }

    // Draw random recipient, redraw if we hit the donor.
    // Kind of dirty, but avoid allocating a new vector without the donor.
    loop {
        let idx = rng.gen_range(0..contemporaries.len());
        let recipient = contemporaries[idx];
        if recipient != donor_species {
            return Some(recipient);
        }
    }
}

/// Finds a transfer recipient with distance-dependent (assortative) weighting.
///
/// The probability of selecting species B as recipient from donor A at time t is:
/// P(B) ∝ exp(-alpha * d(A, B, t))
///
/// where d(A, B, t) = 2 * (t - depth_of_LCA(A, B)) is the phylogenetic distance at time t.
///
/// # Arguments
/// * `lca_depths` - Precomputed matrix where lca_depths[i][j] = depth of LCA(i, j)
/// * `depths` - Time subdivision array
/// * `contemporaneity` - Contemporary species at each time interval
/// * `event_time` - Time of the transfer event
/// * `donor_species` - Index of the donor species
/// * `alpha` - Distance decay parameter (higher = more local transfers)
/// * `rng` - Random number generator
///
/// # Returns
/// Some(recipient_index) or None if donor is alone
fn find_transfer_recipient_assortative<R: Rng>(
    lca_depths: &[Vec<f64>],
    depths: &[f64],
    contemporaneity: &[Vec<usize>],
    event_time: f64,
    donor_species: usize,
    alpha: f64,
    rng: &mut R,
) -> Option<usize> {
    let time_idx = find_time_index(depths, event_time);
    let contemporaries = &contemporaneity[time_idx];

    // If donor is alone, no valid recipient
    if contemporaries.len() <= 1 {
        return None;
    }

    // Compute weights for each potential recipient
    let mut weights: Vec<f64> = Vec::with_capacity(contemporaries.len());
    let mut total_weight = 0.0;

    for &species in contemporaries {
        if species == donor_species {
            weights.push(0.0); // Donor cannot receive from itself
        } else {
            // Distance = 2 * (event_time - LCA_depth)
            let lca_depth = lca_depths[donor_species][species];
            let distance = 2.0 * (event_time - lca_depth);
            let weight = (-alpha * distance).exp();
            weights.push(weight);
            total_weight += weight;
        }
    }

    // No valid recipients (shouldn't happen if len > 1, but be safe)
    if total_weight <= 0.0 {
        return None;
    }

    // Sample from weighted distribution
    let threshold: f64 = rng.gen::<f64>() * total_weight;
    let mut cumulative = 0.0;

    for (i, &weight) in weights.iter().enumerate() {
        cumulative += weight;
        if cumulative >= threshold {
            return Some(contemporaries[i]);
        }
    }
    // panic if we get here
    panic!("Failed to sample transfer recipient in assortative transfer");
}

/// Finalizes simulation by finding the root and handling edge cases.
fn finalize_simulation<'a>(
    species_tree: &'a FlatTree,
    mut state: SimulationState,
    origin_species: usize,
    origin_depth: f64,
) -> (RecTree<'a>, Vec<DTLEvent>) {
    // Handle edge case: gene lost before any events
    // Node name format: {species_idx}_{gene_idx}
    if state.gene_nodes.is_empty() {
        state.gene_nodes.push(FlatNode {
            name: format!("{}_0", origin_species),
            left_child: None,
            right_child: None,
            parent: None,
            depth: Some(origin_depth),
            length: 0.0,
            bd_event: None,
        });
        state.node_mapping.push(origin_species);
        state.event_mapping.push(Event::Loss);
    }

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

/// Simulates a gene tree along a species tree using the DTL model
///
/// This is the main public interface. It computes depths and contemporaneity,
/// then calls the internal simulation function.
///
/// # Arguments
/// * `species_tree` - The species tree (must have depths assigned)
/// * `origin_species` - Species node index where the gene family originates
/// * `lambda_d` - Duplication rate per unit time
/// * `lambda_t` - Transfer rate per unit time
/// * `lambda_l` - Loss rate per unit time
/// * `transfer_alpha` - Optional distance decay for assortative transfers (None = uniform)
/// * `require_extant` - If true, retry until a tree with at least one extant gene is produced
/// * `rng` - Random number generator
///
/// # Returns
/// A tuple containing:
/// - A `RecTree` with the simulated gene tree and mappings to the species tree
/// - A `Vec<DTLEvent>` with all events that occurred during simulation
///
/// # Panics
/// Panics if the species tree doesn't have depths assigned
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
    assert!(lambda_d >= 0.0, "Duplication rate must be non-negative");
    assert!(lambda_t >= 0.0, "Transfer rate must be non-negative");
    assert!(lambda_l >= 0.0, "Loss rate must be non-negative");

    // Compute depths and contemporaneity
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);

    // Compute LCA depths if using assortative transfers
    let lca_depths = transfer_alpha.map(|_| species_tree.precompute_lca_depths());
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
///
/// This is more efficient than calling simulate_dtl multiple times because:
/// - Species tree depths and contemporaneity are computed only once
/// - LCA depths (for assortative transfers) are computed only once
/// - Avoids redundant computation for each simulation
///
/// # Arguments
/// * `species_tree` - The species tree (must have depths assigned)
/// * `origin_species` - Species node index where the gene family originates
/// * `lambda_d` - Duplication rate per unit time
/// * `lambda_t` - Transfer rate per unit time
/// * `lambda_l` - Loss rate per unit time
/// * `transfer_alpha` - Optional distance decay for assortative transfers (None = uniform)
/// * `n_simulations` - Number of gene trees to simulate
/// * `require_extant` - If true, only include trees with at least one extant gene
/// * `rng` - Random number generator
///
/// # Returns
/// A tuple containing:
/// - A `Vec<RecTree>` with all simulated gene trees
/// - A `Vec<Vec<DTLEvent>>` with events for each simulation
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
    assert!(lambda_d >= 0.0, "Duplication rate must be non-negative");
    assert!(lambda_t >= 0.0, "Transfer rate must be non-negative");
    assert!(lambda_l >= 0.0, "Loss rate must be non-negative");

    // Compute depths and contemporaneity once for all simulations
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);

    // Compute LCA depths once if using assortative transfers
    let lca_depths = transfer_alpha.map(|_| species_tree.precompute_lca_depths());
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

// Re-export from io module for backward compatibility
pub use crate::io::save_dtl_events_to_csv as save_events_to_csv;

/// Find the index in the depths array corresponding to a given time
fn find_time_index(depths: &[f64], time: f64) -> usize {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::newick::newick::parse_newick;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn parse_tree(newick: &str) -> FlatTree {
        let mut nodes = parse_newick(newick).unwrap();
        let root = nodes.pop().expect("No tree found");
        root.to_flat_tree()
    }

    #[test]
    fn test_dtl_pure_speciation() {
        // Simple species tree: ((A:1,B:1)AB:1,C:2)root:0;
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // With zero D/T/L rates, should get pure speciation
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 0.0, 0.0, None, false, &mut rng);

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!("Events: S={}, D={}, T={}, L={}, Leaves={}", s, d, t, l, leaves);
        println!("Total event count: {}", events.len());

        // Should have 3 extant genes (one per species leaf)
        let extant = count_extant_genes(&rec_tree);
        assert_eq!(extant, 3, "Should have 3 extant genes with no D/T/L");
    }

    #[test]
    fn test_dtl_with_duplication() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(123);

        // High duplication rate
        let (rec_tree, _events) = simulate_dtl(&species_tree, species_tree.root, 2.0, 0.0, 0.0, None, false, &mut rng);

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!("Events with duplication: S={}, D={}, T={}, L={}, Leaves={}", s, d, t, l, leaves);

        // Should have duplications
        assert!(d > 0 || leaves >= 3, "Should have duplications or at least 3 leaves");
    }

    #[test]
    fn test_dtl_with_loss() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(999);

        // High loss rate
        let (rec_tree, _events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 0.0, 5.0, None, false, &mut rng);

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!("Events with loss: S={}, D={}, T={}, L={}, Leaves={}", s, d, t, l, leaves);

        // May have losses
        let extant = count_extant_genes(&rec_tree);
        println!("Extant genes: {}", extant);
    }

    #[test]
    fn test_dtl_xml_export() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // With mixed events
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 1.0, 0.5, 0.5, None, false, &mut rng);

        // Save events to CSV
        save_events_to_csv(&events, &species_tree, &rec_tree.gene_tree, "test_dtl_events.csv").expect("Failed to write events CSV");
        println!("Events saved to test_dtl_events.csv");

        // Generate XML
        let xml = rec_tree.to_xml();

        // Verify XML has required sections
        assert!(xml.contains("<recPhylo"));
        assert!(xml.contains("<spTree>"));
        assert!(xml.contains("<recGeneTree>"));
        assert!(xml.contains("<branchLength>"));

        // Write to file for inspection
        use std::fs;
        fs::write("test_dtl_output.xml", &xml).expect("Failed to write XML");
        println!("XML output written to test_dtl_output.xml");
    }

    #[test]
    fn test_dtl_with_transfer() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(456);

        // High transfer rate (uniform)
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 2.0, 0.0, None, false, &mut rng);

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!("Events with transfer: S={}, D={}, T={}, L={}, Leaves={}", s, d, t, l, leaves);
        println!("Total event count: {}", events.len());
    }

    #[test]
    fn test_dtl_assortative_transfer() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(456);

        // High transfer rate with assortative selection (alpha=1.0)
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 2.0, 0.0, Some(1.0), false, &mut rng);

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!("Events with assortative transfer: S={}, D={}, T={}, L={}, Leaves={}", s, d, t, l, leaves);
        println!("Total event count: {}", events.len());

        // Verify LCA computation works
        let lca_depths = species_tree.precompute_lca_depths();
        assert!(lca_depths.len() > 0, "LCA depths should be computed");
    }
}
