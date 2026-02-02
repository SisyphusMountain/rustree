// birth-death processes for the generation of trees.
// We start with classical birth-death processes.

use crate::node::{FlatTree, FlatNode};
use rand::Rng;

/// Represents an event in the birth-death process
#[derive(Clone, Debug)]
pub struct TreeEvent {
    /// Time when the event occurred (going backwards from present at 0)
    pub time: f64,
    /// Node ID where the event occurred
    pub node_id: usize,
    /// Type of event: "Speciation", "Extinction", or "Leaf"
    pub event_type: String,
    /// First child node ID (for speciation events)
    pub child1: Option<usize>,
    /// Second child node ID (for speciation events)
    pub child2: Option<usize>,
}

impl TreeEvent {
    /// Convert event to CSV row format, resolving node names from the tree
    pub fn to_csv_row(&self, tree: &FlatTree) -> String {
        format!(
            "{},{},{},{},{}",
            self.time,
            tree.nodes[self.node_id].name,
            self.event_type,
            self.child1.map_or(String::new(), |c| tree.nodes[c].name.clone()),
            self.child2.map_or(String::new(), |c| tree.nodes[c].name.clone())
        )
    }

    /// CSV header for event data
    pub fn csv_header() -> &'static str {
        "time,node_name,event_type,child1_name,child2_name"
    }
}

/// Simulates a tree with a fixed number of extant species using the
/// constant rate birth-death backward process (Stadler 2011)
///
/// This implements the EBDP backward algorithm for the simple case of
/// constant rates (no mass extinctions or rate shifts).
///
/// # Arguments
/// * `n` - Number of extant species (must be > 0)
/// * `lambda` - Speciation/birth rate (must be > 0)
/// * `mu` - Extinction/death rate (must be >= 0 and < lambda)
/// * `rng` - Random number generator
///
/// # Returns
/// A tuple containing:
/// - `FlatTree` with `n` extant species simulated under the birth-death process
/// - `Vec<TreeEvent>` list of events that occurred during the simulation
///
/// Node names are integers: 0 to n-1 are extant species, n onwards are other nodes
///
/// # Panics
/// Panics if lambda <= mu, as the backwards process may grow indefinitely
///
/// # References
/// Stadler, T. (2011). Simulating trees with a fixed number of extant species.
/// Systematic Biology, 60(5), 676-684.
pub fn simulate_bd_tree<R: Rng>(n: usize, lambda: f64, mu: f64, rng: &mut R) -> (FlatTree, Vec<TreeEvent>) {
    assert!(n > 0, "Number of species must be positive");
    assert!(lambda > 0.0, "Speciation rate must be positive");
    assert!(mu >= 0.0, "Extinction rate must be non-negative");
    assert!(lambda > mu, "Speciation rate must be strictly greater than extinction rate");

    // Start at present (time = 0) with n isolated vertices
    let mut time = 0.0;
    let mut num_lineages = n;

    // Track active lineages: (node_index, time_when_lineage_started)
    let mut active_lineages: Vec<(usize, f64)> = Vec::new();

    // Build the tree nodes
    let mut nodes: Vec<FlatNode> = Vec::new();

    // Track events
    let mut events: Vec<TreeEvent> = Vec::new();

    // Create n extant species (leaf nodes) at present
    // Nodes 0 to n-1 are extant species
    for i in 0..n {
        nodes.push(FlatNode {
            name: i.to_string(),
            left_child: None,
            right_child: None,
            parent: None,
            depth: Some(0.0),
            length: 0.0, // Will be set when we add the parent
        });
        active_lineages.push((i, 0.0));

        // Record leaf event
        events.push(TreeEvent {
            time: 0.0,
            node_id: i,
            event_type: "Leaf".to_string(),
            child1: None,
            child2: None,
        });
    }

    // Simulate backwards in time until we reach the origin
    while num_lineages > 0 {
        // Draw waiting time from exponential distribution with rate (λ + μ)N
        // Using inverse transform method: if U ~ Uniform(0,1), then -ln(U)/rate ~ Exp(rate)
        let rate = (lambda + mu) * (num_lineages as f64);
        let uniform_sample: f64 = rng.gen();
        let waiting_time = -uniform_sample.ln() / rate;
        time += waiting_time;

        // Determine event type
        let event_prob: f64 = rng.gen();

        if event_prob < mu / (lambda + mu) {
            // Extinction event (backward in time = lineage appears)
            // Forward in time, this lineage went extinct
            let extinct_idx = nodes.len();
            nodes.push(FlatNode {
                name: extinct_idx.to_string(),
                left_child: None,
                right_child: None,
                parent: None,
                depth: Some(time),
                length: 0.0,
            });
            active_lineages.push((extinct_idx, time));
            num_lineages += 1;

            // Record extinction event
            events.push(TreeEvent {
                time,
                node_id: extinct_idx,
                event_type: "Extinction".to_string(),
                child1: None,
                child2: None,
            });
        } else {
            // Speciation event (backward in time = coalescence)
            if num_lineages == 1 {
                // Reached the origin (root of the tree)
                let (root_idx, root_start) = active_lineages[0];
                nodes[root_idx].length = time - root_start;
                num_lineages = 0;
            } else {
                // Pick two random lineages uniformly and coalesce them
                let idx1 = rng.gen_range(0..active_lineages.len());
                let (child1_idx, child1_start) = active_lineages.remove(idx1);

                let idx2 = rng.gen_range(0..active_lineages.len());
                let (child2_idx, child2_start) = active_lineages.remove(idx2);

                // Create parent node (the speciation event)
                let parent_idx = nodes.len();
                nodes.push(FlatNode {
                    name: parent_idx.to_string(),
                    left_child: Some(child1_idx),
                    right_child: Some(child2_idx),
                    parent: None,
                    depth: Some(time),
                    length: 0.0, // Will be set when we add its parent
                });

                // Update children to point to parent and set branch lengths
                nodes[child1_idx].parent = Some(parent_idx);
                nodes[child1_idx].length = time - child1_start;

                nodes[child2_idx].parent = Some(parent_idx);
                nodes[child2_idx].length = time - child2_start;

                // Parent becomes an active lineage
                active_lineages.push((parent_idx, time));
                num_lineages -= 1;

                // Record speciation event
                events.push(TreeEvent {
                    time,
                    node_id: parent_idx,
                    event_type: "Speciation".to_string(),
                    child1: Some(child1_idx),
                    child2: Some(child2_idx),
                });
            }
        }
    }

    // Find the root (the node with no parent)
    let root_idx = nodes.iter().position(|n| n.parent.is_none()).unwrap();

    let tree = FlatTree {
        nodes,
        root: root_idx,
    };

    (tree, events)
}

/// Generate birth-death events from an existing tree (e.g., parsed from Newick).
///
/// This creates events by treating all internal nodes as speciations and all
/// leaves as extant species. No extinction events are generated since the tree
/// only contains surviving lineages.
///
/// # Arguments
/// * `tree` - A FlatTree with depths assigned
///
/// # Returns
/// A vector of TreeEvent representing the tree's evolutionary history
///
/// # Panics
/// Panics if any node doesn't have a depth assigned
///
/// # Example
/// ```ignore
/// let tree = parse_newick("((A:1,B:1):1,C:2):0;").unwrap();
/// let mut flat_tree = tree.to_flat_tree();
/// flat_tree.assign_depths();
/// let events = generate_events_from_tree(&flat_tree);
/// ```
pub fn generate_events_from_tree(tree: &FlatTree) -> Vec<TreeEvent> {
    let mut events = Vec::with_capacity(tree.nodes.len());

    for (idx, node) in tree.nodes.iter().enumerate() {
        let depth = node.depth.expect("Tree must have depths assigned. Call assign_depths() first.");

        if node.left_child.is_none() && node.right_child.is_none() {
            // Leaf node
            events.push(TreeEvent {
                time: depth,
                node_id: idx,
                event_type: "Leaf".to_string(),
                child1: None,
                child2: None,
            });
        } else {
            // Internal node = speciation
            events.push(TreeEvent {
                time: depth,
                node_id: idx,
                event_type: "Speciation".to_string(),
                child1: node.left_child,
                child2: node.right_child,
            });
        }
    }

    // Sort by time (depth) for consistency
    events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());

    events
}

// Re-export from io module for backward compatibility
pub use crate::io::save_bd_events_to_csv as save_events_to_csv;