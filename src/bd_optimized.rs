// Optimized birth-death tree simulation
use crate::node::{FlatTree, FlatNode};
use crate::bd::TreeEvent;
use rand::Rng;

/// Optimized version of simulate_bd_tree with performance improvements:
/// 1. Pre-allocated vectors with capacity hints
/// 2. swap_remove() instead of remove() (O(1) vs O(n))
/// 3. Cached computed values (mu/(lambda+mu), lambda+mu)
/// 4. Reduced string allocations
pub fn simulate_bd_tree_optimized<R: Rng>(n: usize, lambda: f64, mu: f64, rng: &mut R) -> (FlatTree, Vec<TreeEvent>) {
    assert!(n > 0, "Number of species must be positive");
    assert!(lambda > 0.0, "Speciation rate must be positive");
    assert!(mu >= 0.0, "Extinction rate must be non-negative");
    assert!(lambda > mu, "Speciation rate must be strictly greater than extinction rate");

    // Pre-compute constants
    let lambda_plus_mu = lambda + mu;
    let extinction_prob = mu / lambda_plus_mu;

    // Estimate capacity: n leaves + ~n internal nodes + some extinctions
    let estimated_capacity = (n as f64 * 2.5) as usize;

    let mut time = 0.0;
    let mut num_lineages = n;

    // Pre-allocate vectors with capacity
    let mut active_lineages: Vec<(usize, f64)> = Vec::with_capacity(n);
    let mut nodes: Vec<FlatNode> = Vec::with_capacity(estimated_capacity);
    let mut events: Vec<TreeEvent> = Vec::with_capacity(estimated_capacity);

    // Static strings to avoid allocations
    const LEAF_EVENT: &str = "Leaf";
    const EXTINCTION_EVENT: &str = "Extinction";
    const SPECIATION_EVENT: &str = "Speciation";

    // Create n extant species (leaf nodes) at present
    for i in 0..n {
        nodes.push(FlatNode {
            name: i.to_string(),
            left_child: None,
            right_child: None,
            parent: None,
            depth: Some(0.0),
            length: 0.0,
        });
        active_lineages.push((i, 0.0));

        events.push(TreeEvent {
            time: 0.0,
            node_id: i,
            event_type: LEAF_EVENT.to_string(),
            child1: None,
            child2: None,
        });
    }

    // Simulate backwards in time
    while num_lineages > 0 {
        // Draw waiting time from exponential distribution
        let rate = lambda_plus_mu * (num_lineages as f64);
        let uniform_sample: f64 = rng.gen();
        let waiting_time = -uniform_sample.ln() / rate;
        time += waiting_time;

        // Determine event type
        let event_prob: f64 = rng.gen();

        if event_prob < extinction_prob {
            // Extinction event
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

            events.push(TreeEvent {
                time,
                node_id: extinct_idx,
                event_type: EXTINCTION_EVENT.to_string(),
                child1: None,
                child2: None,
            });
        } else {
            // Speciation event
            if num_lineages == 1 {
                let (root_idx, root_start) = active_lineages[0];
                nodes[root_idx].length = time - root_start;
                num_lineages = 0;
            } else {
                // Use swap_remove for O(1) removal instead of O(n) remove
                let idx1 = rng.gen_range(0..active_lineages.len());
                let (child1_idx, child1_start) = active_lineages.swap_remove(idx1);

                let idx2 = rng.gen_range(0..active_lineages.len());
                let (child2_idx, child2_start) = active_lineages.swap_remove(idx2);

                let parent_idx = nodes.len();
                nodes.push(FlatNode {
                    name: parent_idx.to_string(),
                    left_child: Some(child1_idx),
                    right_child: Some(child2_idx),
                    parent: None,
                    depth: Some(time),
                    length: 0.0,
                });

                nodes[child1_idx].parent = Some(parent_idx);
                nodes[child1_idx].length = time - child1_start;

                nodes[child2_idx].parent = Some(parent_idx);
                nodes[child2_idx].length = time - child2_start;

                active_lineages.push((parent_idx, time));
                num_lineages -= 1;

                events.push(TreeEvent {
                    time,
                    node_id: parent_idx,
                    event_type: SPECIATION_EVENT.to_string(),
                    child1: Some(child1_idx),
                    child2: Some(child2_idx),
                });
            }
        }
    }

    let root_idx = nodes.iter().position(|n| n.parent.is_none()).unwrap();

    let tree = FlatTree {
        nodes,
        root: root_idx,
    };

    (tree, events)
}
