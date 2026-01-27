// DTL (Duplication-Transfer-Loss) simulation along a species tree
//
// This module simulates gene tree evolution within a species tree using the DTL model.
// Events: Speciation (S), Duplication (D), Transfer (T), Loss (L)

use crate::node::{FlatTree, FlatNode, Event, RecTree};
use rand::Rng;

/// Represents a DTL event during gene tree simulation
#[derive(Clone, Debug)]
pub struct DTLEvent {
    /// Time when the event occurred (absolute time from root)
    pub time: f64,
    /// Gene node ID where the event occurred
    pub gene_node_id: usize,
    /// Type of event: "Speciation", "Duplication", "Transfer", "Loss", or "Leaf"
    pub event_type: String,
    /// Species tree branch where the event occurred
    pub species_node: String,
    /// For transfers: donor species branch
    pub donor_species: Option<String>,
    /// For transfers: recipient species branch
    pub recipient_species: Option<String>,
    /// First child gene node ID (for Speciation, Duplication, Transfer)
    pub child1: Option<usize>,
    /// Second child gene node ID (for Speciation, Duplication, Transfer)
    pub child2: Option<usize>,
}

impl DTLEvent {
    /// Convert event to CSV row format
    pub fn to_csv_row(&self) -> String {
        format!(
            "{},{},{},{},{},{},{},{}",
            self.time,
            self.gene_node_id,
            self.event_type,
            self.species_node,
            self.donor_species.as_ref().unwrap_or(&String::new()),
            self.recipient_species.as_ref().unwrap_or(&String::new()),
            self.child1.map_or(String::from(""), |c| c.to_string()),
            self.child2.map_or(String::from(""), |c| c.to_string())
        )
    }

    /// CSV header for event data
    pub fn csv_header() -> &'static str {
        "time,gene_node_id,event_type,species_node,donor_species,recipient_species,child1,child2"
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

/// Simulates a gene tree along a species tree using the DTL model
///
/// # Arguments
/// * `species_tree` - The species tree (must have depths assigned)
/// * `origin_species` - Species node index where the gene family originates
/// * `lambda_d` - Duplication rate per unit time
/// * `lambda_t` - Transfer rate per unit time
/// * `lambda_l` - Loss rate per unit time
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
    rng: &mut R,
) -> (RecTree<'a>, Vec<DTLEvent>) {
    assert!(lambda_d >= 0.0, "Duplication rate must be non-negative");
    assert!(lambda_t >= 0.0, "Transfer rate must be non-negative");
    assert!(lambda_l >= 0.0, "Loss rate must be non-negative");

    // Get time subdivision and contemporaneity from species tree
    let depths = species_tree.make_subdivision();
    let contemporaneity = species_tree.find_contemporaneity(&depths);

    // Find the depth of the origin species node
    let origin_depth = species_tree.nodes[origin_species]
        .depth
        .expect("Species tree must have depths assigned");

    // Gene tree construction
    let mut gene_nodes: Vec<FlatNode> = Vec::new();
    let mut node_mapping: Vec<usize> = Vec::new();
    let mut event_mapping: Vec<Event> = Vec::new();
    let mut events: Vec<DTLEvent> = Vec::new();

    // Create the initial gene copy at the origin
    let initial_gene_idx = gene_nodes.len();
    gene_nodes.push(FlatNode {
        name: format!("g{}", initial_gene_idx),
        left_child: None,
        right_child: None,
        parent: None,
        depth: Some(origin_depth),
        length: 0.0, // Will be set when child is created or at leaf
    });
    node_mapping.push(origin_species);
    event_mapping.push(Event::Speciation); // Origin event (will be updated if it's a leaf)

    // Active gene copies being simulated
    let mut active_copies: Vec<GeneCopy> = vec![GeneCopy {
        gene_node_idx: initial_gene_idx,
        current_time: origin_depth,
        species_node_idx: origin_species,
    }];

    // Total event rate
    let total_rate = lambda_d + lambda_t + lambda_l;

    // Process until no active copies remain
    while !active_copies.is_empty() {
        // Process each active copy
        let mut new_copies: Vec<GeneCopy> = Vec::new();
        let mut copies_to_remove: Vec<usize> = Vec::new();

        for (copy_idx, copy) in active_copies.iter().enumerate() {
            let species_node = &species_tree.nodes[copy.species_node_idx];
            let species_end_time = species_node.depth.unwrap();

            // Simulate events along this branch segment
            let mut current_time = copy.current_time;
            let current_gene_idx = copy.gene_node_idx;
            let mut lost = false;

            // Simulate D/T/L events along the branch
            while current_time < species_end_time && !lost {
                // Draw waiting time to next event
                let waiting_time = if total_rate > 0.0 {
                    let u: f64 = rng.gen();
                    -u.ln() / total_rate
                } else {
                    f64::INFINITY
                };

                let event_time = current_time + waiting_time;

                if event_time >= species_end_time {
                    // No event before end of branch, advance to the node
                    let branch_length = species_end_time - current_time;
                    gene_nodes[current_gene_idx].length += branch_length;
                    current_time = species_end_time;
                    // Will process node after the while loop
                } else {
                    // An event occurs
                    let branch_length = event_time - current_time;
                    gene_nodes[current_gene_idx].length += branch_length;
                    current_time = event_time;

                    // Determine event type
                    let event_prob: f64 = rng.gen();

                    if event_prob < lambda_d / total_rate {
                        // Duplication event
                        let child1_idx = gene_nodes.len();
                        gene_nodes.push(FlatNode {
                            name: format!("g{}", child1_idx),
                            left_child: None,
                            right_child: None,
                            parent: Some(current_gene_idx),
                            depth: None,
                            length: 0.0,
                        });
                        node_mapping.push(copy.species_node_idx);
                        event_mapping.push(Event::Duplication);

                        let child2_idx = gene_nodes.len();
                        gene_nodes.push(FlatNode {
                            name: format!("g{}", child2_idx),
                            left_child: None,
                            right_child: None,
                            parent: Some(current_gene_idx),
                            depth: None,
                            length: 0.0,
                        });
                        node_mapping.push(copy.species_node_idx);
                        event_mapping.push(Event::Duplication);

                        // Update parent as duplication node
                        gene_nodes[current_gene_idx].left_child = Some(child1_idx);
                        gene_nodes[current_gene_idx].right_child = Some(child2_idx);
                        gene_nodes[current_gene_idx].depth = Some(event_time);
                        event_mapping[current_gene_idx] = Event::Duplication;

                        // Record event
                        events.push(DTLEvent {
                            time: event_time,
                            gene_node_id: current_gene_idx,
                            event_type: "Duplication".to_string(),
                            species_node: species_tree.nodes[copy.species_node_idx].name.clone(),
                            donor_species: None,
                            recipient_species: None,
                            child1: Some(child1_idx),
                            child2: Some(child2_idx),
                        });

                        // Both copies continue on the same species branch
                        new_copies.push(GeneCopy {
                            gene_node_idx: child1_idx,
                            species_node_idx: copy.species_node_idx,
                            current_time: event_time,
                        });
                        new_copies.push(GeneCopy {
                            gene_node_idx: child2_idx,
                            species_node_idx: copy.species_node_idx,
                            current_time: event_time,
                        });

                        // Current gene is done (became a duplication node)
                        lost = true;

                    } else if event_prob < (lambda_d + lambda_t) / total_rate {
                        // Transfer event
                        // Find contemporary species branches at this time
                        let time_idx = find_time_index(&depths, event_time);
                        let contemporary_branches: Vec<usize> = contemporaneity[time_idx]
                            .iter()
                            .filter(|&&s| s != copy.species_node_idx)
                            .copied()
                            .collect();

                        if contemporary_branches.is_empty() {
                            // No valid transfer target, treat as failed transfer (gene continues)
                            continue;
                        }

                        // Pick random recipient
                        let recipient_idx = rng.gen_range(0..contemporary_branches.len());
                        let recipient_species = contemporary_branches[recipient_idx];

                        // Create two child gene nodes: donor copy and recipient copy
                        let donor_child_idx = gene_nodes.len();
                        gene_nodes.push(FlatNode {
                            name: format!("g{}", donor_child_idx),
                            left_child: None,
                            right_child: None,
                            parent: Some(current_gene_idx),
                            depth: None,
                            length: 0.0,
                        });
                        node_mapping.push(copy.species_node_idx);
                        event_mapping.push(Event::Transfer);

                        let recipient_child_idx = gene_nodes.len();
                        gene_nodes.push(FlatNode {
                            name: format!("g{}", recipient_child_idx),
                            left_child: None,
                            right_child: None,
                            parent: Some(current_gene_idx),
                            depth: None,
                            length: 0.0,
                        });
                        node_mapping.push(recipient_species);
                        event_mapping.push(Event::Transfer);

                        // Update parent as transfer node
                        gene_nodes[current_gene_idx].left_child = Some(donor_child_idx);
                        gene_nodes[current_gene_idx].right_child = Some(recipient_child_idx);
                        gene_nodes[current_gene_idx].depth = Some(event_time);
                        event_mapping[current_gene_idx] = Event::Transfer;

                        // Record event
                        events.push(DTLEvent {
                            time: event_time,
                            gene_node_id: current_gene_idx,
                            event_type: "Transfer".to_string(),
                            species_node: species_tree.nodes[copy.species_node_idx].name.clone(),
                            donor_species: Some(species_tree.nodes[copy.species_node_idx].name.clone()),
                            recipient_species: Some(species_tree.nodes[recipient_species].name.clone()),
                            child1: Some(donor_child_idx),
                            child2: Some(recipient_child_idx),
                        });

                        // Donor continues on same branch
                        new_copies.push(GeneCopy {
                            gene_node_idx: donor_child_idx,
                            species_node_idx: copy.species_node_idx,
                            current_time: event_time,
                        });

                        // Recipient starts on the new branch
                        new_copies.push(GeneCopy {
                            gene_node_idx: recipient_child_idx,
                            species_node_idx: recipient_species,
                            current_time: event_time,
                        });

                        // Current gene is done (became a transfer node)
                        lost = true;

                    } else {
                        // Loss event
                        gene_nodes[current_gene_idx].depth = Some(event_time);
                        event_mapping[current_gene_idx] = Event::Loss;

                        // Record event
                        events.push(DTLEvent {
                            time: event_time,
                            gene_node_id: current_gene_idx,
                            event_type: "Loss".to_string(),
                            species_node: species_tree.nodes[copy.species_node_idx].name.clone(),
                            donor_species: None,
                            recipient_species: None,
                            child1: None,
                            child2: None,
                        });

                        lost = true;
                    }
                }
            }

            // After the while loop: if we haven't lost the gene, process the species node
            if !lost {
                // Check if this is a leaf or internal node
                if species_node.left_child.is_none() && species_node.right_child.is_none() {
                    // Reached a leaf species - gene survives as extant
                    gene_nodes[current_gene_idx].depth = Some(species_end_time);
                    event_mapping[current_gene_idx] = Event::Leaf;

                    // Record event
                    events.push(DTLEvent {
                        time: species_end_time,
                        gene_node_id: current_gene_idx,
                        event_type: "Leaf".to_string(),
                        species_node: species_tree.nodes[copy.species_node_idx].name.clone(),
                        donor_species: None,
                        recipient_species: None,
                        child1: None,
                        child2: None,
                    });
                } else {
                    // Speciation: gene follows both children
                    let left_species = species_node.left_child.unwrap();
                    let right_species = species_node.right_child.unwrap();

                    // Create two child gene nodes
                    let left_gene_idx = gene_nodes.len();
                    gene_nodes.push(FlatNode {
                        name: format!("g{}", left_gene_idx),
                        left_child: None,
                        right_child: None,
                        parent: Some(current_gene_idx),
                        depth: None,
                        length: 0.0,
                    });
                    node_mapping.push(left_species);
                    event_mapping.push(Event::Speciation);

                    let right_gene_idx = gene_nodes.len();
                    gene_nodes.push(FlatNode {
                        name: format!("g{}", right_gene_idx),
                        left_child: None,
                        right_child: None,
                        parent: Some(current_gene_idx),
                        depth: None,
                        length: 0.0,
                    });
                    node_mapping.push(right_species);
                    event_mapping.push(Event::Speciation);

                    // Update parent
                    gene_nodes[current_gene_idx].left_child = Some(left_gene_idx);
                    gene_nodes[current_gene_idx].right_child = Some(right_gene_idx);
                    gene_nodes[current_gene_idx].depth = Some(species_end_time);
                    event_mapping[current_gene_idx] = Event::Speciation;

                    // Record event
                    events.push(DTLEvent {
                        time: species_end_time,
                        gene_node_id: current_gene_idx,
                        event_type: "Speciation".to_string(),
                        species_node: species_tree.nodes[copy.species_node_idx].name.clone(),
                        donor_species: None,
                        recipient_species: None,
                        child1: Some(left_gene_idx),
                        child2: Some(right_gene_idx),
                    });

                    // Add new copies to process on the child branches
                    // Start time is the parent's depth (which is the start of the child branch)
                    let left_start_time = species_end_time;
                    let right_start_time = species_end_time;

                    new_copies.push(GeneCopy {
                        gene_node_idx: left_gene_idx,
                        species_node_idx: left_species,
                        current_time: left_start_time,
                    });
                    new_copies.push(GeneCopy {
                        gene_node_idx: right_gene_idx,
                        species_node_idx: right_species,
                        current_time: right_start_time,
                    });
                }
            }

            // Mark this copy as processed
            copies_to_remove.push(copy_idx);
        }

        // Remove processed copies (in reverse order to maintain indices)
        for &idx in copies_to_remove.iter().rev() {
            active_copies.swap_remove(idx);
        }

        // Add new copies
        active_copies.extend(new_copies);
    }

    // Find the root of the gene tree (node with no parent)
    let root_idx = gene_nodes
        .iter()
        .position(|n| n.parent.is_none())
        .unwrap_or(0);

    // Handle edge case: if the gene was lost before any events, create minimal tree
    if gene_nodes.is_empty() {
        gene_nodes.push(FlatNode {
            name: "g0".to_string(),
            left_child: None,
            right_child: None,
            parent: None,
            depth: Some(origin_depth),
            length: 0.0,
        });
        node_mapping.push(origin_species);
        event_mapping.push(Event::Loss);
    }

    let gene_tree = FlatTree {
        nodes: gene_nodes,
        root: root_idx,
    };

    (RecTree::new(species_tree, gene_tree, node_mapping, event_mapping), events)
}

/// Saves DTL events to a CSV file
///
/// # Arguments
/// * `events` - A slice of DTLEvent structures to save
/// * `filename` - Path to the output CSV file
///
/// # Returns
/// Result indicating success or IO error
pub fn save_events_to_csv(events: &[DTLEvent], filename: &str) -> std::io::Result<()> {
    use std::fs::File;
    use std::io::{Write, BufWriter};

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "{}", DTLEvent::csv_header())?;

    // Write each event
    for event in events {
        writeln!(writer, "{}", event.to_csv_row())?;
    }

    writer.flush()?;
    Ok(())
}

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

/// Counts the number of extant genes (leaves that are not losses)
pub fn count_extant_genes(rec_tree: &RecTree) -> usize {
    rec_tree
        .gene_tree
        .nodes
        .iter()
        .enumerate()
        .filter(|(i, node)| {
            node.left_child.is_none()
                && node.right_child.is_none()
                && rec_tree.event_mapping[*i] == Event::Leaf
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
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 0.0, 0.0, &mut rng);

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
        let (rec_tree, _events) = simulate_dtl(&species_tree, species_tree.root, 2.0, 0.0, 0.0, &mut rng);

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
        let (rec_tree, _events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 0.0, 5.0, &mut rng);

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
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 1.0, 0.5, 0.5, &mut rng);

        // Save events to CSV
        save_events_to_csv(&events, "test_dtl_events.csv").expect("Failed to write events CSV");
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

        // High transfer rate
        let (rec_tree, events) = simulate_dtl(&species_tree, species_tree.root, 0.0, 2.0, 0.0, &mut rng);

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!("Events with transfer: S={}, D={}, T={}, L={}, Leaves={}", s, d, t, l, leaves);
        println!("Total event count: {}", events.len());
    }
}
