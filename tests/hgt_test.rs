use std::fs;
use std::path::PathBuf;
use pest::Parser;
use rustree::newick::newick::{NewickParser, newick_to_tree, Rule};
use rustree::node::{compare_nodes_topology, FlatTree};
use rustree::surgery::{spr_topology};
use rustree::debug::diffed_flat_tree_table;

fn run_spr_test(donor_name: &str, receiver_name: &str, expected_filename: &str) {
    // Read the original tree file.
    let orig_path = PathBuf::from("./tests/hgt_trees/original_tree.nwk");
    let orig_str = fs::read_to_string(&orig_path)
        .expect("Failed to read original tree file");
    
    // Parse the original tree.
    let mut pairs = NewickParser::parse(Rule::newick, orig_str.trim())
        .expect("Failed to parse original tree");
    let mut nodes = newick_to_tree(pairs.next().unwrap());
    let mut orig_tree = nodes.pop().expect("No tree found in original file");
    
    // Convert to a flat tree.
    let mut flat_tree = orig_tree.to_flat_tree();
    println!("Original tree:");
    let old_table = diffed_flat_tree_table(&flat_tree, None, false);
    
    // Locate donor and receiver nodes by name.
    let donor_index = flat_tree.nodes.iter().position(|node| node.name == donor_name)
        .unwrap_or_else(|| panic!("Donor node '{}' not found", donor_name));
    let receiver_index = flat_tree.nodes.iter().position(|node| node.name == receiver_name)
        .unwrap_or_else(|| panic!("Receiver node '{}' not found", receiver_name));
    
    // Perform the SPR event with a fixed time (e.g., 0.5).
    spr_topology(&mut flat_tree, donor_index, receiver_index);
    
    // For testing purposes, set all branch lengths to 1.0.
    for node in flat_tree.nodes.iter_mut() {
        node.length = 1.0;
    }
    
    // Update the root (if topology has changed).
    let new_root_index = flat_tree.nodes.iter().position(|node| node.parent.is_none())
        .expect("No root found after SPR");
    flat_tree.root = new_root_index;
    println!("New root: {}", flat_tree.nodes[new_root_index].name);
    println!("Modified tree:");
    let _new_table = diffed_flat_tree_table(&flat_tree, Some(&old_table[1..]), true);
    
    // Reconstruct the modified tree.
    let computed_tree = flat_tree.to_node();
    
    // Read and parse the expected tree.
    let expected_path = PathBuf::from(format!("./tests/hgt_trees/{}", expected_filename));
    let expected_str = fs::read_to_string(&expected_path)
        .expect("Failed to read expected tree file");
    let mut pairs = NewickParser::parse(Rule::newick, expected_str.trim())
        .expect("Failed to parse expected tree");
    let mut expected_nodes = newick_to_tree(pairs.next().unwrap());
    let expected_tree = expected_nodes.pop().expect("No expected tree found");
    
    // Compare the computed tree with the expected tree.
    if !compare_nodes_topology(&computed_tree, &expected_tree) {
        println!("Computed tree:");
        println!("{}", computed_tree.to_newick());
        println!("Expected tree:");
        println!("{}", expected_tree.to_newick());
        panic!("SPR tree for donor '{}' -> receiver '{}' does not match expected topology.", donor_name, receiver_name);
    }
}

#[test]
fn test_spr_3_to_2() {
    // For donor "3" and receiver "2", expected topology is in 3_to_2.nwk.
    run_spr_test("3", "2", "3_to_2.nwk");
}

#[test]
fn test_spr_1_to_2() {
    // For donor "1" and receiver "2", expected topology is in 1_to_2.nwk.
    run_spr_test("1", "2", "1_to_2.nwk");
}

#[test]
fn test_spr_4_to_T10() {
    // For donor "4" and receiver "T10", expected topology is in 4_to_T10.nwk.
    run_spr_test("4", "T10", "4_to_T10.nwk");
}