use std::fs;
use std::path::PathBuf;
use pest::Parser;
use rustree::node::{TraversalOrder, compare_nodes};
use rustree::newick::newick::{newick_to_tree, node_to_newick, NewickParser, Rule};
use rustree::sampling::remove_node;



#[test]
fn test_remove_node() {
    // Read the Newick string from the test file
    let newick_path = PathBuf::from("tests/test_tree.nwk");
    let newick_str = fs::read_to_string(newick_path).expect("Failed to read Newick file");

    // Parse the Newick string into a Node tree
    let mut pairs = NewickParser::parse(Rule::newick, &newick_str).expect("Failed to parse Newick string");
    let mut node_tree = newick_to_tree(pairs.next().unwrap());
    let mut root_node = node_tree.pop().unwrap();

    // Convert the Node tree to a FlatTree
    root_node.zero_root_length();
    root_node.assign_depths(0.0);
    let mut flat_tree = root_node.to_flat_tree();

    // Find the index of node "F" and remove it.
    let index_f = flat_tree.iter(TraversalOrder::PreOrder)
                           .position(|node| node.name == "F")
                           .expect("Node 'F' not found");
    remove_node(&mut flat_tree, index_f);

    // Reconstruct the tree and convert it back to Newick (for debug)
    let mut reconstructed_tree = flat_tree.to_node();
    reconstructed_tree.assign_depths(0.0);
    let newick_result = node_to_newick(&reconstructed_tree) + ";";
    println!("Resulting Newick string: {}", newick_result);

    // Now read expected tree from test_tree_F.nwk and compare Nodes
    let expected_path = PathBuf::from("tests/test_tree_F.nwk");
    let expected_newick = fs::read_to_string(expected_path)
                           .expect("Failed to read expected Newick file");
    let mut expected_pairs = NewickParser::parse(Rule::newick, &expected_newick)
                           .expect("Failed to parse expected Newick string");
    let mut expected_tree = newick_to_tree(expected_pairs.next().unwrap()).pop().unwrap();
    expected_tree.assign_depths(0.0);
    // Reconvert both trees to Newick and display them.
    let recon_newick = node_to_newick(&reconstructed_tree) + ";";
    let expected_recon_newick = node_to_newick(&expected_tree) + ";";
    println!("Reconstructed Newick string: {}", recon_newick);
    println!("Expected Newick string: {}", expected_recon_newick);

    // Compare the reconstructed tree and the expected tree (ignoring children permutation)
    if !compare_nodes(&reconstructed_tree, &expected_tree) {
        println!("Reconstructed Node: {:?}", reconstructed_tree);
        println!("Expected Node: {:?}", expected_tree);
        panic!("Reconstructed tree does not match the expected tree.");
    } else {
        // Show both newick strings if the trees match
        println!("Reconstructed tree matches the expected tree.");
    }
}
