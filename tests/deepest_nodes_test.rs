use std::fs;
use std::path::PathBuf;
use pest::Parser;
use rustree::newick::newick::{newick_to_tree, node_to_newick, NewickParser, Rule};
use rustree::node::{FlatNode, FlatTree};
use rustree::extract_extant::find_deepest_nodes;

#[test]
fn test_find_deepest_nodes() {
    // Build a flat tree manually.
    let nodes = vec![
        FlatNode {
            name: "R".to_string(),
            left_child: Some(1),
            right_child: Some(2),
            parent: None,
            depth: Some(0.0),
            length: 0.0,
        },
        FlatNode {
            name: "A".to_string(),
            left_child: None,
            right_child: None,
            parent: Some(0),
            depth: Some(1.0),
            length: 0.1,
        },
        FlatNode {
            name: "B".to_string(),
            left_child: None,
            right_child: None,
            parent: Some(0),
            depth: Some(2.0),
            length: 0.2,
        },
    ];
    let tree = FlatTree { nodes, root: 0 };

    // Retrieve the deepest 1 node.
    let deepest_one = find_deepest_nodes(&tree, 1);
    assert_eq!(deepest_one.len(), 1);
    assert_eq!(tree.nodes[deepest_one[0]].name, "B");

    // Retrieve the deepest 2 nodes; expected order is descending by depth.
    let deepest_two = find_deepest_nodes(&tree, 2);
    assert_eq!(deepest_two.len(), 2);
    assert_eq!(tree.nodes[deepest_two[0]].name, "B");
    assert_eq!(tree.nodes[deepest_two[1]].name, "A");
}

#[test]
fn test_find_deepest_nodes_from_file() {
    // Read the Newick string from test_tree.nwk.
    let tree_path = PathBuf::from("src/test_tree.nwk");
    let newick_str = fs::read_to_string(&tree_path)
        .expect("Failed to read test_tree.nwk")
        .trim()
        .to_string();

    // Parse the Newick string into a Node tree.
    let mut pairs = NewickParser::parse(Rule::newick, &newick_str)
        .expect("Error parsing Newick string");
    let mut node_tree = newick_to_tree(pairs.next().unwrap());
    let mut root_node = node_tree.pop().expect("Failed to pop root node");

    // Prepare the tree by zeroing and assigning depths.
    root_node.zero_root_length();
    root_node.assign_depths(0.0);

    // Convert the Node tree into a FlatTree.
    let flat_tree = root_node.to_flat_tree();

    // Use find_deepest_nodes to extract the deepest nodes.
    let deepest_one = find_deepest_nodes(&flat_tree, 1);
    println!("Deepest one node kept: {}", flat_tree.nodes[deepest_one[0]].name);

    let deepest_two = find_deepest_nodes(&flat_tree, 2);
    let kept_names: Vec<_> = deepest_two.iter().map(|&i| flat_tree.nodes[i].name.clone()).collect();
    println!("Deepest two nodes kept (in order): {:?}", kept_names);

    // Reconvert the FlatTree back to a Node tree.
    let reconstructed_node_tree = flat_tree.to_node();
    let newick_reconstructed = node_to_newick(&reconstructed_node_tree) + ";";
    println!("Reconstructed Newick string: {}", newick_reconstructed);

    // Verification assertions remain.
    assert_eq!(deepest_one.len(), 1);
    assert_eq!(flat_tree.nodes[deepest_one[0]].name, "B");
    assert_eq!(deepest_two.len(), 2);
    assert_eq!(flat_tree.nodes[deepest_two[0]].name, "B");
    assert_eq!(flat_tree.nodes[deepest_two[1]].name, "A");
}
