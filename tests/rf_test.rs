// tests/rf_comparison_table.rs

use std::fs;
use std::path::{Path, PathBuf};
use pest::Parser;
use rustree::node::FlatTree;
use rustree::robinson_foulds::unrooted_robinson_foulds;
use rustree::newick::newick::{NewickParser, Rule, newick_to_tree};

/// Reads a Newick file (which may contain several trees separated by ';')
/// and returns the first tree as a FlatTree.
///
/// This helper mirrors your usual file‑opening logic.
fn flat_tree_from_file<P: AsRef<Path>>(file_path: P) -> FlatTree {
    // Read the contents of the file.
    let tree_str = fs::read_to_string(&file_path)
        .unwrap_or_else(|_| panic!("Failed to read tree file: {:?}", file_path.as_ref()));
    let sanitized = tree_str.trim();
    // Split by semicolon, trim, and append the semicolon back.
    let trees: Vec<String> = sanitized
        .split(';')
        .filter_map(|s| {
            let s = s.trim();
            if s.is_empty() {
                None
            } else {
                Some(format!("{};", s))
            }
        })
        .collect();
    if trees.is_empty() {
        panic!("No tree found in file: {:?}", file_path.as_ref());
    }
    // Use the first valid tree.
    let tree_newick = &trees[0];
    let pairs = NewickParser::parse(Rule::newick, tree_newick)
        .expect("Failed to parse Newick tree");
    let mut node_tree = newick_to_tree(pairs.into_iter().next().unwrap())
        .pop()
        .expect("No tree produced");
    node_tree.to_flat_tree()
}

#[test]
fn test_rf_distance_comparison_table() {
    // Directory containing the Newick trees and the expected RF results.
    let dir_path = "/home/enzo/Documents/git/rustree/tests/RF_trees";
    let dir = Path::new(dir_path);

    // Collect all files with names containing "tree_" and with extension ".newick".
    let mut tree_files: Vec<PathBuf> = fs::read_dir(dir)
        .expect("Failed to read directory")
        .filter_map(|entry| {
            let entry = entry.expect("Error reading directory entry");
            let path = entry.path();
            if path.is_file()
                && path.extension().map(|s| s == "newick").unwrap_or(false)
                && path.to_string_lossy().contains("tree_")
            {
                Some(path)
            } else {
                None
            }
        })
        .collect();
    // Sort the files to ensure a consistent order.
    tree_files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));

    // Parse each file into a FlatTree.
    let flat_trees: Vec<FlatTree> = tree_files
        .iter()
        .map(|path| flat_tree_from_file(path))
        .collect();

    // Print table header.
    println!("+---------+---------+------------+------------+");
    println!("| Tree i  | Tree j  | RF (ape)   | RF (Rust)  |");
    println!("+---------+---------+------------+------------+");

    // Iterate over every unique tree pair.
    for i in 0..flat_trees.len() {
        for j in (i + 1)..flat_trees.len() {
            // Expected RF distance is stored in a file named "rf_tree_XX_YY.txt"
            let expected_file = format!("{}/rf_tree_{:02}_{:02}.txt", dir_path, i + 1, j + 1);
            let expected_rf_str = fs::read_to_string(&expected_file)
                .unwrap_or_else(|_| panic!("Failed to read expected RF file: {}", expected_file));
            let expected_rf: usize = expected_rf_str
                .trim()
                .parse()
                .unwrap_or_else(|_| panic!("Expected RF value in {} is not a number", expected_file));

            // Compute the RF distance using your Rust implementation.
            let computed_rf = unrooted_robinson_foulds(&flat_trees[i], &flat_trees[j]);

            // Print the comparison row.
            println!(
                "| {:>7} | {:>7} | {:>10} | {:>10} |",
                i + 1,
                j + 1,
                expected_rf,
                computed_rf
            );
        }
    }
    println!("+---------+---------+------------+------------+");
}
