// tests/rf_test.rs
//
// Robinson-Foulds distance comparisons against unrooted reference values.
//
// The historical `rf_tree_XX_YY.txt` fixtures in `tests/RF_trees/` were
// generated from rooted `ape::dist.topo` inputs and do not match the fixed
// unrooted RF implementation. The reference values below were re-checked
// against `ape::dist.topo(c(unroot(t1), unroot(t2)), method = "PH85")`.

use rustree::node::FlatTree;
use rustree::parse_newick;
use rustree::robinson_foulds::{robinson_foulds, unrooted_robinson_foulds};
use std::fs;
use std::path::{Path, PathBuf};

fn flat_tree_from_file<P: AsRef<Path>>(file_path: P) -> FlatTree {
    let tree_str = fs::read_to_string(&file_path)
        .unwrap_or_else(|_| panic!("Failed to read tree file: {:?}", file_path.as_ref()));
    let sanitized = tree_str.trim();
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
    let tree_newick = &trees[0];
    let mut nodes = parse_newick(tree_newick).expect("Failed to parse Newick tree");
    let root = nodes.pop().expect("No tree produced");
    root.to_flat_tree()
}

fn flat_tree_by_number(dir_path: &Path, number: usize) -> FlatTree {
    let path = dir_path.join(format!("tree_{:02}.newick", number));
    flat_tree_from_file(path)
}

#[test]
fn test_unrooted_rf_matches_reference_pairs() {
    let dir_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/RF_trees");
    if !dir_path.exists() {
        eprintln!("skipping RF fixture test: missing {}", dir_path.display());
        return;
    }

    let reference_pairs = [
        (1, 2, 12),
        (1, 7, 8),
        (1, 14, 12),
        (1, 71, 12),
        (44, 74, 14),
    ];

    for (left, right, expected_rf) in reference_pairs {
        let tree1 = flat_tree_by_number(&dir_path, left);
        let tree2 = flat_tree_by_number(&dir_path, right);
        let computed_rf = unrooted_robinson_foulds(&tree1, &tree2).unwrap();

        assert_eq!(
            computed_rf, expected_rf,
            "Unrooted RF mismatch for trees {} vs {}: expected {}, got {}",
            left, right, expected_rf, computed_rf
        );
    }
}

#[test]
fn test_unrooted_rf_root_placement_invariant() {
    let mut tree1_nodes = parse_newick("((A:1,B:1):1,(C:1,D:1):1):0;").expect("tree1 parse failed");
    let mut tree2_nodes = parse_newick("(A:1,(B:1,(C:1,D:1):1):1):0;").expect("tree2 parse failed");
    let tree1 = tree1_nodes
        .pop()
        .expect("tree1 missing root")
        .to_flat_tree();
    let tree2 = tree2_nodes
        .pop()
        .expect("tree2 missing root")
        .to_flat_tree();

    assert_eq!(robinson_foulds(&tree1, &tree2).unwrap(), 2);
    assert_eq!(unrooted_robinson_foulds(&tree1, &tree2).unwrap(), 0);
}
