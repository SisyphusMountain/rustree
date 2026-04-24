use crate::error::RustreeError;
use crate::node::FlatTree;
use std::collections::{BTreeSet, HashSet};

/// Computes the set of leaf names in the FlatTree.
/// A leaf is defined as a node with no left or right child.
fn collect_leaves(tree: &FlatTree) -> BTreeSet<String> {
    let mut leaves = BTreeSet::new();
    for node in &tree.nodes {
        if node.left_child.is_none() && node.right_child.is_none() {
            leaves.insert(node.name.clone());
        }
    }
    leaves
}

/// Recursively computes the set of leaf names for the subtree rooted at `index`.
fn subtree_leaves(tree: &FlatTree, index: usize) -> BTreeSet<String> {
    let node = &tree.nodes[index];
    if node.left_child.is_none() && node.right_child.is_none() {
        let mut set = BTreeSet::new();
        set.insert(node.name.clone());
        return set;
    }
    let mut leaves = BTreeSet::new();
    if let Some(left_index) = node.left_child {
        leaves.extend(subtree_leaves(tree, left_index));
    }
    if let Some(right_index) = node.right_child {
        leaves.extend(subtree_leaves(tree, right_index));
    }
    leaves
}

fn encode_partition(leaves: &BTreeSet<String>) -> String {
    leaves.iter().cloned().collect::<Vec<_>>().join(",")
}

fn validate_matching_leaves(
    tree1: &FlatTree,
    tree2: &FlatTree,
) -> Result<BTreeSet<String>, RustreeError> {
    let leaves1 = collect_leaves(tree1);
    let leaves2 = collect_leaves(tree2);
    if leaves1 != leaves2 {
        return Err(RustreeError::Validation(
            "The trees do not have identical leaf labels.".to_string(),
        ));
    }
    Ok(leaves1)
}

/// Extracts rooted clades from a tree.
///
/// Each non-root internal node defines a clade in rooted RF distance.
fn get_rooted_clades(tree: &FlatTree) -> HashSet<String> {
    let mut clades = HashSet::new();

    for (index, node) in tree.nodes.iter().enumerate() {
        if node.parent.is_none() {
            continue;
        }
        if node.left_child.is_none() && node.right_child.is_none() {
            continue;
        }

        let sub_leaves = subtree_leaves(tree, index);
        clades.insert(encode_partition(&sub_leaves));
    }

    clades
}

/// Extracts unrooted bipartitions from a rooted tree.
///
/// Each internal edge defines a bipartition of the leaf set. We represent each
/// bipartition by the smaller side (by cardinality, then lexicographically).
/// This makes the result independent of root placement.
fn get_unrooted_bipartitions(tree: &FlatTree, all_leaves: &BTreeSet<String>) -> HashSet<String> {
    let n = all_leaves.len();
    let mut bipartitions = HashSet::new();

    for (index, node) in tree.nodes.iter().enumerate() {
        if node.parent.is_none() {
            continue;
        } // skip root
        if node.left_child.is_none() && node.right_child.is_none() {
            continue;
        } // skip leaves

        let clade = subtree_leaves(tree, index);
        // Normalize: use the smaller side of the bipartition
        let side = if clade.len() * 2 < n {
            clade
        } else if clade.len() * 2 > n {
            all_leaves.difference(&clade).cloned().collect()
        } else {
            // Equal size: pick lexicographically smaller
            let complement: BTreeSet<String> = all_leaves.difference(&clade).cloned().collect();
            if clade < complement {
                clade
            } else {
                complement
            }
        };

        if side.len() > 1 && side.len() < n - 1 {
            bipartitions.insert(encode_partition(&side));
        }
    }
    bipartitions
}

fn symmetric_difference_size(left: &HashSet<String>, right: &HashSet<String>) -> usize {
    let common: HashSet<_> = left.intersection(right).collect();
    (left.len() - common.len()) + (right.len() - common.len())
}

/// Computes the rooted Robinson-Foulds distance between two trees.
///
pub fn robinson_foulds(tree1: &FlatTree, tree2: &FlatTree) -> Result<usize, RustreeError> {
    validate_matching_leaves(tree1, tree2)?;
    let clades1 = get_rooted_clades(tree1);
    let clades2 = get_rooted_clades(tree2);
    Ok(symmetric_difference_size(&clades1, &clades2))
}

/// Computes the unrooted Robinson-Foulds distance between two trees.
///
/// Compares bipartitions rather than rooted clades, so two trees that differ
/// only in root placement have RF distance 0.
pub fn unrooted_robinson_foulds(tree1: &FlatTree, tree2: &FlatTree) -> Result<usize, RustreeError> {
    let all_leaves = validate_matching_leaves(tree1, tree2)?;
    let bip1 = get_unrooted_bipartitions(tree1, &all_leaves);
    let bip2 = get_unrooted_bipartitions(tree2, &all_leaves);
    Ok(symmetric_difference_size(&bip1, &bip2))
}

/// Backward-compatible alias for the validated unrooted RF implementation.
pub fn true_unrooted_robinson_foulds(
    tree1: &FlatTree,
    tree2: &FlatTree,
) -> Result<usize, RustreeError> {
    unrooted_robinson_foulds(tree1, tree2)
}

#[cfg(test)]
mod tests {
    use super::{robinson_foulds, unrooted_robinson_foulds};
    use crate::newick::parse_newick;

    fn flat_tree(newick: &str) -> crate::node::FlatTree {
        let mut nodes = parse_newick(newick).expect("failed to parse test tree");
        nodes.pop().expect("missing root").to_flat_tree()
    }

    #[test]
    fn test_unrooted_rf_is_root_invariant() {
        let tree1 = flat_tree("((A:1,B:1):1,(C:1,D:1):1):0;");
        let tree2 = flat_tree("(A:1,(B:1,(C:1,D:1):1):1):0;");

        assert_eq!(robinson_foulds(&tree1, &tree2).unwrap(), 2);
        assert_eq!(unrooted_robinson_foulds(&tree1, &tree2).unwrap(), 0);
    }

    #[test]
    fn test_rf_rejects_mismatched_leaf_sets() {
        let tree1 = flat_tree("((A:1,B:1):1,C:2):0;");
        let tree2 = flat_tree("((A:1,B:1):1,D:2):0;");

        let err = unrooted_robinson_foulds(&tree1, &tree2).unwrap_err();
        assert!(err.to_string().contains("identical leaf labels"));
    }
}
