//! Tree sampling utilities - extract induced subtrees from a subset of leaves.

use std::collections::{HashSet, HashMap};
use crate::node::{FlatTree, FlatNode};

/// Status of a node during induced subtree extraction.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum NodeMark {
    /// Node will be kept in the induced tree (sampled leaf or MRCA of sampled leaves)
    Keep,
    /// Node has descendants to keep but will be collapsed (not an MRCA)
    HasDescendant,
    /// Node and all its descendants will be discarded
    Discard,
}

/// Extracts the induced subtree containing only the specified leaves.
///
/// The induced subtree contains:
/// - All leaves in `keep_leaves`
/// - All internal nodes that are MRCAs of kept leaves (nodes where both subtrees have kept descendants)
///
/// Internal nodes with only one subtree containing kept leaves are collapsed
/// (their branch length is added to the descendant).
///
/// # Arguments
/// * `tree` - The original tree
/// * `keep_leaf_indices` - Set of leaf node indices to keep
///
/// # Returns
/// A new `FlatTree` containing only the induced subtree, or `None` if no leaves are kept.
pub fn extract_induced_subtree(tree: &FlatTree, keep_leaf_indices: &HashSet<usize>) -> Option<FlatTree> {
    if keep_leaf_indices.is_empty() {
        return None;
    }

    // Step 1: Mark all nodes (postorder traversal)
    let mut marks = vec![NodeMark::Discard; tree.nodes.len()];
    mark_nodes_postorder(tree, tree.root, keep_leaf_indices, &mut marks);

    // If root is discarded, no valid subtree
    if marks[tree.root] == NodeMark::Discard {
        return None;
    }

    // Step 2: Build the induced tree (preorder traversal)
    let mut new_nodes: Vec<FlatNode> = Vec::new();
    let mut old_to_new: Vec<Option<usize>> = vec![None; tree.nodes.len()];

    build_induced_tree(tree, tree.root, None, 0.0, &marks, &mut new_nodes, &mut old_to_new);

    if new_nodes.is_empty() {
        return None;
    }

    Some(FlatTree {
        nodes: new_nodes,
        root: 0, // Root is always first node added
    })
}

/// Marks nodes using postorder traversal (children before parents).
fn mark_nodes_postorder(
    tree: &FlatTree,
    node_idx: usize,
    keep_leaves: &HashSet<usize>,
    marks: &mut [NodeMark],
) {
    let node = &tree.nodes[node_idx];

    // Process children first (postorder)
    let left_mark = node.left_child.map(|c| {
        mark_nodes_postorder(tree, c, keep_leaves, marks);
        marks[c]
    });
    let right_mark = node.right_child.map(|c| {
        mark_nodes_postorder(tree, c, keep_leaves, marks);
        marks[c]
    });

    // Determine this node's mark
    let is_leaf = node.left_child.is_none() && node.right_child.is_none();

    if is_leaf {
        // Leaf: keep if in the set
        marks[node_idx] = if keep_leaves.contains(&node_idx) {
            NodeMark::Keep
        } else {
            NodeMark::Discard
        };
    } else {
        // Internal node: check children
        let left_has_kept = matches!(left_mark, Some(NodeMark::Keep | NodeMark::HasDescendant));
        let right_has_kept = matches!(right_mark, Some(NodeMark::Keep | NodeMark::HasDescendant));

        marks[node_idx] = match (left_has_kept, right_has_kept) {
            (true, true) => NodeMark::Keep,         // MRCA of kept leaves
            (true, false) | (false, true) => NodeMark::HasDescendant, // Will be collapsed
            (false, false) => NodeMark::Discard,    // Nothing to keep
        };
    }
}

/// Builds the induced tree using preorder traversal.
/// Skips HasDescendant nodes and accumulates their branch lengths.
fn build_induced_tree(
    tree: &FlatTree,
    node_idx: usize,
    new_parent: Option<usize>,
    accumulated_length: f64,
    marks: &[NodeMark],
    new_nodes: &mut Vec<FlatNode>,
    old_to_new: &mut [Option<usize>],
) {
    let node = &tree.nodes[node_idx];
    let mark = marks[node_idx];

    match mark {
        NodeMark::Discard => {
            // Skip entirely
        }
        NodeMark::HasDescendant => {
            // Skip this node, but continue to children with accumulated length
            let new_length = accumulated_length + node.length;
            if let Some(left) = node.left_child {
                build_induced_tree(tree, left, new_parent, new_length, marks, new_nodes, old_to_new);
            }
            if let Some(right) = node.right_child {
                build_induced_tree(tree, right, new_parent, new_length, marks, new_nodes, old_to_new);
            }
        }
        NodeMark::Keep => {
            // Create new node
            let new_idx = new_nodes.len();
            old_to_new[node_idx] = Some(new_idx);

            new_nodes.push(FlatNode {
                name: node.name.clone(),
                left_child: None,  // Will be set when processing children
                right_child: None,
                parent: new_parent,
                depth: node.depth,
                length: accumulated_length + node.length,
                bd_event: node.bd_event,  // Preserve event type from original node
            });

            // Update parent's child pointer
            if let Some(parent_idx) = new_parent {
                if new_nodes[parent_idx].left_child.is_none() {
                    new_nodes[parent_idx].left_child = Some(new_idx);
                } else {
                    new_nodes[parent_idx].right_child = Some(new_idx);
                }
            }

            // Process children
            if let Some(left) = node.left_child {
                build_induced_tree(tree, left, Some(new_idx), 0.0, marks, new_nodes, old_to_new);
            }
            if let Some(right) = node.right_child {
                build_induced_tree(tree, right, Some(new_idx), 0.0, marks, new_nodes, old_to_new);
            }
        }
    }
}

/// Finds leaf indices by their names.
///
/// # Arguments
/// * `tree` - The tree to search
/// * `names` - Names of leaves to find
///
/// # Returns
/// A HashSet of leaf indices matching the given names.
pub fn find_leaf_indices_by_names(tree: &FlatTree, names: &[String]) -> HashSet<usize> {
    let name_set: HashSet<&str> = names.iter().map(|s| s.as_str()).collect();

    tree.nodes
        .iter()
        .enumerate()
        .filter(|(_, node)| {
            node.left_child.is_none()
                && node.right_child.is_none()
                && name_set.contains(node.name.as_str())
        })
        .map(|(i, _)| i)
        .collect()
}

/// Finds all leaf indices in the tree.
pub fn find_all_leaf_indices(tree: &FlatTree) -> Vec<usize> {
    tree.nodes
        .iter()
        .enumerate()
        .filter(|(_, node)| node.left_child.is_none() && node.right_child.is_none())
        .map(|(i, _)| i)
        .collect()
}

/// Extracts the induced subtree keeping only leaves with the given names.
///
/// This is a convenience wrapper around `extract_induced_subtree`.
///
/// # Arguments
/// * `tree` - The original tree
/// * `leaf_names` - Names of leaves to keep
///
/// # Returns
/// A new `FlatTree` containing only the induced subtree, or `None` if no matching leaves found.
pub fn extract_induced_subtree_by_names(tree: &FlatTree, leaf_names: &[String]) -> Option<FlatTree> {
    let keep_indices = find_leaf_indices_by_names(tree, leaf_names);
    extract_induced_subtree(tree, &keep_indices)
}

/// Computes the Lowest Common Ancestor (LCA) of two nodes in a tree.
///
/// # Arguments
/// * `tree` - The tree to search
/// * `node1_idx` - Index of first node
/// * `node2_idx` - Index of second node
///
/// # Returns
/// The index of the LCA node (the deepest node that is an ancestor of both input nodes).
pub fn compute_lca(tree: &FlatTree, node1_idx: usize, node2_idx: usize) -> usize {
    // Build path from node1 to root
    let mut path1 = HashSet::new();
    let mut current = node1_idx;
    path1.insert(current);

    while let Some(parent) = tree.nodes[current].parent {
        path1.insert(parent);
        current = parent;
    }

    // Walk from node2 to root until we find a node in path1
    let mut current = node2_idx;
    if path1.contains(&current) {
        return current;
    }

    while let Some(parent) = tree.nodes[current].parent {
        if path1.contains(&parent) {
            return parent;
        }
        current = parent;
    }

    // If we get here, return the root (should always be in path1)
    tree.root
}

/// Gets all leaf names descended from a given node.
///
/// # Arguments
/// * `tree` - The tree to search
/// * `node_idx` - Index of the node
///
/// # Returns
/// A vector of leaf names descended from this node.
pub fn get_descendant_leaf_names(tree: &FlatTree, node_idx: usize) -> Vec<String> {
    let mut leaves = Vec::new();
    let mut stack = vec![node_idx];

    while let Some(idx) = stack.pop() {
        let node = &tree.nodes[idx];

        if node.left_child.is_none() && node.right_child.is_none() {
            // Leaf node
            leaves.push(node.name.clone());
        } else {
            // Internal node - add children to stack
            if let Some(left) = node.left_child {
                stack.push(left);
            }
            if let Some(right) = node.right_child {
                stack.push(right);
            }
        }
    }

    leaves
}

/// Builds a mapping from all pairs of leaf names to their LCA node index.
///
/// This creates a bidirectional mapping where (leaf1, leaf2) and (leaf2, leaf1)
/// both map to the same LCA.
///
/// # Arguments
/// * `tree` - The tree to analyze
///
/// # Returns
/// A HashMap mapping (leaf_name1, leaf_name2) pairs to LCA node indices.
pub fn build_leaf_pair_lca_map(tree: &FlatTree) -> HashMap<(String, String), usize> {
    let leaf_indices = find_all_leaf_indices(tree);
    let mut lca_map = HashMap::new();

    for i in 0..leaf_indices.len() {
        for j in (i+1)..leaf_indices.len() {
            let leaf1 = &tree.nodes[leaf_indices[i]];
            let leaf2 = &tree.nodes[leaf_indices[j]];
            let lca_idx = compute_lca(tree, leaf_indices[i], leaf_indices[j]);

            // Store both orderings for easy lookup
            lca_map.insert((leaf1.name.clone(), leaf2.name.clone()), lca_idx);
            lca_map.insert((leaf2.name.clone(), leaf1.name.clone()), lca_idx);
        }
    }

    lca_map
}

/// Maps sampled species tree node indices to original species tree node indices.
///
/// Uses LCA-based matching to identify corresponding internal nodes between trees.
///
/// # Arguments
/// * `sampled_tree` - The sampled (subset) tree
/// * `original_tree` - The original (complete) tree
/// * `sampled_lca_map` - LCA map for the sampled tree
/// * `original_lca_map` - LCA map for the original tree
///
/// # Returns
/// A HashMap mapping sampled tree indices to original tree indices.
pub fn build_sampled_to_original_mapping(
    sampled_tree: &FlatTree,
    original_tree: &FlatTree,
    _sampled_lca_map: &HashMap<(String, String), usize>,
    original_lca_map: &HashMap<(String, String), usize>,
) -> HashMap<usize, usize> {
    let mut mapping = HashMap::new();

    // For each node in sampled tree
    for sampled_idx in 0..sampled_tree.nodes.len() {
        let sampled_node = &sampled_tree.nodes[sampled_idx];

        if sampled_node.left_child.is_none() && sampled_node.right_child.is_none() {
            // Leaf node - use name-based lookup
            let original_idx = original_tree.nodes.iter()
                .position(|n| n.name == sampled_node.name)
                .expect("Leaf must exist in original tree");
            mapping.insert(sampled_idx, original_idx);
        } else {
            // Internal node - use LCA-based lookup
            let leaf_names = get_descendant_leaf_names(sampled_tree, sampled_idx);

            if leaf_names.len() >= 2 {
                // Use first two leaves to identify this internal node
                let key = (leaf_names[0].clone(), leaf_names[1].clone());
                let original_idx = original_lca_map.get(&key)
                    .expect("LCA must exist in original tree");
                mapping.insert(sampled_idx, *original_idx);
            }
        }
    }

    mapping
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::newick::newick::parse_newick;

    fn make_tree(newick: &str) -> FlatTree {
        let mut nodes = parse_newick(newick).unwrap();
        nodes.pop().unwrap().to_flat_tree()
    }

    #[test]
    fn test_extract_all_leaves() {
        // Keep all leaves - should get same topology
        let tree = make_tree("((A:1,B:1):1,C:2):0;");
        let leaves: HashSet<usize> = find_all_leaf_indices(&tree).into_iter().collect();

        let induced = extract_induced_subtree(&tree, &leaves).unwrap();
        assert_eq!(induced.nodes.len(), 5); // 3 leaves + 2 internal
    }

    #[test]
    fn test_extract_two_siblings() {
        // Keep A and B (siblings) - should collapse to just (A,B)
        let tree = make_tree("((A:1,B:1):1,C:2):0;");
        let keep = find_leaf_indices_by_names(&tree, &["A".to_string(), "B".to_string()]);

        let induced = extract_induced_subtree(&tree, &keep).unwrap();
        assert_eq!(induced.nodes.len(), 3); // A, B, and their parent
    }

    #[test]
    fn test_extract_distant_leaves() {
        // Keep A and C - should collapse intermediate node
        let tree = make_tree("((A:1,B:1):1,C:2):0;");
        let keep = find_leaf_indices_by_names(&tree, &["A".to_string(), "C".to_string()]);

        let induced = extract_induced_subtree(&tree, &keep).unwrap();
        assert_eq!(induced.nodes.len(), 3); // A, C, and root

        // A should have length 2 (1 + 1 from collapsed node)
        let a_idx = induced.nodes.iter().position(|n| n.name == "A").unwrap();
        assert_eq!(induced.nodes[a_idx].length, 2.0);
    }

    #[test]
    fn test_extract_single_leaf() {
        let tree = make_tree("((A:1,B:1):1,C:2):0;");
        let keep = find_leaf_indices_by_names(&tree, &["A".to_string()]);

        let induced = extract_induced_subtree(&tree, &keep).unwrap();
        assert_eq!(induced.nodes.len(), 1); // Just A
        assert_eq!(induced.nodes[0].name, "A");
        assert_eq!(induced.nodes[0].length, 2.0); // 1 + 1
    }

    #[test]
    fn test_extract_empty() {
        let tree = make_tree("((A:1,B:1):1,C:2):0;");
        let keep: HashSet<usize> = HashSet::new();

        assert!(extract_induced_subtree(&tree, &keep).is_none());
    }
}
