//! Conversion functions between recursive Node and flat FlatTree representations.

use super::{Node, FlatNode, FlatTree, RecTreeOwned};
use std::collections::HashMap;

/// Converts an arborescent tree (a tree where each node owns its children)
/// into a flat structure (a vector of FlatNodes).
///
/// # Arguments
/// * `node` - The root node of the tree.
/// * `flat_tree` - The vector to be filled with flat nodes.
/// * `parent` - The parent index (if any).
///
/// # Returns
/// The index of the node in the flat tree.
pub fn node_to_flat(node: &Node, flat_tree: &mut Vec<FlatNode>, parent: Option<usize>) -> usize {
    let index = flat_tree.len();
    flat_tree.push(FlatNode {
        name: node.name.clone(),
        left_child: None,
        right_child: None,
        parent,
        depth: node.depth,
        length: node.length,
        bd_event: None,
    });

    if let Some(left) = &node.left_child {
        let left_index = node_to_flat(left, flat_tree, Some(index));
        flat_tree[index].left_child = Some(left_index);
    }

    if let Some(right) = &node.right_child {
        let right_index = node_to_flat(right, flat_tree, Some(index));
        flat_tree[index].right_child = Some(right_index);
    }

    index
}

/// Converts a flat tree into a nested `Node` tree structure.
///
/// # Arguments
/// * `flat_tree` - The vector representing the flat tree.
/// * `index` - The index of the node in the flat tree.
///
/// # Returns
/// The corresponding `Node` tree.
pub fn flat_to_node(flat_tree: &[FlatNode], index: usize) -> Option<Node> {
    let flat_node = &flat_tree[index];
    let left_child = flat_node
        .left_child
        .and_then(|i| flat_to_node(flat_tree, i).map(Box::new));
    let right_child = flat_node
        .right_child
        .and_then(|i| flat_to_node(flat_tree, i).map(Box::new));

    Some(Node {
        name: flat_node.name.clone(),
        left_child,
        right_child,
        depth: flat_node.depth,
        length: flat_node.length,
    })
}

// Methods on Node for conversion
impl Node {
    /// Converts a recursive `Node` structure into a `FlatTree`.
    pub fn to_flat_tree(&self) -> FlatTree {
        let mut flat_nodes = Vec::new();
        let root_index = self.node_to_flat_internal(&mut flat_nodes, None);
        FlatTree {
            nodes: flat_nodes,
            root: root_index,
        }
    }

    /// Internal helper method for converting a `Node` to flat nodes.
    fn node_to_flat_internal(&self, flat_nodes: &mut Vec<FlatNode>, parent_index: Option<usize>) -> usize {
        let index = flat_nodes.len();
        flat_nodes.push(FlatNode {
            name: self.name.clone(),
            left_child: None,
            right_child: None,
            parent: parent_index,
            depth: self.depth,
            length: self.length,
            bd_event: None,
        });

        if let Some(ref left_child) = self.left_child {
            let left_index = left_child.node_to_flat_internal(flat_nodes, Some(index));
            flat_nodes[index].left_child = Some(left_index);
        }

        if let Some(ref right_child) = self.right_child {
            let right_index = right_child.node_to_flat_internal(flat_nodes, Some(index));
            flat_nodes[index].right_child = Some(right_index);
        }

        index
    }
}

// Methods on FlatTree for conversion
impl FlatTree {
    /// Converts a `FlatTree` into a recursive `Node` structure.
    pub fn to_node(&self) -> Node {
        self.flat_to_node_internal(self.root)
    }

    /// Find a node by name and return its index.
    pub fn find_node_index(&self, name: &str) -> Option<usize> {
        self.nodes.iter().position(|n| n.name == name)
    }

    /// Internal helper method for converting flat nodes to a `Node`.
    fn flat_to_node_internal(&self, index: usize) -> Node {
        let flat_node = &self.nodes[index];

        let left_child = flat_node.left_child.map(|i| {
            Box::new(self.flat_to_node_internal(i))
        });

        let right_child = flat_node.right_child.map(|i| {
            Box::new(self.flat_to_node_internal(i))
        });

        Node {
            name: flat_node.name.clone(),
            left_child,
            right_child,
            depth: flat_node.depth,
            length: flat_node.length,
        }
    }

    /// Returns the number of nodes in the tree.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Returns true if the tree has no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

/// Helper function to collect node indices in postorder traversal.
fn collect_postorder_indices(tree: &FlatTree, index: usize, result: &mut Vec<usize>) {
    let node = &tree.nodes[index];

    // Visit left subtree
    if let Some(left_idx) = node.left_child {
        collect_postorder_indices(tree, left_idx, result);
    }

    // Visit right subtree
    if let Some(right_idx) = node.right_child {
        collect_postorder_indices(tree, right_idx, result);
    }

    // Visit current node (postorder: left, right, node)
    result.push(index);
}

/// Obtain mapping between two trees with identical topology using postorder traversal.
///
/// This function matches nodes between an original tree and a reconciled tree by comparing
/// their positions in postorder traversal. Trees with the same topology will have corresponding
/// nodes at the same positions in their postorder sequences.
///
/// # Arguments
/// * `original_tree` - The original tree with desired node names
/// * `reconciled_tree` - The reconciled tree (e.g., from ALERax) to be mapped
///
/// # Returns
/// A HashMap where keys are reconciled tree node indices and values are original tree node indices.
/// This maps each node in the reconciled tree to its corresponding node in the original tree.
///
/// # Errors
/// Returns an error if the trees have different topologies (different number of nodes or structure).
///
/// # Example
/// ```rust
/// let mapping = obtain_mapping(&original_tree, &reconciled_tree)?;
/// // For each node in reconciled tree, mapping[idx] gives corresponding original tree node
/// ```
pub fn obtain_mapping(
    original_tree: &FlatTree,
    reconciled_tree: &FlatTree,
) -> Result<HashMap<usize, usize>, String> {
    // Check if trees have the same number of nodes
    if original_tree.nodes.len() != reconciled_tree.nodes.len() {
        return Err(format!(
            "Trees have different number of nodes: original={}, reconciled={}",
            original_tree.nodes.len(),
            reconciled_tree.nodes.len()
        ));
    }

    // Perform postorder traversal on both trees and collect node indices
    let mut original_postorder = Vec::new();
    collect_postorder_indices(original_tree, original_tree.root, &mut original_postorder);

    let mut reconciled_postorder = Vec::new();
    collect_postorder_indices(reconciled_tree, reconciled_tree.root, &mut reconciled_postorder);

    // Verify both traversals have the same length
    if original_postorder.len() != reconciled_postorder.len() {
        return Err("Trees have different structures (postorder lengths differ)".to_string());
    }

    // Create mapping: reconciled_idx -> original_idx
    let mut mapping = HashMap::new();
    for (pos, (&reconciled_idx, &original_idx)) in reconciled_postorder
        .iter()
        .zip(original_postorder.iter())
        .enumerate()
    {
        // Verify structural compatibility by checking if both are leaves or both are internal
        let orig_is_leaf = original_tree.nodes[original_idx].left_child.is_none()
            && original_tree.nodes[original_idx].right_child.is_none();
        let recon_is_leaf = reconciled_tree.nodes[reconciled_idx].left_child.is_none()
            && reconciled_tree.nodes[reconciled_idx].right_child.is_none();

        if orig_is_leaf != recon_is_leaf {
            return Err(format!(
                "Structural mismatch at position {}: original node {} and reconciled node {} differ in leaf/internal status",
                pos, original_idx, reconciled_idx
            ));
        }

        mapping.insert(reconciled_idx, original_idx);
    }

    Ok(mapping)
}

/// Rename a reconciled tree to match the original tree's node names.
///
/// This function renames all nodes in a reconciled tree (e.g., from ALERax) to match the
/// corresponding nodes in the original tree. It uses postorder traversal to match nodes
/// with identical topology.
///
/// The reconciliation data (node_mapping and event_mapping) are indexed by gene tree node
/// indices, so they remain valid after renaming as only names change, not indices.
///
/// # Arguments
/// * `original_tree` - The original tree with desired node names
/// * `reconciled_tree_owned` - Mutable reference to the reconciled tree to rename
///
/// # Returns
/// Ok(()) on success, or an error if trees have incompatible topologies.
///
/// # Example
/// ```rust
/// // After ALERax reconciliation
/// rename_reconciled_tree(&original_tree, &mut reconciled_tree)?;
/// // Now reconciled_tree.gene_tree has the same node names as original_tree
/// ```
pub fn rename_reconciled_tree(
    original_tree: &FlatTree,
    reconciled_tree_owned: &mut RecTreeOwned,
) -> Result<(), String> {
    // Obtain mapping from reconciled to original tree
    let mapping = obtain_mapping(original_tree, &reconciled_tree_owned.gene_tree)?;

    // Rename nodes in the reconciled gene tree
    for (reconciled_idx, &original_idx) in mapping.iter() {
        let original_name = &original_tree.nodes[original_idx].name;
        reconciled_tree_owned.gene_tree.nodes[*reconciled_idx].name = original_name.clone();
    }

    Ok(())
}

#[cfg(test)]
mod mapping_tests {
    use super::*;

    #[test]
    fn test_obtain_mapping_identical_trees() {
        // Create two identical trees
        let mut nodes = vec![
            FlatNode {
                name: "A".to_string(),
                left_child: None,
                right_child: None,
                parent: Some(2),
                depth: Some(1.0),
                length: 1.0,
                bd_event: None,
            },
            FlatNode {
                name: "B".to_string(),
                left_child: None,
                right_child: None,
                parent: Some(2),
                depth: Some(1.0),
                length: 1.0,
                bd_event: None,
            },
            FlatNode {
                name: "Root".to_string(),
                left_child: Some(0),
                right_child: Some(1),
                parent: None,
                depth: Some(0.0),
                length: 0.0,
                bd_event: None,
            },
        ];

        let tree1 = FlatTree {
            nodes: nodes.clone(),
            root: 2,
        };

        // Create second tree with different names
        nodes[0].name = "X".to_string();
        nodes[1].name = "Y".to_string();
        nodes[2].name = "Top".to_string();

        let tree2 = FlatTree { nodes, root: 2 };

        // Obtain mapping
        let mapping = obtain_mapping(&tree1, &tree2).unwrap();

        // Verify mapping
        assert_eq!(mapping.len(), 3);
        assert_eq!(mapping[&0], 0); // Leaf A maps to leaf A
        assert_eq!(mapping[&1], 1); // Leaf B maps to leaf B
        assert_eq!(mapping[&2], 2); // Root maps to Root
    }

    #[test]
    fn test_obtain_mapping_different_sizes() {
        let tree1 = FlatTree {
            nodes: vec![FlatNode {
                name: "A".to_string(),
                left_child: None,
                right_child: None,
                parent: None,
                depth: Some(0.0),
                length: 0.0,
                bd_event: None,
            }],
            root: 0,
        };

        let tree2 = FlatTree {
            nodes: vec![
                FlatNode {
                    name: "X".to_string(),
                    left_child: None,
                    right_child: None,
                    parent: Some(1),
                    depth: Some(1.0),
                    length: 1.0,
                    bd_event: None,
                },
                FlatNode {
                    name: "Root".to_string(),
                    left_child: Some(0),
                    right_child: None,
                    parent: None,
                    depth: Some(0.0),
                    length: 0.0,
                    bd_event: None,
                },
            ],
            root: 1,
        };

        // Should fail due to different number of nodes
        let result = obtain_mapping(&tree1, &tree2);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("different number of nodes"));
    }
}
