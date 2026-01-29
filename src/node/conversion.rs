//! Conversion functions between recursive Node and flat FlatTree representations.

use super::{Node, FlatNode, FlatTree};

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
