// rustree/src/node.rs

use std::ops::{Index, IndexMut};
#[derive(Clone, Debug)]
pub struct FlatNode {
    pub name: String,
    pub left_child: Option<usize>,
    pub right_child: Option<usize>,
    pub parent: Option<usize>,
    pub depth: Option<f64>,
    pub length: f64,
}
#[derive(Clone, Debug)]
/// Flat tree implementation (vector of nodes)
pub struct FlatTree {
    pub nodes: Vec<FlatNode>,
    pub root: usize,
}
#[derive(Clone, Debug)]
pub struct Node {
    pub name: String,
    pub left_child: Option<Box<Node>>,
    pub right_child: Option<Box<Node>>,
    pub depth: Option<f64>,
    pub length: f64,
}
// Implement HasName for &Node
impl<'a> HasName for &'a Node {
    fn name(&self) -> &str {
        &self.name
    }
}

// Implement HasName for &FlatNode
impl<'a> HasName for &'a FlatNode {
    fn name(&self) -> &str {
        &self.name
    }
}




pub enum TraversalOrder {
    PreOrder,
    InOrder,
    PostOrder,
}

pub struct NodeIter<'a> {
    stack: Vec<NodeState<'a>>,
    order: TraversalOrder,
}

enum NodeState<'a> {
    Start(&'a Node),
    Left(&'a Node),
    Right(&'a Node),
    End(&'a Node),
}

impl<'a> Iterator for NodeIter<'a> {
    type Item = &'a Node;
    // Do we need mutable reference?
    fn next(&mut self) -> Option<Self::Item>{
        while let Some(state) = self.stack.pop(){
            match state {
                NodeState::Start(node) => match self.order {
                    TraversalOrder::PreOrder => {
                        // NLR
                        self.stack.push(NodeState::Right(node));
                        self.stack.push(NodeState::Left(node));
                        return Some(node);
                    }
                    TraversalOrder::InOrder => {
                        // LNR
                        self.stack.push(NodeState::Right(node));
                        self.stack.push(NodeState::End(node));
                        self.stack.push(NodeState::Left(node));
                    }
                    TraversalOrder::PostOrder => {
                        // LRN
                        self.stack.push(NodeState::End(node));
                        self.stack.push(NodeState::Right(node));
                        self.stack.push(NodeState::Left(node));
                    }
                }
                NodeState::Left(node) => {
                    if let Some(ref left) = node.left_child {
                        self.stack.push(NodeState::Start(left));
                    }
                }
                NodeState::Right(node) => {
                    if let Some(ref right) = node.right_child {
                        self.stack.push(NodeState::Start(right));
                    }
                }

                NodeState::End(node) => match self.order {
                    TraversalOrder::InOrder => {
                        return Some(node);
                    }
                    TraversalOrder::PostOrder => {
                        return Some(node);
                    }
                    TraversalOrder::PreOrder => {} // No return there
                }
            }
        }
    return None;
}
}
impl Node {
    pub fn iter(&self, order: TraversalOrder) -> NodeIter {
        NodeIter {
            stack: vec![NodeState::Start(self)],
            order,
        }
    }
    /// Converts a recursive `Node` structure into a `FlatTree`.
    ///
    /// # Returns
    /// A `FlatTree` representation of the `Node`.
    pub fn to_flat_tree(&self) -> FlatTree {
        let mut flat_nodes = Vec::new();
        let root_index = self.node_to_flat_internal(&mut flat_nodes, None);
        FlatTree {
            nodes: flat_nodes,
            root: root_index,
        }
    }

    /// Internal helper method for converting a `Node` to flat nodes.
    ///
    /// # Arguments
    /// * `flat_nodes` - The vector to store `FlatNode` instances.
    /// * `parent_index` - The index of the parent node (if any).
    ///
    /// # Returns
    /// The index of the current node in `flat_nodes`.
    fn node_to_flat_internal(&self, flat_nodes: &mut Vec<FlatNode>, parent_index: Option<usize>) -> usize {
        let index = flat_nodes.len();
        flat_nodes.push(FlatNode {
            name: self.name.clone(),
            left_child: None,  // Will be updated after recursive calls
            right_child: None, // Will be updated after recursive calls
            parent: parent_index,
            depth: self.depth,
            length: self.length,
        });

        // Process left child
        if let Some(ref left_child) = self.left_child {
            let left_index = left_child.node_to_flat_internal(flat_nodes, Some(index));
            flat_nodes[index].left_child = Some(left_index);
        }

        // Process right child
        if let Some(ref right_child) = self.right_child {
            let right_index = right_child.node_to_flat_internal(flat_nodes, Some(index));
            flat_nodes[index].right_child = Some(right_index);
        }

        index
    }


}
pub struct FlatTreeIter<'a> {
    tree: &'a FlatTree,
    stack: Vec<FlatTreeState>,
    order: TraversalOrder,
}

enum FlatTreeState {
    Start(usize),
    Left(usize),
    Right(usize),
    End(usize),
}

impl<'a> Iterator for FlatTreeIter<'a> {
    type Item = &'a FlatNode;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(state) = self.stack.pop() {
            match state {
                FlatTreeState::Start(index) => {
                    let node = &self.tree.nodes[index];
                    match self.order {
                        TraversalOrder::PreOrder => {
                            // PreOrder: Visit node, left, right
                            self.stack.push(FlatTreeState::Right(index));
                            self.stack.push(FlatTreeState::Left(index));
                            return Some(node);
                        }
                        TraversalOrder::InOrder => {
                            // InOrder: Left, visit node, right
                            self.stack.push(FlatTreeState::Right(index));
                            self.stack.push(FlatTreeState::End(index));
                            self.stack.push(FlatTreeState::Left(index));
                        }
                        TraversalOrder::PostOrder => {
                            // PostOrder: Left, right, visit node
                            self.stack.push(FlatTreeState::End(index));
                            self.stack.push(FlatTreeState::Right(index));
                            self.stack.push(FlatTreeState::Left(index));
                        }
                    }
                }
                FlatTreeState::Left(index) => {
                    if let Some(left_index) = self.tree.nodes[index].left_child {
                        self.stack.push(FlatTreeState::Start(left_index));
                    }
                }
                FlatTreeState::Right(index) => {
                    if let Some(right_index) = self.tree.nodes[index].right_child {
                        self.stack.push(FlatTreeState::Start(right_index));
                    }
                }
                FlatTreeState::End(index) => {
                    let node = &self.tree.nodes[index];
                    match self.order {
                        TraversalOrder::InOrder | TraversalOrder::PostOrder => {
                            return Some(node);
                        }
                        _ => {}
                    }
                }
            }
        }
        None
    }
}

impl FlatTree {
    pub fn iter(&self, order: TraversalOrder) -> FlatTreeIter {
        FlatTreeIter {
            tree: self,
            stack: vec![FlatTreeState::Start(self.root)],
            order,
        }
    }
    /// Converts a `FlatTree` into a recursive `Node` structure.
    ///
    /// # Returns
    /// A `Node` representation of the `FlatTree`.
    pub fn to_node(&self) -> Node {
        self.flat_to_node_internal(self.root)
    }

    /// Internal helper method for converting flat nodes to a `Node`.
    ///
    /// # Arguments
    /// * `index` - The index of the current node in `self.nodes`.
    ///
    /// # Returns
    /// A `Node` representing the current node.
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
    pub fn len(&self) -> usize {
        self.nodes.len()
    }
}
impl Index<usize> for FlatTree {
    type Output = FlatNode;

    fn index(&self, index: usize) -> &Self::Output {
        &self.nodes[index]
    }
}
impl IndexMut<usize> for FlatTree {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.nodes[index]
    }
}

pub trait HasName {
    fn name(&self) -> &str;
}

impl HasName for Node {
    fn name(&self) -> &str {
        &self.name
    }
}

impl HasName for FlatNode {
    fn name(&self) -> &str {
        &self.name
    }
}


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
        left_child: None,  // Will fill this in a moment
        right_child: None, // Will fill this in a moment
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
/// * `parent_index` - The index of the parent node (if any).
///
/// # Returns
/// The corresponding `Node` tree.
pub fn flat_to_node(
    flat_tree: &[FlatNode],
    index: usize,
) -> Option<Node> {
    let flat_node = &flat_tree[index];
    let left_child = flat_node
        .left_child
        .and_then(|i| flat_to_node(flat_tree, i).map(Box::new));
    let right_child = flat_node
        .right_child
        .and_then(|i| flat_to_node(flat_tree, i,).map(Box::new));

    Some(Node {
        name: flat_node.name.clone(),
        left_child,
        right_child,
        depth: flat_node.depth,
        length: flat_node.length,
    })
}


/// Recursively compares two Node objects (ignoring the order of children).
pub fn compare_nodes(n1: &Node, n2: &Node) -> bool {
    if n1.name != n2.name { return false; }
    if (n1.length - n2.length).abs() > 1e-3 { return false; }
    if n1.depth != n2.depth { return false; }

    // Collect non-None children for each node.
    let mut children1 = Vec::new();
    if let Some(child) = &n1.left_child { children1.push(child); }
    if let Some(child) = &n1.right_child { children1.push(child); }

    let mut children2 = Vec::new();
    if let Some(child) = &n2.left_child { children2.push(child); }
    if let Some(child) = &n2.right_child { children2.push(child); }

    if children1.len() != children2.len() {
        return false;
    }

    match children1.len() {
        0 => true,
        1 => compare_nodes(children1[0], children2[0]),
        2 => {
            (compare_nodes(children1[0], children2[0]) && compare_nodes(children1[1], children2[1]))
            || (compare_nodes(children1[0], children2[1]) && compare_nodes(children1[1], children2[0]))
        },
        _ => true // for trees with more than 2 children (unexpected)
    }
}