//! Tree iteration utilities for recursive and flat trees.

use super::{Node, FlatNode, FlatTree, TraversalOrder};

// ============================================================================
// NodeIter - Iterator for recursive Node trees
// ============================================================================

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

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(state) = self.stack.pop() {
            match state {
                NodeState::Start(node) => match self.order {
                    TraversalOrder::PreOrder => {
                        self.stack.push(NodeState::Right(node));
                        self.stack.push(NodeState::Left(node));
                        return Some(node);
                    }
                    TraversalOrder::InOrder => {
                        self.stack.push(NodeState::Right(node));
                        self.stack.push(NodeState::End(node));
                        self.stack.push(NodeState::Left(node));
                    }
                    TraversalOrder::PostOrder => {
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
                    TraversalOrder::InOrder | TraversalOrder::PostOrder => {
                        return Some(node);
                    }
                    TraversalOrder::PreOrder => {}
                }
            }
        }
        None
    }
}

impl Node {
    pub fn iter(&self, order: TraversalOrder) -> NodeIter<'_> {
        NodeIter {
            stack: vec![NodeState::Start(self)],
            order,
        }
    }
}

// ============================================================================
// FlatTreeIter - Iterator for FlatTree
// ============================================================================

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
                            self.stack.push(FlatTreeState::Right(index));
                            self.stack.push(FlatTreeState::Left(index));
                            return Some(node);
                        }
                        TraversalOrder::InOrder => {
                            self.stack.push(FlatTreeState::Right(index));
                            self.stack.push(FlatTreeState::End(index));
                            self.stack.push(FlatTreeState::Left(index));
                        }
                        TraversalOrder::PostOrder => {
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
    pub fn iter(&self, order: TraversalOrder) -> FlatTreeIter<'_> {
        FlatTreeIter {
            tree: self,
            stack: vec![FlatTreeState::Start(self.root)],
            order,
        }
    }
}
