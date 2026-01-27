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

// SPR operations return Result<T, String> for error messages


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
    pub fn find_node_index(&self, name: &str) -> Option<usize> {
        self.nodes.iter().position(|n| n.name == name)
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

// SPR operations have been moved to the surgery module.
// Use `rustree::surgery::spr_topology(&mut flat_tree, moving_node_index, recipient_idx)` instead.
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

/// Events that can occur during DTL reconciliation
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Event {
    /// Speciation event - gene tree lineage follows species tree split
    Speciation,
    /// Duplication event - gene duplicates within a species
    Duplication,
    /// Horizontal gene transfer event
    Transfer,
    /// Gene loss event
    Loss,
    /// Leaf node (extant gene)
    Leaf,
}

/// Reconciled tree structure for DTL model
///
/// Represents a gene tree reconciled with a species tree, storing the mapping
/// between gene tree nodes and species tree nodes, as well as the evolutionary
/// events (Duplication, Transfer, Loss) at each node.
#[derive(Clone, Debug)]
pub struct RecTree<'a> {
    /// Reference to the species tree
    pub species_tree: &'a FlatTree,
    /// The gene tree (reconciled)
    pub gene_tree: FlatTree,
    /// Maps each gene tree node index to its corresponding species tree node index
    /// node_mapping[i] = species tree node index for gene tree node i
    pub node_mapping: Vec<usize>,
    /// Maps each gene tree node index to its evolutionary event
    /// event_mapping[i] = event type for gene tree node i
    pub event_mapping: Vec<Event>,
}

impl<'a> RecTree<'a> {
    /// Creates a new RecTree
    ///
    /// # Arguments
    /// * `species_tree` - Reference to the species tree
    /// * `gene_tree` - The gene tree
    /// * `node_mapping` - Mapping from gene tree nodes to species tree nodes
    /// * `event_mapping` - Event type for each gene tree node
    ///
    /// # Panics
    /// Panics if the mappings don't match the gene tree size
    pub fn new(
        species_tree: &'a FlatTree,
        gene_tree: FlatTree,
        node_mapping: Vec<usize>,
        event_mapping: Vec<Event>,
    ) -> Self {
        assert_eq!(
            gene_tree.nodes.len(),
            node_mapping.len(),
            "node_mapping must have same length as gene_tree nodes"
        );
        assert_eq!(
            gene_tree.nodes.len(),
            event_mapping.len(),
            "event_mapping must have same length as gene_tree nodes"
        );

        RecTree {
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
        }
    }

    /// Gets the species tree node index for a given gene tree node
    pub fn species_node_for(&self, gene_node_idx: usize) -> usize {
        self.node_mapping[gene_node_idx]
    }

    /// Gets the event type for a given gene tree node
    pub fn event_for(&self, gene_node_idx: usize) -> &Event {
        &self.event_mapping[gene_node_idx]
    }

    /// Gets the gene tree node and its corresponding species node and event
    pub fn get_full_info(&self, gene_node_idx: usize) -> (&FlatNode, usize, &Event) {
        (
            &self.gene_tree.nodes[gene_node_idx],
            self.node_mapping[gene_node_idx],
            &self.event_mapping[gene_node_idx],
        )
    }

    /// Exports the reconciled tree to RecPhyloXML format with branch lengths
    ///
    /// # Returns
    /// A string containing the XML representation of the reconciled gene tree and species tree
    pub fn to_xml(&self) -> String {
        let mut xml = String::new();

        // Header
        xml.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        xml.push_str("<recPhylo\n");
        xml.push_str("\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
        xml.push_str("\txsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\"\n");
        xml.push_str("\txmlns=\"http://www.recg.org\">\n");

        // Species tree section
        xml.push_str("<spTree>\n<phylogeny>\n");
        self.write_species_clade(&mut xml, self.species_tree.root, 0);
        xml.push_str("</phylogeny>\n</spTree>\n");

        // Gene tree section
        xml.push_str("<recGeneTree>\n<phylogeny rooted=\"true\">\n");
        self.write_gene_clade(&mut xml, self.gene_tree.root, 0);
        xml.push_str("</phylogeny>\n</recGeneTree>\n");

        xml.push_str("</recPhylo>\n");
        xml
    }

    /// Helper function to write a species tree clade to XML
    fn write_species_clade(&self, xml: &mut String, node_idx: usize, indent: usize) {
        let node = &self.species_tree.nodes[node_idx];
        let indent_str = "\t".repeat(indent);

        xml.push_str(&format!("{}<clade>\n", indent_str));
        xml.push_str(&format!("{}\t<name>{}</name>\n", indent_str, node.name));
        xml.push_str(&format!("{}\t<branchLength>{}</branchLength>\n", indent_str, node.length));

        // Recursively write children
        if let Some(left_idx) = node.left_child {
            self.write_species_clade(xml, left_idx, indent + 1);
        }
        if let Some(right_idx) = node.right_child {
            self.write_species_clade(xml, right_idx, indent + 1);
        }

        xml.push_str(&format!("{}</clade>\n", indent_str));
    }

    /// Helper function to write a gene tree clade to XML with reconciliation events
    fn write_gene_clade(&self, xml: &mut String, node_idx: usize, indent: usize) {
        let node = &self.gene_tree.nodes[node_idx];
        let species_idx = self.node_mapping[node_idx];
        let species_name = &self.species_tree.nodes[species_idx].name;
        let event = &self.event_mapping[node_idx];
        let indent_str = "\t".repeat(indent);

        xml.push_str(&format!("{}<clade>\n", indent_str));

        // Node name
        let display_name = match event {
            Event::Loss => "loss".to_string(),
            Event::Leaf => node.name.clone(),
            _ => "NULL".to_string(),
        };
        xml.push_str(&format!("{}\t<name>{}</name>\n", indent_str, display_name));
        xml.push_str(&format!("{}\t<branchLength>{}</branchLength>\n", indent_str, node.length));

        // Events section
        xml.push_str(&format!("{}\t<eventsRec>\n", indent_str));

        // Check if this is a transfer recipient (parent is Transfer and species differs from parent)
        let is_transfer_recipient = if let Some(parent_idx) = node.parent {
            if self.event_mapping[parent_idx] == Event::Transfer {
                let parent_species_idx = self.node_mapping[parent_idx];
                species_idx != parent_species_idx
            } else {
                false
            }
        } else {
            false
        };

        // If this is a transfer recipient, add transferBack first
        if is_transfer_recipient {
            xml.push_str(&format!("{}\t\t<transferBack destinationSpecies=\"{}\"/>\n",
                indent_str, species_name));
        }

        // Write the appropriate event tags
        match event {
            Event::Speciation => {
                xml.push_str(&format!("{}\t\t<speciation speciesLocation=\"{}\"/>\n",
                    indent_str, species_name));
            }
            Event::Duplication => {
                xml.push_str(&format!("{}\t\t<duplication speciesLocation=\"{}\"/>\n",
                    indent_str, species_name));
            }
            Event::Transfer => {
                // If this node has children, it's a branching out (donor parent)
                if node.left_child.is_some() && node.right_child.is_some() {
                    xml.push_str(&format!("{}\t\t<branchingOut speciesLocation=\"{}\"/>\n",
                        indent_str, species_name));
                }
                // Recipient nodes will have transferBack added above (is_transfer_recipient check)
            }
            Event::Loss => {
                xml.push_str(&format!("{}\t\t<loss speciesLocation=\"{}\"/>\n",
                    indent_str, species_name));
            }
            Event::Leaf => {
                xml.push_str(&format!("{}\t\t<leaf speciesLocation=\"{}\"/>\n",
                    indent_str, species_name));
            }
        }

        xml.push_str(&format!("{}\t</eventsRec>\n", indent_str));

        // Recursively write children
        if let Some(left_idx) = node.left_child {
            self.write_gene_clade(xml, left_idx, indent + 1);
        }
        if let Some(right_idx) = node.right_child {
            self.write_gene_clade(xml, right_idx, indent + 1);
        }

        xml.push_str(&format!("{}</clade>\n", indent_str));
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


