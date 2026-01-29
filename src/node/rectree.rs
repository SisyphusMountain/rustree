//! Reconciled tree structure for DTL (Duplication-Transfer-Loss) model.

use super::{FlatNode, FlatTree};

/// Events that can occur during DTL reconciliation.
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

/// Reconciled tree structure for DTL model.
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
    pub node_mapping: Vec<usize>,
    /// Maps each gene tree node index to its evolutionary event
    pub event_mapping: Vec<Event>,
}

impl<'a> RecTree<'a> {
    /// Creates a new RecTree.
    ///
    /// # Panics
    /// Panics if the mappings don't match the gene tree size.
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

    /// Gets the species tree node index for a given gene tree node.
    pub fn species_node_for(&self, gene_node_idx: usize) -> usize {
        self.node_mapping[gene_node_idx]
    }

    /// Gets the event type for a given gene tree node.
    pub fn event_for(&self, gene_node_idx: usize) -> &Event {
        &self.event_mapping[gene_node_idx]
    }

    /// Gets the gene tree node and its corresponding species node and event.
    pub fn get_full_info(&self, gene_node_idx: usize) -> (&FlatNode, usize, &Event) {
        (
            &self.gene_tree.nodes[gene_node_idx],
            self.node_mapping[gene_node_idx],
            &self.event_mapping[gene_node_idx],
        )
    }

    /// Exports the reconciled tree to RecPhyloXML format with branch lengths.
    pub fn to_xml(&self) -> String {
        let estimated_size = self.gene_tree.nodes.len() * 200 + 1000;
        let mut xml = String::with_capacity(estimated_size);

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

    /// Helper function to write a species tree clade to XML.
    fn write_species_clade(&self, xml: &mut String, node_idx: usize, indent: usize) {
        let node = &self.species_tree.nodes[node_idx];
        let indent_str = Self::get_indent(indent);

        xml.push_str(indent_str);
        xml.push_str("<clade>\n");
        xml.push_str(indent_str);
        xml.push_str("\t<name>");
        xml.push_str(&node.name);
        xml.push_str("</name>\n");
        xml.push_str(indent_str);
        xml.push_str("\t<branchLength>");
        xml.push_str(&node.length.to_string());
        xml.push_str("</branchLength>\n");

        if let Some(left_idx) = node.left_child {
            self.write_species_clade(xml, left_idx, indent + 1);
        }
        if let Some(right_idx) = node.right_child {
            self.write_species_clade(xml, right_idx, indent + 1);
        }

        xml.push_str(indent_str);
        xml.push_str("</clade>\n");
    }

    /// Helper function to write a gene tree clade to XML with reconciliation events.
    fn write_gene_clade(&self, xml: &mut String, node_idx: usize, indent: usize) {
        let node = &self.gene_tree.nodes[node_idx];
        let species_idx = self.node_mapping[node_idx];
        let species_name = &self.species_tree.nodes[species_idx].name;
        let event = &self.event_mapping[node_idx];
        let indent_str = Self::get_indent(indent);

        xml.push_str(indent_str);
        xml.push_str("<clade>\n");

        // Node name
        xml.push_str(indent_str);
        xml.push_str("\t<name>");
        match event {
            Event::Loss => xml.push_str("loss"),
            Event::Leaf => xml.push_str(&node.name),
            _ => xml.push_str("NULL"),
        }
        xml.push_str("</name>\n");

        xml.push_str(indent_str);
        xml.push_str("\t<branchLength>");
        xml.push_str(&node.length.to_string());
        xml.push_str("</branchLength>\n");

        // Events section
        xml.push_str(indent_str);
        xml.push_str("\t<eventsRec>\n");

        // Check if this is a transfer recipient
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

        if is_transfer_recipient {
            xml.push_str(indent_str);
            xml.push_str("\t\t<transferBack destinationSpecies=\"");
            xml.push_str(species_name);
            xml.push_str("\"/>\n");
        }

        // Write the appropriate event tags
        match event {
            Event::Speciation => {
                xml.push_str(indent_str);
                xml.push_str("\t\t<speciation speciesLocation=\"");
                xml.push_str(species_name);
                xml.push_str("\"/>\n");
            }
            Event::Duplication => {
                xml.push_str(indent_str);
                xml.push_str("\t\t<duplication speciesLocation=\"");
                xml.push_str(species_name);
                xml.push_str("\"/>\n");
            }
            Event::Transfer => {
                if node.left_child.is_some() && node.right_child.is_some() {
                    xml.push_str(indent_str);
                    xml.push_str("\t\t<branchingOut speciesLocation=\"");
                    xml.push_str(species_name);
                    xml.push_str("\"/>\n");
                }
            }
            Event::Loss => {
                xml.push_str(indent_str);
                xml.push_str("\t\t<loss speciesLocation=\"");
                xml.push_str(species_name);
                xml.push_str("\"/>\n");
            }
            Event::Leaf => {
                xml.push_str(indent_str);
                xml.push_str("\t\t<leaf speciesLocation=\"");
                xml.push_str(species_name);
                xml.push_str("\"/>\n");
            }
        }

        xml.push_str(indent_str);
        xml.push_str("\t</eventsRec>\n");

        if let Some(left_idx) = node.left_child {
            self.write_gene_clade(xml, left_idx, indent + 1);
        }
        if let Some(right_idx) = node.right_child {
            self.write_gene_clade(xml, right_idx, indent + 1);
        }

        xml.push_str(indent_str);
        xml.push_str("</clade>\n");
    }

    /// Get indent string for XML formatting.
    fn get_indent(indent: usize) -> &'static str {
        const INDENTS: [&str; 10] = [
            "", "\t", "\t\t", "\t\t\t", "\t\t\t\t", "\t\t\t\t\t",
            "\t\t\t\t\t\t", "\t\t\t\t\t\t\t", "\t\t\t\t\t\t\t\t", "\t\t\t\t\t\t\t\t\t"
        ];
        if indent < INDENTS.len() {
            INDENTS[indent]
        } else {
            // For very deep trees, fall back to dynamic allocation
            // This is a leak but only happens for indent >= 10
            Box::leak(Box::new("\t".repeat(indent)))
        }
    }
}
