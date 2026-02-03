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

/// Owned version of RecTree that owns both the species tree and gene tree.
///
/// This structure is useful when parsing RecPhyloXML files where both trees
/// need to be owned by the same structure. It can be converted to a borrowed
/// RecTree using the `as_rectree()` method.
#[derive(Clone, Debug)]
pub struct RecTreeOwned {
    /// The species tree (owned)
    pub species_tree: FlatTree,
    /// The gene tree (owned)
    pub gene_tree: FlatTree,
    /// Maps each gene tree node index to its corresponding species tree node index
    pub node_mapping: Vec<usize>,
    /// Maps each gene tree node index to its evolutionary event
    pub event_mapping: Vec<Event>,
}

impl RecTreeOwned {
    /// Creates a new RecTreeOwned.
    ///
    /// # Panics
    /// Panics if the mappings don't match the gene tree size.
    pub fn new(
        species_tree: FlatTree,
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

        RecTreeOwned {
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
        }
    }

    /// Parse a RecPhyloXML string and create a RecTreeOwned.
    ///
    /// # Arguments
    /// * `xml_content` - The XML content as a string
    ///
    /// # Returns
    /// A Result containing the RecTreeOwned or an error message
    pub fn from_xml(xml_content: &str) -> Result<Self, String> {
        use super::recphyloxml::parse_recphyloxml;

        let (species_tree, gene_tree, node_mapping, event_mapping) =
            parse_recphyloxml(xml_content).map_err(|e| e.to_string())?;

        Ok(RecTreeOwned::new(
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
        ))
    }

    /// Parse a RecPhyloXML file and create a RecTreeOwned.
    ///
    /// # Arguments
    /// * `filepath` - Path to the RecPhyloXML file
    ///
    /// # Returns
    /// A Result containing the RecTreeOwned or an error message
    pub fn from_xml_file(filepath: &str) -> Result<Self, String> {
        use super::recphyloxml::parse_recphyloxml_file;

        let (species_tree, gene_tree, node_mapping, event_mapping) =
            parse_recphyloxml_file(filepath).map_err(|e| e.to_string())?;

        Ok(RecTreeOwned::new(
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
        ))
    }

    /// Convert to a borrowed RecTree.
    ///
    /// This allows using RecTreeOwned with functions that expect a RecTree.
    pub fn as_rectree(&self) -> RecTree<'_> {
        RecTree {
            species_tree: &self.species_tree,
            gene_tree: self.gene_tree.clone(),
            node_mapping: self.node_mapping.clone(),
            event_mapping: self.event_mapping.clone(),
        }
    }

    /// Exports the reconciled tree to RecPhyloXML format with branch lengths.
    pub fn to_xml(&self) -> String {
        self.as_rectree().to_xml()
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

    /// Parse gene-tree-only RecPhyloXML with a separate species tree.
    ///
    /// This is useful when the species tree is provided separately (e.g., from a Newick file)
    /// and the XML only contains the reconciled gene tree (no `<spTree>` section).
    ///
    /// # Arguments
    /// * `xml_content` - XML string containing only `<recGeneTree>` section
    /// * `species_tree` - Pre-parsed species tree as FlatTree
    ///
    /// # Returns
    /// A Result containing the RecTreeOwned or an error message
    pub fn from_gene_tree_xml(xml_content: &str, species_tree: FlatTree) -> Result<Self, String> {
        use super::recphyloxml::parse_gene_tree_only;

        let (gene_tree, node_mapping, event_mapping) =
            parse_gene_tree_only(xml_content, &species_tree).map_err(|e| e.to_string())?;

        Ok(RecTreeOwned::new(
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
        ))
    }

    /// Parse gene-tree-only RecPhyloXML file with a separate species tree.
    ///
    /// File-based version of `from_gene_tree_xml`.
    ///
    /// # Arguments
    /// * `xml_filepath` - Path to XML file containing only `<recGeneTree>` section
    /// * `species_tree` - Pre-parsed species tree as FlatTree
    ///
    /// # Returns
    /// A Result containing the RecTreeOwned or an error message
    pub fn from_gene_tree_xml_file(
        xml_filepath: &str,
        species_tree: FlatTree,
    ) -> Result<Self, String> {
        use super::recphyloxml::parse_gene_tree_only_file;

        let (gene_tree, node_mapping, event_mapping) =
            parse_gene_tree_only_file(xml_filepath, &species_tree).map_err(|e| e.to_string())?;

        Ok(RecTreeOwned::new(
            species_tree,
            gene_tree,
            node_mapping,
            event_mapping,
        ))
    }

    /// Parse reconciled tree from separate files: Newick species tree + gene tree XML.
    ///
    /// Convenience function that combines parsing a species tree from Newick format
    /// and a gene tree from RecPhyloXML format (without `<spTree>` section).
    ///
    /// # Arguments
    /// * `species_newick_path` - Path to Newick file containing species tree
    /// * `gene_xml_path` - Path to XML file containing only `<recGeneTree>` section
    ///
    /// # Returns
    /// A Result containing the RecTreeOwned or an error message
    ///
    /// # Example
    /// ```no_run
    /// use rustree::RecTreeOwned;
    ///
    /// let rec_tree = RecTreeOwned::from_separate_files(
    ///     "species_tree.nwk",
    ///     "gene_tree_rec.xml"
    /// )?;
    /// # Ok::<(), String>(())
    /// ```
    pub fn from_separate_files(
        species_newick_path: &str,
        gene_xml_path: &str,
    ) -> Result<Self, String> {
        use crate::newick::newick::parse_newick;
        use std::fs;

        // Parse species tree from Newick file
        let species_newick = fs::read_to_string(species_newick_path)
            .map_err(|e| format!("Failed to read species tree file: {}", e))?;

        let mut species_nodes = parse_newick(&species_newick)
            .map_err(|e| format!("Failed to parse species tree Newick: {}", e))?;

        let species_root = species_nodes
            .pop()
            .ok_or_else(|| "No tree found in species Newick file".to_string())?;

        let mut species_tree = species_root.to_flat_tree();
        species_tree.assign_depths();

        // Parse gene tree from XML file
        Self::from_gene_tree_xml_file(gene_xml_path, species_tree)
    }

    /// Sample species tree leaves and filter gene tree accordingly.
    ///
    /// This method samples a subset of species from the species tree and automatically
    /// filters the gene tree to keep only genes that map to the sampled species. The
    /// reconciliation mappings are preserved using an LCA-based approach.
    ///
    /// # Arguments
    /// * `species_leaf_names` - Names of species leaves to keep
    ///
    /// # Returns
    /// A new RecTreeOwned with sampled trees and updated mappings
    ///
    /// # Algorithm
    /// 1. Sample the species tree using the provided leaf names
    /// 2. Build LCA mappings for both original and sampled species trees
    /// 3. Create a mapping from sampled species tree indices to original tree indices
    /// 4. Filter gene tree leaves that map to sampled species
    /// 5. Extract the induced gene tree subtree
    /// 6. Rebuild reconciliation mappings for the sampled gene tree
    ///
    /// # Example
    /// ```no_run
    /// use rustree::RecTreeOwned;
    ///
    /// let rec_tree = RecTreeOwned::from_xml_file("reconciliation.xml")?;
    ///
    /// // Sample only species A, B, and C
    /// let sampled = rec_tree.sample_species_leaves(&[
    ///     "species_A".to_string(),
    ///     "species_B".to_string(),
    ///     "species_C".to_string(),
    /// ])?;
    /// # Ok::<(), String>(())
    /// ```
    pub fn sample_species_leaves(&self, species_leaf_names: &[String]) -> Result<Self, String> {
        use crate::sampling::{
            extract_induced_subtree,
            extract_induced_subtree_by_names,
            build_leaf_pair_lca_map,
            build_sampled_to_original_mapping,
        };
        use std::collections::{HashSet, HashMap};

        // Step 1: Sample species tree
        let sampled_species_tree = extract_induced_subtree_by_names(
            &self.species_tree,
            species_leaf_names
        ).ok_or_else(|| "Failed to sample species tree".to_string())?;

        // Step 2: Build LCA mappings
        let original_species_lca_map = build_leaf_pair_lca_map(&self.species_tree);
        let sampled_species_lca_map = build_leaf_pair_lca_map(&sampled_species_tree);

        // Step 3: Create sampled→original index mapping
        let _species_idx_mapping = build_sampled_to_original_mapping(
            &sampled_species_tree,
            &self.species_tree,
            &sampled_species_lca_map,
            &original_species_lca_map,
        );

        // Step 4: Find gene tree leaves that map to sampled species leaves
        let species_leaf_name_set: HashSet<&str> = species_leaf_names
            .iter()
            .map(|s| s.as_str())
            .collect();

        let gene_leaves_to_keep: HashSet<usize> = self.gene_tree.nodes
            .iter()
            .enumerate()
            .filter(|(idx, node)| {
                // Check if it's a leaf
                node.left_child.is_none() && node.right_child.is_none() &&
                // Check if its species is in sampled set
                {
                    let species_idx = self.node_mapping[*idx];
                    let species_node = &self.species_tree.nodes[species_idx];
                    species_leaf_name_set.contains(species_node.name.as_str())
                }
            })
            .map(|(idx, _)| idx)
            .collect();

        if gene_leaves_to_keep.is_empty() {
            return Err("No gene tree leaves map to the sampled species".to_string());
        }

        // Step 5: Extract gene tree subtree
        let sampled_gene_tree = extract_induced_subtree(&self.gene_tree, &gene_leaves_to_keep)
            .ok_or_else(|| "Failed to sample gene tree".to_string())?;

        // Step 6: Build mapping from old gene indices to new gene indices
        let mut old_to_new_gene_idx: HashMap<usize, usize> = HashMap::new();
        for new_idx in 0..sampled_gene_tree.nodes.len() {
            let new_node = &sampled_gene_tree.nodes[new_idx];
            // Find corresponding node in original tree by name
            if let Some(old_idx) = self.gene_tree.nodes.iter()
                .position(|n| n.name == new_node.name) {
                old_to_new_gene_idx.insert(old_idx, new_idx);
            }
        }

        // Step 7: Rebuild mappings for sampled gene tree
        let mut new_node_mapping = vec![0; sampled_gene_tree.nodes.len()];
        let mut new_event_mapping = vec![Event::Leaf; sampled_gene_tree.nodes.len()];

        for new_gene_idx in 0..sampled_gene_tree.nodes.len() {
            let gene_node = &sampled_gene_tree.nodes[new_gene_idx];

            // Find old gene index by name
            let old_gene_idx = self.gene_tree.nodes.iter()
                .position(|n| n.name == gene_node.name)
                .ok_or_else(|| format!("Gene node {} not found in original tree", gene_node.name))?;

            // Get old species index
            let old_species_idx = self.node_mapping[old_gene_idx];

            // Map to new species index using sampled→original mapping
            // Find the species node by name in the sampled tree
            let species_name = &self.species_tree.nodes[old_species_idx].name;
            let new_species_idx = sampled_species_tree.nodes.iter()
                .position(|n| n.name == *species_name)
                .ok_or_else(|| format!("Species {} not found in sampled tree", species_name))?;

            new_node_mapping[new_gene_idx] = new_species_idx;
            new_event_mapping[new_gene_idx] = self.event_mapping[old_gene_idx].clone();
        }

        Ok(RecTreeOwned::new(
            sampled_species_tree,
            sampled_gene_tree,
            new_node_mapping,
            new_event_mapping,
        ))
    }
}
