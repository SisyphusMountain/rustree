//! Reconciled tree structure for DTL (Duplication-Transfer-Loss) model.

use super::{FlatNode, FlatTree};
use std::io::{self, Write, BufWriter};
use std::fs::File;

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
    /// Reference to the gene tree (reconciled)
    pub gene_tree: &'a FlatTree,
    /// Maps each gene tree node index to its corresponding species tree node index.
    /// `None` means the mapping is unknown (e.g., after pruning/sampling).
    pub node_mapping: &'a [Option<usize>],
    /// Maps each gene tree node index to its evolutionary event
    pub event_mapping: &'a [Event],
}

impl<'a> RecTree<'a> {
    /// Creates a new RecTree.
    ///
    /// # Panics
    /// Panics if the mappings don't match the gene tree size.
    pub fn new(
        species_tree: &'a FlatTree,
        gene_tree: &'a FlatTree,
        node_mapping: &'a [Option<usize>],
        event_mapping: &'a [Event],
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
    /// Returns `None` if the mapping is unknown (e.g., after pruning).
    pub fn species_node_for(&self, gene_node_idx: usize) -> Option<usize> {
        self.node_mapping[gene_node_idx]
    }

    /// Gets the event type for a given gene tree node.
    pub fn event_for(&self, gene_node_idx: usize) -> &Event {
        &self.event_mapping[gene_node_idx]
    }

    /// Gets the gene tree node and its corresponding species node and event.
    pub fn get_full_info(&self, gene_node_idx: usize) -> (&FlatNode, Option<usize>, &Event) {
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
        let species_idx_opt = self.node_mapping[node_idx];
        let species_name = species_idx_opt.map(|idx| self.species_tree.nodes[idx].name.as_str());
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
                // Both must be Some and different for a transfer
                matches!((species_idx_opt, parent_species_idx), (Some(a), Some(b)) if a != b)
            } else {
                false
            }
        } else {
            false
        };

        if is_transfer_recipient {
            if let Some(name) = species_name {
                xml.push_str(indent_str);
                xml.push_str("\t\t<transferBack destinationSpecies=\"");
                xml.push_str(name);
                xml.push_str("\"/>\n");
            }
        }

        // Write the appropriate event tags
        match event {
            Event::Speciation => {
                xml.push_str(indent_str);
                if let Some(name) = species_name {
                    xml.push_str("\t\t<speciation speciesLocation=\"");
                    xml.push_str(name);
                    xml.push_str("\"/>\n");
                } else {
                    xml.push_str("\t\t<speciation/>\n");
                }
            }
            Event::Duplication => {
                xml.push_str(indent_str);
                if let Some(name) = species_name {
                    xml.push_str("\t\t<duplication speciesLocation=\"");
                    xml.push_str(name);
                    xml.push_str("\"/>\n");
                } else {
                    xml.push_str("\t\t<duplication/>\n");
                }
            }
            Event::Transfer => {
                if node.left_child.is_some() && node.right_child.is_some() {
                    xml.push_str(indent_str);
                    if let Some(name) = species_name {
                        xml.push_str("\t\t<branchingOut speciesLocation=\"");
                        xml.push_str(name);
                        xml.push_str("\"/>\n");
                    } else {
                        xml.push_str("\t\t<branchingOut/>\n");
                    }
                }
            }
            Event::Loss => {
                xml.push_str(indent_str);
                if let Some(name) = species_name {
                    xml.push_str("\t\t<loss speciesLocation=\"");
                    xml.push_str(name);
                    xml.push_str("\"/>\n");
                } else {
                    xml.push_str("\t\t<loss/>\n");
                }
            }
            Event::Leaf => {
                xml.push_str(indent_str);
                if let Some(name) = species_name {
                    xml.push_str("\t\t<leaf speciesLocation=\"");
                    xml.push_str(name);
                    xml.push_str("\"/>\n");
                } else {
                    xml.push_str("\t\t<leaf/>\n");
                }
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

    /// Extract structured column data for CSV export and DataFrame creation.
    pub fn to_columns(&self) -> RecTreeColumns {
        let n = self.gene_tree.nodes.len();
        let nodes = &self.gene_tree.nodes;

        let mut cols = RecTreeColumns {
            node_id: Vec::with_capacity(n),
            name: Vec::with_capacity(n),
            parent: Vec::with_capacity(n),
            left_child: Vec::with_capacity(n),
            left_child_name: Vec::with_capacity(n),
            right_child: Vec::with_capacity(n),
            right_child_name: Vec::with_capacity(n),
            length: Vec::with_capacity(n),
            depth: Vec::with_capacity(n),
            species_node: Vec::with_capacity(n),
            species_node_left: Vec::with_capacity(n),
            species_node_right: Vec::with_capacity(n),
            event: Vec::with_capacity(n),
        };

        for (i, node) in nodes.iter().enumerate() {
            cols.node_id.push(i);
            cols.name.push(node.name.clone());
            cols.parent.push(node.parent.map_or(String::new(), |p| p.to_string()));
            cols.left_child.push(node.left_child.map_or(String::new(), |c| c.to_string()));
            cols.left_child_name.push(node.left_child.map_or(String::new(), |c| nodes[c].name.clone()));
            cols.right_child.push(node.right_child.map_or(String::new(), |c| c.to_string()));
            cols.right_child_name.push(node.right_child.map_or(String::new(), |c| nodes[c].name.clone()));
            cols.length.push(node.length);
            cols.depth.push(node.depth.map_or(String::new(), |d| format!("{:.6}", d)));
            cols.species_node.push(match self.node_mapping[i] {
                Some(idx) => self.species_tree.nodes[idx].name.clone(),
                None => String::new(),
            });
            cols.species_node_left.push(node.left_child.map_or(String::new(), |c| {
                match self.node_mapping[c] {
                    Some(idx) => self.species_tree.nodes[idx].name.clone(),
                    None => String::new(),
                }
            }));
            cols.species_node_right.push(node.right_child.map_or(String::new(), |c| {
                match self.node_mapping[c] {
                    Some(idx) => self.species_tree.nodes[idx].name.clone(),
                    None => String::new(),
                }
            }));
            cols.event.push(match self.event_mapping[i] {
                Event::Speciation => "Speciation".to_string(),
                Event::Duplication => "Duplication".to_string(),
                Event::Transfer => "Transfer".to_string(),
                Event::Loss => "Loss".to_string(),
                Event::Leaf => "Leaf".to_string(),
            });
        }

        cols
    }

    /// Export reconciled tree to CSV string (header + rows).
    pub fn to_csv_string(&self) -> String {
        self.to_columns().to_csv_string()
    }

    /// Save reconciled tree to a CSV file.
    pub fn save_csv(&self, filepath: &str) -> io::Result<()> {
        self.to_columns().save_csv(filepath)
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

/// Structured column data from a reconciled tree, for CSV export and DataFrame creation.
///
/// Contains one vector per column. All vectors have the same length (number of gene tree nodes).
pub struct RecTreeColumns {
    pub node_id: Vec<usize>,
    pub name: Vec<String>,
    pub parent: Vec<String>,
    pub left_child: Vec<String>,
    pub left_child_name: Vec<String>,
    pub right_child: Vec<String>,
    pub right_child_name: Vec<String>,
    pub length: Vec<f64>,
    pub depth: Vec<String>,
    pub species_node: Vec<String>,
    pub species_node_left: Vec<String>,
    pub species_node_right: Vec<String>,
    pub event: Vec<String>,
}

impl RecTreeColumns {
    /// CSV header line.
    pub fn csv_header() -> &'static str {
        "node_id,name,parent,left_child,left_child_name,right_child,right_child_name,length,depth,species_node,species_node_left,species_node_right,event"
    }

    /// Convert to a full CSV string (header + rows).
    pub fn to_csv_string(&self) -> String {
        let n = self.node_id.len();
        let mut csv = String::with_capacity(n * 100 + 200);
        csv.push_str(Self::csv_header());
        csv.push('\n');
        for i in 0..n {
            csv.push_str(&format!(
                "{},{},{},{},{},{},{},{:.6},{},{},{},{},{}\n",
                self.node_id[i], self.name[i], self.parent[i],
                self.left_child[i], self.left_child_name[i],
                self.right_child[i], self.right_child_name[i],
                self.length[i], self.depth[i],
                self.species_node[i], self.species_node_left[i], self.species_node_right[i],
                self.event[i],
            ));
        }
        csv
    }

    /// Write to a CSV file.
    pub fn save_csv(&self, filepath: &str) -> io::Result<()> {
        let file = File::create(filepath)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{}", Self::csv_header())?;
        for i in 0..self.node_id.len() {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{:.6},{},{},{},{},{}",
                self.node_id[i], self.name[i], self.parent[i],
                self.left_child[i], self.left_child_name[i],
                self.right_child[i], self.right_child_name[i],
                self.length[i], self.depth[i],
                self.species_node[i], self.species_node_left[i], self.species_node_right[i],
                self.event[i],
            )?;
        }
        writer.flush()
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
    /// Maps each gene tree node index to its corresponding species tree node index.
    /// `None` means the mapping is unknown (e.g., after pruning/sampling).
    pub node_mapping: Vec<Option<usize>>,
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
        node_mapping: Vec<Option<usize>>,
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
            gene_tree: &self.gene_tree,
            node_mapping: &self.node_mapping,
            event_mapping: &self.event_mapping,
        }
    }

    /// Exports the reconciled tree to RecPhyloXML format with branch lengths.
    pub fn to_xml(&self) -> String {
        self.as_rectree().to_xml()
    }

    /// Extract structured column data for CSV export and DataFrame creation.
    pub fn to_columns(&self) -> RecTreeColumns {
        self.as_rectree().to_columns()
    }

    /// Export reconciled tree to CSV string (header + rows).
    pub fn to_csv_string(&self) -> String {
        self.as_rectree().to_csv_string()
    }

    /// Save reconciled tree to a CSV file.
    pub fn save_csv(&self, filepath: &str) -> io::Result<()> {
        self.as_rectree().save_csv(filepath)
    }

    /// Gets the species tree node index for a given gene tree node.
    /// Returns `None` if the mapping is unknown (e.g., after pruning).
    pub fn species_node_for(&self, gene_node_idx: usize) -> Option<usize> {
        self.node_mapping[gene_node_idx]
    }

    /// Gets the event type for a given gene tree node.
    pub fn event_for(&self, gene_node_idx: usize) -> &Event {
        &self.event_mapping[gene_node_idx]
    }

    /// Gets the gene tree node and its corresponding species node and event.
    pub fn get_full_info(&self, gene_node_idx: usize) -> (&FlatNode, Option<usize>, &Event) {
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
        };
        use std::collections::HashSet;

        // Step 1: Sample species tree (get old→new index mapping)
        let (sampled_species_tree, species_old_to_new) = extract_induced_subtree_by_names(
            &self.species_tree,
            species_leaf_names
        ).ok_or_else(|| "Failed to sample species tree".to_string())?;

        // Step 2: Find gene tree leaves that map to sampled species leaves
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
                    if let Some(species_idx) = self.node_mapping[*idx] {
                        let species_node = &self.species_tree.nodes[species_idx];
                        species_leaf_name_set.contains(species_node.name.as_str())
                    } else {
                        false // unmapped genes cannot be in sampled set
                    }
                }
            })
            .map(|(idx, _)| idx)
            .collect();

        if gene_leaves_to_keep.is_empty() {
            return Err("No gene tree leaves map to the sampled species".to_string());
        }

        // Step 3: Extract gene tree subtree (get old→new index mapping)
        let (sampled_gene_tree, gene_old_to_new) = extract_induced_subtree(&self.gene_tree, &gene_leaves_to_keep)
            .ok_or_else(|| "Failed to sample gene tree".to_string())?;

        // Step 4: Build new→old gene index mapping by inverting gene_old_to_new
        let mut gene_new_to_old: Vec<Option<usize>> = vec![None; sampled_gene_tree.nodes.len()];
        for (old_idx, new_idx_opt) in gene_old_to_new.iter().enumerate() {
            if let Some(new_idx) = new_idx_opt {
                gene_new_to_old[*new_idx] = Some(old_idx);
            }
        }

        // Step 5: Rebuild mappings for sampled gene tree using index-based lookup
        let mut new_node_mapping: Vec<Option<usize>> = vec![None; sampled_gene_tree.nodes.len()];
        let mut new_event_mapping = vec![Event::Leaf; sampled_gene_tree.nodes.len()];

        for new_gene_idx in 0..sampled_gene_tree.nodes.len() {
            let old_gene_idx = gene_new_to_old[new_gene_idx]
                .ok_or_else(|| format!("New gene node {} has no original mapping", new_gene_idx))?;

            // Propagate event from original tree
            new_event_mapping[new_gene_idx] = self.event_mapping[old_gene_idx].clone();

            // Map species: old gene → old species → new species (via species_old_to_new)
            if let Some(old_species_idx) = self.node_mapping[old_gene_idx] {
                if let Some(new_species_idx) = species_old_to_new[old_species_idx] {
                    new_node_mapping[new_gene_idx] = Some(new_species_idx);
                }
                // else: species node was collapsed/discarded in sampling → None
            }
            // else: gene was unmapped → stays None
        }

        Ok(RecTreeOwned::new(
            sampled_species_tree,
            sampled_gene_tree,
            new_node_mapping,
            new_event_mapping,
        ))
    }
}
