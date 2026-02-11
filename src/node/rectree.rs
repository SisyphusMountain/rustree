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
