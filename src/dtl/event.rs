// DTL Event type definition

use crate::node::FlatTree;

/// Represents a DTL event during gene tree simulation
#[derive(Clone, Debug)]
pub enum DTLEvent {
    /// Speciation: gene follows both descendant species
    Speciation {
        time: f64,
        gene_id: usize,
        species_id: usize,
        left_child: usize,
        right_child: usize,
    },
    /// Duplication: gene duplicates within a species
    Duplication {
        time: f64,
        gene_id: usize,
        species_id: usize,
        child1: usize,
        child2: usize,
    },
    /// Transfer: gene transfers from one species to another
    Transfer {
        time: f64,
        gene_id: usize,
        from_species: usize,
        to_species: usize,
        donor_child: usize,
        recipient_child: usize,
    },
    /// Loss: gene is lost (lineage terminates)
    Loss {
        time: f64,
        gene_id: usize,
        species_id: usize,
    },
    /// Leaf: gene survives to present in extant species
    Leaf {
        time: f64,
        gene_id: usize,
        species_id: usize,
    },
}

impl DTLEvent {
    /// Convert event to CSV row format, resolving names from trees
    ///
    /// # Arguments
    /// * `species_tree` - The species tree to resolve species names
    /// * `gene_tree` - The gene tree to resolve gene node names
    pub fn to_csv_row(&self, species_tree: &FlatTree, gene_tree: &FlatTree) -> String {
        match self {
            DTLEvent::Speciation { time, gene_id, species_id, left_child, right_child } => {
                format!(
                    "{},{},Speciation,{},,,{},{}",
                    time,
                    gene_tree.nodes[*gene_id].name,
                    species_tree.nodes[*species_id].name,
                    gene_tree.nodes[*left_child].name,
                    gene_tree.nodes[*right_child].name
                )
            }
            DTLEvent::Duplication { time, gene_id, species_id, child1, child2 } => {
                format!(
                    "{},{},Duplication,{},,,{},{}",
                    time,
                    gene_tree.nodes[*gene_id].name,
                    species_tree.nodes[*species_id].name,
                    gene_tree.nodes[*child1].name,
                    gene_tree.nodes[*child2].name
                )
            }
            DTLEvent::Transfer { time, gene_id, from_species, to_species, donor_child, recipient_child } => {
                format!(
                    "{},{},Transfer,{},{},{},{},{}",
                    time,
                    gene_tree.nodes[*gene_id].name,
                    species_tree.nodes[*from_species].name,
                    species_tree.nodes[*from_species].name,
                    species_tree.nodes[*to_species].name,
                    gene_tree.nodes[*donor_child].name,
                    gene_tree.nodes[*recipient_child].name
                )
            }
            DTLEvent::Loss { time, gene_id, species_id } => {
                format!(
                    "{},{},Loss,{},,,,",
                    time,
                    gene_tree.nodes[*gene_id].name,
                    species_tree.nodes[*species_id].name
                )
            }
            DTLEvent::Leaf { time, gene_id, species_id } => {
                format!(
                    "{},{},Leaf,{},,,,",
                    time,
                    gene_tree.nodes[*gene_id].name,
                    species_tree.nodes[*species_id].name
                )
            }
        }
    }

    /// CSV header for event data
    pub fn csv_header() -> &'static str {
        "time,gene_node_name,event_type,species_node,donor_species,recipient_species,child1_name,child2_name"
    }
}
