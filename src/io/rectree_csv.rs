//! CSV export for reconciled trees (RecTree/RecTreeOwned/RecTreeColumns).

use crate::node::{RecTree, RecTreeOwned};
use crate::node::rectree::Event;
use std::io::{self, Write, BufWriter};
use std::fs::File;

// ============================================================================
// RecTreeColumns — structured column data for CSV/DataFrame export
// ============================================================================

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

// ============================================================================
// RecTree CSV methods
// ============================================================================

impl<'a> RecTree<'a> {
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
}

// ============================================================================
// RecTreeOwned CSV methods (delegates to RecTree)
// ============================================================================

impl RecTreeOwned {
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
}
