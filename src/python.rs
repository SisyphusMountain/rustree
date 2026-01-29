//! Python bindings for rustree using PyO3.
//!
//! Provides Python access to:
//! - Birth-death species tree simulation
//! - DTL gene tree simulation
//! - Gene tree sampling (induced subtree extraction)

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::bd::simulate_bd_tree;
use crate::dtl::{simulate_dtl, simulate_dtl_batch};
use crate::node::{FlatTree, Event};
use crate::sampling::{extract_induced_subtree, find_leaf_indices_by_names};
use std::fs;
use std::process::Command;

/// A species tree simulated under the birth-death process.
#[pyclass]
#[derive(Clone)]
pub struct PySpeciesTree {
    tree: FlatTree,
}

#[pymethods]
impl PySpeciesTree {
    /// Convert the species tree to Newick format.
    fn to_newick(&self) -> String {
        self.tree.to_newick() + ";"
    }

    /// Get the number of nodes in the tree.
    fn num_nodes(&self) -> usize {
        self.tree.nodes.len()
    }

    /// Get the number of extant species (leaves).
    fn num_leaves(&self) -> usize {
        self.tree.nodes.iter()
            .filter(|n| n.left_child.is_none() && n.right_child.is_none())
            .count()
    }

    /// Get the total tree height (depth of deepest leaf).
    fn tree_height(&self) -> f64 {
        self.tree.nodes.iter()
            .filter_map(|n| n.depth)
            .fold(0.0, f64::max)
    }

    /// Get the root node index.
    fn root_index(&self) -> usize {
        self.tree.root
    }

    /// Get leaf names.
    fn leaf_names(&self) -> Vec<String> {
        self.tree.nodes.iter()
            .filter(|n| n.left_child.is_none() && n.right_child.is_none())
            .map(|n| n.name.clone())
            .collect()
    }

    /// Save the species tree to a Newick file.
    ///
    /// # Arguments
    /// * `filepath` - Path to save the Newick file
    fn save_newick(&self, filepath: &str) -> PyResult<()> {
        let newick = self.to_newick();
        fs::write(filepath, newick)
            .map_err(|e| PyValueError::new_err(format!("Failed to write Newick file: {}", e)))?;
        Ok(())
    }

    /// Simulate a gene tree along this species tree using the DTL model.
    ///
    /// # Arguments
    /// * `lambda_d` - Duplication rate per unit time
    /// * `lambda_t` - Transfer rate per unit time
    /// * `lambda_l` - Loss rate per unit time
    /// * `transfer_alpha` - Distance decay for assortative transfers (optional, None = uniform)
    /// * `seed` - Random seed for reproducibility (optional)
    ///
    /// # Returns
    /// A PyGeneTree containing the simulated gene tree with its mapping to the species tree.
    #[pyo3(signature = (lambda_d, lambda_t, lambda_l, transfer_alpha=None, seed=None))]
    fn simulate_dtl(&self, lambda_d: f64, lambda_t: f64, lambda_l: f64, transfer_alpha: Option<f64>, seed: Option<u64>) -> PyResult<PyGeneTree> {
        if lambda_d < 0.0 || lambda_t < 0.0 || lambda_l < 0.0 {
            return Err(PyValueError::new_err("Rates must be non-negative"));
        }

        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let origin_species = self.tree.root;
        let (rec_tree, _events) = simulate_dtl(&self.tree, origin_species, lambda_d, lambda_t, lambda_l, transfer_alpha, &mut rng);

        // Extract owned data from RecTree
        let gene_tree = rec_tree.gene_tree;
        let node_mapping = rec_tree.node_mapping;
        let event_mapping = rec_tree.event_mapping;

        Ok(PyGeneTree {
            gene_tree,
            species_tree: self.tree.clone(),
            node_mapping,
            event_mapping,
        })
    }

    /// Simulate multiple gene trees efficiently with shared pre-computed data.
    ///
    /// This is faster than calling simulate_dtl multiple times because depths
    /// and contemporaneity are computed only once.
    ///
    /// # Arguments
    /// * `n` - Number of gene trees to simulate
    /// * `lambda_d` - Duplication rate per unit time
    /// * `lambda_t` - Transfer rate per unit time
    /// * `lambda_l` - Loss rate per unit time
    /// * `transfer_alpha` - Distance decay for assortative transfers (optional, None = uniform)
    /// * `seed` - Random seed for reproducibility (optional)
    ///
    /// # Returns
    /// A list of PyGeneTree objects.
    #[pyo3(signature = (n, lambda_d, lambda_t, lambda_l, transfer_alpha=None, seed=None))]
    fn simulate_dtl_batch(&self, n: usize, lambda_d: f64, lambda_t: f64, lambda_l: f64, transfer_alpha: Option<f64>, seed: Option<u64>) -> PyResult<Vec<PyGeneTree>> {
        if lambda_d < 0.0 || lambda_t < 0.0 || lambda_l < 0.0 {
            return Err(PyValueError::new_err("Rates must be non-negative"));
        }

        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let origin_species = self.tree.root;
        let (rec_trees, _all_events) = simulate_dtl_batch(
            &self.tree,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            n,
            &mut rng,
        );

        let gene_trees: Vec<PyGeneTree> = rec_trees
            .into_iter()
            .map(|rec_tree| PyGeneTree {
                gene_tree: rec_tree.gene_tree,
                species_tree: self.tree.clone(),
                node_mapping: rec_tree.node_mapping,
                event_mapping: rec_tree.event_mapping,
            })
            .collect();

        Ok(gene_trees)
    }
}

/// A gene tree simulated under the DTL model.
#[pyclass]
#[derive(Clone)]
pub struct PyGeneTree {
    gene_tree: FlatTree,
    species_tree: FlatTree,
    node_mapping: Vec<usize>,
    event_mapping: Vec<Event>,
}

#[pymethods]
impl PyGeneTree {
    /// Convert the gene tree to Newick format.
    fn to_newick(&self) -> String {
        self.gene_tree.to_newick() + ";"
    }

    /// Save the gene tree to a Newick file.
    ///
    /// # Arguments
    /// * `filepath` - Path to save the Newick file
    fn save_newick(&self, filepath: &str) -> PyResult<()> {
        let newick = self.to_newick();
        fs::write(filepath, newick)
            .map_err(|e| PyValueError::new_err(format!("Failed to write Newick file: {}", e)))?;
        Ok(())
    }

    /// Get the number of nodes in the gene tree.
    fn num_nodes(&self) -> usize {
        self.gene_tree.nodes.len()
    }

    /// Get the number of extant genes (leaves that survived, not losses).
    fn num_extant(&self) -> usize {
        self.gene_tree.nodes.iter()
            .enumerate()
            .filter(|(i, n)| {
                n.left_child.is_none()
                && n.right_child.is_none()
                && self.event_mapping[*i] == Event::Leaf
            })
            .count()
    }

    /// Get the number of events by type: (speciations, duplications, transfers, losses, leaves).
    fn count_events(&self) -> (usize, usize, usize, usize, usize) {
        let mut speciations = 0;
        let mut duplications = 0;
        let mut transfers = 0;
        let mut losses = 0;
        let mut leaves = 0;

        for event in &self.event_mapping {
            match event {
                Event::Speciation => speciations += 1,
                Event::Duplication => duplications += 1,
                Event::Transfer => transfers += 1,
                Event::Loss => losses += 1,
                Event::Leaf => leaves += 1,
            }
        }

        (speciations, duplications, transfers, losses, leaves)
    }

    /// Get names of extant genes (genes that survived to present).
    fn extant_gene_names(&self) -> Vec<String> {
        self.gene_tree.nodes.iter()
            .enumerate()
            .filter(|(i, n)| {
                n.left_child.is_none()
                && n.right_child.is_none()
                && self.event_mapping[*i] == Event::Leaf
            })
            .map(|(_, n)| n.name.clone())
            .collect()
    }

    /// Sample the gene tree by keeping only genes from extant species.
    ///
    /// Returns a new gene tree containing only the induced subtree of extant genes.
    /// If all genes are lost, returns None (represented as an error in Python).
    fn sample_extant(&self) -> PyResult<PyGeneTree> {
        // Find indices of extant genes (leaves with Event::Leaf)
        let extant_indices: std::collections::HashSet<usize> = self.gene_tree.nodes.iter()
            .enumerate()
            .filter(|(i, n)| {
                n.left_child.is_none()
                && n.right_child.is_none()
                && self.event_mapping[*i] == Event::Leaf
            })
            .map(|(i, _)| i)
            .collect();

        if extant_indices.is_empty() {
            return Err(PyValueError::new_err("No extant genes to sample"));
        }

        // Extract induced subtree
        let sampled_tree = extract_induced_subtree(&self.gene_tree, &extant_indices)
            .ok_or_else(|| PyValueError::new_err("Failed to extract induced subtree"))?;

        // For the sampled tree, we need to rebuild mappings
        // The sampled tree is just the topology - we lose the detailed event mapping
        // since internal nodes may be collapsed
        let num_nodes = sampled_tree.nodes.len();
        let new_event_mapping: Vec<Event> = sampled_tree.nodes.iter()
            .map(|n| {
                if n.left_child.is_none() && n.right_child.is_none() {
                    Event::Leaf
                } else {
                    Event::Speciation // Internal nodes are treated as speciations in sampled tree
                }
            })
            .collect();

        // For node mapping, we can't preserve it accurately after collapsing
        // Just map all nodes to root species as a placeholder
        let new_node_mapping = vec![self.species_tree.root; num_nodes];

        Ok(PyGeneTree {
            gene_tree: sampled_tree,
            species_tree: self.species_tree.clone(),
            node_mapping: new_node_mapping,
            event_mapping: new_event_mapping,
        })
    }

    /// Sample the gene tree by keeping only genes with the specified names.
    ///
    /// # Arguments
    /// * `names` - List of gene names to keep
    ///
    /// # Returns
    /// A new gene tree containing only the induced subtree of the specified genes.
    fn sample_by_names(&self, names: Vec<String>) -> PyResult<PyGeneTree> {
        let keep_indices = find_leaf_indices_by_names(&self.gene_tree, &names);

        if keep_indices.is_empty() {
            return Err(PyValueError::new_err("No matching genes found"));
        }

        let sampled_tree = extract_induced_subtree(&self.gene_tree, &keep_indices)
            .ok_or_else(|| PyValueError::new_err("Failed to extract induced subtree"))?;

        let num_nodes = sampled_tree.nodes.len();
        let new_event_mapping: Vec<Event> = sampled_tree.nodes.iter()
            .map(|n| {
                if n.left_child.is_none() && n.right_child.is_none() {
                    Event::Leaf
                } else {
                    Event::Speciation
                }
            })
            .collect();

        let new_node_mapping = vec![self.species_tree.root; num_nodes];

        Ok(PyGeneTree {
            gene_tree: sampled_tree,
            species_tree: self.species_tree.clone(),
            node_mapping: new_node_mapping,
            event_mapping: new_event_mapping,
        })
    }

    /// Export the reconciled tree to RecPhyloXML format as a string.
    fn to_xml(&self) -> String {
        // Rebuild RecTree temporarily for XML export
        use crate::node::RecTree;
        let rec_tree = RecTree::new(
            &self.species_tree,
            self.gene_tree.clone(),
            self.node_mapping.clone(),
            self.event_mapping.clone(),
        );
        rec_tree.to_xml()
    }

    /// Save the reconciled tree to a RecPhyloXML file.
    ///
    /// # Arguments
    /// * `filepath` - Path to save the XML file
    fn save_xml(&self, filepath: &str) -> PyResult<()> {
        let xml = self.to_xml();
        fs::write(filepath, xml)
            .map_err(|e| PyValueError::new_err(format!("Failed to write XML file: {}", e)))?;
        Ok(())
    }

    /// Generate an SVG visualization of the reconciled tree using thirdkind.
    ///
    /// Requires thirdkind to be installed (`cargo install thirdkind`).
    ///
    /// # Arguments
    /// * `filepath` - Optional path to save the SVG file. If None, returns SVG as string.
    /// * `open_browser` - If True, opens the SVG in a web browser (default: False)
    ///
    /// # Returns
    /// The SVG content as a string.
    #[pyo3(signature = (filepath=None, open_browser=false))]
    fn to_svg(&self, filepath: Option<&str>, open_browser: bool) -> PyResult<String> {
        // Create temp file for XML input
        let temp_dir = std::env::temp_dir();
        let xml_path = temp_dir.join("rustree_temp.recphyloxml");
        let svg_path = temp_dir.join("rustree_temp.svg");

        // Write XML to temp file
        let xml = self.to_xml();
        fs::write(&xml_path, &xml)
            .map_err(|e| PyValueError::new_err(format!("Failed to write temp XML: {}", e)))?;

        // Call thirdkind
        let mut cmd = Command::new("thirdkind");
        cmd.arg("-f").arg(&xml_path)
           .arg("-o").arg(&svg_path);

        if open_browser {
            cmd.arg("-b");
        }

        let output = cmd.output()
            .map_err(|e| PyValueError::new_err(format!(
                "Failed to run thirdkind. Is it installed? (`cargo install thirdkind`)\nError: {}", e
            )))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(PyValueError::new_err(format!("thirdkind failed: {}", stderr)));
        }

        // Read SVG output
        let svg = fs::read_to_string(&svg_path)
            .map_err(|e| PyValueError::new_err(format!("Failed to read SVG output: {}", e)))?;

        // Save to user-specified path if provided
        if let Some(path) = filepath {
            fs::write(path, &svg)
                .map_err(|e| PyValueError::new_err(format!("Failed to write SVG file: {}", e)))?;
        }

        // Cleanup temp files
        let _ = fs::remove_file(&xml_path);
        let _ = fs::remove_file(&svg_path);

        Ok(svg)
    }

    /// Display the reconciled tree visualization in a Jupyter notebook.
    ///
    /// Requires thirdkind to be installed and IPython/Jupyter environment.
    fn display(&self, py: Python) -> PyResult<PyObject> {
        let svg = self.to_svg(None, false)?;

        let ipython_display = py.import("IPython.display")?;
        let svg_class = ipython_display.getattr("SVG")?;
        let display_obj = svg_class.call1((svg,))?;

        Ok(display_obj.into())
    }

    /// Export gene tree data as a pandas DataFrame.
    ///
    /// Returns a DataFrame with columns: node_id, name, parent, left_child, right_child,
    /// length, depth, species_node, event
    ///
    /// # Arguments
    /// * `filepath` - Optional path to save the CSV file
    ///
    /// # Returns
    /// A pandas DataFrame with the gene tree data.
    #[pyo3(signature = (filepath=None))]
    fn to_csv(&self, py: Python, filepath: Option<&str>) -> PyResult<PyObject> {
        // Build CSV data
        let mut rows: Vec<(usize, String, String, String, String, f64, String, String, String)> = Vec::new();

        for (i, node) in self.gene_tree.nodes.iter().enumerate() {
            let parent = node.parent.map_or(String::new(), |p| p.to_string());
            let left = node.left_child.map_or(String::new(), |c| c.to_string());
            let right = node.right_child.map_or(String::new(), |c| c.to_string());
            let depth = node.depth.map_or(String::new(), |d| format!("{:.6}", d));
            let species_node = self.species_tree.nodes[self.node_mapping[i]].name.clone();
            let event = match self.event_mapping[i] {
                Event::Speciation => "Speciation",
                Event::Duplication => "Duplication",
                Event::Transfer => "Transfer",
                Event::Loss => "Loss",
                Event::Leaf => "Leaf",
            };

            rows.push((
                i,
                node.name.clone(),
                parent,
                left,
                right,
                node.length,
                depth,
                species_node,
                event.to_string(),
            ));
        }

        // Save to CSV if filepath provided
        if let Some(path) = filepath {
            let mut csv_content = String::from("node_id,name,parent,left_child,right_child,length,depth,species_node,event\n");
            for row in &rows {
                csv_content.push_str(&format!(
                    "{},{},{},{},{},{:.6},{},{},{}\n",
                    row.0, row.1, row.2, row.3, row.4, row.5, row.6, row.7, row.8
                ));
            }
            fs::write(path, csv_content)
                .map_err(|e| PyValueError::new_err(format!("Failed to write CSV file: {}", e)))?;
        }

        // Create pandas DataFrame
        let pandas = py.import("pandas")?;
        let dict = pyo3::types::PyDict::new(py);

        let node_ids: Vec<usize> = rows.iter().map(|r| r.0).collect();
        let names: Vec<&str> = rows.iter().map(|r| r.1.as_str()).collect();
        let parents: Vec<&str> = rows.iter().map(|r| r.2.as_str()).collect();
        let lefts: Vec<&str> = rows.iter().map(|r| r.3.as_str()).collect();
        let rights: Vec<&str> = rows.iter().map(|r| r.4.as_str()).collect();
        let lengths: Vec<f64> = rows.iter().map(|r| r.5).collect();
        let depths: Vec<&str> = rows.iter().map(|r| r.6.as_str()).collect();
        let species_nodes: Vec<&str> = rows.iter().map(|r| r.7.as_str()).collect();
        let events: Vec<&str> = rows.iter().map(|r| r.8.as_str()).collect();

        dict.set_item("node_id", node_ids)?;
        dict.set_item("name", names)?;
        dict.set_item("parent", parents)?;
        dict.set_item("left_child", lefts)?;
        dict.set_item("right_child", rights)?;
        dict.set_item("length", lengths)?;
        dict.set_item("depth", depths)?;
        dict.set_item("species_node", species_nodes)?;
        dict.set_item("event", events)?;

        let df = pandas.call_method1("DataFrame", (dict,))?;
        Ok(df.into())
    }
}

/// Simulate a birth-death species tree.
///
/// # Arguments
/// * `n` - Number of extant species (must be > 0)
/// * `lambda_` - Speciation/birth rate (must be > 0)
/// * `mu` - Extinction/death rate (must be >= 0 and < lambda)
/// * `seed` - Random seed for reproducibility (optional)
///
/// # Returns
/// A PySpeciesTree with n extant species simulated under the birth-death process.
#[pyfunction]
#[pyo3(signature = (n, lambda_, mu, seed=None))]
fn simulate_species_tree(n: usize, lambda_: f64, mu: f64, seed: Option<u64>) -> PyResult<PySpeciesTree> {
    if n == 0 {
        return Err(PyValueError::new_err("Number of species must be positive"));
    }
    if lambda_ <= 0.0 {
        return Err(PyValueError::new_err("Speciation rate must be positive"));
    }
    if mu < 0.0 {
        return Err(PyValueError::new_err("Extinction rate must be non-negative"));
    }
    if lambda_ <= mu {
        return Err(PyValueError::new_err("Speciation rate must be strictly greater than extinction rate"));
    }

    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let (mut tree, _events) = simulate_bd_tree(n, lambda_, mu, &mut rng);
    tree.assign_depths();

    Ok(PySpeciesTree { tree })
}

/// Parse a Newick string into a species tree.
///
/// # Arguments
/// * `newick_str` - A Newick formatted string representing the tree
///
/// # Returns
/// A PySpeciesTree parsed from the Newick string.
#[pyfunction]
fn parse_species_tree(newick_str: &str) -> PyResult<PySpeciesTree> {
    use crate::newick::newick::parse_newick;

    let mut nodes = parse_newick(newick_str)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse Newick: {}", e)))?;

    let mut root = nodes.pop()
        .ok_or_else(|| PyValueError::new_err("No tree found in Newick string"))?;

    root.zero_root_length();
    root.assign_depths(0.0);
    let tree = root.to_flat_tree();

    Ok(PySpeciesTree { tree })
}

/// Python module for rustree.
#[pymodule]
fn rustree(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(simulate_species_tree, m)?)?;
    m.add_function(wrap_pyfunction!(parse_species_tree, m)?)?;
    m.add_class::<PySpeciesTree>()?;
    m.add_class::<PyGeneTree>()?;
    Ok(())
}
