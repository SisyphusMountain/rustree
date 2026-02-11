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

use crate::bd::simulate_bd_tree_bwd;
use crate::dtl::{simulate_dtl, simulate_dtl_batch, simulate_dtl_per_species, simulate_dtl_per_species_batch};
use crate::node::{FlatTree, Event, RecTree};
use crate::sampling::{extract_induced_subtree, extract_induced_subtree_by_names, find_leaf_indices_by_names};
use std::collections::HashMap;
use std::fs;
use std::process::Command;
use std::sync::Arc;

// ============================================================================
// Helper Functions (Refactored to reduce duplication)
// ============================================================================

/// Validate DTL rates (all must be non-negative)
fn validate_dtl_rates(lambda_d: f64, lambda_t: f64, lambda_l: f64) -> PyResult<()> {
    if lambda_d < 0.0 || lambda_t < 0.0 || lambda_l < 0.0 {
        return Err(PyValueError::new_err("Rates must be non-negative"));
    }
    Ok(())
}

/// Validate replacement_transfer probability parameter
fn validate_replacement_transfer(replacement_transfer: Option<f64>) -> PyResult<()> {
    if let Some(p) = replacement_transfer {
        if !(0.0..=1.0).contains(&p) {
            return Err(PyValueError::new_err(
                "replacement_transfer must be between 0.0 and 1.0"
            ));
        }
    }
    Ok(())
}

/// Initialize RNG from optional seed
fn init_rng(seed: Option<u64>) -> StdRng {
    match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    }
}

/// Check if a node is a leaf (no children)
fn is_leaf(node: &crate::node::FlatNode) -> bool {
    node.left_child.is_none() && node.right_child.is_none()
}

// ============================================================================

/// A species tree simulated under the birth-death process.
#[pyclass]
#[derive(Clone)]
pub struct PySpeciesTree {
    tree: Arc<FlatTree>,
}

#[pymethods]
impl PySpeciesTree {
    /// Convert the species tree to Newick format.
    fn to_newick(&self) -> PyResult<String> {
        let nwk = self.tree.to_newick()
            .map_err(|e| PyValueError::new_err(e))?;
        Ok(nwk + ";")
    }

    /// Get the number of nodes in the tree.
    fn num_nodes(&self) -> usize {
        self.tree.nodes.len()
    }

    /// Get the number of extant species (leaves).
    fn num_leaves(&self) -> usize {
        self.tree.nodes.iter().filter(|n| is_leaf(n)).count()
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
            .filter(|n| is_leaf(n))
            .map(|n| n.name.clone())
            .collect()
    }

    /// Extract an induced subtree keeping only the specified leaf names.
    ///
    /// This method creates a new species tree containing only the specified leaves
    /// and their most recent common ancestors (MRCAs). The resulting tree preserves
    /// the topology and branch lengths among the selected species.
    ///
    /// # Arguments
    /// * `names` - List of leaf names to keep in the subtree
    ///
    /// # Returns
    /// A new PySpeciesTree containing only the specified leaves and their MRCAs.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The names list is empty
    /// - No matching leaves are found in the tree
    /// - Failed to construct the induced subtree
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Simulate or load a species tree
    /// tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.5, seed=42)
    ///
    /// # Get all leaf names
    /// all_leaves = tree.leaf_names()
    /// print(f"Original tree has {len(all_leaves)} leaves")
    ///
    /// # Extract a subtree with only 3 species
    /// selected_species = all_leaves[:3]
    /// subtree = tree.extract_induced_subtree_by_names(selected_species)
    /// print(f"Subtree has {subtree.num_leaves()} leaves")
    /// ```
    fn extract_induced_subtree_by_names(&self, names: Vec<String>) -> PyResult<PySpeciesTree> {
        if names.is_empty() {
            return Err(PyValueError::new_err("Leaf names list cannot be empty"));
        }

        let (induced_tree, _) = extract_induced_subtree_by_names(&self.tree, &names)
            .ok_or_else(|| PyValueError::new_err(
                "Failed to extract induced subtree (no matching leaves found or subtree extraction failed)"
            ))?;

        Ok(PySpeciesTree {
            tree: Arc::new(induced_tree),
        })
    }

    /// Save the species tree to a Newick file.
    ///
    /// # Arguments
    /// * `filepath` - Path to save the Newick file
    fn save_newick(&self, filepath: &str) -> PyResult<()> {
        let newick = self.to_newick()?;
        fs::write(filepath, newick)
            .map_err(|e| PyValueError::new_err(format!("Failed to write Newick file: {}", e)))?;
        Ok(())
    }

    /// Generate an SVG visualization of the species tree using thirdkind.
    ///
    /// Requires thirdkind to be installed (`cargo install thirdkind`).
    ///
    /// # Arguments
    /// * `filepath` - Optional path to save the SVG file. If None, returns SVG as string.
    /// * `open_browser` - If True, opens the SVG in a web browser (default: False)
    /// * `internal_names` - If True, display internal node names (default: True)
    ///
    /// # Returns
    /// The SVG content as a string.
    #[pyo3(signature = (filepath=None, open_browser=false, internal_names=true))]
    fn to_svg(&self, filepath: Option<&str>, open_browser: bool, internal_names: bool) -> PyResult<String> {
        let temp_dir = std::env::temp_dir();
        let nwk_path = temp_dir.join("rustree_species_temp.nwk");
        let svg_path = temp_dir.join("rustree_species_temp.svg");

        // Write Newick to temp file
        let newick = self.to_newick()?;
        fs::write(&nwk_path, &newick)
            .map_err(|e| PyValueError::new_err(format!("Failed to write temp Newick: {}", e)))?;

        // Call thirdkind
        let mut cmd = Command::new("thirdkind");
        cmd.arg("--input-file").arg(&nwk_path)
           .arg("-o").arg(&svg_path);

        if internal_names {
            cmd.arg("-i");
        }

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

        let svg = fs::read_to_string(&svg_path)
            .map_err(|e| PyValueError::new_err(format!("Failed to read SVG output: {}", e)))?;

        if let Some(path) = filepath {
            fs::write(path, &svg)
                .map_err(|e| PyValueError::new_err(format!("Failed to write SVG file: {}", e)))?;
        }

        let _ = fs::remove_file(&nwk_path);
        let _ = fs::remove_file(&svg_path);

        Ok(svg)
    }

    /// Display the species tree visualization in a Jupyter notebook.
    ///
    /// Requires thirdkind to be installed and IPython/Jupyter environment.
    ///
    /// # Arguments
    /// * `internal_names` - If True, display internal node names (default: True)
    #[pyo3(signature = (internal_names=true))]
    fn display(&self, py: Python, internal_names: bool) -> PyResult<PyObject> {
        let svg = self.to_svg(None, false, internal_names)?;

        let ipython_display = py.import("IPython.display")?;
        let svg_class = ipython_display.getattr("SVG")?;
        let display_obj = svg_class.call1((svg,))?;

        Ok(display_obj.into())
    }

    /// Save birth-death events to a CSV file.
    ///
    /// For trees from simulate_species_tree, saves the actual simulation events.
    /// For trees from parse_species_tree, generates events by treating internal nodes
    /// as speciations and leaves as extant species (no extinction events).
    ///
    /// # Arguments
    /// * `filepath` - Path to save the CSV file
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Simulate a tree
    /// tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.3, seed=42)
    ///
    /// # Save birth-death events to CSV
    /// tree.save_bd_events_csv("events.csv")
    /// ```
    #[pyo3(signature = (filepath, eps=None))]
    fn save_bd_events_csv(&self, filepath: &str, eps: Option<f64>) -> PyResult<()> {
        use crate::bd::{generate_events_from_tree, generate_events_with_extinction};
        use crate::io::save_bd_events_to_csv;

        let events = match eps {
            Some(eps) => generate_events_with_extinction(&self.tree, eps),
            None => generate_events_from_tree(&self.tree),
        }.map_err(|e| PyValueError::new_err(e))?;

        save_bd_events_to_csv(&events, &self.tree, filepath)
            .map_err(|e| PyValueError::new_err(format!("Failed to write CSV file: {}", e)))?;

        Ok(())
    }

    /// Get birth-death events as a dictionary.
    ///
    /// Returns a dictionary with keys:
    /// - 'time': List of event times (floats)
    /// - 'node_name': List of node names (strings)
    /// - 'event_type': List of event types ('Speciation', 'Extinction', or 'Leaf')
    /// - 'child1_name': List of first child names (strings, empty for non-speciation events)
    /// - 'child2_name': List of second child names (strings, empty for non-speciation events)
    ///
    /// For trees from simulate_species_tree, returns the actual simulation events.
    /// For trees from parse_species_tree, generates events from tree structure.
    ///
    /// # Returns
    /// A dictionary containing lists of event data, suitable for pandas DataFrame creation.
    ///
    /// # Example
    /// ```python
    /// import rustree
    /// import pandas as pd
    ///
    /// # Simulate a tree
    /// tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.3, seed=42)
    ///
    /// # Get events as dict
    /// events = tree.get_bd_events()
    ///
    /// # Convert to pandas DataFrame
    /// df = pd.DataFrame(events)
    /// print(df.head())
    /// ```
    #[pyo3(signature = (eps=None))]
    fn get_bd_events(&self, py: Python, eps: Option<f64>) -> PyResult<PyObject> {
        use crate::bd::{generate_events_from_tree, generate_events_with_extinction};
        use pyo3::types::PyDict;

        let events = match eps {
            Some(eps) => generate_events_with_extinction(&self.tree, eps),
            None => generate_events_from_tree(&self.tree),
        }.map_err(|e| PyValueError::new_err(e))?;

        // Build dictionary with lists
        let dict = PyDict::new(py);

        let times: Vec<f64> = events.iter().map(|e| e.time).collect();
        let node_names: Vec<String> = events.iter()
            .map(|e| self.tree.nodes[e.node_id].name.clone())
            .collect();
        let event_types: Vec<&str> = events.iter()
            .map(|e| e.event_type.as_str())
            .collect();
        let child1_names: Vec<String> = events.iter()
            .map(|e| e.child1.map_or(String::new(), |c| self.tree.nodes[c].name.clone()))
            .collect();
        let child2_names: Vec<String> = events.iter()
            .map(|e| e.child2.map_or(String::new(), |c| self.tree.nodes[c].name.clone()))
            .collect();

        dict.set_item("time", times)?;
        dict.set_item("node_name", node_names)?;
        dict.set_item("event_type", event_types)?;
        dict.set_item("child1_name", child1_names)?;
        dict.set_item("child2_name", child2_names)?;

        Ok(dict.into())
    }

    /// Get Lineages Through Time (LTT) data for plotting.
    ///
    /// Returns a dictionary with:
    /// - `times`: Time points (going backward from present at 0)
    /// - `lineages`: Number of lineages alive at each time point
    ///
    /// This data can be directly plotted to visualize lineage accumulation/extinction
    /// over the tree's history.
    ///
    /// # Example
    /// ```python
    /// import rustree
    /// import matplotlib.pyplot as plt
    ///
    /// # Simulate or parse a species tree
    /// tree = rustree.simulate_species_tree(100, lambda_val=1.0, mu=0.5, seed=42)
    ///
    /// # Get LTT data
    /// ltt = tree.get_ltt_data()
    ///
    /// # Plot
    /// plt.plot(ltt['times'], ltt['lineages'])
    /// plt.xlabel('Time (past to present)')
    /// plt.ylabel('Number of lineages')
    /// plt.title('Lineages Through Time')
    /// plt.show()
    /// ```
    #[pyo3(signature = (eps=None))]
    fn get_ltt_data(&self, py: Python, eps: Option<f64>) -> PyResult<PyObject> {
        use crate::bd::{generate_events_from_tree, generate_events_with_extinction};
        use pyo3::types::PyDict;

        // Get all events
        let mut events = match eps {
            Some(eps) => generate_events_with_extinction(&self.tree, eps),
            None => generate_events_from_tree(&self.tree),
        }.map_err(|e| PyValueError::new_err(e))?;

        // Sort events by time (ascending, from present to past)
        events.sort_by(|a, b| a.time.total_cmp(&b.time));

        // Build LTT data
        let mut times = Vec::new();
        let mut lineages = Vec::new();

        // Start at present (time = 0) with n extant species
        let n_extant = self.tree.nodes.iter()
            .filter(|n| n.left_child.is_none() && n.right_child.is_none())
            .count();

        let mut current_lineages = n_extant;
        times.push(0.0);
        lineages.push(current_lineages);

        // Process events backward in time
        for event in events.iter() {
            match event.event_type.as_str() {
                "Leaf" => {
                    // Leaf events at present don't change count (already counted)
                    continue;
                }
                "Extinction" => {
                    // Going backward: an extinction means a lineage appears
                    current_lineages += 1;
                }
                "Speciation" => {
                    // Going backward: a speciation means two lineages merge into one
                    current_lineages -= 1;
                }
                _ => continue,
            }
            times.push(event.time);
            lineages.push(current_lineages);
        }

        // Create dictionary
        let dict = PyDict::new(py);
        dict.set_item("times", times)?;
        dict.set_item("lineages", lineages)?;

        Ok(dict.into())
    }

    /// Plot Lineages Through Time (LTT) using matplotlib.
    ///
    /// Creates and displays an LTT plot showing how the number of lineages
    /// changes over time.
    ///
    /// # Arguments
    /// * `filepath` - Optional path to save the plot. If None, displays interactively.
    /// * `title` - Optional plot title (default: "Lineages Through Time")
    /// * `xlabel` - Optional x-axis label (default: "Time (past to present)")
    /// * `ylabel` - Optional y-axis label (default: "Number of lineages")
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// tree = rustree.simulate_species_tree(100, lambda_val=1.0, mu=0.5, seed=42)
    ///
    /// # Display plot
    /// tree.plot_ltt()
    ///
    /// # Save to file
    /// tree.plot_ltt(filepath="ltt_plot.png", title="My LTT Plot")
    /// ```
    #[pyo3(signature = (filepath=None, title=None, xlabel=None, ylabel=None))]
    fn plot_ltt(
        &self,
        py: Python,
        filepath: Option<&str>,
        title: Option<&str>,
        xlabel: Option<&str>,
        ylabel: Option<&str>,
    ) -> PyResult<()> {
        // Get LTT data
        let ltt_data = self.get_ltt_data(py, None)?;

        // Import matplotlib
        let plt = py.import("matplotlib.pyplot")?;

        // Extract data from dictionary
        let dict = ltt_data.downcast::<pyo3::types::PyDict>(py)?;
        let times = dict.get_item("times")
            .ok_or_else(|| PyValueError::new_err("Missing 'times' in LTT data"))?;
        let lineages = dict.get_item("lineages")
            .ok_or_else(|| PyValueError::new_err("Missing 'lineages' in LTT data"))?;

        // Create plot
        plt.call_method1("plot", (times, lineages))?;
        plt.call_method1("xlabel", (xlabel.unwrap_or("Time (past to present)"),))?;
        plt.call_method1("ylabel", (ylabel.unwrap_or("Number of lineages)"),))?;
        plt.call_method1("title", (title.unwrap_or("Lineages Through Time"),))?;
        plt.call_method0("grid")?;

        // Save or show
        if let Some(path) = filepath {
            plt.call_method1("savefig", (path,))?;
        } else {
            plt.call_method0("show")?;
        }

        Ok(())
    }

    /// Simulate a gene tree along this species tree using the DTL model.
    ///
    /// # Arguments
    /// * `lambda_d` - Duplication rate per unit time
    /// * `lambda_t` - Transfer rate per unit time
    /// * `lambda_l` - Loss rate per unit time
    /// * `transfer_alpha` - Distance decay for assortative transfers (optional, None = uniform)
    /// * `require_extant` - If true, retry until a tree with extant genes is produced (default false)
    /// * `seed` - Random seed for reproducibility (optional)
    ///
    /// # Returns
    /// A PyGeneTree containing the simulated gene tree with its mapping to the species tree.
    #[pyo3(signature = (lambda_d, lambda_t, lambda_l, transfer_alpha=None, replacement_transfer=None, require_extant=false, seed=None))]
    fn simulate_dtl(&self, lambda_d: f64, lambda_t: f64, lambda_l: f64, transfer_alpha: Option<f64>, replacement_transfer: Option<f64>, require_extant: bool, seed: Option<u64>) -> PyResult<PyGeneTree> {
        validate_dtl_rates(lambda_d, lambda_t, lambda_l)?;
        validate_replacement_transfer(replacement_transfer)?;
        let mut rng = init_rng(seed);
        let origin_species = self.tree.root;
        let (rec_tree, _events) = simulate_dtl(&self.tree, origin_species, lambda_d, lambda_t, lambda_l, transfer_alpha, replacement_transfer, require_extant, &mut rng)
            .map_err(|e| PyValueError::new_err(e))?;

        // Extract owned data from RecTree
        let gene_tree = rec_tree.gene_tree;
        let node_mapping = rec_tree.node_mapping;
        let event_mapping = rec_tree.event_mapping;

        Ok(PyGeneTree {
            gene_tree,
            species_tree: Arc::clone(&self.tree),
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
    /// * `require_extant` - If true, only include trees with extant genes (default false)
    /// * `seed` - Random seed for reproducibility (optional)
    ///
    /// # Returns
    /// A list of PyGeneTree objects.
    #[pyo3(signature = (n, lambda_d, lambda_t, lambda_l, transfer_alpha=None, replacement_transfer=None, require_extant=false, seed=None))]
    fn simulate_dtl_batch(&self, n: usize, lambda_d: f64, lambda_t: f64, lambda_l: f64, transfer_alpha: Option<f64>, replacement_transfer: Option<f64>, require_extant: bool, seed: Option<u64>) -> PyResult<Vec<PyGeneTree>> {
        validate_dtl_rates(lambda_d, lambda_t, lambda_l)?;
        validate_replacement_transfer(replacement_transfer)?;
        let mut rng = init_rng(seed);

        let origin_species = self.tree.root;
        let (rec_trees, _all_events) = simulate_dtl_batch(
            &self.tree,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
            n,
            require_extant,
            &mut rng,
        ).map_err(|e| PyValueError::new_err(e))?;

        let gene_trees: Vec<PyGeneTree> = rec_trees
            .into_iter()
            .map(|rec_tree| PyGeneTree {
                gene_tree: rec_tree.gene_tree,
                species_tree: Arc::clone(&self.tree),
                node_mapping: rec_tree.node_mapping,
                event_mapping: rec_tree.event_mapping,
            })
            .collect();

        Ok(gene_trees)
    }

    /// Simulate a gene tree using the Zombi-style per-species DTL model.
    ///
    /// In this model, the event rate is proportional to the number of SPECIES
    /// with gene copies, NOT the number of gene copies themselves. This means
    /// duplications don't increase the event rate.
    ///
    /// When an event occurs:
    /// 1. A random species (among those with genes) is chosen
    /// 2. A random gene copy within that species is affected
    ///
    /// This is the model used by Zombi and is appropriate when DTL events are
    /// driven by species-level processes rather than gene-copy-level processes.
    ///
    /// # Arguments
    /// * `lambda_d` - Duplication rate per species per unit time
    /// * `lambda_t` - Transfer rate per species per unit time
    /// * `lambda_l` - Loss rate per species per unit time
    /// * `transfer_alpha` - Distance decay for assortative transfers (optional, None = uniform)
    /// * `require_extant` - If true, retry until at least one extant gene exists (default false)
    /// * `seed` - Random seed for reproducibility (optional)
    ///
    /// # Returns
    /// A PyGeneTree containing the simulated gene tree with its mapping to the species tree.
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Simulate species tree
    /// species_tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)
    ///
    /// # Compare per-gene-copy model vs per-species model
    /// gene_tree_per_copy = species_tree.simulate_dtl(0.5, 0.2, 0.3, seed=123)
    /// gene_tree_per_species = species_tree.simulate_dtl_per_species(0.5, 0.2, 0.3, seed=123)
    ///
    /// print(f"Per-copy model: {gene_tree_per_copy.num_leaves()} leaves")
    /// print(f"Per-species model: {gene_tree_per_species.num_leaves()} leaves")
    /// # Per-species typically has fewer leaves since duplications don't increase rate
    /// ```
    #[pyo3(signature = (lambda_d, lambda_t, lambda_l, transfer_alpha=None, replacement_transfer=None, require_extant=false, seed=None))]
    fn simulate_dtl_per_species(&self, lambda_d: f64, lambda_t: f64, lambda_l: f64, transfer_alpha: Option<f64>, replacement_transfer: Option<f64>, require_extant: bool, seed: Option<u64>) -> PyResult<PyGeneTree> {
        validate_dtl_rates(lambda_d, lambda_t, lambda_l)?;
        validate_replacement_transfer(replacement_transfer)?;
        let mut rng = init_rng(seed);
        let origin_species = self.tree.root;
        let (rec_tree, _events) = simulate_dtl_per_species(&self.tree, origin_species, lambda_d, lambda_t, lambda_l, transfer_alpha, replacement_transfer, require_extant, &mut rng)
            .map_err(|e| PyValueError::new_err(e))?;

        // Extract owned data from RecTree
        let gene_tree = rec_tree.gene_tree;
        let node_mapping = rec_tree.node_mapping;
        let event_mapping = rec_tree.event_mapping;

        Ok(PyGeneTree {
            gene_tree,
            species_tree: Arc::clone(&self.tree),
            node_mapping,
            event_mapping,
        })
    }

    /// Simulate multiple gene trees using the Zombi-style per-species DTL model.
    ///
    /// This is faster than calling simulate_dtl_per_species multiple times because
    /// species tree events and LCA depths are computed only once.
    ///
    /// # Arguments
    /// * `n` - Number of gene trees to simulate
    /// * `lambda_d` - Duplication rate per species per unit time
    /// * `lambda_t` - Transfer rate per species per unit time
    /// * `lambda_l` - Loss rate per species per unit time
    /// * `transfer_alpha` - Distance decay for assortative transfers (optional, None = uniform)
    /// * `require_extant` - If true, only include trees with extant genes (default false)
    /// * `seed` - Random seed for reproducibility (optional)
    ///
    /// # Returns
    /// A list of PyGeneTree objects.
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Simulate species tree
    /// species_tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)
    ///
    /// # Simulate 100 gene trees efficiently
    /// gene_trees = species_tree.simulate_dtl_per_species_batch(
    ///     n=100,
    ///     lambda_d=0.5,
    ///     lambda_t=0.2,
    ///     lambda_l=0.3,
    ///     seed=123
    /// )
    /// print(f"Generated {len(gene_trees)} gene trees")
    /// ```
    #[pyo3(signature = (n, lambda_d, lambda_t, lambda_l, transfer_alpha=None, replacement_transfer=None, require_extant=false, seed=None))]
    fn simulate_dtl_per_species_batch(&self, n: usize, lambda_d: f64, lambda_t: f64, lambda_l: f64, transfer_alpha: Option<f64>, replacement_transfer: Option<f64>, require_extant: bool, seed: Option<u64>) -> PyResult<Vec<PyGeneTree>> {
        validate_dtl_rates(lambda_d, lambda_t, lambda_l)?;
        validate_replacement_transfer(replacement_transfer)?;
        let mut rng = init_rng(seed);
        let origin_species = self.tree.root;
        let (rec_trees, _all_events) = simulate_dtl_per_species_batch(
            &self.tree,
            origin_species,
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
            n,
            require_extant,
            &mut rng,
        ).map_err(|e| PyValueError::new_err(e))?;

        let gene_trees: Vec<PyGeneTree> = rec_trees
            .into_iter()
            .map(|rec_tree| PyGeneTree {
                gene_tree: rec_tree.gene_tree,
                species_tree: Arc::clone(&self.tree),
                node_mapping: rec_tree.node_mapping,
                event_mapping: rec_tree.event_mapping,
            })
            .collect();

        Ok(gene_trees)
    }

    /// Compute all pairwise distances between nodes in the species tree.
    ///
    /// # Arguments
    /// * `distance_type` - Type of distance: "topological" (number of edges) or "metric" (sum of branch lengths)
    /// * `leaves_only` - If true, only compute distances between leaf nodes (default true)
    ///
    /// # Returns
    /// A pandas DataFrame with columns: node1, node2, distance
    ///
    /// # Example
    /// ```python
    /// import rustree
    /// tree = rustree.simulate_species_tree(5, 1.0, 0.5, seed=42)
    /// df = tree.pairwise_distances("metric", leaves_only=True)
    /// print(df)
    /// ```
    #[pyo3(signature = (distance_type, leaves_only=true))]
    fn pairwise_distances(&self, py: Python, distance_type: &str, leaves_only: bool) -> PyResult<PyObject> {
        use crate::metric_functions::DistanceType;

        let dist_type = match distance_type.to_lowercase().as_str() {
            "topological" | "topo" => DistanceType::Topological,
            "metric" | "branch" | "patristic" => DistanceType::Metric,
            _ => return Err(PyValueError::new_err(format!(
                "Invalid distance_type '{}'. Use 'topological' or 'metric'.",
                distance_type
            ))),
        };

        let distances = self.tree.pairwise_distances(dist_type, leaves_only)
            .map_err(|e| PyValueError::new_err(format!("Failed to compute pairwise distances: {}", e)))?;

        let node1: Vec<&str> = distances.iter().map(|d| d.node1).collect();
        let node2: Vec<&str> = distances.iter().map(|d| d.node2).collect();
        let dist: Vec<f64> = distances.iter().map(|d| d.distance).collect();

        // Create pandas DataFrame
        let pandas = py.import("pandas")?;
        let dict = pyo3::types::PyDict::new(py);

        dict.set_item("node1", node1)?;
        dict.set_item("node2", node2)?;
        dict.set_item("distance", dist)?;

        let df = pandas.call_method1("DataFrame", (dict,))?;
        Ok(df.into())
    }

    /// Save pairwise distances between nodes to a CSV file.
    ///
    /// Computes all pairwise distances and writes them to a CSV file with the format:
    /// node1,node2,distance
    ///
    /// # Arguments
    /// * `filepath` - Path to save the CSV file
    /// * `distance_type` - Type of distance: "topological" (number of edges) or "metric" (sum of branch lengths)
    /// * `leaves_only` - If true, only compute distances between leaf nodes (default: true)
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Simulate a species tree
    /// tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.5, seed=42)
    ///
    /// # Save topological distances between all leaves
    /// tree.save_pairwise_distances_csv("distances.csv", "topological", leaves_only=True)
    ///
    /// # Save metric distances between all nodes
    /// tree.save_pairwise_distances_csv("all_distances.csv", "metric", leaves_only=False)
    /// ```
    #[pyo3(signature = (filepath, distance_type, leaves_only=true))]
    fn save_pairwise_distances_csv(&self, filepath: &str, distance_type: &str, leaves_only: bool) -> PyResult<()> {
        use crate::metric_functions::{DistanceType, PairwiseDistance};

        // Parse distance type
        let dist_type = match distance_type.to_lowercase().as_str() {
            "topological" | "topo" => DistanceType::Topological,
            "metric" | "branch" | "patristic" => DistanceType::Metric,
            _ => return Err(PyValueError::new_err(format!(
                "Invalid distance_type '{}'. Use 'topological' or 'metric'.",
                distance_type
            ))),
        };

        // Compute pairwise distances
        let distances = self.tree.pairwise_distances(dist_type, leaves_only)
            .map_err(|e| PyValueError::new_err(format!("Failed to compute pairwise distances: {}", e)))?;

        // Write to CSV
        let mut csv = String::from(PairwiseDistance::csv_header());
        csv.push('\n');
        for d in &distances {
            csv.push_str(&d.to_csv_row());
            csv.push('\n');
        }

        fs::write(filepath, csv)
            .map_err(|e| PyValueError::new_err(format!("Failed to write CSV file: {}", e)))?;

        Ok(())
    }
}

/// A gene tree simulated under the DTL model.
#[pyclass]
#[derive(Clone)]
pub struct PyGeneTree {
    gene_tree: FlatTree,
    species_tree: Arc<FlatTree>,
    node_mapping: Vec<Option<usize>>,
    event_mapping: Vec<Event>,
}

#[pymethods]
impl PyGeneTree {
    /// Convert the gene tree to Newick format.
    fn to_newick(&self) -> PyResult<String> {
        let nwk = self.gene_tree.to_newick()
            .map_err(|e| PyValueError::new_err(e))?;
        Ok(nwk + ";")
    }

    /// Save the gene tree to a Newick file.
    ///
    /// # Arguments
    /// * `filepath` - Path to save the Newick file
    fn save_newick(&self, filepath: &str) -> PyResult<()> {
        let newick = self.to_newick()?;
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

        // Extract induced subtree with old→new index mapping
        let (sampled_tree, old_to_new) = extract_induced_subtree(&self.gene_tree, &extant_indices)
            .ok_or_else(|| PyValueError::new_err("Failed to extract induced subtree"))?;

        // Build new→old mapping by inverting old_to_new
        let mut new_to_old: Vec<Option<usize>> = vec![None; sampled_tree.nodes.len()];
        for (old_idx, new_idx_opt) in old_to_new.iter().enumerate() {
            if let Some(new_idx) = new_idx_opt {
                new_to_old[*new_idx] = Some(old_idx);
            }
        }

        // Rebuild event mapping using index-based lookup
        let new_event_mapping: Vec<Event> = new_to_old.iter()
            .enumerate()
            .map(|(new_idx, old_idx_opt)| {
                if let Some(old_idx) = old_idx_opt {
                    self.event_mapping[*old_idx].clone()
                } else {
                    // Collapsed node — shouldn't happen since extract_induced_subtree
                    // only keeps nodes marked Keep
                    if sampled_tree.nodes[new_idx].left_child.is_none() {
                        Event::Leaf
                    } else {
                        Event::Speciation
                    }
                }
            })
            .collect();

        // Rebuild node mapping: leaves get species from gene name, internal nodes get None
        let new_node_mapping: Vec<Option<usize>> = sampled_tree.nodes.iter()
            .map(|n| {
                if n.left_child.is_none() && n.right_child.is_none() {
                    if let Some(pos) = n.name.rfind('_') {
                        let species_name = &n.name[..pos];
                        self.species_tree.nodes.iter()
                            .position(|sn| sn.name == species_name)
                    } else {
                        None
                    }
                } else {
                    None // Internal node: no valid mapping after sampling
                }
            })
            .collect();

        Ok(PyGeneTree {
            gene_tree: sampled_tree,
            species_tree: Arc::clone(&self.species_tree),
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

        let (sampled_tree, old_to_new) = extract_induced_subtree(&self.gene_tree, &keep_indices)
            .ok_or_else(|| PyValueError::new_err("Failed to extract induced subtree"))?;

        // Build new→old mapping by inverting old_to_new
        let mut new_to_old: Vec<Option<usize>> = vec![None; sampled_tree.nodes.len()];
        for (old_idx, new_idx_opt) in old_to_new.iter().enumerate() {
            if let Some(new_idx) = new_idx_opt {
                new_to_old[*new_idx] = Some(old_idx);
            }
        }

        // Rebuild event mapping using index-based lookup
        let new_event_mapping: Vec<Event> = new_to_old.iter()
            .enumerate()
            .map(|(new_idx, old_idx_opt)| {
                if let Some(old_idx) = old_idx_opt {
                    self.event_mapping[*old_idx].clone()
                } else {
                    if sampled_tree.nodes[new_idx].left_child.is_none() {
                        Event::Leaf
                    } else {
                        Event::Speciation
                    }
                }
            })
            .collect();

        // Rebuild node mapping: leaves get species from gene name, internal nodes get None
        let new_node_mapping: Vec<Option<usize>> = sampled_tree.nodes.iter()
            .map(|n| {
                if n.left_child.is_none() && n.right_child.is_none() {
                    if let Some(pos) = n.name.rfind('_') {
                        let species_name = &n.name[..pos];
                        self.species_tree.nodes.iter()
                            .position(|sn| sn.name == species_name)
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();

        Ok(PyGeneTree {
            gene_tree: sampled_tree,
            species_tree: Arc::clone(&self.species_tree),
            node_mapping: new_node_mapping,
            event_mapping: new_event_mapping,
        })
    }

    /// Sample species tree leaves and filter gene tree accordingly.
    ///
    /// This method samples a subset of species from the species tree and automatically
    /// filters the gene tree to keep only genes that map to the sampled species. The
    /// reconciliation mappings are preserved using an LCA-based approach.
    ///
    /// # Arguments
    /// * `species_leaf_names` - List of species leaf names to keep
    ///
    /// # Returns
    /// A new PyGeneTree with sampled species and gene trees, with preserved reconciliation mappings.
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Load a reconciled tree
    /// gt = rustree.parse_recphyloxml("reconciliation.xml")
    ///
    /// # Sample only species A, B, and C
    /// sampled_gt = gt.sample_species_leaves(["species_A", "species_B", "species_C"])
    ///
    /// print(f"Original: {gt.num_nodes()} gene nodes, {len(gt.to_newick())} chars")
    /// print(f"Sampled: {sampled_gt.num_nodes()} gene nodes, {len(sampled_gt.to_newick())} chars")
    /// ```
    fn sample_species_leaves(&self, species_leaf_names: Vec<String>) -> PyResult<PyGeneTree> {
        use crate::node::RecTreeOwned;

        // Build a RecTreeOwned from self (deep-clone species tree out of Arc for owned use)
        let rec_tree_owned = RecTreeOwned::new(
            (*self.species_tree).clone(),
            self.gene_tree.clone(),
            self.node_mapping.clone(),
            self.event_mapping.clone(),
        );

        // Sample using the RecTreeOwned method
        let sampled_rec_tree = rec_tree_owned
            .sample_species_leaves(&species_leaf_names)
            .map_err(|e| PyValueError::new_err(format!("Failed to sample species leaves: {}", e)))?;

        // Convert back to PyGeneTree (wrap species tree in Arc for sharing)
        Ok(PyGeneTree {
            gene_tree: sampled_rec_tree.gene_tree,
            species_tree: Arc::new(sampled_rec_tree.species_tree),
            node_mapping: sampled_rec_tree.node_mapping,
            event_mapping: sampled_rec_tree.event_mapping,
        })
    }

    /// Export the reconciled tree to RecPhyloXML format as a string.
    fn to_xml(&self) -> String {
        // Rebuild RecTree temporarily for XML export
        use crate::node::RecTree;
        let rec_tree = RecTree::new(
            &*self.species_tree,
            &self.gene_tree,
            &self.node_mapping,
            &self.event_mapping,
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
    /// Returns a DataFrame with columns: node_id, name, parent, left_child,
    /// left_child_name, right_child, right_child_name, length, depth,
    /// species_node, species_node_left, species_node_right, event
    ///
    /// # Arguments
    /// * `filepath` - Optional path to save the CSV file
    ///
    /// # Returns
    /// A pandas DataFrame with the gene tree data.
    #[pyo3(signature = (filepath=None))]
    fn to_csv(&self, py: Python, filepath: Option<&str>) -> PyResult<PyObject> {
        let rec_tree = RecTree::new(
            &self.species_tree,
            &self.gene_tree,
            &self.node_mapping,
            &self.event_mapping,
        );
        let cols = rec_tree.to_columns();

        if let Some(path) = filepath {
            cols.save_csv(path)
                .map_err(|e| PyValueError::new_err(format!("Failed to write CSV: {}", e)))?;
        }

        let pandas = py.import("pandas")?;
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("node_id", cols.node_id)?;
        dict.set_item("name", cols.name)?;
        dict.set_item("parent", cols.parent)?;
        dict.set_item("left_child", cols.left_child)?;
        dict.set_item("left_child_name", cols.left_child_name)?;
        dict.set_item("right_child", cols.right_child)?;
        dict.set_item("right_child_name", cols.right_child_name)?;
        dict.set_item("length", cols.length)?;
        dict.set_item("depth", cols.depth)?;
        dict.set_item("species_node", cols.species_node)?;
        dict.set_item("species_node_left", cols.species_node_left)?;
        dict.set_item("species_node_right", cols.species_node_right)?;
        dict.set_item("event", cols.event)?;

        let df = pandas.call_method1("DataFrame", (dict,))?;
        Ok(df.into())
    }

    /// Compute all pairwise distances between nodes in the gene tree.
    ///
    /// # Arguments
    /// * `distance_type` - Type of distance: "topological" (number of edges) or "metric" (sum of branch lengths)
    /// * `leaves_only` - If true, only compute distances between leaf nodes (default true)
    ///
    /// # Returns
    /// A pandas DataFrame with columns: node1, node2, distance
    ///
    /// # Example
    /// ```python
    /// import rustree
    /// species_tree = rustree.simulate_species_tree(5, 1.0, 0.5, seed=42)
    /// gene_tree = species_tree.simulate_dtl(0.5, 0.2, 0.3, seed=123)
    /// df = gene_tree.pairwise_distances("metric", leaves_only=True)
    /// print(df)
    /// ```
    #[pyo3(signature = (distance_type, leaves_only=true))]
    fn pairwise_distances(&self, py: Python, distance_type: &str, leaves_only: bool) -> PyResult<PyObject> {
        use crate::metric_functions::DistanceType;

        let dist_type = match distance_type.to_lowercase().as_str() {
            "topological" | "topo" => DistanceType::Topological,
            "metric" | "branch" | "patristic" => DistanceType::Metric,
            _ => return Err(PyValueError::new_err(format!(
                "Invalid distance_type '{}'. Use 'topological' or 'metric'.",
                distance_type
            ))),
        };

        let distances = self.gene_tree.pairwise_distances(dist_type, leaves_only)
            .map_err(|e| PyValueError::new_err(format!("Failed to compute pairwise distances: {}", e)))?;

        let node1: Vec<&str> = distances.iter().map(|d| d.node1).collect();
        let node2: Vec<&str> = distances.iter().map(|d| d.node2).collect();
        let dist: Vec<f64> = distances.iter().map(|d| d.distance).collect();

        // Create pandas DataFrame
        let pandas = py.import("pandas")?;
        let dict = pyo3::types::PyDict::new(py);

        dict.set_item("node1", node1)?;
        dict.set_item("node2", node2)?;
        dict.set_item("distance", dist)?;

        let df = pandas.call_method1("DataFrame", (dict,))?;
        Ok(df.into())
    }

    /// Save pairwise distances between nodes in the gene tree to a CSV file.
    ///
    /// Computes all pairwise distances and writes them to a CSV file with the format:
    /// node1,node2,distance
    ///
    /// # Arguments
    /// * `filepath` - Path to save the CSV file
    /// * `distance_type` - Type of distance: "topological" (number of edges) or "metric" (sum of branch lengths)
    /// * `leaves_only` - If true, only compute distances between leaf nodes (default: true)
    ///
    /// # Example
    /// ```python
    /// import rustree
    ///
    /// # Simulate a species tree and gene tree
    /// species_tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.5, seed=42)
    /// gene_tree = species_tree.simulate_dtl(0.5, 0.2, 0.3, seed=123)
    ///
    /// # Save topological distances between all leaves
    /// gene_tree.save_pairwise_distances_csv("gene_distances.csv", "topological", leaves_only=True)
    ///
    /// # Save metric distances between all nodes
    /// gene_tree.save_pairwise_distances_csv("all_gene_distances.csv", "metric", leaves_only=False)
    /// ```
    #[pyo3(signature = (filepath, distance_type, leaves_only=true))]
    fn save_pairwise_distances_csv(&self, filepath: &str, distance_type: &str, leaves_only: bool) -> PyResult<()> {
        use crate::metric_functions::{DistanceType, PairwiseDistance};

        // Parse distance type
        let dist_type = match distance_type.to_lowercase().as_str() {
            "topological" | "topo" => DistanceType::Topological,
            "metric" | "branch" | "patristic" => DistanceType::Metric,
            _ => return Err(PyValueError::new_err(format!(
                "Invalid distance_type '{}'. Use 'topological' or 'metric'.",
                distance_type
            ))),
        };

        // Compute pairwise distances
        let distances = self.gene_tree.pairwise_distances(dist_type, leaves_only)
            .map_err(|e| PyValueError::new_err(format!("Failed to compute pairwise distances: {}", e)))?;

        // Write to CSV
        let mut csv = String::from(PairwiseDistance::csv_header());
        csv.push('\n');
        for d in &distances {
            csv.push_str(&d.to_csv_row());
            csv.push('\n');
        }

        fs::write(filepath, csv)
            .map_err(|e| PyValueError::new_err(format!("Failed to write CSV file: {}", e)))?;

        Ok(())
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

    let (mut tree, _events) = simulate_bd_tree_bwd(n, lambda_, mu, &mut rng);
    tree.assign_depths();

    Ok(PySpeciesTree { tree: Arc::new(tree) })
}

/// Parse a Newick string or file into a species tree.
///
/// # Arguments
/// * `newick_str` - A Newick formatted string or path to a Newick file
///
/// # Returns
/// A PySpeciesTree parsed from the Newick string.
///
/// # Example
/// ```python
/// # From Newick string
/// tree = rustree.parse_species_tree("(A:1.0,B:1.0):0.5;")
///
/// # From file
/// tree = rustree.parse_species_tree("species.nwk")
/// ```
#[pyfunction]
fn parse_species_tree(newick_str: &str) -> PyResult<PySpeciesTree> {
    use crate::newick::newick::parse_newick;
    use std::path::Path;

    // Check if input is a file path
    let newick_content = if Path::new(newick_str).exists() {
        fs::read_to_string(newick_str)
            .map_err(|e| PyValueError::new_err(format!("Failed to read file '{}': {}", newick_str, e)))?
    } else {
        newick_str.to_string()
    };

    let mut nodes = parse_newick(&newick_content)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse Newick: {}", e)))?;

    let mut root = nodes.pop()
        .ok_or_else(|| PyValueError::new_err("No tree found in Newick string"))?;

    root.assign_depths(0.0);
    let tree = root.to_flat_tree();

    Ok(PySpeciesTree { tree: Arc::new(tree) })
}

/// Parse a RecPhyloXML file into a gene tree with reconciliation information.
///
/// # Arguments
/// * `filepath` - Path to the RecPhyloXML file (e.g., from ALERax)
///
/// # Returns
/// A PyGeneTree parsed from the RecPhyloXML file, containing both the species tree
/// and reconciled gene tree with event mappings.
#[pyfunction]
fn parse_recphyloxml(filepath: &str) -> PyResult<PyGeneTree> {
    use crate::node::RecTreeOwned;

    let rec_tree_owned = RecTreeOwned::from_xml_file(filepath)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse RecPhyloXML: {}", e)))?;

    Ok(PyGeneTree {
        gene_tree: rec_tree_owned.gene_tree,
        species_tree: Arc::new(rec_tree_owned.species_tree),
        node_mapping: rec_tree_owned.node_mapping,
        event_mapping: rec_tree_owned.event_mapping,
    })
}

/// Event counts from reconciliation analysis
#[pyclass]
#[derive(Clone)]
pub struct PyEventCounts {
    /// Number of speciation events
    #[pyo3(get)]
    pub speciations: f64,

    /// Number of speciation loss events
    #[pyo3(get)]
    pub speciation_losses: f64,

    /// Number of duplication events
    #[pyo3(get)]
    pub duplications: f64,

    /// Number of duplication loss events
    #[pyo3(get)]
    pub duplication_losses: f64,

    /// Number of transfer events
    #[pyo3(get)]
    pub transfers: f64,

    /// Number of transfer loss events
    #[pyo3(get)]
    pub transfer_losses: f64,

    /// Number of loss events
    #[pyo3(get)]
    pub losses: f64,

    /// Number of leaf nodes
    #[pyo3(get)]
    pub leaves: f64,
}

#[pymethods]
impl PyEventCounts {
    fn __repr__(&self) -> String {
        format!(
            "EventCounts(S={:.1}, D={:.1}, T={:.1}, L={:.1})",
            self.speciations,
            self.duplications,
            self.transfers,
            self.losses
        )
    }
}

/// Summary statistics for reconciliation analysis
#[pyclass]
#[derive(Clone)]
pub struct PyReconciliationStatistics {
    /// Mean event counts across all samples
    #[pyo3(get)]
    pub mean_event_counts: PyEventCounts,

    /// Mean transfers between species pairs
    /// Returns a dict of dicts: {source_species: {dest_species: mean_count}}
    #[pyo3(get)]
    pub mean_transfers: HashMap<String, HashMap<String, f64>>,

    /// Mean events per species node
    /// Returns a dict: {species_name: EventCounts}
    #[pyo3(get)]
    pub events_per_species: HashMap<String, PyEventCounts>,
}

#[pymethods]
impl PyReconciliationStatistics {
    fn __repr__(&self) -> String {
        let total_transfers: f64 = self.mean_transfers.values()
            .flat_map(|dests| dests.values())
            .sum();
        format!(
            "ReconciliationStatistics(mean_events={}, mean_transfers={:.1}, species_count={})",
            self.mean_event_counts.__repr__(),
            total_transfers,
            self.events_per_species.len()
        )
    }

    /// Returns a pandas DataFrame of mean transfers between species pairs.
    ///
    /// Columns: source, destination, mean_count
    fn transfers_df(&self, py: Python) -> PyResult<PyObject> {
        let pandas = py.import("pandas")?;
        let dict = pyo3::types::PyDict::new(py);

        let mut sources: Vec<String> = Vec::new();
        let mut destinations: Vec<String> = Vec::new();
        let mut counts: Vec<f64> = Vec::new();

        for (source, dests) in &self.mean_transfers {
            for (dest, count) in dests {
                sources.push(source.clone());
                destinations.push(dest.clone());
                counts.push(*count);
            }
        }

        dict.set_item("source", sources)?;
        dict.set_item("destination", destinations)?;
        dict.set_item("mean_count", counts)?;

        let df = pandas.call_method1("DataFrame", (dict,))?;
        Ok(df.into())
    }

    /// Returns a pandas DataFrame of mean events per species.
    ///
    /// Columns: species, speciations, speciation_losses, duplications,
    ///          duplication_losses, transfers, transfer_losses, losses, leaves
    fn events_df(&self, py: Python) -> PyResult<PyObject> {
        let pandas = py.import("pandas")?;
        let dict = pyo3::types::PyDict::new(py);

        let mut species_names: Vec<String> = Vec::new();
        let mut speciations: Vec<f64> = Vec::new();
        let mut speciation_losses: Vec<f64> = Vec::new();
        let mut duplications: Vec<f64> = Vec::new();
        let mut duplication_losses: Vec<f64> = Vec::new();
        let mut transfers: Vec<f64> = Vec::new();
        let mut transfer_losses: Vec<f64> = Vec::new();
        let mut losses: Vec<f64> = Vec::new();
        let mut leaves: Vec<f64> = Vec::new();

        // Sort by species name for consistent output
        let mut sorted_species: Vec<(&String, &PyEventCounts)> =
            self.events_per_species.iter().collect();
        sorted_species.sort_by(|(a, _), (b, _)| a.cmp(b));

        for (species, counts) in sorted_species {
            species_names.push(species.clone());
            speciations.push(counts.speciations);
            speciation_losses.push(counts.speciation_losses);
            duplications.push(counts.duplications);
            duplication_losses.push(counts.duplication_losses);
            transfers.push(counts.transfers);
            transfer_losses.push(counts.transfer_losses);
            losses.push(counts.losses);
            leaves.push(counts.leaves);
        }

        dict.set_item("species", species_names)?;
        dict.set_item("speciations", speciations)?;
        dict.set_item("speciation_losses", speciation_losses)?;
        dict.set_item("duplications", duplications)?;
        dict.set_item("duplication_losses", duplication_losses)?;
        dict.set_item("transfers", transfers)?;
        dict.set_item("transfer_losses", transfer_losses)?;
        dict.set_item("losses", losses)?;
        dict.set_item("leaves", leaves)?;

        let df = pandas.call_method1("DataFrame", (dict,))?;
        Ok(df.into())
    }
}

/// Result of ALERax reconciliation for a gene family.
///
/// Contains all reconciliation samples, estimated evolutionary rates,
/// log-likelihood, and summary statistics.
#[pyclass]
pub struct PyAleRaxResult {
    /// All reconciliation samples (typically 100)
    #[pyo3(get)]
    gene_trees: Vec<PyGeneTree>,

    /// Estimated duplication rate
    #[pyo3(get)]
    duplication_rate: f64,

    /// Estimated loss rate
    #[pyo3(get)]
    loss_rate: f64,

    /// Estimated transfer rate
    #[pyo3(get)]
    transfer_rate: f64,

    /// Log-likelihood of the reconciliation
    #[pyo3(get)]
    likelihood: f64,

    /// Summary statistics across all reconciliation samples
    #[pyo3(get)]
    statistics: PyReconciliationStatistics,
}

#[pymethods]
impl PyAleRaxResult {
    /// Get a summary string of the reconciliation result.
    fn __repr__(&self) -> String {
        format!(
            "PyAleRaxResult(samples={}, D={:.4}, L={:.4}, T={:.4}, logL={:.2})",
            self.gene_trees.len(),
            self.duplication_rate,
            self.loss_rate,
            self.transfer_rate,
            self.likelihood
        )
    }
}

/// Reconcile gene trees with species tree using ALERax.
///
/// Calls the external ALERax tool to perform phylogenetic reconciliation,
/// inferring duplication, transfer, and loss events. Returns reconciliation
/// samples that can be analyzed, visualized, and exported.
///
/// # Arguments
/// * `species_tree` - Species tree as PySpeciesTree, Newick string, or file path
/// * `gene_trees` - Gene trees as:
///   - Single Newick string or file path (creates "family_0")
///   - List of Newick strings/paths (creates "family_0", "family_1", ...)
///   - Dict mapping family names to Newick strings/paths
/// * `output_dir` - Optional output directory (default: temporary directory)
/// * `num_samples` - Number of reconciliation samples per family (default: 100)
/// * `model` - Model parametrization: "PER-FAMILY" or "GLOBAL" (default: "PER-FAMILY")
/// * `seed` - Random seed for reproducibility (default: None)
/// * `keep_output` - Whether to preserve ALERax output files (default: False)
/// * `alerax_path` - Path to alerax executable (default: "alerax")
///
/// # Returns
/// Dictionary mapping family names to PyAleRaxResult objects containing
/// reconciled gene trees and estimated evolutionary rates.
///
/// # Example
/// ```python
/// import rustree
///
/// species_tree = rustree.parse_species_tree("species.nwk")
/// results = rustree.reconcile_with_alerax(species_tree, "gene.nwk", seed=42)
///
/// result = results["family_0"]
/// print(f"D={result.duplication_rate}, T={result.transfer_rate}")
///
/// # Access reconciliation samples
/// best_tree = result.gene_trees[0]
/// best_tree.to_svg("reconciled.svg")
/// ```
#[pyfunction]
#[pyo3(signature = (species_tree, gene_trees, output_dir=None, num_samples=100, model="PER-FAMILY".to_string(), seed=None, keep_output=false, alerax_path="alerax".to_string()))]
fn reconcile_with_alerax(
    py: Python,
    species_tree: PyObject,
    gene_trees: PyObject,
    output_dir: Option<String>,
    num_samples: usize,
    model: String,
    seed: Option<u64>,
    keep_output: bool,
    alerax_path: String,
) -> PyResult<HashMap<String, PyAleRaxResult>> {
    use crate::alerax::{run_alerax, AleRaxConfig, GeneFamily, ModelType, validate_inputs};
    use std::fs;
    use std::path::PathBuf;
    use tempfile::TempDir;

    // Parse model type
    let model_type = match model.to_uppercase().as_str() {
        "PER-FAMILY" | "PER_FAMILY" => ModelType::PerFamily,
        "GLOBAL" => ModelType::Global,
        _ => return Err(PyValueError::new_err(format!(
            "Invalid model type '{}'. Must be 'PER-FAMILY' or 'GLOBAL'", model
        ))),
    };

    // Parse species tree
    let species_tree_obj = species_tree.extract::<PySpeciesTree>(py)
        .or_else(|_| {
            // Try as string (Newick or file path)
            let tree_str = species_tree.extract::<String>(py)?;

            // Check if it's a file path
            if PathBuf::from(&tree_str).exists() {
                // Read file and parse
                let content = fs::read_to_string(&tree_str)
                    .map_err(|e| PyValueError::new_err(format!("Failed to read species tree file: {}", e)))?;
                parse_species_tree(&content)
            } else {
                // Parse as Newick string
                parse_species_tree(&tree_str)
            }
        })?;

    let species_tree_newick = species_tree_obj.tree.to_newick()
        .map_err(|e| PyValueError::new_err(e))?;

    // Parse gene trees
    let mut gene_tree_list: Vec<(String, String)> = Vec::new();

    // Try as dict first
    if let Ok(dict) = gene_trees.extract::<HashMap<String, String>>(py) {
        for (name, newick_or_path) in dict {
            let newick = if PathBuf::from(&newick_or_path).exists() {
                fs::read_to_string(&newick_or_path)
                    .map_err(|e| PyValueError::new_err(format!(
                        "Failed to read gene tree file '{}': {}", newick_or_path, e
                    )))?
            } else {
                newick_or_path
            };
            gene_tree_list.push((name, newick));
        }
    } else if let Ok(list) = gene_trees.extract::<Vec<String>>(py) {
        // List of Newick strings or file paths
        for (idx, newick_or_path) in list.iter().enumerate() {
            let (name, newick) = if PathBuf::from(newick_or_path).exists() {
                let path = PathBuf::from(newick_or_path);
                let name = path.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or(&format!("family_{}", idx))
                    .to_string();
                let newick = fs::read_to_string(&path)
                    .map_err(|e| PyValueError::new_err(format!(
                        "Failed to read gene tree file '{}': {}", newick_or_path, e
                    )))?;
                (name, newick)
            } else {
                (format!("family_{}", idx), newick_or_path.clone())
            };
            gene_tree_list.push((name, newick));
        }
    } else if let Ok(single) = gene_trees.extract::<String>(py) {
        // Single Newick string or file path
        let (name, newick) = if PathBuf::from(&single).exists() {
            let path = PathBuf::from(&single);
            let name = path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("family_0")
                .to_string();
            let newick = fs::read_to_string(&path)
                .map_err(|e| PyValueError::new_err(format!(
                    "Failed to read gene tree file '{}': {}", single, e
                )))?;
            (name, newick)
        } else {
            ("family_0".to_string(), single)
        };
        gene_tree_list.push((name, newick));
    } else {
        return Err(PyValueError::new_err(
            "gene_trees must be a Newick string, file path, list of strings/paths, or dict mapping names to strings/paths"
        ));
    }

    if gene_tree_list.is_empty() {
        return Err(PyValueError::new_err("No gene trees provided"));
    }

    // Validate inputs before running ALERax
    validate_inputs(&species_tree_obj.tree, &gene_tree_list)
        .map_err(|e| PyValueError::new_err(format!("Input validation failed: {}", e)))?;

    // Create output directory
    let (output_path, _temp_dir) = if let Some(dir) = output_dir {
        (PathBuf::from(dir), None)
    } else {
        let temp_dir = TempDir::new()
            .map_err(|e| PyValueError::new_err(format!("Failed to create temp directory: {}", e)))?;
        let path = temp_dir.path().to_path_buf();
        (path, Some(temp_dir))
    };

    // Create input directory for ALERax files
    let input_dir = output_path.join("input");
    fs::create_dir_all(&input_dir)
        .map_err(|e| PyValueError::new_err(format!("Failed to create input directory: {}", e)))?;

    let species_tree_path = input_dir.join("species_tree.newick");
    let families_file_path = input_dir.join("families.txt");

    // Prepare gene families
    let families: Vec<GeneFamily> = gene_tree_list.iter().map(|(name, newick)| {
        GeneFamily {
            name: name.clone(),
            gene_tree_path: input_dir.join(format!("{}.newick", name)),
            gene_tree_newick: newick.clone(),
        }
    }).collect();

    let alerax_output_dir = output_path.join("alerax_output");

    // Create ALERax config
    let config = AleRaxConfig {
        species_tree_path,
        families,
        families_file_path,
        output_dir: alerax_output_dir,
        num_samples,
        model_parametrization: model_type,
        seed,
        alerax_path,
    };

    // Run ALERax
    let results = run_alerax(config, &species_tree_newick)
        .map_err(|e| PyValueError::new_err(format!("ALERax execution failed: {}", e)))?;

    // Convert to Python types
    let mut py_results = HashMap::new();
    for (family_name, result) in results {
        // Wrap species tree in Arc so all gene trees in this family share the same allocation.
        // Take the species tree from the first reconciled tree; all trees in a family share
        // the same species tree.
        let mut rec_trees_iter = result.reconciled_trees.into_iter();
        let first_rec_tree = rec_trees_iter.next();
        let (shared_species_tree, first_py_gene_tree) = match first_rec_tree {
            Some(rt) => {
                let arc_st = Arc::new(rt.species_tree);
                let py_gt = PyGeneTree {
                    gene_tree: rt.gene_tree,
                    species_tree: Arc::clone(&arc_st),
                    node_mapping: rt.node_mapping,
                    event_mapping: rt.event_mapping,
                };
                (arc_st, Some(py_gt))
            }
            None => {
                // No reconciled trees — produce empty result
                (Arc::new(species_tree_obj.tree.as_ref().clone()), None)
            }
        };
        let mut py_gene_trees: Vec<PyGeneTree> = first_py_gene_tree.into_iter().collect();
        py_gene_trees.extend(rec_trees_iter.map(|rec_tree| PyGeneTree {
            gene_tree: rec_tree.gene_tree,
            species_tree: Arc::clone(&shared_species_tree),
            node_mapping: rec_tree.node_mapping,
            event_mapping: rec_tree.event_mapping,
        }));

        // Convert statistics
        let mean_event_counts = PyEventCounts {
            speciations: result.statistics.mean_event_counts.speciations,
            speciation_losses: result.statistics.mean_event_counts.speciation_losses,
            duplications: result.statistics.mean_event_counts.duplications,
            duplication_losses: result.statistics.mean_event_counts.duplication_losses,
            transfers: result.statistics.mean_event_counts.transfers,
            transfer_losses: result.statistics.mean_event_counts.transfer_losses,
            losses: result.statistics.mean_event_counts.losses,
            leaves: result.statistics.mean_event_counts.leaves,
        };

        let events_per_species: HashMap<String, PyEventCounts> = result.statistics.events_per_species
            .into_iter()
            .map(|(species, counts)| {
                (species, PyEventCounts {
                    speciations: counts.speciations,
                    speciation_losses: counts.speciation_losses,
                    duplications: counts.duplications,
                    duplication_losses: counts.duplication_losses,
                    transfers: counts.transfers,
                    transfer_losses: counts.transfer_losses,
                    losses: counts.losses,
                    leaves: counts.leaves,
                })
            })
            .collect();

        let statistics = PyReconciliationStatistics {
            mean_event_counts,
            mean_transfers: result.statistics.mean_transfers,
            events_per_species,
        };

        py_results.insert(family_name, PyAleRaxResult {
            gene_trees: py_gene_trees,
            duplication_rate: result.duplication_rate,
            loss_rate: result.loss_rate,
            transfer_rate: result.transfer_rate,
            likelihood: result.likelihood,
            statistics,
        });
    }

    // Clean up temp directory if not keeping output
    if !keep_output {
        drop(_temp_dir);
    }

    Ok(py_results)
}

/// Python module for rustree.
#[pymodule]
fn rustree(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(simulate_species_tree, m)?)?;
    m.add_function(wrap_pyfunction!(parse_species_tree, m)?)?;
    m.add_function(wrap_pyfunction!(parse_recphyloxml, m)?)?;
    m.add_function(wrap_pyfunction!(reconcile_with_alerax, m)?)?;
    m.add_class::<PySpeciesTree>()?;
    m.add_class::<PyGeneTree>()?;
    m.add_class::<PyEventCounts>()?;
    m.add_class::<PyReconciliationStatistics>()?;
    m.add_class::<PyAleRaxResult>()?;
    Ok(())
}
