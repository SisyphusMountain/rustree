//! ML/training tensor construction for Python bindings.
//!
//! Contains training sample creation, tensor building, GCN normalization,
//! batch collation, and inference batch construction.

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::sync::Arc;

use crate::node::{Event, RecTree, remap_gene_tree_indices};
use crate::sampling::extract_induced_subtree;
use crate::bd::simulate_bd_tree_bwd;
use crate::sampling::extract_extant_subtree;
use crate::dtl::simulate_dtl_per_species;

use super::{PyGeneTree, validate_dtl_rates};


/// Create a training sample from species tree, pruned gene tree, and reconciliation XML.
///
/// This replaces the Python `create_sample()` function with a fast Rust implementation.
/// The Zombi-format XML (containing `<recGeneTree>` without `<spTree>`) is parsed for
/// event annotations. The pruned gene tree from Newick provides the observable topology.
///
/// # Arguments
/// * `sp_newick_path` - Path to species tree Newick file
/// * `g_newick_path` - Path to pruned gene tree Newick file
/// * `xml_path` - Path to reconciliation XML file (Zombi format)
///
/// # Returns
/// A Python dict with the same keys as the Python `create_sample()`:
/// species_names, g_root_name, gene_names, g_neighbors, g_leaves_names,
/// sp_children, sp_parents, true_states, true_events, true_root,
/// nb_sp_leaves, nb_g_leaves, sp_tree_path, g_tree_path, xml_g_tree_path
#[pyfunction]
pub fn create_training_sample(
    py: Python,
    sp_newick_path: &str,
    g_newick_path: &str,
    xml_path: &str,
) -> PyResult<PyObject> {
    use crate::newick::newick::parse_newick;
    use crate::node::TraversalOrder;
    use crate::io::recphyloxml::parse_xml_gene_annotations_file;
    use pyo3::types::{PyDict, PyList};

    // 1. Parse species tree
    let sp_newick = fs::read_to_string(sp_newick_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to read species tree: {}", e)))?;
    let mut sp_nodes = parse_newick(&sp_newick)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse species Newick: {}", e)))?;
    let sp_root = sp_nodes.pop()
        .ok_or_else(|| PyValueError::new_err("No tree found in species Newick"))?;
    let sp_tree = sp_root.to_flat_tree();

    // 2. Parse pruned gene tree
    let g_newick = fs::read_to_string(g_newick_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to read gene tree: {}", e)))?;
    let mut g_nodes = parse_newick(&g_newick)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse gene Newick: {}", e)))?;
    let g_root = g_nodes.pop()
        .ok_or_else(|| PyValueError::new_err("No tree found in gene Newick"))?;
    let g_tree = g_root.to_flat_tree();

    // 3. Parse XML annotations
    let annotations = parse_xml_gene_annotations_file(xml_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse XML: {}", e)))?;

    // 4. Build species names (preorder)
    let sp_preorder: Vec<usize> = sp_tree.iter_indices(TraversalOrder::PreOrder).collect();
    let species_names: Vec<&str> = sp_preorder.iter()
        .map(|&i| sp_tree.nodes[i].name.as_str())
        .collect();
    let sp_name_set: std::collections::HashSet<&str> = species_names.iter().copied().collect();

    // 5. Build species children/parents dicts
    let sp_children = PyDict::new(py);
    let sp_parents = PyDict::new(py);
    for &idx in &sp_preorder {
        let node = &sp_tree.nodes[idx];
        let name = &node.name;
        let mut children_list = Vec::new();
        if let Some(left) = node.left_child {
            children_list.push(sp_tree.nodes[left].name.as_str());
        }
        if let Some(right) = node.right_child {
            children_list.push(sp_tree.nodes[right].name.as_str());
        }
        if !children_list.is_empty() {
            sp_children.set_item(name, children_list)?;
        }
        if let Some(parent) = node.parent {
            sp_parents.set_item(name, vec![sp_tree.nodes[parent].name.as_str()])?;
        }
    }

    // 6. Build gene names (preorder), separate root
    let g_preorder: Vec<usize> = g_tree.iter_indices(TraversalOrder::PreOrder).collect();
    let g_root_name = &g_tree.nodes[g_tree.root].name;
    let gene_names: Vec<&str> = g_preorder.iter()
        .filter(|&&i| i != g_tree.root)
        .map(|&i| g_tree.nodes[i].name.as_str())
        .collect();

    // 7. Build gene leaves
    let g_leaves_names: Vec<&str> = g_preorder.iter()
        .filter(|&&i| {
            let n = &g_tree.nodes[i];
            n.left_child.is_none() && n.right_child.is_none()
        })
        .map(|&i| g_tree.nodes[i].name.as_str())
        .collect();
    let nb_g_leaves = g_leaves_names.len();

    // 8. Build gene undirected adjacency (including root)
    let g_neighbors = PyDict::new(py);
    for &idx in &g_preorder {
        let node = &g_tree.nodes[idx];
        let mut neighbors = Vec::new();
        if let Some(parent) = node.parent {
            neighbors.push(g_tree.nodes[parent].name.as_str());
        }
        if let Some(left) = node.left_child {
            neighbors.push(g_tree.nodes[left].name.as_str());
        }
        if let Some(right) = node.right_child {
            neighbors.push(g_tree.nodes[right].name.as_str());
        }
        g_neighbors.set_item(&node.name, neighbors)?;
    }

    // 9. Build true_root (children of gene tree root)
    let g_root_node = &g_tree.nodes[g_tree.root];
    let mut root_children = Vec::new();
    if let Some(left) = g_root_node.left_child {
        root_children.push(g_tree.nodes[left].name.as_str());
    }
    if let Some(right) = g_root_node.right_child {
        root_children.push(g_tree.nodes[right].name.as_str());
    }

    // 10. Build true_states and true_events from XML annotations
    // Event mapping: speciation=0, duplication=1, T=2, P=3
    let true_states = PyDict::new(py);
    let true_events = PyDict::new(py);
    let all_g_names: Vec<&str> = g_preorder.iter()
        .map(|&i| g_tree.nodes[i].name.as_str())
        .collect();
    for &g_name in &all_g_names {
        let (sp_loc, event) = annotations.get(g_name)
            .ok_or_else(|| PyValueError::new_err(
                format!("Gene node '{}' not found in XML annotations", g_name)
            ))?;
        if !sp_name_set.contains(sp_loc.as_str()) {
            return Err(PyValueError::new_err(
                format!("Gene node '{}' maps to unknown species '{}'", g_name, sp_loc)
            ));
        }
        true_states.set_item(g_name, sp_loc)?;
        let event_idx: i32 = match event {
            Event::Speciation => 0,
            Event::Duplication => 1,
            Event::Transfer => 2,
            Event::Leaf => 3,
            Event::Loss => 3,  // shouldn't appear in pruned tree, but map safely
        };
        true_events.set_item(g_name, event_idx)?;
    }

    // 11. Count species leaves
    let nb_sp_leaves = sp_preorder.iter()
        .filter(|&&i| {
            let n = &sp_tree.nodes[i];
            n.left_child.is_none() && n.right_child.is_none()
        })
        .count();

    // 12. Build result dict
    let result = PyDict::new(py);
    result.set_item("species_names", PyList::new(py, &species_names)?)?;
    result.set_item("g_root_name", g_root_name)?;
    result.set_item("_gene_root_name", g_root_name)?;
    result.set_item("gene_names", PyList::new(py, &gene_names)?)?;
    result.set_item("g_neighbors", g_neighbors)?;
    result.set_item("g_leaves_names", PyList::new(py, &g_leaves_names)?)?;
    result.set_item("sp_children", sp_children)?;
    result.set_item("sp_parents", sp_parents)?;
    result.set_item("true_states", true_states)?;
    result.set_item("true_events", true_events)?;
    result.set_item("true_root", PyList::new(py, &root_children)?)?;
    result.set_item("nb_sp_leaves", nb_sp_leaves)?;
    result.set_item("nb_g_leaves", nb_g_leaves)?;
    result.set_item("sp_tree_path", sp_newick_path)?;
    result.set_item("g_tree_path", g_newick_path)?;
    result.set_item("xml_g_tree_path", xml_path)?;

    Ok(result.into())
}

/// Create a training sample dict directly from simulated tree objects.
///
/// Equivalent to `create_training_sample` but takes a `PyGeneTree` (after
/// `sample_extant()`) instead of file paths.  All reconciliation information
/// (node mapping, event types) is read from the internal `RecTree`, so no
/// intermediate files or pandas DataFrames are needed.
///
/// # Arguments
/// * `gene_tree` – a `PyGeneTree` returned by `sample_extant()`.
///
/// # Returns
/// A Python dict with the same keys as `create_training_sample` (minus the
/// three `*_path` keys which are set to empty strings).
#[cfg(feature = "python")]
#[pyfunction]
pub fn create_training_sample_from_sim(
    py: Python,
    gene_tree: &PyGeneTree,
) -> PyResult<PyObject> {
    use crate::node::TraversalOrder;
    use pyo3::types::{PyDict, PyList};

    let rt = &gene_tree.rec_tree;
    let sp_tree = &*rt.species_tree;
    let g_tree = &rt.gene_tree;

    // 1. Species names (preorder)
    let sp_preorder: Vec<usize> = sp_tree.iter_indices(TraversalOrder::PreOrder).collect();
    let species_names: Vec<&str> = sp_preorder.iter()
        .map(|&i| sp_tree.nodes[i].name.as_str())
        .collect();
    let sp_name_set: std::collections::HashSet<&str> = species_names.iter().copied().collect();

    // 2. Species children / parents
    let sp_children = PyDict::new(py);
    let sp_parents = PyDict::new(py);
    for &idx in &sp_preorder {
        let node = &sp_tree.nodes[idx];
        let mut children_list = Vec::new();
        if let Some(left) = node.left_child {
            children_list.push(sp_tree.nodes[left].name.as_str());
        }
        if let Some(right) = node.right_child {
            children_list.push(sp_tree.nodes[right].name.as_str());
        }
        if !children_list.is_empty() {
            sp_children.set_item(&node.name, children_list)?;
        }
        if let Some(parent) = node.parent {
            sp_parents.set_item(&node.name, vec![sp_tree.nodes[parent].name.as_str()])?;
        }
    }

    // 3. Gene names (preorder), root separate
    let g_preorder: Vec<usize> = g_tree.iter_indices(TraversalOrder::PreOrder).collect();
    let g_root_name = &g_tree.nodes[g_tree.root].name;
    let gene_names: Vec<&str> = g_preorder.iter()
        .filter(|&&i| i != g_tree.root)
        .map(|&i| g_tree.nodes[i].name.as_str())
        .collect();

    // 4. Gene leaves
    let g_leaves_names: Vec<&str> = g_preorder.iter()
        .filter(|&&i| {
            let n = &g_tree.nodes[i];
            n.left_child.is_none() && n.right_child.is_none()
        })
        .map(|&i| g_tree.nodes[i].name.as_str())
        .collect();
    let nb_g_leaves = g_leaves_names.len();

    // 5. Gene undirected adjacency
    let g_neighbors = PyDict::new(py);
    for &idx in &g_preorder {
        let node = &g_tree.nodes[idx];
        let mut neighbors = Vec::new();
        if let Some(parent) = node.parent {
            neighbors.push(g_tree.nodes[parent].name.as_str());
        }
        if let Some(left) = node.left_child {
            neighbors.push(g_tree.nodes[left].name.as_str());
        }
        if let Some(right) = node.right_child {
            neighbors.push(g_tree.nodes[right].name.as_str());
        }
        g_neighbors.set_item(&node.name, neighbors)?;
    }

    // 6. true_root (children of gene root)
    let g_root_node = &g_tree.nodes[g_tree.root];
    let mut root_children = Vec::new();
    if let Some(left) = g_root_node.left_child {
        root_children.push(g_tree.nodes[left].name.as_str());
    }
    if let Some(right) = g_root_node.right_child {
        root_children.push(g_tree.nodes[right].name.as_str());
    }

    // 7. true_states and true_events from RecTree mappings
    let true_states = PyDict::new(py);
    let true_events = PyDict::new(py);
    for &idx in &g_preorder {
        let g_name = g_tree.nodes[idx].name.as_str();
        // Species mapping
        let sp_idx = rt.node_mapping[idx]
            .ok_or_else(|| PyValueError::new_err(
                format!("Gene node '{}' has no species mapping", g_name)
            ))?;
        let sp_name = sp_tree.nodes[sp_idx].name.as_str();
        if !sp_name_set.contains(sp_name) {
            return Err(PyValueError::new_err(
                format!("Gene node '{}' maps to unknown species '{}'", g_name, sp_name)
            ));
        }
        true_states.set_item(g_name, sp_name)?;
        // Event mapping
        let event_idx: i32 = match &rt.event_mapping[idx] {
            Event::Speciation => 0,
            Event::Duplication => 1,
            Event::Transfer => 2,
            Event::Leaf => 3,
            Event::Loss => 3,
        };
        true_events.set_item(g_name, event_idx)?;
    }

    // 8. Count species leaves
    let nb_sp_leaves = sp_preorder.iter()
        .filter(|&&i| {
            let n = &sp_tree.nodes[i];
            n.left_child.is_none() && n.right_child.is_none()
        })
        .count();

    // 9. Build result dict
    let result = PyDict::new(py);
    result.set_item("species_names", PyList::new(py, &species_names)?)?;
    result.set_item("g_root_name", g_root_name)?;
    result.set_item("_gene_root_name", g_root_name)?;
    result.set_item("gene_names", PyList::new(py, &gene_names)?)?;
    result.set_item("g_neighbors", g_neighbors)?;
    result.set_item("g_leaves_names", PyList::new(py, &g_leaves_names)?)?;
    result.set_item("sp_children", sp_children)?;
    result.set_item("sp_parents", sp_parents)?;
    result.set_item("true_states", true_states)?;
    result.set_item("true_events", true_events)?;
    result.set_item("true_root", PyList::new(py, &root_children)?)?;
    result.set_item("nb_sp_leaves", nb_sp_leaves)?;
    result.set_item("nb_g_leaves", nb_g_leaves)?;
    result.set_item("sp_tree_path", "")?;
    result.set_item("g_tree_path", "")?;
    result.set_item("xml_g_tree_path", "")?;

    Ok(result.into())
}

/// Build all tensors needed for a training sample in one Rust call.
///
/// Replaces the Python functions: unroot_tree, sample_gene_coloring,
/// _build_species_tensors, _build_gene_edges, _build_gene_event_types,
/// and the masking/assembly logic in get_sample.
///
/// Returns a dict of numpy arrays ready to be wrapped into HeteroData.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (base_sample, seed, force_mask_last_added, predict_root_position, sample_order = "random"))]
pub fn build_training_tensors(
    py: Python,
    base_sample: &Bound<'_, pyo3::types::PyDict>,
    seed: u64,
    force_mask_last_added: bool,
    predict_root_position: bool,
    sample_order: &str,
) -> PyResult<PyObject> {
    use numpy::PyArray1;
    use pyo3::types::PyDict;
    use rand::seq::SliceRandom;
    use rand::Rng;
    use std::collections::{HashMap, HashSet, VecDeque};

    // ---- Extract fields from base_sample dict ----
    let g_root_name: String = base_sample
        .get_item("g_root_name")?.ok_or_else(|| PyValueError::new_err("missing g_root_name"))?
        .extract::<String>()?;
    let nb_g_leaves: usize = base_sample
        .get_item("nb_g_leaves")?.ok_or_else(|| PyValueError::new_err("missing nb_g_leaves"))?
        .extract::<usize>()?;
    let species_names: Vec<String> = base_sample
        .get_item("species_names")?.ok_or_else(|| PyValueError::new_err("missing species_names"))?
        .extract::<Vec<String>>()?;
    let g_leaves_names: Vec<String> = base_sample
        .get_item("g_leaves_names")?.ok_or_else(|| PyValueError::new_err("missing g_leaves_names"))?
        .extract::<Vec<String>>()?;
    let true_root: Vec<String> = base_sample
        .get_item("true_root")?.ok_or_else(|| PyValueError::new_err("missing true_root"))?
        .extract::<Vec<String>>()?;

    // g_neighbors: Dict[str, List[str]]
    let g_neighbors: HashMap<String, Vec<String>> = base_sample
        .get_item("g_neighbors")?.ok_or_else(|| PyValueError::new_err("missing g_neighbors"))?
        .extract()?;

    // sp_children: Dict[str, List[str]]
    let sp_children: HashMap<String, Vec<String>> = base_sample
        .get_item("sp_children")?.ok_or_else(|| PyValueError::new_err("missing sp_children"))?
        .extract()?;

    // true_states: Dict[str, str]
    let true_states: HashMap<String, String> = base_sample
        .get_item("true_states")?.ok_or_else(|| PyValueError::new_err("missing true_states"))?
        .extract()?;

    // true_events: Dict[str, int]
    let true_events: HashMap<String, i32> = base_sample
        .get_item("true_events")?.ok_or_else(|| PyValueError::new_err("missing true_events"))?
        .extract()?;

    // ---- RNG setup ----
    let mut rng = StdRng::seed_from_u64(seed);

    // ---- Determine sampling parameters ----
    let max_to_sample = if !force_mask_last_added || predict_root_position {
        nb_g_leaves.saturating_sub(2)
    } else {
        nb_g_leaves.saturating_sub(1)
    };

    let nb_sampled = if predict_root_position {
        max_to_sample
    } else {
        rng.gen_range(0..=max_to_sample)
    };

    // ---- Unroot tree ----
    let mut g_neighbors_unrooted: HashMap<String, Vec<String>> = HashMap::new();
    for (k, v) in &g_neighbors {
        if k != &g_root_name {
            let filtered: Vec<String> = v.iter()
                .filter(|n| *n != &g_root_name)
                .cloned()
                .collect();
            g_neighbors_unrooted.insert(k.clone(), filtered);
        }
    }
    // Connect root's two children
    if let Some(root_nbrs) = g_neighbors.get(&g_root_name) {
        if root_nbrs.len() == 2 {
            let c1 = &root_nbrs[0];
            let c2 = &root_nbrs[1];
            g_neighbors_unrooted.entry(c1.clone()).or_default().push(c2.clone());
            g_neighbors_unrooted.entry(c2.clone()).or_default().push(c1.clone());
        }
    }

    // ---- Frontier sampling / node selection ----
    let include_root_node: bool;
    let selected: HashSet<String>;
    let frontier_names: Vec<String>;
    let last_added_name: Option<String>;

    if predict_root_position {
        include_root_node = false;
        selected = g_neighbors_unrooted.keys().cloned().collect();
        frontier_names = Vec::new();
        last_added_name = None;
    } else if nb_sampled == max_to_sample && !force_mask_last_added {
        include_root_node = true;
        selected = g_neighbors.keys()
            .filter(|k| *k != &g_root_name)
            .cloned().collect();
        frontier_names = vec![g_root_name.clone()];
        last_added_name = None;
    } else if nb_sampled == max_to_sample && force_mask_last_added {
        include_root_node = true;
        selected = g_neighbors.keys().cloned().collect();
        frontier_names = Vec::new();
        last_added_name = Some(g_root_name.clone());
    } else {
        include_root_node = false;
        // ---- Inline frontier sampling (sample_gene_coloring) ----
        let node_list: Vec<String> = g_neighbors_unrooted.keys().cloned().collect();
        let n = node_list.len();
        let name_to_idx: HashMap<&str, usize> = node_list.iter()
            .enumerate()
            .map(|(i, name)| (name.as_str(), i))
            .collect();

        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (u, vs) in &g_neighbors_unrooted {
            let ui = name_to_idx[u.as_str()];
            for v in vs {
                let vi = name_to_idx[v.as_str()];
                if ui != vi {
                    adj[ui].push(vi);
                }
            }
        }
        // Deduplicate adjacency
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }

        let degrees: Vec<usize> = adj.iter().map(|a| a.len()).collect();
        let leaves_idx: Vec<usize> = degrees.iter().enumerate()
            .filter(|(_, &d)| d == 1)
            .map(|(i, _)| i)
            .collect();
        let internal_idx: Vec<usize> = degrees.iter().enumerate()
            .filter(|(_, &d)| d == 3)
            .map(|(i, _)| i)
            .collect();

        // Leaf distance BFS for bottom_up ordering
        let leaf_dist: Option<Vec<i32>> = if sample_order == "bottom_up" {
            let mut dist = vec![-1i32; n];
            let mut queue = VecDeque::new();
            for &li in &leaves_idx {
                dist[li] = 0;
                queue.push_back(li);
            }
            while let Some(u) = queue.pop_front() {
                for &v in &adj[u] {
                    if dist[v] == -1 {
                        dist[v] = dist[u] + 1;
                        queue.push_back(v);
                    }
                }
            }
            Some(dist)
        } else {
            None
        };

        let mut mapped = vec![false; n];
        for &li in &leaves_idx {
            mapped[li] = true;
        }
        let internal_set: HashSet<usize> = internal_idx.iter().copied().collect();

        let mut mapped_cnt = vec![0u32; n];
        for &ii in &internal_idx {
            let cnt = adj[ii].iter().filter(|&&nb| mapped[nb]).count() as u32;
            mapped_cnt[ii] = cnt;
        }

        let mut frontier_idx: Vec<usize> = internal_idx.iter()
            .filter(|&&ii| mapped_cnt[ii] >= 2)
            .copied()
            .collect();
        let mut added_order: Vec<usize> = Vec::new();
        let mut last_added_i: Option<usize> = None;

        while added_order.len() < nb_sampled && !frontier_idx.is_empty() {
            let pick = if let Some(ref dist) = leaf_dist {
                let min_d = frontier_idx.iter().map(|&v| dist[v]).min().unwrap();
                let closest: Vec<usize> = frontier_idx.iter().enumerate()
                    .filter(|(_, &v)| dist[v] == min_d)
                    .map(|(i, _)| i)
                    .collect();
                *closest.choose(&mut rng).unwrap()
            } else {
                rng.gen_range(0..frontier_idx.len())
            };
            let v = frontier_idx.swap_remove(pick);
            mapped[v] = true;
            added_order.push(v);
            last_added_i = Some(v);

            for &u in &adj[v] {
                if !mapped[u] {
                    mapped_cnt[u] += 1;
                    if internal_set.contains(&u) && mapped_cnt[u] == 2 {
                        frontier_idx.push(u);
                    }
                }
            }
        }

        let mut sel = HashSet::new();
        for &li in &leaves_idx {
            sel.insert(node_list[li].clone());
        }
        for &ai in &added_order {
            sel.insert(node_list[ai].clone());
        }
        selected = sel;
        frontier_names = frontier_idx.iter().map(|&i| node_list[i].clone()).collect();
        last_added_name = last_added_i.map(|i| node_list[i].clone());
    };

    // ---- Species tensors ----
    let _n_sp = species_names.len();
    let mut sp_name_to_id: HashMap<&str, i32> = HashMap::new();
    let mut sp_name_to_idx: HashMap<&str, usize> = HashMap::new();
    for (i, name) in species_names.iter().enumerate() {
        sp_name_to_id.insert(name.as_str(), (i + 1) as i32);
        sp_name_to_idx.insert(name.as_str(), i);
    }
    let x_sp: Vec<i32> = species_names.iter()
        .map(|name: &String| sp_name_to_id[name.as_str()])
        .collect();

    let mut sp_child_src: Vec<i32> = Vec::new();
    let mut sp_child_dst: Vec<i32> = Vec::new();
    let mut sp_parent_src: Vec<i32> = Vec::new();
    let mut sp_parent_dst: Vec<i32> = Vec::new();
    for (p, children) in &sp_children {
        let pi = sp_name_to_idx[p.as_str()] as i32;
        for c in children {
            let ci = sp_name_to_idx[c.as_str()] as i32;
            sp_child_src.push(pi);
            sp_child_dst.push(ci);
            sp_parent_src.push(ci);
            sp_parent_dst.push(pi);
        }
    }

    // ---- Gene node ordering ----
    let active_neighbors = if include_root_node { &g_neighbors } else { &g_neighbors_unrooted };
    let mut g_names: Vec<String> = active_neighbors.keys().cloned().collect();
    g_names.sort();
    let n_gene = g_names.len();
    let g_name_to_idx: HashMap<&str, usize> = g_names.iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();

    // ---- Gene feature tensors ----
    let leaf_set: HashSet<&str> = g_leaves_names.iter().map(|s: &String| s.as_str()).collect();
    let frontier_set: HashSet<&str> = frontier_names.iter().map(|s| s.as_str()).collect();

    let mut x_gene = vec![0i32; n_gene];
    let mut g_true_sp = vec![0i32; n_gene];
    let mut event_true = vec![0i32; n_gene];
    let mut event_input = vec![0i32; n_gene];
    let mut frontier_mask = vec![0u8; n_gene]; // bool as u8
    let mut is_leaf = vec![0u8; n_gene];
    let mut mask_label_node = vec![0u8; n_gene];

    for (i, name) in g_names.iter().enumerate() {
        // Species mapping
        let sp_name = true_states.get(name)
            .ok_or_else(|| PyValueError::new_err(format!("Gene node '{}' not in true_states", name)))?;
        let sp_id = sp_name_to_id.get(sp_name.as_str())
            .ok_or_else(|| PyValueError::new_err(format!("Species '{}' not found", sp_name)))?;
        g_true_sp[i] = *sp_id;

        // Input features: masked for unselected nodes
        if selected.contains(name) {
            x_gene[i] = *sp_id;
        }
        // else x_gene[i] stays 0

        // Event types
        if leaf_set.contains(name.as_str()) {
            event_true[i] = 3; // P/leaf
        } else {
            event_true[i] = *true_events.get(name)
                .ok_or_else(|| PyValueError::new_err(format!("Gene node '{}' not in true_events", name)))?;
        }

        // Event input: masked for unselected nodes
        if selected.contains(name) {
            event_input[i] = event_true[i];
        } else {
            event_input[i] = 4; // mask
        }

        // Frontier mask
        if frontier_set.contains(name.as_str()) {
            frontier_mask[i] = 1;
        }

        // Is leaf
        if leaf_set.contains(name.as_str()) {
            is_leaf[i] = 1;
        }
    }

    // Handle last_added masking
    let last_idx: i32;
    let mask_last_added_event: bool;
    if let Some(ref la) = last_added_name {
        if let Some(&idx) = g_name_to_idx.get(la.as_str()) {
            last_idx = idx as i32;
            if force_mask_last_added {
                event_input[idx] = 4; // mask
                mask_last_added_event = true;
                mask_label_node[idx] = 1;
            } else {
                mask_last_added_event = false;
            }
        } else {
            last_idx = -1;
            mask_last_added_event = false;
        }
    } else {
        last_idx = -1;
        mask_last_added_event = false;
    }

    // ---- Gene edges ----
    let mut undirected_edges: HashSet<(i32, i32)> = HashSet::new();
    for (u, vs) in active_neighbors {
        let ui = g_name_to_idx[u.as_str()] as i32;
        for v in vs {
            if let Some(&vi_usize) = g_name_to_idx.get(v.as_str()) {
                let vi = vi_usize as i32;
                if ui != vi {
                    let e = if ui < vi { (ui, vi) } else { (vi, ui) };
                    undirected_edges.insert(e);
                }
            }
        }
    }

    let mut directed: Vec<(i32, i32)> = undirected_edges.into_iter().collect();
    directed.sort();

    let n_edges = directed.len();
    let mut g_edge_src = Vec::with_capacity(n_edges * 2);
    let mut g_edge_dst = Vec::with_capacity(n_edges * 2);
    let mut g_dir_src = Vec::with_capacity(n_edges);
    let mut g_dir_dst = Vec::with_capacity(n_edges);
    for &(a, b) in &directed {
        // Directed (canonical)
        g_dir_src.push(a);
        g_dir_dst.push(b);
        // Undirected (both directions)
        g_edge_src.push(a);
        g_edge_dst.push(b);
        g_edge_src.push(b);
        g_edge_dst.push(a);
    }

    // Root edge target
    let root_pair_0 = g_name_to_idx.get(true_root[0].as_str()).map(|&v| v as i32).unwrap_or(-1);
    let root_pair_1 = g_name_to_idx.get(true_root[1].as_str()).map(|&v| v as i32).unwrap_or(-1);
    let root_edge_target: i32 = if !include_root_node {
        let pair_fwd = (root_pair_0.min(root_pair_1), root_pair_0.max(root_pair_1));
        directed.iter().position(|&e| e == pair_fwd)
            .map(|i| i as i32)
            .unwrap_or(-1)
    } else {
        -1
    };

    // ---- Build result dict of numpy arrays ----
    // Use from_slice / into_pyarray so numpy fully owns the memory (no Rust
    // allocator callback on dealloc). This prevents cross-allocator corruption
    // when DataLoader workers tear down.
    let result = PyDict::new(py);
    result.set_item("x_sp", PyArray1::from_slice(py, &x_sp))?;
    result.set_item("x_gene", PyArray1::from_slice(py, &x_gene))?;
    result.set_item("g_true_sp", PyArray1::from_slice(py, &g_true_sp))?;
    result.set_item("event_true", PyArray1::from_slice(py, &event_true))?;
    result.set_item("event_input", PyArray1::from_slice(py, &event_input))?;
    result.set_item("frontier_mask", PyArray1::from_slice(py, &frontier_mask))?;
    result.set_item("is_leaf", PyArray1::from_slice(py, &is_leaf))?;
    result.set_item("mask_label_node", PyArray1::from_slice(py, &mask_label_node))?;
    result.set_item("last_added_index", last_idx)?;
    result.set_item("mask_last_added_event", mask_last_added_event)?;

    // Edge indices as 2×E arrays — build via ndarray then convert to numpy.
    let edge_2d = |src: &[i32], dst: &[i32], name: &str| -> PyResult<PyObject> {
        use numpy::ndarray::Array2;
        use numpy::ToPyArray;
        let n = src.len();
        let mut flat = Vec::with_capacity(2 * n);
        flat.extend_from_slice(src);
        flat.extend_from_slice(dst);
        let arr = Array2::from_shape_vec([2, n], flat)
            .map_err(|e| PyValueError::new_err(format!("{}: {}", name, e)))?;
        Ok(arr.to_pyarray(py).into_any().unbind())
    };
    result.set_item("sp_child_edge", edge_2d(&sp_child_src, &sp_child_dst, "sp_child_edge")?)?;
    result.set_item("sp_parent_edge", edge_2d(&sp_parent_src, &sp_parent_dst, "sp_parent_edge")?)?;
    result.set_item("g_edge", edge_2d(&g_edge_src, &g_edge_dst, "g_edge")?)?;
    result.set_item("g_dir_edge", edge_2d(&g_dir_src, &g_dir_dst, "g_dir_edge")?)?;
    result.set_item("root_edge_target", root_edge_target)?;

    Ok(result.into())
}

// ============================================================================
// Batch OTF generation: internal helpers and build_otf_batch
// ============================================================================

/// Internal representation of a training sample (mirrors Python base_sample dict).
struct RustBaseSample {
    g_root_name: String,
    nb_g_leaves: usize,
    species_names: Vec<String>,
    g_leaves_names: Vec<String>,
    true_root: Vec<String>,
    g_neighbors: HashMap<String, Vec<String>>,
    sp_children: HashMap<String, Vec<String>>,
    true_states: HashMap<String, String>,
    true_events: HashMap<String, i32>,
}

/// Per-sample task tensors (before collation).
struct TaskTensors {
    x_gene: Vec<i32>,
    g_true_sp: Vec<i32>,   // 1-based species IDs
    event_true: Vec<i32>,
    event_input: Vec<i32>,
    frontier_mask: Vec<u8>,
    is_leaf: Vec<u8>,
    mask_label_node: Vec<u8>,
    g_edge_src: Vec<i32>,
    g_edge_dst: Vec<i32>,
    g_dir_src: Vec<i32>,
    g_dir_dst: Vec<i32>,
    root_edge_target: i32,
}

/// GCN normalization: add self-loops, compute D^{-1/2} A D^{-1/2} edge weights.
///
/// Matches PyG's `gcn_norm(edge_index, None, num_nodes, improved=False,
/// add_self_loops=True, flow="source_to_target")`.
///
/// Returns (new_src, new_dst, edge_weights) with self-loops included.
fn gcn_norm_internal(src: &[i32], dst: &[i32], num_nodes: usize) -> (Vec<i32>, Vec<i32>, Vec<f32>) {
    let n_orig = src.len();
    let n_total = n_orig + num_nodes; // original edges + self-loops

    let mut new_src = Vec::with_capacity(n_total);
    let mut new_dst = Vec::with_capacity(n_total);
    new_src.extend_from_slice(src);
    new_dst.extend_from_slice(dst);
    for i in 0..num_nodes {
        new_src.push(i as i32);
        new_dst.push(i as i32);
    }

    // Degree = number of edges targeting each node (flow="source_to_target")
    let mut deg = vec![0u32; num_nodes];
    for &d in &new_dst {
        deg[d as usize] += 1;
    }

    // D^{-1/2}
    let inv_sqrt: Vec<f32> = deg.iter().map(|&d| {
        if d == 0 { 0.0f32 } else { (d as f32).powf(-0.5) }
    }).collect();

    // Edge weights: inv_sqrt[src] * inv_sqrt[dst]
    let weights: Vec<f32> = new_src.iter().zip(new_dst.iter())
        .map(|(&s, &d)| inv_sqrt[s as usize] * inv_sqrt[d as usize])
        .collect();

    (new_src, new_dst, weights)
}

/// Compute varlen attention metadata from sizes array.
/// Returns (prefix_sums_with_leading_zero, max_size).
fn compute_varlen_metadata(sizes: &[i32]) -> (Vec<i32>, i32) {
    let mut cu = Vec::with_capacity(sizes.len() + 1);
    cu.push(0);
    let mut cumsum: i32 = 0;
    let mut max_val: i32 = 0;
    for &s in sizes {
        cumsum += s;
        cu.push(cumsum);
        if s > max_val { max_val = s; }
    }
    (cu, max_val)
}

/// Collated tensors for one task across a full batch.
struct CollatedTask {
    gene_x: Vec<i32>,
    gene_y: Vec<i32>,
    event: Vec<i32>,
    event_true: Vec<i32>,
    frontier_mask: Vec<u8>,
    is_leaf: Vec<u8>,
    mask_label_node: Vec<u8>,
    // GCN-normalized gene edges (includes self-loops)
    gene_ei_src: Vec<i32>,
    gene_ei_dst: Vec<i32>,
    gene_ew: Vec<f32>,
    g_dir_src: Vec<i32>,
    g_dir_dst: Vec<i32>,
    true_root_index: Vec<i32>,
    g_sizes: Vec<i32>,
    g_ptr: Vec<i32>,
    g_batch: Vec<i32>,
    root_edges_ptr: Vec<i32>,
    // Varlen attention metadata
    cu_g: Vec<i32>,
    max_g: i32,
}

/// Extract base sample data from a RecTree (no Python boundary).
fn extract_base_sample_internal(rt: &RecTree) -> Result<RustBaseSample, String> {
    use crate::node::TraversalOrder;

    let sp_tree = &*rt.species_tree;
    let g_tree = &rt.gene_tree;

    let sp_preorder: Vec<usize> = sp_tree.iter_indices(TraversalOrder::PreOrder).collect();
    let species_names: Vec<String> = sp_preorder.iter()
        .map(|&i| sp_tree.nodes[i].name.clone())
        .collect();
    let sp_name_set: std::collections::HashSet<&str> = species_names.iter().map(|s| s.as_str()).collect();

    let mut sp_children: HashMap<String, Vec<String>> = HashMap::new();
    for &idx in &sp_preorder {
        let node = &sp_tree.nodes[idx];
        let mut children = Vec::new();
        if let Some(left) = node.left_child {
            children.push(sp_tree.nodes[left].name.clone());
        }
        if let Some(right) = node.right_child {
            children.push(sp_tree.nodes[right].name.clone());
        }
        if !children.is_empty() {
            sp_children.insert(node.name.clone(), children);
        }
    }

    let g_preorder: Vec<usize> = g_tree.iter_indices(TraversalOrder::PreOrder).collect();
    let g_root_name = g_tree.nodes[g_tree.root].name.clone();

    let g_leaves_names: Vec<String> = g_preorder.iter()
        .filter(|&&i| g_tree.nodes[i].left_child.is_none() && g_tree.nodes[i].right_child.is_none())
        .map(|&i| g_tree.nodes[i].name.clone())
        .collect();
    let nb_g_leaves = g_leaves_names.len();

    let mut g_neighbors: HashMap<String, Vec<String>> = HashMap::new();
    for &idx in &g_preorder {
        let node = &g_tree.nodes[idx];
        let mut neighbors = Vec::new();
        if let Some(parent) = node.parent {
            neighbors.push(g_tree.nodes[parent].name.clone());
        }
        if let Some(left) = node.left_child {
            neighbors.push(g_tree.nodes[left].name.clone());
        }
        if let Some(right) = node.right_child {
            neighbors.push(g_tree.nodes[right].name.clone());
        }
        g_neighbors.insert(node.name.clone(), neighbors);
    }

    let g_root_node = &g_tree.nodes[g_tree.root];
    let mut true_root = Vec::new();
    if let Some(left) = g_root_node.left_child {
        true_root.push(g_tree.nodes[left].name.clone());
    }
    if let Some(right) = g_root_node.right_child {
        true_root.push(g_tree.nodes[right].name.clone());
    }

    let mut true_states: HashMap<String, String> = HashMap::new();
    let mut true_events: HashMap<String, i32> = HashMap::new();
    for &idx in &g_preorder {
        let g_name = g_tree.nodes[idx].name.clone();
        let sp_idx = rt.node_mapping[idx]
            .ok_or_else(|| format!("Gene node '{}' has no species mapping", g_name))?;
        let sp_name = sp_tree.nodes[sp_idx].name.clone();
        if !sp_name_set.contains(sp_name.as_str()) {
            return Err(format!("Gene node '{}' maps to unknown species '{}'", g_name, sp_name));
        }
        true_states.insert(g_name.clone(), sp_name);
        let event_idx: i32 = match &rt.event_mapping[idx] {
            Event::Speciation => 0,
            Event::Duplication => 1,
            Event::Transfer => 2,
            Event::Leaf => 3,
            Event::Loss => 3,
        };
        true_events.insert(g_name, event_idx);
    }

    Ok(RustBaseSample {
        g_root_name, nb_g_leaves, species_names, g_leaves_names,
        true_root, g_neighbors, sp_children, true_states, true_events,
    })
}

/// Build task tensors for one sample (internal; mirrors build_training_tensors logic).
/// sp_name_to_id maps species name (owned) → 1-based ID.
fn build_task_tensors_internal(
    sample: &RustBaseSample,
    sp_name_to_id: &HashMap<String, i32>,
    rng: &mut StdRng,
    force_mask_last_added: bool,
    predict_root_position: bool,
    sample_order: &str,
) -> Result<TaskTensors, String> {
    use rand::Rng;
    use rand::seq::SliceRandom;
    use std::collections::{HashSet, VecDeque};

    let g_root_name = &sample.g_root_name;
    let nb_g_leaves = sample.nb_g_leaves;
    let g_neighbors = &sample.g_neighbors;
    let g_leaves_names = &sample.g_leaves_names;
    let true_root = &sample.true_root;
    let true_states = &sample.true_states;
    let true_events = &sample.true_events;

    let max_to_sample = if !force_mask_last_added || predict_root_position {
        nb_g_leaves.saturating_sub(2)
    } else {
        nb_g_leaves.saturating_sub(1)
    };
    let nb_sampled = if predict_root_position {
        max_to_sample
    } else {
        rng.gen_range(0..=max_to_sample)
    };

    // Unroot tree: remove root node and reconnect its children
    let mut g_neighbors_unrooted: HashMap<String, Vec<String>> = HashMap::new();
    for (k, v) in g_neighbors {
        if k != g_root_name {
            let filtered: Vec<String> = v.iter()
                .filter(|n| n.as_str() != g_root_name.as_str())
                .cloned()
                .collect();
            g_neighbors_unrooted.insert(k.clone(), filtered);
        }
    }
    if let Some(root_nbrs) = g_neighbors.get(g_root_name) {
        // Connect all former root children to each other (handles bifurcating
        // and multifurcating roots alike).
        for i in 0..root_nbrs.len() {
            for j in (i + 1)..root_nbrs.len() {
                let ci = root_nbrs[i].clone();
                let cj = root_nbrs[j].clone();
                g_neighbors_unrooted.entry(ci.clone()).or_default().push(cj.clone());
                g_neighbors_unrooted.entry(cj.clone()).or_default().push(ci.clone());
            }
        }
    }
    // Safety: ensure every node referenced as a neighbor is also a key.
    // Collect missing keys first to avoid borrowing issues.
    let all_referenced: Vec<String> = g_neighbors_unrooted.values()
        .flat_map(|vs| vs.iter().cloned())
        .collect();
    for name in all_referenced {
        g_neighbors_unrooted.entry(name).or_default();
    }

    // Node selection
    let include_root_node: bool;
    let selected: HashSet<String>;
    let frontier_names: Vec<String>;
    let last_added_name: Option<String>;

    if predict_root_position {
        include_root_node = false;
        selected = g_neighbors_unrooted.keys().cloned().collect();
        frontier_names = Vec::new();
        last_added_name = None;
    } else if nb_sampled == max_to_sample && !force_mask_last_added {
        include_root_node = true;
        selected = g_neighbors.keys().filter(|k| k.as_str() != g_root_name.as_str()).cloned().collect();
        frontier_names = vec![g_root_name.clone()];
        last_added_name = None;
    } else if nb_sampled == max_to_sample && force_mask_last_added {
        include_root_node = true;
        selected = g_neighbors.keys().cloned().collect();
        frontier_names = Vec::new();
        last_added_name = Some(g_root_name.clone());
    } else {
        include_root_node = false;
        let node_list: Vec<String> = g_neighbors_unrooted.keys().cloned().collect();
        let n = node_list.len();
        let name_to_idx: HashMap<&str, usize> = node_list.iter().enumerate()
            .map(|(i, name)| (name.as_str(), i)).collect();

        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (u, vs) in &g_neighbors_unrooted {
            let ui = name_to_idx[u.as_str()];
            for v in vs {
                let vi = match name_to_idx.get(v.as_str()) {
                    Some(&idx) => idx,
                    None => {
                        eprintln!(
                            "[rustree] WARNING: neighbor '{}' of node '{}' not found \
                             in g_neighbors_unrooted keys — skipping edge",
                            v, u
                        );
                        continue;
                    }
                };
                if ui != vi { adj[ui].push(vi); }
            }
        }
        for a in &mut adj { a.sort_unstable(); a.dedup(); }

        let degrees: Vec<usize> = adj.iter().map(|a| a.len()).collect();
        let leaves_idx: Vec<usize> = degrees.iter().enumerate()
            .filter(|(_, &d)| d == 1).map(|(i, _)| i).collect();
        let internal_idx: Vec<usize> = degrees.iter().enumerate()
            .filter(|(_, &d)| d == 3).map(|(i, _)| i).collect();

        let leaf_dist: Option<Vec<i32>> = if sample_order == "bottom_up" {
            let mut dist = vec![-1i32; n];
            let mut queue = VecDeque::new();
            for &li in &leaves_idx { dist[li] = 0; queue.push_back(li); }
            while let Some(u) = queue.pop_front() {
                for &v in &adj[u] {
                    if dist[v] == -1 { dist[v] = dist[u] + 1; queue.push_back(v); }
                }
            }
            Some(dist)
        } else { None };

        let mut mapped = vec![false; n];
        for &li in &leaves_idx { mapped[li] = true; }
        let internal_set: HashSet<usize> = internal_idx.iter().copied().collect();

        let mut mapped_cnt = vec![0u32; n];
        for &ii in &internal_idx {
            mapped_cnt[ii] = adj[ii].iter().filter(|&&nb| mapped[nb]).count() as u32;
        }

        let mut frontier_idx: Vec<usize> = internal_idx.iter()
            .filter(|&&ii| mapped_cnt[ii] >= 2).copied().collect();
        let mut added_order: Vec<usize> = Vec::new();
        let mut last_added_i: Option<usize> = None;

        while added_order.len() < nb_sampled && !frontier_idx.is_empty() {
            let pick = if let Some(ref dist) = leaf_dist {
                let min_d = frontier_idx.iter().map(|&v| dist[v]).min().unwrap();
                let closest: Vec<usize> = frontier_idx.iter().enumerate()
                    .filter(|(_, &v)| dist[v] == min_d).map(|(i, _)| i).collect();
                *closest.choose(rng).unwrap()
            } else {
                rng.gen_range(0..frontier_idx.len())
            };
            let v = frontier_idx.swap_remove(pick);
            mapped[v] = true;
            added_order.push(v);
            last_added_i = Some(v);
            for &u in &adj[v] {
                if !mapped[u] {
                    mapped_cnt[u] += 1;
                    if internal_set.contains(&u) && mapped_cnt[u] == 2 {
                        frontier_idx.push(u);
                    }
                }
            }
        }

        let mut sel = HashSet::new();
        for &li in &leaves_idx { sel.insert(node_list[li].clone()); }
        for &ai in &added_order { sel.insert(node_list[ai].clone()); }
        selected = sel;
        frontier_names = frontier_idx.iter().map(|&i| node_list[i].clone()).collect();
        last_added_name = last_added_i.map(|i| node_list[i].clone());
    }

    let active_neighbors = if include_root_node { g_neighbors } else { &g_neighbors_unrooted };
    let mut g_names: Vec<String> = active_neighbors.keys().cloned().collect();
    g_names.sort();
    let n_gene = g_names.len();
    let g_name_to_idx: HashMap<&str, usize> = g_names.iter().enumerate()
        .map(|(i, name)| (name.as_str(), i)).collect();

    let leaf_set: HashSet<&str> = g_leaves_names.iter().map(|s| s.as_str()).collect();
    let frontier_set: HashSet<&str> = frontier_names.iter().map(|s| s.as_str()).collect();

    let mut x_gene = vec![0i32; n_gene];
    let mut g_true_sp = vec![0i32; n_gene];
    let mut event_true = vec![0i32; n_gene];
    let mut event_input = vec![0i32; n_gene];
    let mut frontier_mask = vec![0u8; n_gene];
    let mut is_leaf = vec![0u8; n_gene];
    let mut mask_label_node = vec![0u8; n_gene];

    for (i, name) in g_names.iter().enumerate() {
        let sp_name = true_states.get(name)
            .ok_or_else(|| format!("Gene node '{}' not in true_states", name))?;
        let sp_id = sp_name_to_id.get(sp_name.as_str())
            .ok_or_else(|| format!("Species '{}' not found in sp_name_to_id", sp_name))?;
        g_true_sp[i] = *sp_id;

        if selected.contains(name) { x_gene[i] = *sp_id; }

        if leaf_set.contains(name.as_str()) {
            event_true[i] = 3;
        } else {
            event_true[i] = *true_events.get(name)
                .ok_or_else(|| format!("Gene node '{}' not in true_events", name))?;
        }

        event_input[i] = if selected.contains(name) { event_true[i] } else { 4 };

        if frontier_set.contains(name.as_str()) { frontier_mask[i] = 1; }
        if leaf_set.contains(name.as_str()) { is_leaf[i] = 1; }
    }

    if let Some(ref la) = last_added_name {
        if let Some(&idx) = g_name_to_idx.get(la.as_str()) {
            if force_mask_last_added {
                event_input[idx] = 4;
                mask_label_node[idx] = 1;
            }
        }
    }

    let mut undirected_edges: HashSet<(i32, i32)> = HashSet::new();
    for (u, vs) in active_neighbors {
        let ui = g_name_to_idx[u.as_str()] as i32;
        for v in vs {
            if let Some(&vi_usize) = g_name_to_idx.get(v.as_str()) {
                let vi = vi_usize as i32;
                if ui != vi {
                    let e = if ui < vi { (ui, vi) } else { (vi, ui) };
                    undirected_edges.insert(e);
                }
            }
        }
    }

    let mut directed: Vec<(i32, i32)> = undirected_edges.into_iter().collect();
    directed.sort();
    let n_edges = directed.len();
    let mut g_edge_src = Vec::with_capacity(n_edges * 2);
    let mut g_edge_dst = Vec::with_capacity(n_edges * 2);
    let mut g_dir_src = Vec::with_capacity(n_edges);
    let mut g_dir_dst = Vec::with_capacity(n_edges);
    for &(a, b) in &directed {
        g_dir_src.push(a); g_dir_dst.push(b);
        g_edge_src.push(a); g_edge_dst.push(b);
        g_edge_src.push(b); g_edge_dst.push(a);
    }

    let root_pair_0 = true_root.first()
        .and_then(|s| g_name_to_idx.get(s.as_str())).map(|&v| v as i32).unwrap_or(-1);
    let root_pair_1 = true_root.get(1)
        .and_then(|s| g_name_to_idx.get(s.as_str())).map(|&v| v as i32).unwrap_or(-1);
    let root_edge_target: i32 = if !include_root_node && root_pair_0 >= 0 && root_pair_1 >= 0 {
        let pair_fwd = (root_pair_0.min(root_pair_1), root_pair_0.max(root_pair_1));
        directed.iter().position(|&e| e == pair_fwd).map(|i| i as i32).unwrap_or(-1)
    } else {
        -1
    };

    Ok(TaskTensors {
        x_gene, g_true_sp, event_true, event_input,
        frontier_mask, is_leaf, mask_label_node,
        g_edge_src, g_edge_dst, g_dir_src, g_dir_dst, root_edge_target,
    })
}

/// Collate per-sample task tensors into a single batch.
///
/// Also applies GCN normalization to gene edges and computes varlen attention metadata.
fn collate_task_tensors(tensors: &[TaskTensors]) -> CollatedTask {
    let mut gene_x = Vec::new();
    let mut gene_y = Vec::new();
    let mut event = Vec::new();
    let mut event_true = Vec::new();
    let mut frontier_mask = Vec::new();
    let mut is_leaf = Vec::new();
    let mut mask_label_node = Vec::new();
    let mut g_edge_src_raw = Vec::new();
    let mut g_edge_dst_raw = Vec::new();
    let mut g_dir_src = Vec::new();
    let mut g_dir_dst = Vec::new();
    let mut true_root_index = Vec::new();
    let mut g_sizes: Vec<i32> = Vec::new();
    let mut g_ptr: Vec<i32> = vec![0];
    let mut g_batch: Vec<i32> = Vec::new();
    let mut root_edges_ptr: Vec<i32> = vec![0];

    let mut cum_g: i32 = 0;
    let mut cum_re: i32 = 0;

    for (s_idx, t) in tensors.iter().enumerate() {
        let ng = t.x_gene.len() as i32;
        let ndir = t.g_dir_src.len() as i32;

        gene_x.extend_from_slice(&t.x_gene);
        gene_y.extend(t.g_true_sp.iter().map(|&v| v - 1));
        event.extend_from_slice(&t.event_input);
        event_true.extend_from_slice(&t.event_true);
        frontier_mask.extend_from_slice(&t.frontier_mask);
        is_leaf.extend_from_slice(&t.is_leaf);
        mask_label_node.extend_from_slice(&t.mask_label_node);

        for &s in &t.g_edge_src { g_edge_src_raw.push(s + cum_g); }
        for &d in &t.g_edge_dst { g_edge_dst_raw.push(d + cum_g); }
        for &s in &t.g_dir_src { g_dir_src.push(s + cum_g); }
        for &d in &t.g_dir_dst { g_dir_dst.push(d + cum_g); }

        true_root_index.push(t.root_edge_target);
        g_sizes.push(ng);
        for _ in 0..ng { g_batch.push(s_idx as i32); }
        cum_g += ng;
        g_ptr.push(cum_g);
        cum_re += ndir;
        root_edges_ptr.push(cum_re);
    }

    // Apply GCN normalization to gene edges
    let total_gene_nodes = cum_g as usize;
    let (gene_ei_src, gene_ei_dst, gene_ew) =
        gcn_norm_internal(&g_edge_src_raw, &g_edge_dst_raw, total_gene_nodes);

    // Compute varlen attention metadata
    let (cu_g, max_g) = compute_varlen_metadata(&g_sizes);

    CollatedTask {
        gene_x, gene_y, event, event_true,
        frontier_mask, is_leaf, mask_label_node,
        gene_ei_src, gene_ei_dst, gene_ew,
        g_dir_src, g_dir_dst,
        true_root_index, g_sizes, g_ptr, g_batch, root_edges_ptr,
        cu_g, max_g,
    }
}

/// Build a complete batched dict of numpy arrays for on-the-fly training.
///
/// Simulates one species tree and n_gene_trees gene trees, builds all three
/// task tensors (map / evt / root), collates them, and returns flat numpy arrays
/// ready to be passed to batch_to_heterodata() in Python.
///
/// Parameters
/// ----------
/// n_sp : int          Target number of extant species.
/// lambda_birth : float  Speciation rate (must be > mu_death).
/// mu_death : float    Extinction rate (must be >= 0).
/// sp_seed : int       Seed for species tree simulation.
/// n_gene_trees : int  Number of gene trees per batch.
/// gt_seeds : list[int]  Per-gene-tree simulation seeds (length n_gene_trees).
/// lambda_d/t/l : float  DTL rates for gene tree simulation.
/// enable_event : bool   Include evt_* keys in returned dict.
/// enable_root : bool    Include root_* keys in returned dict.
/// sample_order : str    "random" or "bottom_up".
/// map/evt/root_coloring_seeds : list[int]  Per-sample RNG seeds for tensor construction.
/// min_gene_leaves : int   Minimum extant gene leaves required.
/// max_gene_nodes : int    Maximum total gene nodes (0 = no limit).
/// max_retries : int       Retry attempts per gene tree slot on simulation failure.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (n_sp, lambda_birth, mu_death, sp_seed, n_gene_trees, gt_seeds, lambda_d, lambda_t, lambda_l, enable_event, enable_root, sample_order, map_coloring_seeds, evt_coloring_seeds, root_coloring_seeds, min_gene_leaves, max_gene_nodes, max_retries))]
pub fn build_otf_batch(
    py: Python,
    n_sp: usize,
    lambda_birth: f64,
    mu_death: f64,
    sp_seed: u64,
    n_gene_trees: usize,
    gt_seeds: Vec<u64>,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    enable_event: bool,
    enable_root: bool,
    sample_order: &str,
    map_coloring_seeds: Vec<u64>,
    evt_coloring_seeds: Vec<u64>,
    root_coloring_seeds: Vec<u64>,
    min_gene_leaves: usize,
    max_gene_nodes: usize,
    max_retries: usize,
) -> PyResult<PyObject> {
    use numpy::PyArray1;
    use pyo3::types::PyDict;
    use crate::node::TraversalOrder;

    // Validate
    if gt_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "gt_seeds length {} != n_gene_trees {}", gt_seeds.len(), n_gene_trees)));
    }
    if map_coloring_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "map_coloring_seeds length {} != n_gene_trees {}", map_coloring_seeds.len(), n_gene_trees)));
    }
    if enable_event && evt_coloring_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "evt_coloring_seeds length {} != n_gene_trees {}", evt_coloring_seeds.len(), n_gene_trees)));
    }
    if enable_root && root_coloring_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "root_coloring_seeds length {} != n_gene_trees {}", root_coloring_seeds.len(), n_gene_trees)));
    }
    validate_dtl_rates(lambda_d, lambda_t, lambda_l)?;

    // ---- GIL-free section ------------------------------------------------
    // All simulation, tensor building, and collation is pure Rust with no
    // Python object interaction.  Releasing the GIL here lets other Python
    // threads (prefetcher producers, GPU transfer thread) run concurrently.
    let sample_order_owned = sample_order.to_owned();
    let result_or_err = py.allow_threads(move || -> Result<_, String> {
        let sample_order: &str = &sample_order_owned;

        // 1. Simulate species tree
        let mut sp_rng = StdRng::seed_from_u64(sp_seed);
        let (mut sp_tree_raw, _) = simulate_bd_tree_bwd(n_sp, lambda_birth, mu_death, &mut sp_rng)?;
        sp_tree_raw.assign_depths();
        let sp_tree_raw_arc = Arc::new(sp_tree_raw);
        let (extant_sp, _) = extract_extant_subtree(&sp_tree_raw_arc)
            .ok_or_else(|| "build_otf_batch: failed to extract extant species subtree".to_string())?;
        let extant_sp_arc = Arc::new(extant_sp);
        let n_sp_nodes = extant_sp_arc.nodes.len();

        // 2. Simulate and filter gene trees (one per slot, individual seeds)
        //    Parallelised with rayon: each slot is independent (own seed, own RNG).
        let species_identity: Vec<Option<usize>> = (0..n_sp_nodes).map(|i| Some(i)).collect();

        let valid_gene_trees: Vec<RecTree> = (0..n_gene_trees).into_par_iter()
            .map(|i| {
                for attempt in 0..=(max_retries as u64) {
                    let gt_seed = gt_seeds[i].wrapping_add(attempt);
                    let mut gt_rng = StdRng::seed_from_u64(gt_seed);

                    let (mut rec_tree, events) = match simulate_dtl_per_species(
                        &extant_sp_arc, extant_sp_arc.root,
                        lambda_d, lambda_t, lambda_l,
                        None, None, false,
                        &mut gt_rng,
                    ) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };
                    rec_tree.species_tree = Arc::clone(&extant_sp_arc);
                    rec_tree.dtl_events = Some(events);

                    // Extract extant gene leaves
                    let extant_gene_indices: std::collections::HashSet<usize> = rec_tree.gene_tree.nodes
                        .iter().enumerate()
                        .filter(|(idx, n)| {
                            n.left_child.is_none() && n.right_child.is_none()
                                && rec_tree.event_mapping[*idx] == Event::Leaf
                        })
                        .map(|(idx, _)| idx)
                        .collect();

                    if extant_gene_indices.is_empty() { continue; }
                    if extant_gene_indices.len() < min_gene_leaves { continue; }

                    let (extant_gene_tree, gene_old_to_new) = match extract_induced_subtree(
                        &rec_tree.gene_tree, &extant_gene_indices,
                    ) {
                        Some(r) => r,
                        None => continue,
                    };

                    if max_gene_nodes > 0 && extant_gene_tree.nodes.len() > max_gene_nodes { continue; }

                    let (new_node_mapping, new_event_mapping) = match remap_gene_tree_indices(
                        &extant_gene_tree, &gene_old_to_new,
                        &rec_tree.node_mapping, &rec_tree.event_mapping,
                        &species_identity,
                    ) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };

                    return Ok(RecTree::new(
                        Arc::clone(&extant_sp_arc),
                        extant_gene_tree,
                        new_node_mapping,
                        new_event_mapping,
                    ));
                }
                Err(format!(
                    "build_otf_batch: gene tree slot {} failed after {} retries", i, max_retries))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // 3. Species tree tensors (computed once, replicated B times during collation)
        let sp_preorder: Vec<usize> = extant_sp_arc.iter_indices(TraversalOrder::PreOrder).collect();
        let species_names: Vec<String> = sp_preorder.iter()
            .map(|&i| extant_sp_arc.nodes[i].name.clone()).collect();

        let mut sp_name_to_id: HashMap<String, i32> = HashMap::new();
        let mut sp_name_to_idx_map: HashMap<String, usize> = HashMap::new();
        for (i, name) in species_names.iter().enumerate() {
            sp_name_to_id.insert(name.clone(), (i + 1) as i32);
            sp_name_to_idx_map.insert(name.clone(), i);
        }

        let x_sp_single: Vec<i32> = species_names.iter()
            .map(|name| sp_name_to_id[name]).collect();

        let mut sp_child_src_s: Vec<i32> = Vec::new();
        let mut sp_child_dst_s: Vec<i32> = Vec::new();
        let mut sp_parent_src_s: Vec<i32> = Vec::new();
        let mut sp_parent_dst_s: Vec<i32> = Vec::new();
        for &idx in &sp_preorder {
            let node = &extant_sp_arc.nodes[idx];
            let pi = sp_name_to_idx_map[&node.name] as i32;
            for child_opt in [node.left_child, node.right_child] {
                if let Some(child) = child_opt {
                    let ci = sp_name_to_idx_map[&extant_sp_arc.nodes[child].name] as i32;
                    sp_child_src_s.push(pi); sp_child_dst_s.push(ci);
                    sp_parent_src_s.push(ci); sp_parent_dst_s.push(pi);
                }
            }
        }

        // 4. Extract base samples from gene trees (parallel)
        let base_samples: Vec<RustBaseSample> = valid_gene_trees.par_iter()
            .map(|rt| extract_base_sample_internal(rt))
            .collect::<Result<Vec<_>, _>>()?;

        // 5. Build task tensors (parallel — each sample has its own seed/RNG)
        let map_tensors: Vec<TaskTensors> = base_samples.par_iter().enumerate()
            .map(|(i, s)| {
                let mut rng = StdRng::seed_from_u64(map_coloring_seeds[i]);
                build_task_tensors_internal(s, &sp_name_to_id, &mut rng, false, false, sample_order)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let evt_tensors: Vec<TaskTensors> = if enable_event {
            base_samples.par_iter().enumerate()
                .map(|(i, s)| {
                    let mut rng = StdRng::seed_from_u64(evt_coloring_seeds[i]);
                    build_task_tensors_internal(s, &sp_name_to_id, &mut rng, true, false, sample_order)
                })
                .collect::<Result<Vec<_>, _>>()?
        } else { Vec::new() };

        let root_tensors: Vec<TaskTensors> = if enable_root {
            base_samples.par_iter().enumerate()
                .map(|(i, s)| {
                    let mut rng = StdRng::seed_from_u64(root_coloring_seeds[i]);
                    build_task_tensors_internal(s, &sp_name_to_id, &mut rng, false, true, sample_order)
                })
                .collect::<Result<Vec<_>, _>>()?
        } else { Vec::new() };

        // 6. Compute species GCN norm once on single-tree edges, then replicate B times
        let (sp_child_norm_src_s, sp_child_norm_dst_s, sp_child_ew_s) =
            gcn_norm_internal(&sp_child_src_s, &sp_child_dst_s, n_sp_nodes);
        let (sp_parent_norm_src_s, sp_parent_norm_dst_s, sp_parent_ew_s) =
            gcn_norm_internal(&sp_parent_src_s, &sp_parent_dst_s, n_sp_nodes);

        let b = n_gene_trees;
        let n_sp_child_e = sp_child_norm_src_s.len(); // includes self-loops
        let n_sp_parent_e = sp_parent_norm_src_s.len();
        let mut sp_x: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_child_src: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_dst: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_ew: Vec<f32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_parent_src: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_dst: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_ew: Vec<f32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_batch: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_sizes: Vec<i32> = Vec::with_capacity(b);
        let mut sp_ptr: Vec<i32> = Vec::with_capacity(b + 1);
        sp_ptr.push(0);

        for s_idx in 0..b {
            let offset = (s_idx * n_sp_nodes) as i32;
            sp_x.extend_from_slice(&x_sp_single);
            for (&s, &d) in sp_child_norm_src_s.iter().zip(sp_child_norm_dst_s.iter()) {
                sp_child_src.push(s + offset); sp_child_dst.push(d + offset);
            }
            sp_child_ew.extend_from_slice(&sp_child_ew_s);
            for (&s, &d) in sp_parent_norm_src_s.iter().zip(sp_parent_norm_dst_s.iter()) {
                sp_parent_src.push(s + offset); sp_parent_dst.push(d + offset);
            }
            sp_parent_ew.extend_from_slice(&sp_parent_ew_s);
            for _ in 0..n_sp_nodes { sp_batch.push(s_idx as i32); }
            sp_sizes.push(n_sp_nodes as i32);
            sp_ptr.push((s_idx + 1) as i32 * n_sp_nodes as i32);
        }

        // Compute species varlen attention metadata
        let (cu_sp, max_sp) = compute_varlen_metadata(&sp_sizes);

        // 7. Collate task tensors (includes gene GCN norm + varlen metadata)
        let map_c = collate_task_tensors(&map_tensors);
        let evt_c = if enable_event { Some(collate_task_tensors(&evt_tensors)) } else { None };
        let root_c = if enable_root { Some(collate_task_tensors(&root_tensors)) } else { None };

        Ok((sp_x, sp_child_src, sp_child_dst, sp_child_ew,
            sp_parent_src, sp_parent_dst, sp_parent_ew,
            sp_ptr, sp_sizes, sp_batch, cu_sp, max_sp,
            map_c, evt_c, root_c))
    });
    // ---- End GIL-free section --------------------------------------------

    let (sp_x, sp_child_src, sp_child_dst, sp_child_ew,
         sp_parent_src, sp_parent_dst, sp_parent_ew,
         sp_ptr, sp_sizes, sp_batch, cu_sp, max_sp,
         map_c, evt_c, root_c) =
        result_or_err.map_err(|e| PyValueError::new_err(e))?;

    // 8. Build result dict of numpy arrays
    let result = PyDict::new(py);

    let edge_2d = |src: &[i32], dst: &[i32], name: &str| -> PyResult<PyObject> {
        use numpy::ndarray::Array2;
        use numpy::ToPyArray;
        let n = src.len();
        let mut flat = Vec::with_capacity(2 * n);
        flat.extend_from_slice(src);
        flat.extend_from_slice(dst);
        let arr = Array2::from_shape_vec([2, n], flat)
            .map_err(|e| PyValueError::new_err(format!("{}: {}", name, e)))?;
        Ok(arr.to_pyarray(py).into_any().unbind())
    };

    // Species (GCN-normalized edges with self-loops)
    result.set_item("sp_x", PyArray1::from_slice(py, &sp_x))?;
    result.set_item("sp_child_ei", edge_2d(&sp_child_src, &sp_child_dst, "sp_child_ei")?)?;
    result.set_item("sp_child_ew", PyArray1::from_slice(py, &sp_child_ew))?;
    result.set_item("sp_parent_ei", edge_2d(&sp_parent_src, &sp_parent_dst, "sp_parent_ei")?)?;
    result.set_item("sp_parent_ew", PyArray1::from_slice(py, &sp_parent_ew))?;
    result.set_item("sp_ptr", PyArray1::from_slice(py, &sp_ptr))?;
    result.set_item("sp_sizes", PyArray1::from_slice(py, &sp_sizes))?;
    result.set_item("sp_batch", PyArray1::from_slice(py, &sp_batch))?;
    result.set_item("cu_sp", PyArray1::from_slice(py, &cu_sp))?;
    result.set_item("max_sp", max_sp)?;

    // Helper: write one collated task into result dict under a given prefix
    let write_task = |c: &CollatedTask, prefix: &str| -> PyResult<()> {
        result.set_item(format!("{}_gene_x", prefix), PyArray1::from_slice(py, &c.gene_x))?;
        result.set_item(format!("{}_gene_y", prefix), PyArray1::from_slice(py, &c.gene_y))?;
        result.set_item(format!("{}_event", prefix), PyArray1::from_slice(py, &c.event))?;
        result.set_item(format!("{}_event_true", prefix), PyArray1::from_slice(py, &c.event_true))?;
        result.set_item(format!("{}_is_leaf", prefix), PyArray1::from_slice(py, &c.is_leaf))?;
        result.set_item(format!("{}_frontier_mask", prefix), PyArray1::from_slice(py, &c.frontier_mask))?;
        result.set_item(format!("{}_mask_label_node", prefix), PyArray1::from_slice(py, &c.mask_label_node))?;
        result.set_item(format!("{}_gene_ei", prefix), edge_2d(&c.gene_ei_src, &c.gene_ei_dst, &format!("{}_gene_ei", prefix))?)?;
        result.set_item(format!("{}_gene_ew", prefix), PyArray1::from_slice(py, &c.gene_ew))?;
        result.set_item(format!("{}_g_dir_edge", prefix), edge_2d(&c.g_dir_src, &c.g_dir_dst, &format!("{}_g_dir_edge", prefix))?)?;
        result.set_item(format!("{}_true_root_index", prefix), PyArray1::from_slice(py, &c.true_root_index))?;
        result.set_item(format!("{}_g_ptr", prefix), PyArray1::from_slice(py, &c.g_ptr))?;
        result.set_item(format!("{}_g_sizes", prefix), PyArray1::from_slice(py, &c.g_sizes))?;
        result.set_item(format!("{}_g_batch", prefix), PyArray1::from_slice(py, &c.g_batch))?;
        result.set_item(format!("{}_root_edges_ptr", prefix), PyArray1::from_slice(py, &c.root_edges_ptr))?;
        result.set_item(format!("{}_cu_g", prefix), PyArray1::from_slice(py, &c.cu_g))?;
        result.set_item(format!("{}_max_g", prefix), c.max_g)?;
        Ok(())
    };

    write_task(&map_c, "map")?;
    if let Some(ref c) = evt_c { write_task(c, "evt")?; }
    if let Some(ref c) = root_c { write_task(c, "root")?; }

    Ok(result.into())
}

/// Compute GCN normalization on edge indices.
///
/// Adds self-loops and computes D^{-1/2} A D^{-1/2} edge weights.
/// Matches PyG's `gcn_norm(edge_index, None, num_nodes, improved=False,
/// add_self_loops=True, flow="source_to_target")`.
///
/// Parameters
/// ----------
/// edge_index : numpy array of shape (2, E), int32
/// num_nodes : int
///
/// Returns (edge_index_with_self_loops, edge_weights) as numpy arrays.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (edge_index, num_nodes))]
pub fn compute_gcn_norm(py: Python, edge_index: &Bound<'_, numpy::PyArray2<i32>>, num_nodes: usize) -> PyResult<(PyObject, PyObject)> {
    use numpy::{PyArray1, PyArrayMethods};
    use numpy::ndarray::Array2;
    use numpy::ToPyArray;

    let ei = unsafe { edge_index.as_array() };
    let src: Vec<i32> = ei.row(0).to_vec();
    let dst: Vec<i32> = ei.row(1).to_vec();

    let (new_src, new_dst, weights) = gcn_norm_internal(&src, &dst, num_nodes);

    let n = new_src.len();
    let mut flat = Vec::with_capacity(2 * n);
    flat.extend_from_slice(&new_src);
    flat.extend_from_slice(&new_dst);
    let arr = Array2::from_shape_vec([2, n], flat)
        .map_err(|e| PyValueError::new_err(format!("compute_gcn_norm: {}", e)))?;

    let ei_out = arr.to_pyarray(py).into_any().unbind();
    let ew_out = PyArray1::from_slice(py, &weights).into_any().unbind();
    Ok((ei_out, ew_out))
}

/// Build a batched dict of numpy arrays for inference from a base_sample dict.
///
/// Replaces the Python get_sample() + collate() path. Builds task tensors,
/// collates n_copies, applies GCN norm and computes varlen metadata.
/// Returns a flat numpy dict in the same format as build_otf_batch().
///
/// Only mapping task tensors are returned (inference doesn't need separate evt/root batches).
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (base_sample, n_copies, seed, sample_order))]
pub fn build_inference_batch(
    py: Python,
    base_sample: &Bound<'_, pyo3::types::PyDict>,
    n_copies: usize,
    seed: u64,
    sample_order: &str,
) -> PyResult<PyObject> {
    use numpy::PyArray1;
    use pyo3::types::PyDict;

    // Parse the base_sample dict into RustBaseSample
    let base = parse_base_sample_dict(base_sample)?;

    // Build species name → 1-based ID mapping
    let mut sp_name_to_id: HashMap<String, i32> = HashMap::new();
    for (i, name) in base.species_names.iter().enumerate() {
        sp_name_to_id.insert(name.clone(), (i + 1) as i32);
    }

    let n_sp_nodes = base.species_names.len();
    let sample_order_owned = sample_order.to_owned();

    let result_or_err = py.allow_threads(move || -> Result<_, String> {
        let sample_order: &str = &sample_order_owned;

        // Build species tensors (single copy)
        let x_sp_single: Vec<i32> = base.species_names.iter()
            .map(|name| sp_name_to_id[name]).collect();

        let mut sp_child_src_s: Vec<i32> = Vec::new();
        let mut sp_child_dst_s: Vec<i32> = Vec::new();
        let mut sp_parent_src_s: Vec<i32> = Vec::new();
        let mut sp_parent_dst_s: Vec<i32> = Vec::new();

        // Build species name to index mapping
        let mut sp_name_to_idx: HashMap<String, usize> = HashMap::new();
        for (i, name) in base.species_names.iter().enumerate() {
            sp_name_to_idx.insert(name.clone(), i);
        }

        for (parent_name, children) in &base.sp_children {
            if let Some(&pi) = sp_name_to_idx.get(parent_name) {
                for child_name in children {
                    if let Some(&ci) = sp_name_to_idx.get(child_name) {
                        sp_child_src_s.push(pi as i32);
                        sp_child_dst_s.push(ci as i32);
                        sp_parent_src_s.push(ci as i32);
                        sp_parent_dst_s.push(pi as i32);
                    }
                }
            }
        }

        // GCN norm on species edges (single copy)
        let (sp_child_norm_src_s, sp_child_norm_dst_s, sp_child_ew_s) =
            gcn_norm_internal(&sp_child_src_s, &sp_child_dst_s, n_sp_nodes);
        let (sp_parent_norm_src_s, sp_parent_norm_dst_s, sp_parent_ew_s) =
            gcn_norm_internal(&sp_parent_src_s, &sp_parent_dst_s, n_sp_nodes);

        // Build mapping task tensors for each copy
        let mut map_tensors: Vec<TaskTensors> = Vec::with_capacity(n_copies);
        for copy_idx in 0..n_copies {
            let copy_seed = seed.wrapping_add(copy_idx as u64 * 1000003);
            let mut rng = StdRng::seed_from_u64(copy_seed);
            let t = build_task_tensors_internal(&base, &sp_name_to_id, &mut rng, false, false, sample_order)?;
            map_tensors.push(t);
        }

        // Collate
        let map_c = collate_task_tensors(&map_tensors);

        // Replicate species B times
        let b = n_copies;
        let n_sp_child_e = sp_child_norm_src_s.len();
        let n_sp_parent_e = sp_parent_norm_src_s.len();
        let mut sp_x: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_child_src: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_dst: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_ew: Vec<f32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_parent_src: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_dst: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_ew: Vec<f32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_batch: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_sizes: Vec<i32> = Vec::with_capacity(b);
        let mut sp_ptr: Vec<i32> = Vec::with_capacity(b + 1);
        sp_ptr.push(0);

        for s_idx in 0..b {
            let offset = (s_idx * n_sp_nodes) as i32;
            sp_x.extend_from_slice(&x_sp_single);
            for (&s, &d) in sp_child_norm_src_s.iter().zip(sp_child_norm_dst_s.iter()) {
                sp_child_src.push(s + offset); sp_child_dst.push(d + offset);
            }
            sp_child_ew.extend_from_slice(&sp_child_ew_s);
            for (&s, &d) in sp_parent_norm_src_s.iter().zip(sp_parent_norm_dst_s.iter()) {
                sp_parent_src.push(s + offset); sp_parent_dst.push(d + offset);
            }
            sp_parent_ew.extend_from_slice(&sp_parent_ew_s);
            for _ in 0..n_sp_nodes { sp_batch.push(s_idx as i32); }
            sp_sizes.push(n_sp_nodes as i32);
            sp_ptr.push((s_idx + 1) as i32 * n_sp_nodes as i32);
        }

        let (cu_sp, max_sp) = compute_varlen_metadata(&sp_sizes);

        Ok((sp_x, sp_child_src, sp_child_dst, sp_child_ew,
            sp_parent_src, sp_parent_dst, sp_parent_ew,
            sp_ptr, sp_sizes, sp_batch, cu_sp, max_sp, map_c))
    });

    let (sp_x, sp_child_src, sp_child_dst, sp_child_ew,
         sp_parent_src, sp_parent_dst, sp_parent_ew,
         sp_ptr, sp_sizes, sp_batch, cu_sp, max_sp, map_c) =
        result_or_err.map_err(|e| PyValueError::new_err(e))?;

    let result = PyDict::new(py);

    let edge_2d = |src: &[i32], dst: &[i32], name: &str| -> PyResult<PyObject> {
        use numpy::ndarray::Array2;
        use numpy::ToPyArray;
        let n = src.len();
        let mut flat = Vec::with_capacity(2 * n);
        flat.extend_from_slice(src);
        flat.extend_from_slice(dst);
        let arr = Array2::from_shape_vec([2, n], flat)
            .map_err(|e| PyValueError::new_err(format!("{}: {}", name, e)))?;
        Ok(arr.to_pyarray(py).into_any().unbind())
    };

    // Species
    result.set_item("sp_x", PyArray1::from_slice(py, &sp_x))?;
    result.set_item("sp_child_ei", edge_2d(&sp_child_src, &sp_child_dst, "sp_child_ei")?)?;
    result.set_item("sp_child_ew", PyArray1::from_slice(py, &sp_child_ew))?;
    result.set_item("sp_parent_ei", edge_2d(&sp_parent_src, &sp_parent_dst, "sp_parent_ei")?)?;
    result.set_item("sp_parent_ew", PyArray1::from_slice(py, &sp_parent_ew))?;
    result.set_item("sp_ptr", PyArray1::from_slice(py, &sp_ptr))?;
    result.set_item("sp_sizes", PyArray1::from_slice(py, &sp_sizes))?;
    result.set_item("sp_batch", PyArray1::from_slice(py, &sp_batch))?;
    result.set_item("cu_sp", PyArray1::from_slice(py, &cu_sp))?;
    result.set_item("max_sp", max_sp)?;

    // Mapping task (only task needed for inference)
    let c = &map_c;
    result.set_item("map_gene_x", PyArray1::from_slice(py, &c.gene_x))?;
    result.set_item("map_gene_y", PyArray1::from_slice(py, &c.gene_y))?;
    result.set_item("map_event", PyArray1::from_slice(py, &c.event))?;
    result.set_item("map_event_true", PyArray1::from_slice(py, &c.event_true))?;
    result.set_item("map_is_leaf", PyArray1::from_slice(py, &c.is_leaf))?;
    result.set_item("map_frontier_mask", PyArray1::from_slice(py, &c.frontier_mask))?;
    result.set_item("map_mask_label_node", PyArray1::from_slice(py, &c.mask_label_node))?;
    result.set_item("map_gene_ei", edge_2d(&c.gene_ei_src, &c.gene_ei_dst, "map_gene_ei")?)?;
    result.set_item("map_gene_ew", PyArray1::from_slice(py, &c.gene_ew))?;
    result.set_item("map_g_dir_edge", edge_2d(&c.g_dir_src, &c.g_dir_dst, "map_g_dir_edge")?)?;
    result.set_item("map_true_root_index", PyArray1::from_slice(py, &c.true_root_index))?;
    result.set_item("map_g_ptr", PyArray1::from_slice(py, &c.g_ptr))?;
    result.set_item("map_g_sizes", PyArray1::from_slice(py, &c.g_sizes))?;
    result.set_item("map_g_batch", PyArray1::from_slice(py, &c.g_batch))?;
    result.set_item("map_root_edges_ptr", PyArray1::from_slice(py, &c.root_edges_ptr))?;
    result.set_item("map_cu_g", PyArray1::from_slice(py, &c.cu_g))?;
    result.set_item("map_max_g", c.max_g)?;

    Ok(result.into())
}

/// Parse a Python base_sample dict into a RustBaseSample.
fn parse_base_sample_dict(d: &Bound<'_, pyo3::types::PyDict>) -> PyResult<RustBaseSample> {
    let species_names: Vec<String> = d.get_item("species_names")?
        .ok_or_else(|| PyValueError::new_err("missing species_names"))?.extract()?;
    let g_root_name: String = d.get_item("g_root_name")?
        .ok_or_else(|| PyValueError::new_err("missing g_root_name"))?.extract()?;
    let g_leaves_names: Vec<String> = d.get_item("g_leaves_names")?
        .ok_or_else(|| PyValueError::new_err("missing g_leaves_names"))?.extract()?;
    let nb_g_leaves = g_leaves_names.len();
    let true_root: Vec<String> = d.get_item("true_root")?
        .ok_or_else(|| PyValueError::new_err("missing true_root"))?.extract()?;
    let g_neighbors: HashMap<String, Vec<String>> = d.get_item("g_neighbors")?
        .ok_or_else(|| PyValueError::new_err("missing g_neighbors"))?.extract()?;
    let sp_children: HashMap<String, Vec<String>> = d.get_item("sp_children")?
        .ok_or_else(|| PyValueError::new_err("missing sp_children"))?.extract()?;
    let true_states: HashMap<String, String> = d.get_item("true_states")?
        .ok_or_else(|| PyValueError::new_err("missing true_states"))?.extract()?;
    let true_events: HashMap<String, i32> = d.get_item("true_events")?
        .ok_or_else(|| PyValueError::new_err("missing true_events"))?.extract()?;

    Ok(RustBaseSample {
        g_root_name, nb_g_leaves, species_names, g_leaves_names,
        true_root, g_neighbors, sp_children, true_states, true_events,
    })
}

/// Build a PyGeneTree from a species tree, gene tree, and per-node reconciliation annotations.
///
/// This allows constructing a fully annotated reconciled gene tree without going through
/// RecPhyloXML, e.g. from model predictions.
///
/// # Arguments
/// * `sp_newick` - Species tree in Newick format (string, not file path)
/// * `g_newick` - Gene tree in Newick format (string, not file path)
/// * `node_species` - Dict mapping gene node name → species node name
/// * `node_events` - Dict mapping gene node name → event code (0=Speciation, 1=Duplication, 2=Transfer, 3=Leaf)
///
/// # Returns
/// A PyGeneTree with the given reconciliation annotations.
#[pyfunction]
#[pyo3(signature = (sp_newick, g_newick, node_species, node_events))]
pub fn from_reconciliation(
    sp_newick: &str,
    g_newick: &str,
    node_species: &Bound<'_, pyo3::types::PyDict>,
    node_events: &Bound<'_, pyo3::types::PyDict>,
) -> PyResult<PyGeneTree> {
    use crate::newick::newick::parse_newick;
    use crate::node::rectree::{Event, RecTree};

    // Parse species tree
    let mut sp_nodes = parse_newick(sp_newick)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse species Newick: {}", e)))?;
    let sp_root = sp_nodes.pop()
        .ok_or_else(|| PyValueError::new_err("No tree found in species Newick"))?;
    let mut sp_tree = sp_root.to_flat_tree();
    sp_tree.assign_depths();

    // Build species name → index map
    let sp_name_to_idx: std::collections::HashMap<String, usize> = sp_tree
        .nodes
        .iter()
        .enumerate()
        .map(|(idx, node)| (node.name.clone(), idx))
        .collect();

    // Parse gene tree
    let mut g_nodes = parse_newick(g_newick)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse gene Newick: {}", e)))?;
    let g_root = g_nodes.pop()
        .ok_or_else(|| PyValueError::new_err("No tree found in gene Newick"))?;
    let mut g_tree = g_root.to_flat_tree();
    g_tree.assign_depths();

    let n_gene_nodes = g_tree.nodes.len();

    // Build node_mapping and event_mapping
    let mut node_mapping: Vec<Option<usize>> = vec![None; n_gene_nodes];
    let mut event_mapping: Vec<Event> = vec![Event::Speciation; n_gene_nodes];

    for (idx, node) in g_tree.nodes.iter().enumerate() {
        let name = &node.name;

        // Species mapping
        if node_species.contains(name.as_str())? {
            let sp_name_obj = node_species.get_item(name.as_str())?
                .ok_or_else(|| PyValueError::new_err(format!("Key '{}' disappeared from node_species", name)))?;
            let sp_name: String = sp_name_obj.extract()?;
            if let Some(&sp_idx) = sp_name_to_idx.get(&sp_name) {
                node_mapping[idx] = Some(sp_idx);
            } else {
                return Err(PyValueError::new_err(format!(
                    "Species '{}' (mapped from gene node '{}') not found in species tree",
                    sp_name, name
                )));
            }
        }
        // If not in dict, mapping stays None (unknown)

        // Event mapping
        if node_events.contains(name.as_str())? {
            let evt_obj = node_events.get_item(name.as_str())?
                .ok_or_else(|| PyValueError::new_err(format!("Key '{}' disappeared from node_events", name)))?;
            let evt_code: i32 = evt_obj.extract()?;
            event_mapping[idx] = match evt_code {
                0 => Event::Speciation,
                1 => Event::Duplication,
                2 => Event::Transfer,
                3 => Event::Leaf,
                _ => Event::Speciation, // fallback
            };
        } else {
            // Default: leaf if no children, speciation otherwise
            event_mapping[idx] = if node.left_child.is_none() && node.right_child.is_none() {
                Event::Leaf
            } else {
                Event::Speciation
            };
        }
    }

    let rec_tree = RecTree::new_owned(sp_tree, g_tree, node_mapping, event_mapping);
    Ok(PyGeneTree { rec_tree })
}
