//! R bindings for rustree using extendr.
//!
//! Provides R access to:
//! - Newick tree parsing
//! - Birth-death species tree simulation
//! - DTL gene tree simulation
//! - Export to Newick, XML, CSV
//! - SVG visualization via thirdkind

mod conversions;
use conversions::*;

use extendr_api::prelude::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::fs;
use std::process::Command;

use crate::bd::{simulate_bd_tree_bwd, generate_events_from_tree};
use crate::dtl::{simulate_dtl, simulate_dtl_batch, simulate_dtl_per_species, simulate_dtl_per_species_batch};
use crate::dtl::{simulate_dtl_iter, simulate_dtl_per_species_iter};
use crate::induced_transfers::induced_transfers;
use crate::node::{FlatTree, Event, RecTree, GeneForest, remap_gene_tree_indices};
use crate::sampling::{extract_induced_subtree, extract_induced_subtree_by_names};

/// Simulate a birth-death species tree.
///
/// @param n Number of extant species
/// @param lambda Speciation rate (birth rate)
/// @param mu Extinction rate (death rate), must be less than lambda
/// @param seed Optional random seed for reproducibility
/// @return A list containing the tree data
/// @export
#[extendr]
fn simulate_species_tree_r(n: i32, lambda: f64, mu: f64, seed: Robj) -> Result<List> {
    if n <= 0 {
        return Err("Number of species must be positive".into());
    }
    if lambda <= 0.0 {
        return Err("Speciation rate must be positive".into());
    }
    if mu < 0.0 {
        return Err("Extinction rate must be non-negative".into());
    }
    if lambda <= mu {
        return Err("Speciation rate must be greater than extinction rate".into());
    }

    let mut rng = if seed.is_null() || seed.is_na() {
        StdRng::from_entropy()
    } else {
        match seed.as_integer() {
            Some(s) => StdRng::seed_from_u64(s as u64),
            None => return Err("seed must be an integer".into()),
        }
    };

    let (mut tree, events) = simulate_bd_tree_bwd(n as usize, lambda, mu, &mut rng)
        .map_err(|e| Error::Other(e))?;
    tree.assign_depths();

    // Serialize tree and events to R list
    let names: Vec<String> = tree.nodes.iter().map(|n| n.name.clone()).collect();
    let parents: Vec<Rint> = tree.nodes.iter()
        .map(|n| n.parent.map(|p| Rint::from(p as i32)).unwrap_or(Rint::na()))
        .collect();
    let left_children: Vec<Rint> = tree.nodes.iter()
        .map(|n| n.left_child.map(|c| Rint::from(c as i32)).unwrap_or(Rint::na()))
        .collect();
    let right_children: Vec<Rint> = tree.nodes.iter()
        .map(|n| n.right_child.map(|c| Rint::from(c as i32)).unwrap_or(Rint::na()))
        .collect();
    let lengths: Vec<f64> = tree.nodes.iter().map(|n| n.length).collect();
    let depths: Vec<Rfloat> = tree.nodes.iter()
        .map(|n| n.depth.map(|d| Rfloat::from(d)).unwrap_or(Rfloat::na()))
        .collect();
    // bd_event: convert Option<BDEvent> to character vector for R
    let bd_events: Vec<Rstr> = tree.nodes.iter()
        .map(|n| match n.bd_event {
            Some(e) => Rstr::from(e.as_str()),
            None => Rstr::na(),
        })
        .collect();

    let events_list = bd_events_to_rlist(&events);

    Ok(list!(
        name = names,
        parent = parents,
        left_child = left_children,
        right_child = right_children,
        length = lengths,
        depth = depths,
        root = tree.root as i32,
        bd_event = bd_events,
        events = events_list
    ))
}

/// Parse a Newick string into a species tree.
///
/// @param newick_str A Newick formatted string
/// @return A list containing the tree data
/// @export
#[extendr]
fn parse_newick_r(newick_str: &str) -> Result<List> {
    use crate::newick::parse_newick;

    let mut nodes = parse_newick(newick_str)
        .map_err(|e| format!("Failed to parse Newick: {}", e))?;

    let mut root = nodes.pop()
        .ok_or("No tree found in Newick string")?;

    root.zero_root_length();
    root.assign_depths(0.0);
    let tree = root.to_flat_tree();

    Ok(flattree_to_rlist(&tree))
}

/// Convert a species tree to Newick format.
///
/// @param tree_list A tree list from simulate_species_tree_r or parse_newick_r
/// @return Newick string representation
/// @export
#[extendr]
fn tree_to_newick_r(tree_list: List) -> Result<String> {
    let tree = rlist_to_flattree(&tree_list)?;
    let nwk = tree.to_newick().map_err(|e| Error::Other(e))?;
    Ok(nwk + ";")
}

/// Get the number of leaves in a tree.
///
/// @param tree_list A tree list
/// @return Number of leaves
/// @export
#[extendr]
fn tree_num_leaves_r(tree_list: List) -> Result<i32> {
    let tree = rlist_to_flattree(&tree_list)?;
    let count = tree.nodes.iter()
        .filter(|n| n.left_child.is_none() && n.right_child.is_none())
        .count();
    Ok(count as i32)
}

/// Get leaf names from a tree.
///
/// @param tree_list A tree list
/// @return Character vector of leaf names
/// @export
#[extendr]
fn tree_leaf_names_r(tree_list: List) -> Result<Vec<String>> {
    let tree = rlist_to_flattree(&tree_list)?;
    let names: Vec<String> = tree.nodes.iter()
        .filter(|n| n.left_child.is_none() && n.right_child.is_none())
        .map(|n| n.name.clone())
        .collect();
    Ok(names)
}

/// Simulate a single DTL gene tree.
///
/// @param species_tree_list Species tree from simulate_species_tree_r
/// @param lambda_d Duplication rate
/// @param lambda_t Transfer rate
/// @param lambda_l Loss rate
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, retry until a tree with extant genes is produced (default FALSE)
/// @param seed Optional random seed
/// @return A list containing the gene tree data with reconciliation info
/// @export
#[extendr]
fn simulate_dtl_r(
    species_tree_list: List,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
) -> Result<List> {
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = rlist_to_flattree(&species_tree_list)?;

    let mut rng = if seed.is_null() || seed.is_na() {
        StdRng::from_entropy()
    } else {
        match seed.as_integer() {
            Some(s) => StdRng::seed_from_u64(s as u64),
            None => return Err("seed must be an integer".into()),
        }
    };

    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    let (rec_tree, events) = simulate_dtl(&species_tree, origin_species, lambda_d, lambda_t, lambda_l, alpha, repl, require_extant, &mut rng)
        .map_err(|e| Error::Other(e))?;

    let result = rectree_to_rlist(&species_tree, &rec_tree.gene_tree, &rec_tree.node_mapping, &rec_tree.event_mapping)?;
    result.set_attrib("dtl_events", dtl_events_to_rlist(&events, &species_tree))?;
    Ok(result)
}

/// Simulate a batch of DTL gene trees efficiently.
///
/// @param species_tree_list Species tree from simulate_species_tree_r
/// @param n Number of gene trees to simulate
/// @param lambda_d Duplication rate
/// @param lambda_t Transfer rate
/// @param lambda_l Loss rate
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, only include trees with extant genes (default FALSE)
/// @param seed Optional random seed
/// @return A list of gene tree lists
/// @export
#[extendr]
fn simulate_dtl_batch_r(
    species_tree_list: List,
    n: i32,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
) -> Result<List> {
    if n <= 0 {
        return Err("Number of trees must be positive".into());
    }
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = rlist_to_flattree(&species_tree_list)?;

    let mut rng = if seed.is_null() || seed.is_na() {
        StdRng::from_entropy()
    } else {
        match seed.as_integer() {
            Some(s) => StdRng::seed_from_u64(s as u64),
            None => return Err("seed must be an integer".into()),
        }
    };

    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    let (rec_trees, all_events) = simulate_dtl_batch(
        &species_tree,
        origin_species,
        lambda_d,
        lambda_t,
        lambda_l,
        alpha,
        repl,
        n as usize,
        require_extant,
        &mut rng,
    ).map_err(|e| Error::Other(e))?;

    // Convert to list of lists, attaching DTL events to each gene tree
    let gene_tree_lists: Vec<List> = rec_trees
        .iter()
        .zip(all_events.iter())
        .map(|(rec_tree, events)| {
            let gt_list = rectree_to_rlist(&species_tree, &rec_tree.gene_tree, &rec_tree.node_mapping, &rec_tree.event_mapping)?;
            gt_list.set_attrib("dtl_events", dtl_events_to_rlist(events, &species_tree))?;
            Ok(gt_list)
        })
        .collect::<Result<Vec<List>>>()?;

    Ok(List::from_values(gene_tree_lists))
}

/// Simulate a single DTL gene tree using the Zombi-style per-species model.
///
/// In this model, the event rate is proportional to the number of SPECIES
/// with gene copies, NOT the number of gene copies themselves. This means
/// duplications don't increase the event rate.
///
/// @param species_tree_list Species tree from simulate_species_tree_r
/// @param lambda_d Duplication rate per species per unit time
/// @param lambda_t Transfer rate per species per unit time
/// @param lambda_l Loss rate per species per unit time
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, retry until a tree with extant genes is produced (default FALSE)
/// @param seed Optional random seed
/// @return A list containing the gene tree data with reconciliation info
/// @export
#[extendr]
fn simulate_dtl_per_species_r(
    species_tree_list: List,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
) -> Result<List> {
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = rlist_to_flattree(&species_tree_list)?;

    let mut rng = if seed.is_null() || seed.is_na() {
        StdRng::from_entropy()
    } else {
        match seed.as_integer() {
            Some(s) => StdRng::seed_from_u64(s as u64),
            None => return Err("seed must be an integer".into()),
        }
    };

    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    let (rec_tree, events) = simulate_dtl_per_species(&species_tree, origin_species, lambda_d, lambda_t, lambda_l, alpha, repl, require_extant, &mut rng)
        .map_err(|e| Error::Other(e))?;

    let result = rectree_to_rlist(&species_tree, &rec_tree.gene_tree, &rec_tree.node_mapping, &rec_tree.event_mapping)?;
    result.set_attrib("dtl_events", dtl_events_to_rlist(&events, &species_tree))?;
    Ok(result)
}

/// Simulate a batch of DTL gene trees using the Zombi-style per-species model.
///
/// This is faster than calling simulate_dtl_per_species_r multiple times because
/// species tree events and LCA depths are computed only once.
///
/// @param species_tree_list Species tree from simulate_species_tree_r
/// @param n Number of gene trees to simulate
/// @param lambda_d Duplication rate per species per unit time
/// @param lambda_t Transfer rate per species per unit time
/// @param lambda_l Loss rate per species per unit time
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, only include trees with extant genes (default FALSE)
/// @param seed Optional random seed
/// @return A list of gene tree lists
/// @export
#[extendr]
fn simulate_dtl_per_species_batch_r(
    species_tree_list: List,
    n: i32,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
) -> Result<List> {
    if n <= 0 {
        return Err("Number of trees must be positive".into());
    }
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = rlist_to_flattree(&species_tree_list)?;

    let mut rng = if seed.is_null() || seed.is_na() {
        StdRng::from_entropy()
    } else {
        match seed.as_integer() {
            Some(s) => StdRng::seed_from_u64(s as u64),
            None => return Err("seed must be an integer".into()),
        }
    };

    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    let (rec_trees, all_events) = simulate_dtl_per_species_batch(
        &species_tree,
        origin_species,
        lambda_d,
        lambda_t,
        lambda_l,
        alpha,
        repl,
        n as usize,
        require_extant,
        &mut rng,
    ).map_err(|e| Error::Other(e))?;

    // Convert to list of lists, attaching DTL events to each gene tree
    let gene_tree_lists: Vec<List> = rec_trees
        .iter()
        .zip(all_events.iter())
        .map(|(rec_tree, events)| {
            let gt_list = rectree_to_rlist(&species_tree, &rec_tree.gene_tree, &rec_tree.node_mapping, &rec_tree.event_mapping)?;
            gt_list.set_attrib("dtl_events", dtl_events_to_rlist(events, &species_tree))?;
            Ok(gt_list)
        })
        .collect::<Result<Vec<List>>>()?;

    Ok(List::from_values(gene_tree_lists))
}

/// Get the number of extant genes in a gene tree.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @return Number of extant genes (non-loss leaves)
/// @export
#[extendr]
fn gene_tree_num_extant_r(gene_tree_list: List) -> Result<i32> {
    let events: Vec<String> = gene_tree_list.dollar("event")?.as_str_vector()
        .ok_or("Failed to get event column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let left_children: Vec<i32> = gene_tree_list.dollar("left_child")?.as_integer_vector()
        .ok_or("Failed to get left_child column")?;

    let count = events.iter()
        .zip(left_children.iter())
        .filter(|(e, lc)| *e == "Leaf" && lc.is_na())
        .count();

    Ok(count as i32)
}

/// Convert a gene tree to Newick format.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @return Newick string representation
/// @export
#[extendr]
fn gene_tree_to_newick_r(gene_tree_list: List) -> Result<String> {
    let (gene_tree, _, _, _) = rlist_to_genetree(&gene_tree_list)?;
    let nwk = gene_tree.to_newick().map_err(|e| Error::Other(e))?;
    Ok(nwk + ";")
}

/// Export a gene tree to RecPhyloXML format.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @return RecPhyloXML string
/// @export
#[extendr]
fn gene_tree_to_xml_r(gene_tree_list: List) -> Result<String> {
    let (gene_tree, species_tree, node_mapping, event_mapping) = rlist_to_genetree(&gene_tree_list)?;

    let rec_tree = RecTree::new_owned(species_tree, gene_tree, node_mapping, event_mapping);
    Ok(rec_tree.to_xml())
}

/// Save a tree to a Newick file.
///
/// @param tree_list A tree list (species or gene tree)
/// @param filepath Path to save the file
/// @export
#[extendr]
fn save_newick_r(tree_list: List, filepath: &str) -> Result<()> {
    // Check for species_node field to determine if this is a gene tree
    // Gene trees have species_node, species trees do not
    let is_gene_tree = tree_list.dollar("species_node")
        .map(|robj| !robj.is_null())
        .unwrap_or(false);

    let newick = if is_gene_tree {
        gene_tree_to_newick_r(tree_list)?
    } else {
        tree_to_newick_r(tree_list)?
    };

    fs::write(filepath, newick)
        .map_err(|e| format!("Failed to write file: {}", e))?;
    Ok(())
}

/// Save a gene tree to RecPhyloXML format.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @param filepath Path to save the XML file
/// @export
#[extendr]
fn save_xml_r(gene_tree_list: List, filepath: &str) -> Result<()> {
    let xml = gene_tree_to_xml_r(gene_tree_list)?;
    fs::write(filepath, xml)
        .map_err(|e| format!("Failed to write file: {}", e))?;
    Ok(())
}

/// Save a gene tree to CSV format.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @param filepath Path to save the CSV file
/// @export
#[extendr]
fn save_csv_r(gene_tree_list: List, filepath: &str) -> Result<()> {
    let (gene_tree, species_tree, node_mapping, event_mapping) = rlist_to_genetree(&gene_tree_list)?;
    let rec_tree = RecTree::new_owned(species_tree, gene_tree, node_mapping, event_mapping);
    rec_tree.save_csv(filepath)
        .map_err(|e| extendr_api::Error::Other(format!("Failed to write CSV: {}", e)))?;
    Ok(())
}

/// Save birth-death events from a species tree to CSV format.
///
/// For trees from simulate_species_tree, saves the actual simulation events.
/// For trees from parse_newick, generates events by treating internal nodes
/// as speciations and leaves as extant species (no extinction events).
///
/// @param species_tree_list A species tree list from simulate_species_tree_r or parse_newick_r
/// @param filepath Path to save the CSV file
/// @export
#[extendr]
fn save_bd_events_csv_r(species_tree_list: List, filepath: &str) -> Result<()> {
    // Get the tree structure first (needed in both cases)
    let tree = rlist_to_flattree(&species_tree_list)?;

    // Try to get events from the tree (exists for simulated trees)
    let events = match species_tree_list.dollar("events") {
        Ok(events_robj) if !events_robj.is_null() => {
            // Tree has events (from simulate_species_tree)
            let events_list: List = events_robj
                .try_into()
                .map_err(|_| "Failed to parse events list")?;
            rlist_to_bd_events(&events_list)?
        }
        _ => {
            // No events (from parse_newick) - generate them from tree structure
            generate_events_from_tree(&tree).map_err(|e| Error::Other(e))?
        }
    };

    // Call the existing Rust function
    use crate::io::save_bd_events_to_csv;
    save_bd_events_to_csv(&events, &tree, filepath)
        .map_err(|e| format!("Failed to write CSV file: {}", e))?;

    Ok(())
}

/// Generate SVG visualization using thirdkind.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @param filepath Optional path to save the SVG (NULL to return string only)
/// @param open_browser Whether to open in browser (default FALSE)
/// @return SVG content as string
/// @export
#[extendr]
fn gene_tree_to_svg_r(gene_tree_list: List, filepath: Robj, open_browser: bool) -> Result<String> {
    let xml = gene_tree_to_xml_r(gene_tree_list)?;

    let temp_dir = std::env::temp_dir();
    let xml_path = temp_dir.join("rustree_temp.recphyloxml");
    let svg_path = temp_dir.join("rustree_temp.svg");

    fs::write(&xml_path, &xml)
        .map_err(|e| format!("Failed to write temp XML: {}", e))?;

    let mut cmd = Command::new("thirdkind");
    cmd.arg("-f").arg(&xml_path)
       .arg("-o").arg(&svg_path);

    if open_browser {
        cmd.arg("-b");
    }

    let output = cmd.output()
        .map_err(|e| format!("Failed to run thirdkind. Is it installed? (`cargo install thirdkind`)\nError: {}", e))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("thirdkind failed: {}", stderr).into());
    }

    let svg = fs::read_to_string(&svg_path)
        .map_err(|e| format!("Failed to read SVG output: {}", e))?;

    // Save to user-specified path if provided
    if !filepath.is_null() {
        let path = filepath.as_str().ok_or("filepath must be a string")?;
        fs::write(path, &svg)
            .map_err(|e| format!("Failed to write SVG file: {}", e))?;
    }

    // Cleanup
    let _ = fs::remove_file(&xml_path);
    let _ = fs::remove_file(&svg_path);

    Ok(svg)
}

/// Sample extant genes from a gene tree.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r
/// @return A new gene tree containing only extant genes
/// @export
#[extendr]
fn sample_extant_r(gene_tree_list: List) -> Result<List> {
    let (gene_tree, species_tree, node_mapping, event_mapping) = rlist_to_genetree(&gene_tree_list)?;

    // Find extant gene indices
    let extant_indices: std::collections::HashSet<usize> = gene_tree.nodes.iter()
        .enumerate()
        .filter(|(i, n)| {
            n.left_child.is_none()
            && n.right_child.is_none()
            && event_mapping[*i] == Event::Leaf
        })
        .map(|(i, _)| i)
        .collect();

    if extant_indices.is_empty() {
        return Err("No extant genes to sample".into());
    }

    let (sampled_tree, gene_old_to_new) = extract_induced_subtree(&gene_tree, &extant_indices)
        .ok_or("Failed to extract induced subtree")?;

    // Species tree is unchanged, use identity mapping
    let species_identity: Vec<Option<usize>> = (0..species_tree.nodes.len())
        .map(Some)
        .collect();

    let (new_node_mapping, new_event_mapping) = remap_gene_tree_indices(
        &sampled_tree, &gene_old_to_new,
        &node_mapping, &event_mapping,
        &species_identity,
    ).map_err(|e| format!("Remap failed: {}", e))?;

    rectree_to_rlist(&species_tree, &sampled_tree, &new_node_mapping, &new_event_mapping)
}

/// Sample leaves and filter all trees accordingly.
///
/// This function samples a subset of species from the species tree and automatically
/// filters the gene tree to keep only genes that map to the sampled species. The
/// reconciliation mappings are preserved using an LCA-based approach.
///
/// @param gene_tree_list A gene tree list from simulate_dtl_r or parse_recphyloxml_r
/// @param species_leaf_names Character vector of species leaf names to keep
/// @return A new gene tree with sampled species and gene trees
/// @export
#[extendr]
fn sample_leaves_r(gene_tree_list: List, species_leaf_names: Robj) -> Result<List> {
    // Convert species_leaf_names to Vec<String>
    let species_names: Vec<String> = if species_leaf_names.is_string() {
        species_leaf_names.as_str_vector()
            .ok_or("Failed to convert species_leaf_names to string vector")?
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        return Err("species_leaf_names must be a character vector".into());
    };

    // Parse gene tree list
    let (gene_tree, species_tree, node_mapping, event_mapping) = rlist_to_genetree(&gene_tree_list)?;

    // Build a single-tree GeneForest and sample
    let forest = GeneForest::with_gene_trees(species_tree, vec![(gene_tree, node_mapping, event_mapping)]);
    let sampled_forest = forest.sample_leaves(&species_names)
        .map_err(|e| format!("Failed to sample leaves: {}", e))?;
    let sampled = sampled_forest.get(0)
        .ok_or("No gene trees in sampled result")?;

    // Convert back to R list
    rectree_to_rlist(
        &sampled.species_tree,
        &sampled.gene_tree,
        &sampled.node_mapping,
        &sampled.event_mapping,
    )
}

/// Extract an induced subtree keeping only specified leaves.
///
/// This function works on both species trees and gene trees.
/// For species trees, returns a pruned tree with only the specified leaves.
/// For gene trees, returns the pruned gene tree (reconciliation info is not preserved).
///
/// @param tree_list A tree list (species tree or gene tree)
/// @param leaf_names Character vector of leaf names to keep
/// @return A new tree containing only the specified leaves and their MRCAs
/// @export
#[extendr]
fn extract_induced_subtree_by_names_r(tree_list: List, leaf_names: Robj) -> Result<List> {
    // Convert leaf_names to Vec<String>
    let names: Vec<String> = leaf_names.as_str_vector()
        .ok_or("leaf_names must be a character vector")?
        .iter()
        .map(|s| s.to_string())
        .collect();

    if names.is_empty() {
        return Err("leaf_names cannot be empty".into());
    }

    // Extract the tree structure (works for both species and gene trees)
    let tree = rlist_to_flattree(&tree_list)?;

    // Extract the induced subtree
    let (induced_tree, _) = extract_induced_subtree_by_names(&tree, &names)
        .ok_or("Failed to extract induced subtree (no matching leaves found)")?;

    // Return as a simple tree (no reconciliation info)
    Ok(flattree_to_rlist(&induced_tree))
}

// Conversion helpers are in conversions.rs

/// Compute all pairwise distances between nodes in a tree.
///
/// @param tree_list A tree list from simulate_species_tree_r or parse_newick_r
/// @param distance_type Type of distance: "topological" (number of edges) or "metric" (sum of branch lengths)
/// @param leaves_only If TRUE, only compute distances between leaf nodes (default TRUE)
/// @return A list with three vectors: node1, node2, distance (suitable for conversion to data.frame)
/// @export
#[extendr]
fn pairwise_distances_r(tree_list: List, distance_type: &str, leaves_only: bool) -> Result<List> {
    let tree = rlist_to_flattree(&tree_list)?;

    let dist_type = crate::bindings_common::parse_distance_type(distance_type)
        .map_err(|e| Error::Other(e))?;

    let distances = tree.pairwise_distances(dist_type, leaves_only)
        .map_err(|e| format!("Failed to compute pairwise distances: {}", e))?;

    let node1: Vec<String> = distances.iter().map(|d| d.node1.to_string()).collect();
    let node2: Vec<String> = distances.iter().map(|d| d.node2.to_string()).collect();
    let dist: Vec<f64> = distances.iter().map(|d| d.distance).collect();

    Ok(list!(
        node1 = node1,
        node2 = node2,
        distance = dist
    ))
}

/// Save pairwise distances to a CSV file.
///
/// @param tree_list A tree list from simulate_species_tree_r or parse_newick_r
/// @param filepath Path to save the CSV file
/// @param distance_type Type of distance: "topological" or "metric"
/// @param leaves_only If TRUE, only compute distances between leaf nodes (default TRUE)
/// @export
#[extendr]
fn save_pairwise_distances_csv_r(tree_list: List, filepath: &str, distance_type: &str, leaves_only: bool) -> Result<()> {
    use crate::metric_functions::PairwiseDistance;

    let tree = rlist_to_flattree(&tree_list)?;

    let dist_type = crate::bindings_common::parse_distance_type(distance_type)
        .map_err(|e| Error::Other(e))?;

    let distances = tree.pairwise_distances(dist_type, leaves_only)
        .map_err(|e| format!("Failed to compute pairwise distances: {}", e))?;

    // Write to CSV
    let mut csv = String::from(PairwiseDistance::csv_header());
    csv.push('\n');
    for d in &distances {
        csv.push_str(&d.to_csv_row());
        csv.push('\n');
    }

    fs::write(filepath, csv)
        .map_err(|e| format!("Failed to write CSV file: {}", e))?;

    Ok(())
}

/// Parse a RecPhyloXML file into a gene tree with reconciliation information.
///
/// @param filepath Path to the RecPhyloXML file (e.g., from ALERax)
/// @return A gene tree list with reconciliation information
/// @export
#[extendr]
fn parse_recphyloxml_r(filepath: &str) -> Result<List> {
    let rec_tree = RecTree::from_xml_file(filepath)
        .map_err(|e| format!("Failed to parse RecPhyloXML: {}", e))?;

    // Convert to R list using the rectree_to_rlist function
    rectree_to_rlist(
        &rec_tree.species_tree,
        &rec_tree.gene_tree,
        &rec_tree.node_mapping,
        &rec_tree.event_mapping,
    )
}

// ---------------------------------------------------------------------------
// Streaming / lazy-iterator DTL functions
// ---------------------------------------------------------------------------

/// Helper: parse a Newick string into a FlatTree with depths assigned.
fn parse_newick_to_flattree(newick: &str) -> Result<FlatTree> {
    use crate::newick::parse_newick;

    let mut nodes = parse_newick(newick)
        .map_err(|e| format!("Failed to parse Newick: {}", e))?;
    let mut root = nodes.pop()
        .ok_or("No tree found in Newick string")?;
    root.zero_root_length();
    root.assign_depths(0.0);
    let mut tree = root.to_flat_tree();
    tree.assign_depths();
    Ok(tree)
}

/// Helper: create an StdRng from an optional R seed value.
fn make_rng(seed: &Robj) -> Result<StdRng> {
    if seed.is_null() || seed.is_na() {
        Ok(StdRng::from_entropy())
    } else {
        match seed.as_integer() {
            Some(s) => Ok(StdRng::seed_from_u64(s as u64)),
            None => Err("seed must be an integer".into()),
        }
    }
}

/// Helper: extract optional transfer_alpha from an Robj.
fn extract_alpha(transfer_alpha: &Robj) -> Result<Option<f64>> {
    if transfer_alpha.is_null() || transfer_alpha.is_na() {
        Ok(None)
    } else {
        Ok(Some(transfer_alpha.as_real().ok_or("transfer_alpha must be a number")?))
    }
}

/// Helper: extract optional replacement_transfer from an Robj.
fn extract_replacement(replacement_transfer: &Robj) -> Result<Option<f64>> {
    if replacement_transfer.is_null() || replacement_transfer.is_na() {
        Ok(None)
    } else {
        Ok(Some(replacement_transfer.as_real().ok_or("replacement_transfer must be a number")?))
    }
}

/// Stream DTL gene trees (per-gene model) to RecPhyloXML files.
///
/// Simulates `n` gene trees one at a time using the per-gene-copy DTL model
/// and writes each tree as a RecPhyloXML file in `output_dir`. Files are
/// named `gene_0000.xml`, `gene_0001.xml`, etc.
///
/// This is memory-efficient for large batches because only one tree is held
/// in memory at a time.
///
/// @param newick Species tree in Newick format
/// @param n Number of gene trees to simulate
/// @param lambda_d Duplication rate
/// @param lambda_t Transfer rate
/// @param lambda_l Loss rate
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, retry until a tree with extant genes is produced
/// @param seed Optional random seed for reproducibility
/// @param output_dir Directory to write XML files to (created if it doesn't exist)
/// @return The output directory path
/// @export
#[extendr]
fn simulate_dtl_stream_xml_r(
    newick: &str,
    n: i32,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
    output_dir: &str,
) -> Result<String> {
    if n <= 0 {
        return Err("Number of trees must be positive".into());
    }
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = parse_newick_to_flattree(newick)?;
    let mut rng = make_rng(&seed)?;
    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    simulate_dtl_iter(
        &species_tree, origin_species,
        lambda_d, lambda_t, lambda_l,
        alpha, repl,
        n as usize, require_extant, &mut rng,
    )
    .map_err(|e| Error::Other(e))?
    .save_xml(output_dir)
    .map_err(|e| Error::Other(e))?;

    Ok(output_dir.to_string())
}

/// Stream DTL gene trees (per-gene model) to Newick files.
///
/// Simulates `n` gene trees one at a time using the per-gene-copy DTL model
/// and writes each gene tree in Newick format to `output_dir`. Files are
/// named `gene_0000.nwk`, `gene_0001.nwk`, etc.
///
/// Note: Newick format does not preserve reconciliation information (species
/// mapping, events). Use `simulate_dtl_stream_xml_r` for full reconciled trees.
///
/// @param newick Species tree in Newick format
/// @param n Number of gene trees to simulate
/// @param lambda_d Duplication rate
/// @param lambda_t Transfer rate
/// @param lambda_l Loss rate
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, retry until a tree with extant genes is produced
/// @param seed Optional random seed for reproducibility
/// @param output_dir Directory to write Newick files to (created if it doesn't exist)
/// @return The output directory path
/// @export
#[extendr]
fn simulate_dtl_stream_newick_r(
    newick: &str,
    n: i32,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
    output_dir: &str,
) -> Result<String> {
    if n <= 0 {
        return Err("Number of trees must be positive".into());
    }
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = parse_newick_to_flattree(newick)?;
    let mut rng = make_rng(&seed)?;
    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    simulate_dtl_iter(
        &species_tree, origin_species,
        lambda_d, lambda_t, lambda_l,
        alpha, repl,
        n as usize, require_extant, &mut rng,
    )
    .map_err(|e| Error::Other(e))?
    .save_newick(output_dir)
    .map_err(|e| Error::Other(e))?;

    Ok(output_dir.to_string())
}

/// Stream DTL gene trees (per-species/Zombi model) to RecPhyloXML files.
///
/// Simulates `n` gene trees one at a time using the Zombi-style per-species
/// DTL model and writes each tree as a RecPhyloXML file in `output_dir`.
/// Files are named `gene_0000.xml`, `gene_0001.xml`, etc.
///
/// In the per-species model, the event rate is proportional to the number of
/// SPECIES with gene copies, NOT the number of gene copies themselves.
///
/// @param newick Species tree in Newick format
/// @param n Number of gene trees to simulate
/// @param lambda_d Duplication rate per species per unit time
/// @param lambda_t Transfer rate per species per unit time
/// @param lambda_l Loss rate per species per unit time
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, retry until a tree with extant genes is produced
/// @param seed Optional random seed for reproducibility
/// @param output_dir Directory to write XML files to (created if it doesn't exist)
/// @return The output directory path
/// @export
#[extendr]
fn simulate_dtl_per_species_stream_xml_r(
    newick: &str,
    n: i32,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
    output_dir: &str,
) -> Result<String> {
    if n <= 0 {
        return Err("Number of trees must be positive".into());
    }
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = parse_newick_to_flattree(newick)?;
    let mut rng = make_rng(&seed)?;
    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    simulate_dtl_per_species_iter(
        &species_tree, origin_species,
        lambda_d, lambda_t, lambda_l,
        alpha, repl,
        n as usize, require_extant, &mut rng,
    )
    .map_err(|e| Error::Other(e))?
    .save_xml(output_dir)
    .map_err(|e| Error::Other(e))?;

    Ok(output_dir.to_string())
}

/// Stream DTL gene trees (per-species/Zombi model) to Newick files.
///
/// Simulates `n` gene trees one at a time using the Zombi-style per-species
/// DTL model and writes each gene tree in Newick format to `output_dir`.
/// Files are named `gene_0000.nwk`, `gene_0001.nwk`, etc.
///
/// Note: Newick format does not preserve reconciliation information (species
/// mapping, events). Use `simulate_dtl_per_species_stream_xml_r` for full
/// reconciled trees.
///
/// @param newick Species tree in Newick format
/// @param n Number of gene trees to simulate
/// @param lambda_d Duplication rate per species per unit time
/// @param lambda_t Transfer rate per species per unit time
/// @param lambda_l Loss rate per species per unit time
/// @param transfer_alpha Optional distance decay for assortative transfers (NULL = uniform)
/// @param replacement_transfer Optional probability of replacement transfer (NULL = no replacement)
/// @param require_extant If TRUE, retry until a tree with extant genes is produced
/// @param seed Optional random seed for reproducibility
/// @param output_dir Directory to write Newick files to (created if it doesn't exist)
/// @return The output directory path
/// @export
#[extendr]
fn simulate_dtl_per_species_stream_newick_r(
    newick: &str,
    n: i32,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    transfer_alpha: Robj,
    replacement_transfer: Robj,
    require_extant: bool,
    seed: Robj,
    output_dir: &str,
) -> Result<String> {
    if n <= 0 {
        return Err("Number of trees must be positive".into());
    }
    crate::bindings_common::validate_dtl_rates(lambda_d, lambda_t, lambda_l)
        .map_err(|e| Error::Other(e))?;

    let species_tree = parse_newick_to_flattree(newick)?;
    let mut rng = make_rng(&seed)?;
    let alpha = extract_alpha(&transfer_alpha)?;
    let repl = extract_replacement(&replacement_transfer)?;

    let origin_species = species_tree.root;
    simulate_dtl_per_species_iter(
        &species_tree, origin_species,
        lambda_d, lambda_t, lambda_l,
        alpha, repl,
        n as usize, require_extant, &mut rng,
    )
    .map_err(|e| Error::Other(e))?
    .save_newick(output_dir)
    .map_err(|e| Error::Other(e))?;

    Ok(output_dir.to_string())
}

// DTL event conversion helpers are in conversions.rs

/// Compute induced transfers by projecting transfer events onto a sampled subtree.
///
/// Given a species tree, a set of sampled leaf names, and DTL events from
/// simulation, this function projects each transfer's donor and recipient
/// species onto the sampled subtree.
///
/// @param species_tree_list The complete species tree
/// @param sampled_leaf_names Character vector of species leaf names to keep
/// @param dtl_events_list DTL events list (from attr(gene_tree, "dtl_events"))
/// @return A data.frame with columns: time, gene_id, from_complete, to_complete,
///         from_sampled, to_sampled
/// @export
#[extendr]
fn induced_transfers_r(
    species_tree_list: List,
    sampled_leaf_names: Robj,
    dtl_events_list: List,
) -> Result<List> {
    let species_tree = rlist_to_flattree(&species_tree_list)?;

    let names: Vec<String> = sampled_leaf_names.as_str_vector()
        .ok_or("sampled_leaf_names must be a character vector")?
        .iter().map(|s| s.to_string()).collect();

    if names.is_empty() {
        return Err("sampled_leaf_names cannot be empty".into());
    }

    let events = rlist_to_dtl_events(&dtl_events_list, &species_tree)?;

    let result = induced_transfers(&species_tree, &names, &events);

    // Build the sampled tree to resolve sampled indices to names
    let (sampled_tree, _) = extract_induced_subtree_by_names(&species_tree, &names)
        .ok_or("Failed to extract sampled subtree")?;

    let n = result.len();
    let mut times: Vec<f64> = Vec::with_capacity(n);
    let mut gene_ids: Vec<i32> = Vec::with_capacity(n);
    let mut from_complete: Vec<String> = Vec::with_capacity(n);
    let mut to_complete: Vec<String> = Vec::with_capacity(n);
    let mut from_sampled: Vec<Rstr> = Vec::with_capacity(n);
    let mut to_sampled: Vec<Rstr> = Vec::with_capacity(n);

    for it in &result {
        times.push(it.time);
        gene_ids.push(it.gene_id as i32);
        from_complete.push(species_tree.nodes[it.from_species_complete].name.clone());
        to_complete.push(species_tree.nodes[it.to_species_complete].name.clone());
        from_sampled.push(match it.from_species_sampled {
            Some(idx) => Rstr::from(sampled_tree.nodes[idx].name.clone()),
            None => Rstr::na(),
        });
        to_sampled.push(match it.to_species_sampled {
            Some(idx) => Rstr::from(sampled_tree.nodes[idx].name.clone()),
            None => Rstr::na(),
        });
    }

    Ok(list!(
        time = times,
        gene_id = gene_ids,
        from_complete = from_complete,
        to_complete = to_complete,
        from_sampled = from_sampled,
        to_sampled = to_sampled
    ))
}

/// Name unnamed internal nodes as `internal0`, `internal1`, etc.
///
/// Assigns names to internal nodes (nodes with children) that have empty names.
/// Errors if any existing node name starts with "internal".
///
/// @param tree_list A tree list from parse_newick or simulate_species_tree
/// @return The tree with internal nodes named
/// @export
#[extendr]
fn name_internal_nodes_r(tree_list: List) -> Result<List> {
    let mut tree = rlist_to_flattree(&tree_list)?;
    if tree.nodes.iter().any(|n| n.name.starts_with("internal")) {
        return Err(Error::Other(
            "Cannot auto-name internal nodes: at least one node already has a name starting with \"internal\"".to_string()
        ));
    }
    tree.name_internal_nodes();
    Ok(flattree_to_rlist(&tree))
}

// Macro to generate exports
extendr_module! {
    mod rustree;
    fn simulate_species_tree_r;
    fn parse_newick_r;
    fn name_internal_nodes_r;
    fn tree_to_newick_r;
    fn tree_num_leaves_r;
    fn tree_leaf_names_r;
    fn simulate_dtl_r;
    fn simulate_dtl_batch_r;
    fn simulate_dtl_per_species_r;
    fn simulate_dtl_per_species_batch_r;
    fn gene_tree_num_extant_r;
    fn gene_tree_to_newick_r;
    fn gene_tree_to_xml_r;
    fn save_newick_r;
    fn save_xml_r;
    fn save_csv_r;
    fn save_bd_events_csv_r;
    fn gene_tree_to_svg_r;
    fn sample_extant_r;
    fn sample_leaves_r;
    fn extract_induced_subtree_by_names_r;
    fn pairwise_distances_r;
    fn save_pairwise_distances_csv_r;
    fn parse_recphyloxml_r;
    fn induced_transfers_r;
    fn simulate_dtl_stream_xml_r;
    fn simulate_dtl_stream_newick_r;
    fn simulate_dtl_per_species_stream_xml_r;
    fn simulate_dtl_per_species_stream_newick_r;
}