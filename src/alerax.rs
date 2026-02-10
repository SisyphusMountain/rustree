/// ALERax integration module for gene tree reconciliation
///
/// This module provides functionality to call the external ALERax tool
/// for reconciling gene trees with species trees and parsing the results.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use crate::node::{FlatTree, RecTreeOwned};

/// Configuration for ALERax reconciliation
pub struct AleRaxConfig {
    pub species_tree_path: PathBuf,
    pub families: Vec<GeneFamily>,
    pub families_file_path: PathBuf,
    pub output_dir: PathBuf,
    pub num_samples: usize,
    pub model_parametrization: ModelType,
    pub seed: Option<u64>,
    pub alerax_path: String,
}

/// Model parametrization type for ALERax
#[derive(Clone, Copy)]
pub enum ModelType {
    PerFamily,
    Global,
}

impl ModelType {
    fn as_str(&self) -> &str {
        match self {
            ModelType::PerFamily => "PER-FAMILY",
            ModelType::Global => "GLOBAL",
        }
    }
}

/// Gene family definition
pub struct GeneFamily {
    pub name: String,
    pub gene_tree_path: PathBuf,
    pub gene_tree_newick: String,
}

/// Summary statistics for reconciliation events
#[derive(Clone, Debug)]
pub struct ReconciliationStatistics {
    /// Mean event counts across all samples
    pub mean_event_counts: EventCounts,
    /// Mean transfers between species pairs (source -> destination -> count)
    pub mean_transfers: HashMap<String, HashMap<String, f64>>,
    /// Mean events per species node (species -> event type -> count)
    pub events_per_species: HashMap<String, EventCounts>,
}

/// Event counts for reconciliation
#[derive(Clone, Debug, Default)]
pub struct EventCounts {
    pub speciations: f64,
    pub speciation_losses: f64,
    pub duplications: f64,
    pub duplication_losses: f64,
    pub transfers: f64,
    pub transfer_losses: f64,
    pub losses: f64,
    pub leaves: f64,
}

/// Result of ALERax reconciliation for one family
pub struct AleRaxFamilyResult {
    pub family_name: String,
    pub reconciled_trees: Vec<RecTreeOwned>,
    pub duplication_rate: f64,
    pub loss_rate: f64,
    pub transfer_rate: f64,
    pub likelihood: f64,
    pub statistics: ReconciliationStatistics,
}

/// Check if ALERax is installed and accessible
fn check_alerax_installed(alerax_path: &str) -> Result<(), String> {
    // ALERax doesn't support --version, so just try to run it with --help
    // If it fails to execute at all, we know it's not installed
    Command::new(alerax_path)
        .arg("--help")
        .output()
        .map_err(|_| format!(
            "ALERax not found at '{}'. Please install ALERax:
  conda install -c bioconda alerax
  or download from: https://github.com/BenoitMorel/AleRax",
            alerax_path
        ))?;
    Ok(())
}

/// Extract species name from gene name (format: species_geneid)
fn extract_species_from_gene_name(gene_name: &str) -> String {
    gene_name.split('_')
        .next()
        .unwrap_or(gene_name)
        .to_string()
}

/// Validate that gene tree leaves map to species tree leaves
pub fn validate_inputs(
    species_tree: &FlatTree,
    gene_trees: &[(String, String)],  // (family_name, newick_string)
) -> Result<(), String> {
    // Check species tree is valid
    if species_tree.nodes.is_empty() {
        return Err("Species tree is empty".to_string());
    }

    let species_leaf_names: HashSet<String> = species_tree.nodes.iter()
        .filter(|n| n.left_child.is_none() && n.right_child.is_none())
        .map(|n| n.name.clone())
        .collect();

    if species_leaf_names.is_empty() {
        return Err("Species tree has no leaves".to_string());
    }

    // Check each gene tree
    for (family_name, newick) in gene_trees {
        // Parse gene tree to check leaves
        let mut nodes = crate::newick::newick::parse_newick(newick)
            .map_err(|e| format!("Failed to parse gene tree '{}': {}", family_name, e))?;

        let root = nodes.pop()
            .ok_or_else(|| format!("Gene tree for family '{}' is empty", family_name))?;

        let gene_tree = root.to_flat_tree();

        // Extract gene leaf names and map to species
        let gene_leaf_species: HashSet<String> = gene_tree.nodes.iter()
            .filter(|n| n.left_child.is_none() && n.right_child.is_none())
            .map(|n| extract_species_from_gene_name(&n.name))
            .collect();

        // Check if gene leaves are subset of species leaves
        let unmapped: Vec<&String> = gene_leaf_species.iter()
            .filter(|name| !species_leaf_names.contains(*name))
            .collect();

        if !unmapped.is_empty() {
            return Err(format!(
                "Gene tree '{}' has leaves not in species tree: {:?}",
                family_name, unmapped
            ));
        }

        // Check minimum coverage (ALERax requires at least 4 species)
        if gene_leaf_species.len() < 4 {
            return Err(format!(
                "Gene tree '{}' covers only {} species (minimum 4 required for reconciliation)",
                family_name, gene_leaf_species.len()
            ));
        }
    }

    Ok(())
}

/// Create families.txt file in ALERax format
fn create_families_file(families: &[GeneFamily]) -> String {
    let mut content = String::from("[FAMILIES]\n");
    for family in families {
        content.push_str(&format!("- {}\n", family.name));
        content.push_str(&format!("gene_tree = {}\n", family.gene_tree_path.display()));
    }
    content
}

/// Prepare ALERax input files
fn prepare_alerax_input(config: &AleRaxConfig, species_tree_newick: &str) -> Result<(), String> {
    // Write species tree to file (ensure it ends with semicolon)
    let species_newick = if species_tree_newick.trim().ends_with(';') {
        species_tree_newick.to_string()
    } else {
        format!("{};", species_tree_newick.trim())
    };

    fs::write(&config.species_tree_path, species_newick)
        .map_err(|e| format!("Failed to write species tree file: {}", e))?;

    // Write each gene tree to file (ensure they end with semicolon)
    for family in &config.families {
        let gene_newick = if family.gene_tree_newick.trim().ends_with(';') {
            family.gene_tree_newick.clone()
        } else {
            format!("{};", family.gene_tree_newick.trim())
        };

        fs::write(&family.gene_tree_path, gene_newick)
            .map_err(|e| format!("Failed to write gene tree for family '{}': {}", family.name, e))?;
    }

    // Create families.txt file
    let families_txt = create_families_file(&config.families);
    fs::write(&config.families_file_path, families_txt)
        .map_err(|e| format!("Failed to write families file: {}", e))?;

    Ok(())
}

/// Build ALERax command
fn build_alerax_command(config: &AleRaxConfig) -> Command {
    let mut cmd = Command::new(&config.alerax_path);

    cmd.arg("-s").arg(&config.species_tree_path);
    cmd.arg("-f").arg(&config.families_file_path);
    cmd.arg("-p").arg(&config.output_dir);
    cmd.arg("--model-parametrization").arg(config.model_parametrization.as_str());
    cmd.arg("--gene-tree-samples").arg(config.num_samples.to_string());

    if let Some(seed) = config.seed {
        cmd.arg("--seed").arg(seed.to_string());
    }

    cmd
}

/// Parse ALERax rates file (D L T)
fn parse_rates_file(path: &Path) -> Result<(f64, f64, f64), String> {
    let content = fs::read_to_string(path)
        .map_err(|e| format!("Failed to read rates file: {}", e))?;

    let lines: Vec<&str> = content.lines().collect();
    if lines.len() < 2 {
        return Err("Invalid rates file format (expected at least 2 lines)".to_string());
    }

    // Line 1: "D L T"
    // Line 2: "0.098 1e-10 0.066"
    let values: Vec<&str> = lines[1].split_whitespace().collect();
    if values.len() != 3 {
        return Err(format!("Expected 3 rate values, found {}", values.len()));
    }

    let d = values[0].parse()
        .map_err(|_| format!("Failed to parse duplication rate: {}", values[0]))?;
    let l = values[1].parse()
        .map_err(|_| format!("Failed to parse loss rate: {}", values[1]))?;
    let t = values[2].parse()
        .map_err(|_| format!("Failed to parse transfer rate: {}", values[2]))?;

    Ok((d, l, t))
}

/// Parse ALERax per-family likelihoods file
fn parse_likelihoods_file(path: &Path) -> Result<HashMap<String, f64>, String> {
    let content = fs::read_to_string(path)
        .map_err(|e| format!("Failed to read likelihoods file: {}", e))?;

    let mut likelihoods = HashMap::new();
    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 {
            let family_name = parts[0].to_string();
            let likelihood = parts[1].parse()
                .map_err(|_| format!("Failed to parse likelihood for family {}", family_name))?;
            likelihoods.insert(family_name, likelihood);
        }
    }

    Ok(likelihoods)
}

/// Parse ALERax eventCounts file for a single sample
fn parse_event_counts_file(path: &Path) -> Result<EventCounts, String> {
    let content = fs::read_to_string(path)
        .map_err(|e| format!("Failed to read eventCounts file: {}", e))?;

    let mut counts = EventCounts::default();

    for line in content.lines() {
        let parts: Vec<&str> = line.split(':').collect();
        if parts.len() == 2 {
            let event_type = parts[0].trim();
            let count: f64 = parts[1].trim().parse()
                .map_err(|_| format!("Failed to parse count for {}", event_type))?;

            match event_type {
                "S" => counts.speciations = count,
                "SL" => counts.speciation_losses = count,
                "D" => counts.duplications = count,
                "DL" => counts.duplication_losses = count,
                "T" => counts.transfers = count,
                "TL" => counts.transfer_losses = count,
                "L" => counts.losses = count,
                "Leaf" => counts.leaves = count,
                _ => {} // Ignore unknown event types
            }
        }
    }

    Ok(counts)
}

/// Parse ALERax transfers file for a single sample
fn parse_transfers_file(path: &Path) -> Result<HashMap<String, HashMap<String, f64>>, String> {
    let content = fs::read_to_string(path)
        .map_err(|e| format!("Failed to read transfers file: {}", e))?;

    let mut transfers: HashMap<String, HashMap<String, f64>> = HashMap::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 3 {
            let source = parts[0].to_string();
            let destination = parts[1].to_string();
            let count: f64 = parts[2].parse()
                .map_err(|_| format!("Failed to parse transfer count in line: {}", line))?;

            transfers.entry(source)
                .or_insert_with(HashMap::new)
                .insert(destination, count);
        }
    }

    Ok(transfers)
}

/// Compute events per species from a reconciled tree
fn compute_events_per_species(rec_tree: &RecTreeOwned) -> HashMap<String, EventCounts> {
    use crate::node::Event;

    let mut species_events: HashMap<String, EventCounts> = HashMap::new();

    // Traverse all nodes in the reconciled tree
    for i in 0..rec_tree.gene_tree.nodes.len() {
        // Get the species node this gene node maps to
        let species_idx = rec_tree.node_mapping[i];
        let species_node = &rec_tree.species_tree.nodes[species_idx];
        let species_name = species_node.name.clone();

        // Get or create event counts for this species
        let counts = species_events.entry(species_name).or_insert_with(EventCounts::default);

        // Increment count based on event type
        match rec_tree.event_mapping[i] {
            Event::Speciation => counts.speciations += 1.0,
            Event::Duplication => counts.duplications += 1.0,
            Event::Transfer => counts.transfers += 1.0,
            Event::Loss => counts.losses += 1.0,
            Event::Leaf => counts.leaves += 1.0,
        }
    }

    species_events
}

/// Aggregate statistics across multiple samples
fn aggregate_statistics(
    output_dir: &Path,
    family_name: &str,
    num_samples: usize,
    reconciled_trees: &[RecTreeOwned],
) -> Result<ReconciliationStatistics, String> {
    let reconciliation_dir = output_dir.join("reconciliations").join("all");

    // Aggregate event counts across samples
    let mut total_counts = EventCounts::default();
    for sample_idx in 0..num_samples {
        let event_counts_path = reconciliation_dir.join(format!("{}_eventCounts_{}.txt", family_name, sample_idx));
        if !event_counts_path.exists() {
            continue;
        }

        let counts = parse_event_counts_file(&event_counts_path)?;
        total_counts.speciations += counts.speciations;
        total_counts.speciation_losses += counts.speciation_losses;
        total_counts.duplications += counts.duplications;
        total_counts.duplication_losses += counts.duplication_losses;
        total_counts.transfers += counts.transfers;
        total_counts.transfer_losses += counts.transfer_losses;
        total_counts.losses += counts.losses;
        total_counts.leaves += counts.leaves;
    }

    // Compute mean event counts
    let n_samples = reconciled_trees.len() as f64;
    let mean_event_counts = EventCounts {
        speciations: total_counts.speciations / n_samples,
        speciation_losses: total_counts.speciation_losses / n_samples,
        duplications: total_counts.duplications / n_samples,
        duplication_losses: total_counts.duplication_losses / n_samples,
        transfers: total_counts.transfers / n_samples,
        transfer_losses: total_counts.transfer_losses / n_samples,
        losses: total_counts.losses / n_samples,
        leaves: total_counts.leaves / n_samples,
    };

    // Aggregate transfers across samples
    let mut transfer_totals: HashMap<String, HashMap<String, f64>> = HashMap::new();
    for sample_idx in 0..num_samples {
        let transfers_path = reconciliation_dir.join(format!("{}_transfers_{}.txt", family_name, sample_idx));
        if !transfers_path.exists() {
            continue;
        }

        let transfers = parse_transfers_file(&transfers_path)?;
        for (source, destinations) in transfers {
            for (dest, count) in destinations {
                *transfer_totals.entry(source.clone())
                    .or_insert_with(HashMap::new)
                    .entry(dest.clone())
                    .or_insert(0.0) += count;
            }
        }
    }

    // Compute mean transfers
    let mut mean_transfers: HashMap<String, HashMap<String, f64>> = HashMap::new();
    for (source, destinations) in transfer_totals {
        for (dest, total) in destinations {
            mean_transfers.entry(source.clone())
                .or_insert_with(HashMap::new)
                .insert(dest, total / n_samples);
        }
    }

    // Compute events per species across all samples
    let mut species_event_totals: HashMap<String, EventCounts> = HashMap::new();
    for rec_tree in reconciled_trees {
        let species_events = compute_events_per_species(rec_tree);
        for (species, counts) in species_events {
            let total = species_event_totals.entry(species).or_insert_with(EventCounts::default);
            total.speciations += counts.speciations;
            total.speciation_losses += counts.speciation_losses;
            total.duplications += counts.duplications;
            total.duplication_losses += counts.duplication_losses;
            total.transfers += counts.transfers;
            total.transfer_losses += counts.transfer_losses;
            total.losses += counts.losses;
            total.leaves += counts.leaves;
        }
    }

    // Compute mean events per species
    let mut events_per_species: HashMap<String, EventCounts> = HashMap::new();
    for (species, totals) in species_event_totals {
        events_per_species.insert(species, EventCounts {
            speciations: totals.speciations / n_samples,
            speciation_losses: totals.speciation_losses / n_samples,
            duplications: totals.duplications / n_samples,
            duplication_losses: totals.duplication_losses / n_samples,
            transfers: totals.transfers / n_samples,
            transfer_losses: totals.transfer_losses / n_samples,
            losses: totals.losses / n_samples,
            leaves: totals.leaves / n_samples,
        });
    }

    Ok(ReconciliationStatistics {
        mean_event_counts,
        mean_transfers,
        events_per_species,
    })
}

/// Parse ALERax reconciliation samples
fn parse_alerax_results(
    output_dir: &Path,
    families: &[GeneFamily],
    num_samples: usize,
) -> Result<HashMap<String, AleRaxFamilyResult>, String> {
    let mut results = HashMap::new();

    // Parse likelihoods file once
    let likelihoods_path = output_dir.join("per_fam_likelihoods.txt");
    let likelihoods = parse_likelihoods_file(&likelihoods_path)?;

    for family in families {
        // Parse rates
        let rates_path = output_dir
            .join("model_parameters")
            .join(format!("{}_rates.txt", family.name));
        let (d, l, t) = parse_rates_file(&rates_path)?;

        // Parse reconciliation samples
        let reconciliation_dir = output_dir.join("reconciliations").join("all");
        let mut reconciled_trees = Vec::new();

        for sample_idx in 0..num_samples {
            let xml_path = reconciliation_dir
                .join(format!("{}_sample_{}.xml", family.name, sample_idx));

            // Skip if sample doesn't exist (ALERax might generate fewer samples)
            if !xml_path.exists() {
                continue;
            }

            // Use existing RecTreeOwned::from_xml_file to parse
            let rec_tree = RecTreeOwned::from_xml_file(xml_path.to_str().unwrap())
                .map_err(|e| format!("Failed to parse RecPhyloXML for family '{}' sample {}: {}",
                    family.name, sample_idx, e))?;
            reconciled_trees.push(rec_tree);
        }

        if reconciled_trees.is_empty() {
            return Err(format!("No reconciliation samples found for family '{}'", family.name));
        }

        let likelihood = likelihoods.get(&family.name)
            .copied()
            .ok_or_else(|| format!("No likelihood found for family '{}'", family.name))?;

        // Compute summary statistics
        let statistics = aggregate_statistics(output_dir, &family.name, num_samples, &reconciled_trees)?;

        results.insert(family.name.clone(), AleRaxFamilyResult {
            family_name: family.name.clone(),
            reconciled_trees,
            duplication_rate: d,
            loss_rate: l,
            transfer_rate: t,
            likelihood,
            statistics,
        });
    }

    Ok(results)
}

/// Run ALERax and return reconciliation results
pub fn run_alerax(
    config: AleRaxConfig,
    species_tree_newick: &str,
) -> Result<HashMap<String, AleRaxFamilyResult>, String> {
    // Check if ALERax is installed
    check_alerax_installed(&config.alerax_path)?;

    // Prepare input files
    prepare_alerax_input(&config, species_tree_newick)?;

    // Build and execute command
    let mut cmd = build_alerax_command(&config);

    // Show command for debugging
    let cmd_str = format!("{:?}", cmd);
    eprintln!("Running ALERax command: {}", cmd_str);

    let output = cmd.output()
        .map_err(|e| format!("Failed to execute ALERax: {}", e))?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    if !output.status.success() {
        return Err(format!(
            "ALERax failed with exit code {:?}:\nStderr:\n{}\nStdout:\n{}",
            output.status.code(), stderr, stdout
        ));
    }

    // Check if output directory was created
    if !config.output_dir.exists() {
        return Err(format!(
            "ALERax completed but output directory does not exist: {}\nStdout:\n{}\nStderr:\n{}",
            config.output_dir.display(), stdout, stderr
        ));
    }

    // Parse results
    parse_alerax_results(&config.output_dir, &config.families, config.num_samples)
        .map_err(|e| format!(
            "{}.\n\nALERax output directory: {}\nStdout:\n{}\nStderr:\n{}",
            e, config.output_dir.display(), stdout, stderr
        ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_species_from_gene_name() {
        assert_eq!(extract_species_from_gene_name("n5_20"), "n5");
        assert_eq!(extract_species_from_gene_name("species_A_gene123"), "species");
        assert_eq!(extract_species_from_gene_name("single"), "single");
    }

    #[test]
    fn test_create_families_file() {
        let families = vec![
            GeneFamily {
                name: "fam1".to_string(),
                gene_tree_path: PathBuf::from("/tmp/tree1.nwk"),
                gene_tree_newick: "".to_string(),
            },
            GeneFamily {
                name: "fam2".to_string(),
                gene_tree_path: PathBuf::from("/tmp/tree2.nwk"),
                gene_tree_newick: "".to_string(),
            },
        ];

        let content = create_families_file(&families);
        assert!(content.contains("[FAMILIES]"));
        assert!(content.contains("- fam1"));
        assert!(content.contains("gene_tree = /tmp/tree1.nwk"));
        assert!(content.contains("- fam2"));
        assert!(content.contains("gene_tree = /tmp/tree2.nwk"));
    }

    #[test]
    fn test_model_type_as_str() {
        assert_eq!(ModelType::PerFamily.as_str(), "PER-FAMILY");
        assert_eq!(ModelType::Global.as_str(), "GLOBAL");
    }
}
