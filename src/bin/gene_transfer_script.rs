use rustree::newick::newick::{NewickParser, newick_to_tree, Rule};
use rustree::node::FlatTree;
use rustree::sampling::spr;
use pest::Parser;
use std::env;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::Path;
use csv::Writer;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use regex::Regex;
use std::collections::HashMap;
use rand::distributions::WeightedIndex; // <-- already imported
use rand::seq::SliceRandom;

/// Computes the interval intensity by summing transfer rates of all species present
/// in each time interval (given by the contemporaneity vector).
///
/// # Arguments
/// * `contemporaneity` - A slice of vectors where each inner vector contains node indices for an interval.
/// * `transfer_rates` - A slice of f64 values with the transfer rate for each species (indexed by node index).
///
/// # Returns
/// A vector of intensities for each interval.
fn compute_interval_intensity(contemporaneity: &[Vec<usize>], transfer_rates: &[f64]) -> Vec<f64> {
    contemporaneity.iter().map(|species_indices| {
        species_indices.iter().map(|&i| transfer_rates.get(i).cloned().unwrap_or(0.0)).sum()
    }).collect()
}

/// Computes the cumulative distribution function (CDF) from a vector of intervals and a vector of transfer intensities.
/// The CDF is normalized to sum to 1.
/// 
/// # Arguments
/// * `intervals` - A vector of time intervals.
/// * `transfer_intensity` - A vector of transfer intensities for each interval.
/// These transfer intensities are the sum of transfer rates for all species present in the interval.
/// 
/// # Returns
/// A vector of CDF values.
fn make_cdf(intervals: Vec<f64>, transfer_intensity: Vec<f64>) -> Vec<f64> {
    let n = intervals.len();
    let mut cdf: Vec<f64> = Vec::new();
    cdf.push(intervals[0] * transfer_intensity[0]);
    for i in 1..n {
        cdf.push(cdf[i - 1] + intervals[i] * transfer_intensity[i]);
    }
    let total_value = cdf[n - 1];
    for i in 0..n {
        cdf[i] = cdf[i] / total_value;
    }
    cdf
}

fn choose_from_cdf(cdf: &Vec<f64>, depths: &Vec<f64>, rng: &mut StdRng) -> (f64, usize) {
    let r: f64 = rng.gen_range(0.0..1.0);
    let index = match cdf.binary_search_by(|&probe| probe.partial_cmp(&r).unwrap_or(std::cmp::Ordering::Less)) {
        Ok(index) => index,
        Err(index) => {
            if index == cdf.len() {
                panic!("Random value exceeds CDF range. CDF: {:?}", cdf);
            } else {
                index
            }
        }
    };
    let depth = (r - cdf[index - 1]) / (cdf[index] - cdf[index - 1]) * (depths[index] - depths[index - 1]) + depths[index - 1];
    (depth, index)
}

/// Builds a WeightedIndex from species indices using their transfer_rates.
/// This avoids recomputing the weights each time we sample.
fn build_weighted_index(vec: &Vec<usize>, transfer_rates: &[f64]) -> WeightedIndex<f64> {
    let weights: Vec<f64> = vec.iter().map(|&i| transfer_rates[i]).collect();
    WeightedIndex::new(&weights).unwrap()
}


fn random_pair_with_weights(
    contemporaneous_species: &[usize],
    weighted: &WeightedIndex<f64>,
    rng: &mut StdRng,
) -> Option<(usize, usize)> {
    let n = contemporaneous_species.len();
    if n < 2 {
        return None;
    }
    let donor_idx = rng.sample(weighted);
    let donor = contemporaneous_species[donor_idx];
    // Efficiently choose a receiver by sampling an index in range 0..(n-1)
    let rand_idx = rng.gen_range(0..(n - 1));
    // Adjust index to skip the donor index
    let receiver_idx = if rand_idx >= donor_idx { rand_idx + 1 } else { rand_idx };
    let receiver = contemporaneous_species[receiver_idx];
    Some((donor, receiver))
}

fn generate_transfers(
    n_transfers: usize,
    contemporaneity: &Vec<Vec<usize>>,
    weighted_indices: &Vec<Option<WeightedIndex<f64>>>,
    cdf: &Vec<f64>,
    depths: &Vec<f64>,
    rng: &mut StdRng,
    use_weighed_sampling: bool,
) -> Vec<(usize, usize, f64)> {
    let mut transfers: Vec<(usize, usize, f64)> = Vec::new();
    for _ in 0..n_transfers {
        let (depth, index) = choose_from_cdf(cdf, depths, rng);
        let contemporaneous_species = &contemporaneity[index];
        let weighted = weighted_indices[index].as_ref().expect("No weighted index available");
        // To ensure we get the same results as before, we need to use the old function when not using weighted sampling.
        if use_weighed_sampling {
            if let Some((donor, receiver)) = random_pair_with_weights(contemporaneous_species, weighted, rng) {
                transfers.push((donor, receiver, depth));
            }
        } else {
            if let Some((donor, receiver)) = random_pair(contemporaneous_species, rng) {
                transfers.push((donor, receiver, depth));
            }
        }
    }
    transfers.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    transfers
}

/// Writes transfer events to a CSV file.
///
/// # Arguments
///
/// * `flat_tree` - The flat representation of the tree.
/// * `transfers` - A vector of transfers to write.
/// * `filename` - The path to the output CSV file.
///
/// # Returns
///
/// An `io::Result<()>` indicating success or failure.
fn make_transfers_csv(
    flat_tree: &mut FlatTree,
    transfers: &Vec<(usize, usize, f64)>,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = Writer::from_path(filename)?;

    // Writing the header
    writer.write_record(&["Donor", "Recipient", "Depth"])?;

    // Write transfer data
    for transfer in transfers.iter() {
        let donor_name = &flat_tree[transfer.0].name;
        let recipient_name = &flat_tree[transfer.1].name;
        let depth = transfer.2;

        // Writing the data rows
        writer.write_record(&[donor_name, recipient_name, &depth.to_string()])?;
    }

    writer.flush()?;
    Ok(())
}


fn one_gene_sim_to_string(
    flat_tree: &mut FlatTree,
    n_transfers: usize,
    contemporaneity: &Vec<Vec<usize>>,
    cdf: &Vec<f64>,
    depths: &Vec<f64>,
    gene_index: u32,
    output_dir: &str,
    rng: &mut StdRng,
    weighted_indices: &Vec<Option<WeightedIndex<f64>>>,
    use_weighed_sampling: bool,
) -> Result<String, io::Error> {
    // Ensure the output directory exists
    fs::create_dir_all(output_dir)?;

    // Paths for transfers and genes directories
    let transfers_dir = Path::new(output_dir).join("transfers");
    let genes_dir = Path::new(output_dir).join("genes");

    // Ensure the transfers and genes directories exist
    fs::create_dir_all(&transfers_dir)?;
    fs::create_dir_all(&genes_dir)?;

    // Create transfers
    let transfers = generate_transfers(n_transfers,
        contemporaneity,
        weighted_indices,
        cdf,
        depths,
        rng,
        use_weighed_sampling);

    // Export them to a CSV file in the transfers directory
    let transfer_filename = transfers_dir.join(format!("transfers_{}.csv", gene_index));
    make_transfers_csv(flat_tree, &transfers, transfer_filename.to_str().unwrap())
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    
    // Create the gene tree using SPR moves
    for transfer in transfers {
        spr(flat_tree, transfer.0, transfer.1, transfer.2);
    }

    // Convert flat tree to node tree
    let mut reconstructed_tree = flat_tree.to_node();

    // Update lengths based on depths
    let root_depth = reconstructed_tree
        .depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    reconstructed_tree.depths_to_lengths(root_depth);

    // Convert the tree to Newick format
    let reconstructed_newick = (&reconstructed_tree).to_newick() + ";";

    // Save the gene tree as a Newick string in the genes directory
    let gene_filename = genes_dir.join(format!("gene_{}.nwk", gene_index));
    let mut gene_file = File::create(gene_filename)?;
    gene_file.write_all(reconstructed_newick.as_bytes())?;

    Ok(reconstructed_newick)
}

fn main() {
    // Updated usage: now 5 or 6 arguments are accepted.
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 && args.len() != 6 {
        eprintln!(
            "Usage: {} <path_to_nwk_file> <output_directory> <path_to_transfers_file> <rng_seed> [transfer_rate_file]",
            args[0]
        );
        return;
    }

    let nwk_file_path = &args[1];
    let output_dir = &args[2];
    let transfers_file_path = &args[3];
    let rng_seed_str = &args[4];

    // Convert rng_seed_str to u64
    let seed = rng_seed_str.parse::<u64>().unwrap_or_else(|e| {
        eprintln!("Error parsing RNG seed: {}", e);
        std::process::exit(1);
    });
    let mut rng = StdRng::seed_from_u64(seed);

    // Ensure the output directory exists
    fs::create_dir_all(output_dir).expect("Failed to create output directory");


    let content = fs::read_to_string(nwk_file_path).expect("Failed to read the file");
    let re = Regex::new(r"\s+").unwrap(); // Matches any whitespace characters
    let sanitized_content = re.replace_all(&content, "");
    let trees: Vec<String> = sanitized_content
        .split(';')
        .filter_map(|tree| {
            let trimmed = tree.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string() + ";")
            }
        })
        .collect();

    let transfers_content = fs::read_to_string(transfers_file_path).expect("Failed to read transfers file");
    let n_transfers: Vec<usize> = transfers_content
        .trim()
        .split(',')
        .filter_map(|n| n.parse().ok())
        .collect();

    // Normally, this loop only processes a single tree.
    for (index, tree) in trees.iter().enumerate() {
        let tree_output_dir = Path::new(output_dir).join(format!("tree_{}", index));
        fs::create_dir_all(&tree_output_dir).expect("Failed to create tree output directory");

        let pairs = NewickParser::parse(Rule::newick, &tree).unwrap_or_else(|e| panic!("Failed to parse Newick tree."));
        for pair in pairs {
            let mut vec_trees = newick_to_tree(pair);
            let mut node_tree = vec_trees.pop().unwrap();
            node_tree.zero_root_length();
            node_tree.assign_depths(0.0);

            let mut flat_tree = node_tree.to_flat_tree();
            let depths = flat_tree.make_subdivision();
            let intervals = flat_tree.make_intervals();
            let contemporaneity = flat_tree.find_contemporaneity(&depths);
            
            let transfer_rates: Vec<f64> = if args.len() == 6 {
                // If a transfer rate file is provided, read it into a map.
                let rate_file = &args[5];
                let mut rate_map = HashMap::new();
                let mut rdr = csv::ReaderBuilder::new().has_headers(false)
                                  .from_path(rate_file)
                                  .expect("Failed to read transfer rate file");
                for result in rdr.records() {
                    let record = result.expect("CSV error");
                    if record.len() >= 2 {
                        rate_map.insert(record[0].trim().to_string(), record[1].trim().parse::<f64>().unwrap_or(1.0));
                    }
                }
                flat_tree.nodes.iter()
                    .map(|node| rate_map.get(&node.name).cloned().unwrap_or(1.0))
                    .collect()
            } else {
                // Else, use a uniform rate of 1.0 for every node.
                vec![1.0; flat_tree.nodes.len()]
            };

            // Default transfer rates are all 1.0
            let intensities = compute_interval_intensity(&contemporaneity, &transfer_rates);
            let cdf = make_cdf(intervals, intensities);
            
            // Precompute weighted indices once for the whole tree.
            let weighted_indices: Vec<Option<WeightedIndex<f64>>> = contemporaneity
                .iter()
                .map(|species_vec| {
                    if species_vec.len() >= 2 {
                        let weights: Vec<f64> = species_vec.iter().map(|&i| transfer_rates[i]).collect();
                        Some(WeightedIndex::new(&weights).unwrap())
                    } else {
                        None
                    }
                })
                .collect();
            
            



            let mut copied_flat_tree = flat_tree.clone();
            // Updated: pass weighted_indices to create_many_genes so generate_transfers can be called repeatedly.
            create_many_genes(
                &mut copied_flat_tree,
                n_transfers.clone(),
                &contemporaneity,
                &weighted_indices, // <-- new parameter
                &cdf,
                &depths,
                &tree_output_dir.to_string_lossy().to_string(),
                &transfer_rates,
                &mut rng,
            );
        }
    }
}