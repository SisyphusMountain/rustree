use rustree::sampling::{remove_node, find_deepest_nodes};
use rustree::node::TraversalOrder;
use rustree::newick::newick::{node_to_newick, NewickParser, newick_to_tree, Rule};
use std::env;
use std::path::Path;
use std::fs::{self, File};
use std::io::{self, Write};
use pest::Parser;


fn main()-> Result<(), Box<dyn std::error::Error>> {
    // Read the arguments
    let args: Vec<String> = env::args().collect();
    // This script takes the n most recent nodes, samples them from a tree, and returns the sampled tree.
    // If we know the species tree has n extant nodes, we can sample the n most recent nodes to get the extant species tree.
    // Ensure the correct number of arguments are provided
    if args.len() < 4 || args.len() > 5 {
        eprintln!(
            "Usage: {} <species_tree_path> <n_extant_nodes> <output_dir> [--verbose]",
            args[0]
        );
        eprintln!("Received arguments: {:?}", args);
        panic!("Error with the input arguments! See error above.");
    }

    let verbose = args.len() == 5 && args[4] == "--verbose";

    let species_tree_path = &args[1];
    let n_extant = match args[2].parse::<usize>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!(
                "Error: n_extant_nodes must be an integer. Received: {}",
                args[2]
            );
            eprintln!("All arguments: {:?}", args);
            return Ok(()); // Changed from bare "return;"
        }
    };
    let output_dir = &args[3];





    let output_path = Path::new(output_dir);
    if !output_path.exists() {
        fs::create_dir_all(output_path)?;
    }

    // Open the species tree and convert it to a flat tree.
    let species_tree_str = fs::read_to_string(species_tree_path)?;
    let species_tree_str = species_tree_str.trim();

    // Parse the Newick string, mapping the pest error to io::Error
    let mut pairs = NewickParser::parse(Rule::newick, species_tree_str)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

    let mut node_tree = newick_to_tree(
        pairs
            .next()
            .expect("Error converting the Newick file")
    )
    .pop()
    .expect("Error: no tree found");

    // Assign depths
    node_tree.zero_root_length();
    node_tree.assign_depths(0.0);

    // Convert to FlatTree
    let mut flat_tree = node_tree.to_flat_tree();

    // Print species depths (leaves only) sorted in descending order if verbose.
    if verbose {
        let mut species_depths: Vec<_> = flat_tree.iter(TraversalOrder::PreOrder)
            .filter(|node| node.left_child.is_none() && node.right_child.is_none())
            .map(|node| (node.name.clone(), node.depth.unwrap_or(0.0)))
            .collect();
        species_depths.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("Error comparing depths"));
        println!("Species and depths (sorted descending):");
        for (name, depth) in species_depths {
            println!("  Species: {}, Depth: {}", name, depth);
        }
    }

    // Sample the leaves.
    let extant_leaves = find_deepest_nodes(&flat_tree, n_extant);

    // Construct the species tree by removing the unsampled leaves (those that are extinct)
    let leaves: Vec<usize> = flat_tree
        .iter(TraversalOrder::PreOrder)
        .enumerate()
        .filter(|(_, node)| node.left_child.is_none() && node.right_child.is_none())
        .map(|(i, _)| i)
        .collect();

    let leaves_to_be_removed: Vec<usize> = leaves
        .iter()
        .filter(|leaf| !extant_leaves.contains(leaf))
        .cloned()
        .collect();

    for index in &leaves_to_be_removed {
        remove_node(&mut flat_tree, *index);
    }


    // Convert the flat tree back to a Node tree.
    let mut reconstructed_tree = flat_tree.to_node();
    
    // Update lengths based on depths.
    let root_depth = reconstructed_tree
        .depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    reconstructed_tree.depths_to_lengths(root_depth);

    // Convert the tree to Newick format.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the species tree to the output directory.
    let species_filename = Path::new(output_dir).join("extant_species_tree.nwk");
    let mut species_file = File::create(species_filename)?;
    species_file.write_all(reconstructed_newick.as_bytes())?;

    // Return the Newick string and the lists of sampled and removed leaf names.
    let sampled_leaves_names: Vec<String> = extant_leaves
        .iter()
        .map(|i| flat_tree[*i].name.clone())
        .collect();
    let leaves_to_be_removed_names: Vec<String> = leaves_to_be_removed
        .iter()
        .map(|i| flat_tree[*i].name.clone())
        .collect();

    if verbose {
        println!("\n=== Species Tree Sampling Summary ===");
        println!("Number of species sampled: {}", sampled_leaves_names.len());
        println!("Sampled species: {:?}", sampled_leaves_names);
        println!("Number of species removed: {}", leaves_to_be_removed_names.len());
        println!("Removed species: {:?}", leaves_to_be_removed_names);
        println!("Resulting Newick tree: {}", reconstructed_newick);
    }

    let result: Result<(String, Vec<String>, Vec<String>), Box<dyn std::error::Error>> = Ok((reconstructed_newick, sampled_leaves_names, leaves_to_be_removed_names));
    match result {
        Ok((_, sampled_names, removed_names)) => {
            // You can use sampled_names and removed_names if needed
            println!("Sampled Leaves: {:?}", sampled_names);
            println!("Removed Leaves: {:?}", removed_names);
        }
        Err(e) => {
            eprintln!("Error during species tree sampling: {}", e);
            eprintln!("Species Tree Path: {}", species_tree_path);
            eprintln!("Number of Sampled Nodes: {}", n_extant);
            eprintln!("Output Directory: {}", output_dir);
        }
    };


    Ok(())
}
