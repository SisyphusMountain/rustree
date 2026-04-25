use rand::rngs::StdRng;
use rand::SeedableRng;
use rustree::bd::simulate_bd_tree_bwd;
use rustree::bindings_common::extract_extant_gene_tree;
use rustree::dtl::{count_extant_genes, simulate_dtl_batch};
use std::collections::HashSet;
use std::env;
use std::hint::black_box;
use std::time::Instant;

struct BenchmarkCase {
    species_leaves: usize,
    gene_trees: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    seed: u64,
}

fn main() {
    let cases = match parse_args() {
        Ok(cases) => cases,
        Err(err) => {
            eprintln!("{err}");
            std::process::exit(2);
        }
    };

    println!("Unrooted topology hash benchmark");
    println!("================================");
    println!("Gene trees are benchmarked on the extant induced subtree.");
    println!("DTL rates start conservative and scale only after tree sizes stay bounded.");
    println!();

    for case in cases {
        run_case(case);
    }
}

fn default_cases() -> Vec<BenchmarkCase> {
    vec![
        BenchmarkCase {
            species_leaves: 20,
            gene_trees: 200,
            lambda_d: 0.05,
            lambda_t: 0.02,
            lambda_l: 0.05,
            seed: 7,
        },
        BenchmarkCase {
            species_leaves: 40,
            gene_trees: 500,
            lambda_d: 0.10,
            lambda_t: 0.05,
            lambda_l: 0.10,
            seed: 11,
        },
        BenchmarkCase {
            species_leaves: 80,
            gene_trees: 500,
            lambda_d: 0.08,
            lambda_t: 0.04,
            lambda_l: 0.08,
            seed: 19,
        },
    ]
}

fn parse_args() -> Result<Vec<BenchmarkCase>, String> {
    let mut species_leaves = None;
    let mut gene_trees = None;
    let mut lambda_d = None;
    let mut lambda_t = None;
    let mut lambda_l = None;
    let mut seed = 42_u64;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        let value = args
            .next()
            .ok_or_else(|| format!("Missing value for argument '{arg}'"))?;
        match arg.as_str() {
            "--species" => {
                species_leaves = Some(
                    value
                        .parse::<usize>()
                        .map_err(|e| format!("Invalid --species value '{value}': {e}"))?,
                );
            }
            "--gene-trees" => {
                gene_trees = Some(
                    value
                        .parse::<usize>()
                        .map_err(|e| format!("Invalid --gene-trees value '{value}': {e}"))?,
                );
            }
            "--lambda-d" => {
                lambda_d = Some(
                    value
                        .parse::<f64>()
                        .map_err(|e| format!("Invalid --lambda-d value '{value}': {e}"))?,
                );
            }
            "--lambda-t" => {
                lambda_t = Some(
                    value
                        .parse::<f64>()
                        .map_err(|e| format!("Invalid --lambda-t value '{value}': {e}"))?,
                );
            }
            "--lambda-l" => {
                lambda_l = Some(
                    value
                        .parse::<f64>()
                        .map_err(|e| format!("Invalid --lambda-l value '{value}': {e}"))?,
                );
            }
            "--seed" => {
                seed = value
                    .parse::<u64>()
                    .map_err(|e| format!("Invalid --seed value '{value}': {e}"))?;
            }
            other => return Err(format!("Unknown argument '{other}'")),
        }
    }

    match (species_leaves, gene_trees, lambda_d, lambda_t, lambda_l) {
        (None, None, None, None, None) => Ok(default_cases()),
        (Some(species_leaves), Some(gene_trees), Some(lambda_d), Some(lambda_t), Some(lambda_l)) => {
            Ok(vec![BenchmarkCase {
                species_leaves,
                gene_trees,
                lambda_d,
                lambda_t,
                lambda_l,
                seed,
            }])
        }
        _ => Err(
            "Provide either no arguments for the default suite, or all of: --species --gene-trees --lambda-d --lambda-t --lambda-l [--seed]"
                .to_string(),
        ),
    }
}

fn run_case(case: BenchmarkCase) {
    println!(
        "Case: species_leaves={}, gene_trees={}, λ_d={:.2}, λ_t={:.2}, λ_l={:.2}",
        case.species_leaves, case.gene_trees, case.lambda_d, case.lambda_t, case.lambda_l
    );

    let mut species_rng = StdRng::seed_from_u64(case.seed);
    let species_start = Instant::now();
    let (mut species_tree, _) =
        simulate_bd_tree_bwd(case.species_leaves, 1.0, 0.3, &mut species_rng).unwrap();
    species_tree.assign_depths();
    let species_elapsed = species_start.elapsed();

    let species_labeled = species_tree.unrooted_topology_hash().unwrap();
    let species_unlabeled = species_tree.unrooted_shape_hash().unwrap();

    let mut gene_rng = StdRng::seed_from_u64(case.seed + 1);
    let gene_start = Instant::now();
    let (rec_trees, _) = simulate_dtl_batch(
        &species_tree,
        species_tree.root,
        case.lambda_d,
        case.lambda_t,
        case.lambda_l,
        None,
        None,
        case.gene_trees,
        false,
        &mut gene_rng,
    )
    .unwrap();
    let gene_elapsed = gene_start.elapsed();

    let mut total_extant = 0usize;
    let mut dropped = 0usize;
    let mut extant_trees = Vec::with_capacity(rec_trees.len());
    for rec_tree in &rec_trees {
        let extant = count_extant_genes(rec_tree);
        total_extant += extant;
        match extract_extant_gene_tree(rec_tree) {
            Ok(tree) => extant_trees.push(tree),
            Err(_) => dropped += 1,
        }
    }

    if extant_trees.is_empty() {
        println!("  all simulated gene trees lost every extant lineage; skipping hash timing");
        println!();
        return;
    }

    let total_nodes: usize = extant_trees.iter().map(|tree| tree.nodes.len()).sum();

    let labeled_start = Instant::now();
    let labeled_hashes: Vec<u64> = extant_trees
        .iter()
        .map(|tree| black_box(tree.unrooted_topology_hash().unwrap()))
        .collect();
    let labeled_elapsed = labeled_start.elapsed();

    let unlabeled_start = Instant::now();
    let unlabeled_hashes: Vec<u64> = extant_trees
        .iter()
        .map(|tree| black_box(tree.unrooted_shape_hash().unwrap()))
        .collect();
    let unlabeled_elapsed = unlabeled_start.elapsed();

    let labeled_unique: HashSet<u64> = labeled_hashes.iter().copied().collect();
    let unlabeled_unique: HashSet<u64> = unlabeled_hashes.iter().copied().collect();

    println!(
        "  species generation: {:>8.3} ms",
        species_elapsed.as_secs_f64() * 1_000.0
    );
    println!(
        "  species hashes: labeled={} unlabeled={}",
        species_labeled, species_unlabeled
    );
    println!(
        "  gene simulation:   {:>8.3} ms",
        gene_elapsed.as_secs_f64() * 1_000.0
    );
    println!(
        "  extant trees: {:>4} kept / {:>4} total | dropped={} | avg extant leaves: {:>6.2} | avg nodes: {:>6.2}",
        extant_trees.len(),
        rec_trees.len(),
        dropped,
        total_extant as f64 / rec_trees.len() as f64,
        total_nodes as f64 / extant_trees.len() as f64
    );
    println!(
        "  labeled hash:   {:>8.3} ms total | {:>8.3} us/tree | unique={}",
        labeled_elapsed.as_secs_f64() * 1_000.0,
        labeled_elapsed.as_secs_f64() * 1_000_000.0 / extant_trees.len() as f64,
        labeled_unique.len()
    );
    println!(
        "                  {:>8.1} trees/sec",
        extant_trees.len() as f64 / labeled_elapsed.as_secs_f64()
    );
    println!(
        "  unlabeled hash: {:>8.3} ms total | {:>8.3} us/tree | unique={}",
        unlabeled_elapsed.as_secs_f64() * 1_000.0,
        unlabeled_elapsed.as_secs_f64() * 1_000_000.0 / extant_trees.len() as f64,
        unlabeled_unique.len()
    );
    println!(
        "                  {:>8.1} trees/sec",
        extant_trees.len() as f64 / unlabeled_elapsed.as_secs_f64()
    );
    println!();
}
