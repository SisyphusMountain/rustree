use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rustree::bd::simulate_bd_tree_bwd;
use rustree::dtl::{
    count_extant_genes, simulate_dtl, simulate_dtl_batch, simulate_dtl_batch_with_branch_rates,
    simulate_dtl_per_species_with_branch_rates, simulate_dtl_with_branch_rates, BranchDTLRates,
    DTLEvent,
};
use rustree::{parse_newick, Event, FlatTree, RecTree};
use std::collections::BTreeSet;
use std::time::Instant;

fn branch_test_tree() -> FlatTree {
    let mut nodes = parse_newick("((A:1.0,B:1.0)AB:1.0,C:2.0)root:0.0;").unwrap();
    let mut tree = nodes.pop().expect("test tree should parse").to_flat_tree();
    tree.assign_depths();
    tree
}

fn species_index(tree: &FlatTree, name: &str) -> usize {
    tree.nodes
        .iter()
        .position(|node| node.name == name)
        .unwrap_or_else(|| panic!("species {name} not found"))
}

fn expanded_branch_rates(
    tree: &FlatTree,
    origin_species: usize,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
) -> BranchDTLRates {
    let n = tree.nodes.len();
    let mut origination_probability = vec![0.0; n];
    origination_probability[origin_species] = 1.0;

    BranchDTLRates::new(
        vec![lambda_d; n],
        vec![lambda_t; n],
        vec![lambda_l; n],
        origination_probability,
    )
    .unwrap()
}

fn local_loss_branch_rates(
    tree: &FlatTree,
    origin_species: usize,
    loss_species: usize,
    lambda_l: f64,
) -> BranchDTLRates {
    let n = tree.nodes.len();
    let mut origination_probability = vec![0.0; n];
    origination_probability[origin_species] = 1.0;

    let mut losses = vec![0.0; n];
    losses[loss_species] = lambda_l;

    BranchDTLRates::new(vec![0.0; n], vec![0.0; n], losses, origination_probability).unwrap()
}

#[derive(Debug, PartialEq)]
struct RecTreeSignature {
    nodes: Vec<NodeSignature>,
    node_mapping: Vec<Option<usize>>,
    event_mapping: Vec<String>,
}

#[derive(Debug, PartialEq)]
struct NodeSignature {
    name: String,
    left_child: Option<usize>,
    right_child: Option<usize>,
    parent: Option<usize>,
    depth: Option<u64>,
    length: u64,
}

fn rec_tree_signature(rec_tree: &RecTree) -> RecTreeSignature {
    RecTreeSignature {
        nodes: rec_tree
            .gene_tree
            .nodes
            .iter()
            .map(|node| NodeSignature {
                name: node.name.clone(),
                left_child: node.left_child,
                right_child: node.right_child,
                parent: node.parent,
                depth: node.depth.map(f64::to_bits),
                length: node.length.to_bits(),
            })
            .collect(),
        node_mapping: rec_tree.node_mapping.clone(),
        event_mapping: rec_tree
            .event_mapping
            .iter()
            .map(|event| format!("{event:?}"))
            .collect(),
    }
}

fn dtl_event_signatures(events: &[DTLEvent]) -> Vec<String> {
    events
        .iter()
        .map(|event| match event {
            DTLEvent::Speciation {
                time,
                gene_id,
                species_id,
                left_child,
                right_child,
            } => format!(
                "S:{}:{gene_id}:{species_id}:{left_child}:{right_child}",
                time.to_bits()
            ),
            DTLEvent::Duplication {
                time,
                gene_id,
                species_id,
                child1,
                child2,
            } => format!(
                "D:{}:{gene_id}:{species_id}:{child1}:{child2}",
                time.to_bits()
            ),
            DTLEvent::Transfer {
                time,
                gene_id,
                species_id,
                from_species,
                to_species,
                donor_child,
                recipient_child,
            } => format!(
                "T:{}:{gene_id}:{species_id}:{from_species}:{to_species}:{donor_child}:{recipient_child}",
                time.to_bits()
            ),
            DTLEvent::Loss {
                time,
                gene_id,
                species_id,
            } => format!("L:{}:{gene_id}:{species_id}", time.to_bits()),
            DTLEvent::Leaf {
                time,
                gene_id,
                species_id,
            } => format!("F:{}:{gene_id}:{species_id}", time.to_bits()),
        })
        .collect()
}

fn extant_leaf_species(rec_tree: &RecTree) -> BTreeSet<String> {
    rec_tree
        .gene_tree
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(idx, node)| {
            if node.left_child.is_none()
                && node.right_child.is_none()
                && rec_tree.event_mapping[idx] == Event::Leaf
            {
                rec_tree.node_mapping[idx]
                    .map(|species| rec_tree.species_tree.nodes[species].name.clone())
            } else {
                None
            }
        })
        .collect()
}

fn assert_single_loss_on(rec_tree: &RecTree, events: &[DTLEvent], expected_species: usize) {
    assert_eq!(
        count_extant_genes(rec_tree),
        0,
        "forced local loss should leave no extant genes"
    );
    assert_eq!(events.len(), 1, "expected one DTL event before speciation");

    match &events[0] {
        DTLEvent::Loss { species_id, .. } => assert_eq!(*species_id, expected_species),
        other => panic!("expected one loss event, got {other:?}"),
    }
}

#[test]
fn scalar_expanded_branch_rates_match_scalar_loss_before_speciation() {
    let tree = branch_test_tree();
    let origin_species = species_index(&tree, "AB");
    let loss_rate = 1_000_000.0;

    let mut scalar_rng = StdRng::seed_from_u64(101);
    let _: f64 = scalar_rng.gen();
    let (scalar_tree, scalar_events) = simulate_dtl(
        &tree,
        origin_species,
        0.0,
        0.0,
        loss_rate,
        None,
        None,
        false,
        &mut scalar_rng,
    )
    .unwrap();

    let branch_rates = expanded_branch_rates(&tree, origin_species, 0.0, 0.0, loss_rate);
    let mut branch_rng = StdRng::seed_from_u64(101);
    let (branch_tree, branch_events) =
        simulate_dtl_with_branch_rates(&tree, branch_rates, None, None, false, &mut branch_rng)
            .unwrap();

    assert_eq!(
        rec_tree_signature(&scalar_tree),
        rec_tree_signature(&branch_tree)
    );
    assert_eq!(
        dtl_event_signatures(&scalar_events),
        dtl_event_signatures(&branch_events)
    );
    assert_single_loss_on(&branch_tree, &branch_events, origin_species);
}

#[test]
fn categorical_origination_forces_non_root_branch() {
    let tree = branch_test_tree();
    let origin_species = species_index(&tree, "AB");
    let branch_rates = expanded_branch_rates(&tree, origin_species, 0.0, 0.0, 0.0);

    let mut rng = StdRng::seed_from_u64(202);
    let (rec_tree, _events) =
        simulate_dtl_with_branch_rates(&tree, branch_rates, None, None, false, &mut rng).unwrap();

    let expected = BTreeSet::from(["A".to_string(), "B".to_string()]);
    assert_eq!(extant_leaf_species(&rec_tree), expected);
    assert_eq!(count_extant_genes(&rec_tree), 2);
}

#[test]
fn per_gene_branch_rates_use_local_rates_on_occupied_branch() {
    let tree = branch_test_tree();
    let origin_species = species_index(&tree, "AB");
    let branch_rates = local_loss_branch_rates(&tree, origin_species, origin_species, 1_000_000.0);

    let mut rng = StdRng::seed_from_u64(303);
    let (rec_tree, events) =
        simulate_dtl_with_branch_rates(&tree, branch_rates, None, None, false, &mut rng).unwrap();

    assert_single_loss_on(&rec_tree, &events, origin_species);
}

#[test]
fn per_species_branch_rates_use_local_rates_on_occupied_branch() {
    let tree = branch_test_tree();
    let origin_species = species_index(&tree, "AB");
    let branch_rates = local_loss_branch_rates(&tree, origin_species, origin_species, 1_000_000.0);

    let mut rng = StdRng::seed_from_u64(404);
    let (rec_tree, events) = simulate_dtl_per_species_with_branch_rates(
        &tree,
        branch_rates,
        None,
        None,
        false,
        &mut rng,
    )
    .unwrap();

    assert_single_loss_on(&rec_tree, &events, origin_species);
}

#[test]
#[ignore = "timing smoke benchmark; run explicitly"]
fn scalar_vs_branch_rate_full_tables_smoke_benchmark() {
    let mut tree_rng = StdRng::seed_from_u64(505);
    let (mut tree, _) = simulate_bd_tree_bwd(20, 1.0, 0.0, &mut tree_rng).unwrap();
    tree.assign_depths();

    let root = tree.root;
    let (lambda_d, lambda_t, lambda_l) = (0.05, 0.02, 0.02);
    let batch_size = 10usize;
    let reps = 20usize;

    let mut scalar_rng = StdRng::seed_from_u64(606);
    let scalar_start = Instant::now();
    for _ in 0..reps {
        simulate_dtl_batch(
            &tree,
            root,
            lambda_d,
            lambda_t,
            lambda_l,
            None,
            None,
            batch_size,
            false,
            &mut scalar_rng,
        )
        .unwrap();
    }
    let scalar_elapsed = scalar_start.elapsed();

    let branch_rates = expanded_branch_rates(&tree, root, lambda_d, lambda_t, lambda_l);
    let mut branch_rng = StdRng::seed_from_u64(606);
    let branch_start = Instant::now();
    for _ in 0..reps {
        simulate_dtl_batch_with_branch_rates(
            &tree,
            branch_rates.clone(),
            None,
            None,
            batch_size,
            false,
            &mut branch_rng,
        )
        .unwrap();
    }
    let branch_elapsed = branch_start.elapsed();

    println!(
        "scalar_batch_total={:?} branch_rate_full_batch_total={:?} reps={} batch_size={} species_nodes={}",
        scalar_elapsed,
        branch_elapsed,
        reps,
        batch_size,
        tree.nodes.len()
    );

    assert!(scalar_elapsed.as_nanos() > 0);
    assert!(branch_elapsed.as_nanos() > 0);
}
