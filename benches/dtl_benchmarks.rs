use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rustree::bd::simulate_bd_tree_bwd;
use rustree::dtl::{
    simulate_dtl, simulate_dtl_batch, simulate_dtl_batch_with_branch_rates,
    simulate_dtl_per_species, BranchDTLRates,
};
use rustree::FlatTree;
use std::time::Duration;

/// Create a species tree of size n with deterministic seed.
fn make_species_tree(n: usize, seed: u64) -> FlatTree {
    let mut rng = StdRng::seed_from_u64(seed);
    let (mut tree, _) = simulate_bd_tree_bwd(n, 1.0, 0.3, &mut rng).unwrap();
    tree.assign_depths();
    tree
}

/// Create a full branch-rate table equivalent to scalar DTL rates.
fn expanded_branch_rates(
    species_tree: &FlatTree,
    origin_species: usize,
    d: f64,
    t: f64,
    l: f64,
) -> BranchDTLRates {
    let n = species_tree.nodes.len();
    let mut origination_probability = vec![0.0; n];
    origination_probability[origin_species] = 1.0;
    BranchDTLRates::new(vec![d; n], vec![t; n], vec![l; n], origination_probability).unwrap()
}

/// Pilot 10 per-gene DTL simulations to estimate average events per run.
fn pilot_dtl_events(species_tree: &FlatTree, d: f64, t: f64, l: f64) -> u64 {
    let root = species_tree.root;
    let mut rng = StdRng::seed_from_u64(42);
    let total: usize = (0..10)
        .map(|_| {
            simulate_dtl(species_tree, root, d, t, l, None, None, false, &mut rng)
                .unwrap()
                .1
                .len()
        })
        .sum();
    (total / 10).max(1) as u64
}

/// Pilot 10 per-species DTL simulations to estimate average events per run.
fn pilot_dtl_per_species_events(species_tree: &FlatTree, d: f64, t: f64, l: f64) -> u64 {
    let root = species_tree.root;
    let mut rng = StdRng::seed_from_u64(42);
    let total: usize = (0..10)
        .map(|_| {
            simulate_dtl_per_species(species_tree, root, d, t, l, None, None, false, &mut rng)
                .unwrap()
                .1
                .len()
        })
        .sum();
    (total / 10).max(1) as u64
}

fn dtl_pergene_scaling(c: &mut Criterion) {
    let (d, t, l) = (0.5, 0.2, 0.3);
    for &n in &[10, 50, 100, 200] {
        let species_tree = make_species_tree(n, 99);
        let avg_events = pilot_dtl_events(&species_tree, d, t, l);
        let mut group = c.benchmark_group("dtl_pergene_scaling");
        group.throughput(Throughput::Elements(avg_events));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            let root = species_tree.root;
            let mut rng = StdRng::seed_from_u64(123);
            b.iter(|| {
                simulate_dtl(&species_tree, root, d, t, l, None, None, false, &mut rng).unwrap()
            });
        });
        group.finish();
    }
}

fn dtl_pergene_rate_sweep(c: &mut Criterion) {
    let species_tree = make_species_tree(20, 99);
    let root = species_tree.root;
    let base = (1.0, 0.5, 0.3);
    for &mult in &[0.1, 0.5, 1.0, 2.0] {
        let (d, t, l) = (base.0 * mult, base.1 * mult, base.2 * mult);
        let avg_events = pilot_dtl_events(&species_tree, d, t, l);
        let mut group = c.benchmark_group("dtl_pergene_rate_sweep");
        group.throughput(Throughput::Elements(avg_events));
        let label = format!("mult_{}", mult);
        group.bench_with_input(BenchmarkId::new("rate", &label), &mult, |b, _| {
            let mut rng = StdRng::seed_from_u64(123);
            b.iter(|| {
                simulate_dtl(&species_tree, root, d, t, l, None, None, false, &mut rng).unwrap()
            });
        });
        group.finish();
    }
}

fn dtl_perspecies_scaling(c: &mut Criterion) {
    let (d, t, l) = (0.5, 0.2, 0.3);
    for &n in &[10, 50, 100, 200] {
        let species_tree = make_species_tree(n, 99);
        let avg_events = pilot_dtl_per_species_events(&species_tree, d, t, l);
        let mut group = c.benchmark_group("dtl_perspecies_scaling");
        group.throughput(Throughput::Elements(avg_events));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            let root = species_tree.root;
            let mut rng = StdRng::seed_from_u64(123);
            b.iter(|| {
                simulate_dtl_per_species(&species_tree, root, d, t, l, None, None, false, &mut rng)
                    .unwrap()
            });
        });
        group.finish();
    }
}

fn dtl_perspecies_rate_sweep(c: &mut Criterion) {
    let species_tree = make_species_tree(20, 99);
    let root = species_tree.root;
    let base = (1.0, 0.5, 0.3);
    for &mult in &[0.1, 0.5, 1.0, 2.0] {
        let (d, t, l) = (base.0 * mult, base.1 * mult, base.2 * mult);
        let avg_events = pilot_dtl_per_species_events(&species_tree, d, t, l);
        let mut group = c.benchmark_group("dtl_perspecies_rate_sweep");
        group.throughput(Throughput::Elements(avg_events));
        let label = format!("mult_{}", mult);
        group.bench_with_input(BenchmarkId::new("rate", &label), &mult, |b, _| {
            let mut rng = StdRng::seed_from_u64(123);
            b.iter(|| {
                simulate_dtl_per_species(&species_tree, root, d, t, l, None, None, false, &mut rng)
                    .unwrap()
            });
        });
        group.finish();
    }
}

fn dtl_batch_pergene(c: &mut Criterion) {
    let species_tree = make_species_tree(50, 99);
    let root = species_tree.root;
    let (d, t, l) = (0.2, 0.2, 0.1);
    for &batch_size in &[100usize, 1000] {
        // Estimate total events for the batch
        let avg_per_sim = pilot_dtl_events(&species_tree, d, t, l);
        let total_events = avg_per_sim * batch_size as u64;
        let mut group = c.benchmark_group("dtl_batch_pergene");
        group.throughput(Throughput::Elements(total_events));
        group.bench_with_input(
            BenchmarkId::from_parameter(batch_size),
            &batch_size,
            |b, &bs| {
                let mut rng = StdRng::seed_from_u64(123);
                b.iter(|| {
                    simulate_dtl_batch(
                        &species_tree,
                        root,
                        d,
                        t,
                        l,
                        None,
                        None,
                        bs,
                        false,
                        &mut rng,
                    )
                    .unwrap()
                });
            },
        );
        group.finish();
    }
}

fn dtl_scalar_vs_branch_rate_full_tables(c: &mut Criterion) {
    let species_tree = make_species_tree(20, 99);
    let root = species_tree.root;
    let (d, t, l) = (0.05, 0.02, 0.02);
    let batch_size = 10usize;
    let avg_per_sim = pilot_dtl_events(&species_tree, d, t, l);
    let total_events = avg_per_sim * batch_size as u64;
    let branch_rates = expanded_branch_rates(&species_tree, root, d, t, l);

    let mut group = c.benchmark_group("dtl_scalar_vs_branch_rate_full_tables");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_millis(500));
    group.throughput(Throughput::Elements(total_events));

    group.bench_function("scalar_batch", |b| {
        let mut rng = StdRng::seed_from_u64(123);
        b.iter(|| {
            simulate_dtl_batch(
                &species_tree,
                root,
                d,
                t,
                l,
                None,
                None,
                batch_size,
                false,
                &mut rng,
            )
            .unwrap()
        });
    });

    group.bench_function("branch_rate_full_batch", |b| {
        let mut rng = StdRng::seed_from_u64(123);
        b.iter(|| {
            simulate_dtl_batch_with_branch_rates(
                &species_tree,
                branch_rates.clone(),
                None,
                None,
                batch_size,
                false,
                &mut rng,
            )
            .unwrap()
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    dtl_pergene_scaling,
    dtl_pergene_rate_sweep,
    dtl_perspecies_scaling,
    dtl_perspecies_rate_sweep,
    dtl_batch_pergene,
    dtl_scalar_vs_branch_rate_full_tables
);
criterion_main!(benches);
