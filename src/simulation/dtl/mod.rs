// DTL (Duplication-Transfer-Loss) simulation along a species tree
//
// This module simulates gene tree evolution within a species tree using the DTL model.
// Events: Speciation (S), Duplication (D), Transfer (T), Loss (L)

mod event;
pub(crate) mod gillespie;
mod per_gene;
mod per_species;
mod state;
pub(crate) mod stream;
mod utils;

// Re-export main types and functions
use crate::error::RustreeError;

pub use event::DTLEvent;
pub use per_gene::{
    simulate_dtl, simulate_dtl_batch, simulate_dtl_batch_with_branch_rates, simulate_dtl_iter,
    simulate_dtl_iter_with_branch_rates, simulate_dtl_iter_with_config,
    simulate_dtl_with_branch_rates,
};
pub use per_species::{
    simulate_dtl_per_species, simulate_dtl_per_species_batch,
    simulate_dtl_per_species_batch_with_branch_rates, simulate_dtl_per_species_iter,
    simulate_dtl_per_species_iter_with_branch_rates, simulate_dtl_per_species_iter_with_config,
    simulate_dtl_per_species_with_branch_rates,
};
pub use stream::DtlSimIter;
pub(crate) use utils::prepare_simulation;
pub use utils::{count_events, count_extant_genes};

const ORIGINATION_SUM_TOLERANCE: f64 = 1e-8;

/// Full per-branch DTL rates and origination probabilities.
///
/// Branches are identified by species-tree node index. The value at index `i`
/// applies while a gene copy resides on species branch/node `i`; the root index
/// represents the root stem branch.
#[derive(Clone, Debug)]
pub struct BranchDTLRates {
    /// Duplication rate per species-tree branch/node index
    pub lambda_d: Vec<f64>,
    /// Transfer rate per species-tree branch/node index
    pub lambda_t: Vec<f64>,
    /// Loss rate per species-tree branch/node index
    pub lambda_l: Vec<f64>,
    /// Precomputed total DTL event rate per species-tree branch/node index
    pub lambda_total: Vec<f64>,
    /// Categorical origination probability per species-tree branch/node index
    pub origination_probability: Vec<f64>,
    /// Precomputed cumulative origination probability per species-tree branch/node index
    pub origination_cdf: Vec<f64>,
}

impl BranchDTLRates {
    /// Create a validated full branch-rate table.
    pub fn new(
        lambda_d: Vec<f64>,
        lambda_t: Vec<f64>,
        lambda_l: Vec<f64>,
        origination_probability: Vec<f64>,
    ) -> Result<Self, RustreeError> {
        let (lambda_total, origination_cdf) = Self::validate_and_precompute(
            &lambda_d,
            &lambda_t,
            &lambda_l,
            &origination_probability,
        )?;

        Ok(Self {
            lambda_d,
            lambda_t,
            lambda_l,
            lambda_total,
            origination_probability,
            origination_cdf,
        })
    }

    /// Validate vector shape and numeric values.
    pub fn validate(&self) -> Result<(), RustreeError> {
        let (lambda_total, origination_cdf) = Self::validate_and_precompute(
            &self.lambda_d,
            &self.lambda_t,
            &self.lambda_l,
            &self.origination_probability,
        )?;

        Self::validate_precomputed_vector("lambda_total", &self.lambda_total, &lambda_total)?;
        Self::validate_precomputed_vector(
            "origination_cdf",
            &self.origination_cdf,
            &origination_cdf,
        )?;

        Ok(())
    }

    fn validate_and_precompute(
        lambda_d: &[f64],
        lambda_t: &[f64],
        lambda_l: &[f64],
        origination_probability: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), RustreeError> {
        let n = lambda_d.len();
        if n == 0 {
            return Err(RustreeError::Validation(
                "branch DTL rates cannot be empty".to_string(),
            ));
        }
        if lambda_t.len() != n || lambda_l.len() != n {
            return Err(RustreeError::Validation(format!(
                "branch DTL rate vectors must have the same length: lambda_d={}, lambda_t={}, lambda_l={}",
                lambda_d.len(),
                lambda_t.len(),
                lambda_l.len()
            )));
        }
        if origination_probability.len() != n {
            return Err(RustreeError::Validation(format!(
                "origination_probability length ({}) must match branch rate length ({n})",
                origination_probability.len()
            )));
        }

        for (idx, &value) in lambda_d.iter().enumerate() {
            validate_rate(&format!("Duplication rate at branch {idx}"), value)?;
        }
        for (idx, &value) in lambda_t.iter().enumerate() {
            validate_rate(&format!("Transfer rate at branch {idx}"), value)?;
        }
        for (idx, &value) in lambda_l.iter().enumerate() {
            validate_rate(&format!("Loss rate at branch {idx}"), value)?;
        }

        let mut sum = 0.0;
        let mut origination_cdf = Vec::with_capacity(n);
        let mut last_positive = None;
        for (idx, &p) in origination_probability.iter().enumerate() {
            if p < 0.0 || !p.is_finite() {
                return Err(RustreeError::Validation(format!(
                    "origination_probability at branch {idx} must be non-negative and finite, got {p}"
                )));
            }
            sum += p;
            origination_cdf.push(sum);
            if p > 0.0 {
                last_positive = Some(idx);
            }
        }
        if (sum - 1.0).abs() > ORIGINATION_SUM_TOLERANCE {
            return Err(RustreeError::Validation(format!(
                "origination_probability must sum to 1.0 (within {ORIGINATION_SUM_TOLERANCE}), got {sum}"
            )));
        }
        if sum < 1.0 {
            if let Some(last_positive) = last_positive {
                for cdf in &mut origination_cdf[last_positive..] {
                    *cdf = 1.0;
                }
            }
        }

        let mut lambda_total = Vec::with_capacity(n);
        for idx in 0..n {
            let total = lambda_d[idx] + lambda_t[idx] + lambda_l[idx];
            if !total.is_finite() {
                return Err(RustreeError::Validation(format!(
                    "total DTL rate at branch {idx} must be finite, got {total}"
                )));
            }
            lambda_total.push(total);
        }

        Ok((lambda_total, origination_cdf))
    }

    fn validate_precomputed_vector(
        name: &str,
        actual: &[f64],
        expected: &[f64],
    ) -> Result<(), RustreeError> {
        if actual.len() != expected.len() {
            return Err(RustreeError::Validation(format!(
                "{name} length ({}) must match branch rate length ({})",
                actual.len(),
                expected.len()
            )));
        }

        for (idx, (&actual, &expected)) in actual.iter().zip(expected.iter()).enumerate() {
            let tolerance = f64::EPSILON * actual.abs().max(expected.abs()).max(1.0) * 8.0;
            if !actual.is_finite() || (actual - expected).abs() > tolerance {
                return Err(RustreeError::Validation(format!(
                    "precomputed {name} at branch {idx} is out of date: got {actual}, expected {expected}"
                )));
            }
        }

        Ok(())
    }

    /// Validate this table against a concrete species tree size.
    pub fn validate_for_tree(&self, node_count: usize) -> Result<(), RustreeError> {
        self.validate()?;
        if self.lambda_d.len() != node_count {
            return Err(RustreeError::Validation(format!(
                "branch DTL tables must have one value per species-tree node/branch: got {}, expected {}",
                self.lambda_d.len(),
                node_count
            )));
        }
        Ok(())
    }

    #[inline]
    pub(crate) fn total_rate(&self, species_idx: usize) -> f64 {
        self.lambda_total[species_idx]
    }

    #[inline]
    pub(crate) fn event_rates(&self, species_idx: usize) -> (f64, f64, f64) {
        debug_assert!(
            (self.lambda_total[species_idx]
                - (self.lambda_d[species_idx]
                    + self.lambda_t[species_idx]
                    + self.lambda_l[species_idx]))
                .abs()
                <= f64::EPSILON * self.lambda_total[species_idx].abs().max(1.0) * 8.0
        );
        (
            self.lambda_d[species_idx],
            self.lambda_t[species_idx],
            self.lambda_l[species_idx],
        )
    }

    pub(crate) fn sample_origin<R: rand::Rng>(&self, rng: &mut R) -> usize {
        let threshold = rng.gen::<f64>();
        let idx = self
            .origination_cdf
            .partition_point(|&cdf| threshold >= cdf);
        idx.min(self.origination_cdf.len().saturating_sub(1))
    }
}

/// Configuration for DTL (Duplication-Transfer-Loss) simulation.
#[derive(Clone, Debug)]
pub struct DTLConfig {
    /// Duplication rate
    pub lambda_d: f64,
    /// Transfer rate
    pub lambda_t: f64,
    /// Loss rate
    pub lambda_l: f64,
    /// Assortative transfer parameter (None = uniform transfer)
    pub transfer_alpha: Option<f64>,
    /// Replacement transfer probability (None = additive only)
    pub replacement_transfer: Option<f64>,
    /// Optional full per-branch DTL rates and origination probabilities
    pub branch_rates: Option<BranchDTLRates>,
}

impl DTLConfig {
    /// Create a validated DTL simulation configuration.
    pub fn new(
        lambda_d: f64,
        lambda_t: f64,
        lambda_l: f64,
        transfer_alpha: Option<f64>,
        replacement_transfer: Option<f64>,
    ) -> Result<Self, RustreeError> {
        let config = Self {
            lambda_d,
            lambda_t,
            lambda_l,
            transfer_alpha,
            replacement_transfer,
            branch_rates: None,
        };
        config.validate()?;
        Ok(config)
    }

    /// Create a validated DTL configuration backed by full per-branch rates.
    pub fn with_branch_rates(
        branch_rates: BranchDTLRates,
        transfer_alpha: Option<f64>,
        replacement_transfer: Option<f64>,
    ) -> Result<Self, RustreeError> {
        let config = Self {
            lambda_d: 0.0,
            lambda_t: 0.0,
            lambda_l: 0.0,
            transfer_alpha,
            replacement_transfer,
            branch_rates: Some(branch_rates),
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate all numeric configuration values.
    pub fn validate(&self) -> Result<(), RustreeError> {
        validate_rate("Duplication", self.lambda_d)?;
        validate_rate("Transfer", self.lambda_t)?;
        validate_rate("Loss", self.lambda_l)?;

        if let Some(alpha) = self.transfer_alpha {
            if !alpha.is_finite() {
                return Err(RustreeError::Validation(format!(
                    "transfer_alpha must be finite, got {alpha}"
                )));
            }
        }

        if let Some(p) = self.replacement_transfer {
            if !(0.0..=1.0).contains(&p) {
                return Err(RustreeError::Validation(format!(
                    "replacement_transfer must be finite and between 0.0 and 1.0, got {p}"
                )));
            }
        }

        if let Some(branch_rates) = &self.branch_rates {
            branch_rates.validate()?;
        }

        Ok(())
    }

    /// Validate branch-table shape against a concrete species tree.
    pub fn validate_for_tree(&self, node_count: usize) -> Result<(), RustreeError> {
        self.validate()?;
        if let Some(branch_rates) = &self.branch_rates {
            branch_rates.validate_for_tree(node_count)?;
        }
        Ok(())
    }

    #[inline]
    pub(crate) fn branch_total_rate(&self, species_idx: usize) -> f64 {
        self.branch_rates
            .as_ref()
            .map_or(self.lambda_d + self.lambda_t + self.lambda_l, |rates| {
                rates.total_rate(species_idx)
            })
    }

    #[inline]
    pub(crate) fn branch_event_rates(&self, species_idx: usize) -> (f64, f64, f64) {
        self.branch_rates
            .as_ref()
            .map_or((self.lambda_d, self.lambda_t, self.lambda_l), |rates| {
                rates.event_rates(species_idx)
            })
    }

    pub(crate) fn sample_origin<R: rand::Rng>(&self, default_origin: usize, rng: &mut R) -> usize {
        self.branch_rates
            .as_ref()
            .map_or(default_origin, |rates| rates.sample_origin(rng))
    }

    pub(crate) fn uses_branch_rates(&self) -> bool {
        self.branch_rates.is_some()
    }
}

fn validate_rate(name: &str, value: f64) -> Result<(), RustreeError> {
    if value < 0.0 || !value.is_finite() {
        return Err(RustreeError::Validation(format!(
            "{name} rate must be non-negative and finite, got {value}"
        )));
    }
    Ok(())
}

// Re-export from io module for backward compatibility
pub use crate::io::save_dtl_events_to_csv as save_events_to_csv;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::newick::parse_newick;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn parse_tree(newick: &str) -> crate::node::FlatTree {
        let mut nodes = parse_newick(newick).unwrap();
        let root = nodes.pop().expect("No tree found");
        root.to_flat_tree()
    }

    fn zero_branch_rates_with_origin(
        species_tree: &crate::node::FlatTree,
        origin_species: usize,
    ) -> BranchDTLRates {
        let n = species_tree.nodes.len();
        let mut origination_probability = vec![0.0; n];
        origination_probability[origin_species] = 1.0;
        BranchDTLRates::new(
            vec![0.0; n],
            vec![0.0; n],
            vec![0.0; n],
            origination_probability,
        )
        .unwrap()
    }

    #[test]
    fn test_branch_dtl_rates_validation() {
        assert!(BranchDTLRates::new(vec![0.0], vec![0.0], vec![0.0], vec![1.0]).is_ok());
        assert!(BranchDTLRates::new(vec![0.0], vec![0.0, 0.0], vec![0.0], vec![1.0]).is_err());
        assert!(BranchDTLRates::new(vec![0.0], vec![0.0], vec![-0.1], vec![1.0]).is_err());
        assert!(BranchDTLRates::new(
            vec![0.0, 0.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
            vec![0.5, 0.25]
        )
        .is_err());
    }

    #[test]
    fn test_branch_dtl_origin_can_be_non_root() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();
        let origin_species = species_tree
            .nodes
            .iter()
            .position(|node| node.name == "AB")
            .unwrap();

        let branch_rates = zero_branch_rates_with_origin(&species_tree, origin_species);
        let mut rng = StdRng::seed_from_u64(42);
        let (rec_tree, _events) = simulate_dtl_with_branch_rates(
            &species_tree,
            branch_rates,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let extant_leaf_species: std::collections::BTreeSet<_> = rec_tree
            .gene_tree
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| {
                if node.left_child.is_none()
                    && node.right_child.is_none()
                    && rec_tree.event_mapping[idx] == crate::node::Event::Leaf
                {
                    rec_tree.node_mapping[idx]
                        .map(|sp| rec_tree.species_tree.nodes[sp].name.clone())
                } else {
                    None
                }
            })
            .collect();

        assert_eq!(
            extant_leaf_species,
            ["A".to_string(), "B".to_string()].into_iter().collect()
        );
    }

    #[test]
    fn test_branch_dtl_rates_on_unoccupied_branches_do_not_fire() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();
        let origin_species = species_tree
            .nodes
            .iter()
            .position(|node| node.name == "AB")
            .unwrap();
        let c_species = species_tree
            .nodes
            .iter()
            .position(|node| node.name == "C")
            .unwrap();

        let n = species_tree.nodes.len();
        let mut lambda_l = vec![0.0; n];
        lambda_l[c_species] = 1_000.0;
        let mut origination_probability = vec![0.0; n];
        origination_probability[origin_species] = 1.0;
        let branch_rates = BranchDTLRates::new(
            vec![0.0; n],
            vec![0.0; n],
            lambda_l,
            origination_probability,
        )
        .unwrap();

        let mut rng = StdRng::seed_from_u64(42);
        let (rec_tree, _events) = simulate_dtl_with_branch_rates(
            &species_tree,
            branch_rates,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        assert_eq!(count_extant_genes(&rec_tree), 2);
    }

    #[test]
    fn test_dtl_pure_speciation() {
        // Simple species tree: ((A:1,B:1)AB:1,C:2)root:0;
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // With zero D/T/L rates, should get pure speciation
        let (rec_tree, events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            0.0,
            0.0,
            0.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Events: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );
        println!("Total event count: {}", events.len());

        // Should have 3 extant genes (one per species leaf)
        let extant = count_extant_genes(&rec_tree);
        assert_eq!(extant, 3, "Should have 3 extant genes with no D/T/L");
    }

    #[test]
    fn test_dtl_with_duplication() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(123);

        // High duplication rate
        let (rec_tree, _events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            2.0,
            0.0,
            0.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Events with duplication: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );

        // Should have duplications
        assert!(
            d > 0 || leaves >= 3,
            "Should have duplications or at least 3 leaves"
        );
    }

    #[test]
    fn test_dtl_with_loss() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(999);

        // High loss rate
        let (rec_tree, _events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            0.0,
            0.0,
            5.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Events with loss: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );

        // May have losses
        let extant = count_extant_genes(&rec_tree);
        println!("Extant genes: {}", extant);
    }

    #[test]
    fn test_dtl_xml_export() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // With mixed events
        let (rec_tree, events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            1.0,
            0.5,
            0.5,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        // Save events to CSV in a temp directory
        let tmp = tempfile::tempdir().unwrap();
        let csv_path = tmp.path().join("test_dtl_events.csv");
        save_events_to_csv(
            &events,
            &species_tree,
            &rec_tree.gene_tree,
            csv_path.to_str().unwrap(),
        )
        .expect("Failed to write events CSV");

        // Generate XML
        let xml = rec_tree.to_xml();

        // Verify XML has required sections
        assert!(xml.contains("<recPhylo"));
        assert!(xml.contains("<spTree>"));
        assert!(xml.contains("<recGeneTree>"));
        assert!(xml.contains("<branchLength>"));
    }

    #[test]
    fn test_dtl_with_transfer() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(456);

        // High transfer rate (uniform)
        let (rec_tree, events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            0.0,
            2.0,
            0.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Events with transfer: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );
        println!("Total event count: {}", events.len());
    }

    #[test]
    fn test_dtl_assortative_transfer() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(456);

        // High transfer rate with assortative selection (alpha=1.0)
        let (rec_tree, events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            0.0,
            2.0,
            0.0,
            Some(1.0),
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Events with assortative transfer: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );
        println!("Total event count: {}", events.len());

        // Verify LCA computation works
        let lca_depths = species_tree
            .precompute_lca_depths()
            .expect("Failed to precompute LCA depths");
        assert!(!lca_depths.is_empty(), "LCA depths should be computed");
    }

    // Tests for per-species (Zombi-style) simulation

    #[test]
    fn test_per_species_pure_speciation() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // With zero D/T/L rates, should get pure speciation
        let (rec_tree, events) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            0.0,
            0.0,
            0.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Per-species pure speciation: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );
        println!("Total event count: {}", events.len());

        // Should have 3 extant genes (one per species leaf)
        let extant = count_extant_genes(&rec_tree);
        assert_eq!(extant, 3, "Should have 3 extant genes with no D/T/L");
    }

    #[test]
    fn test_per_species_with_duplication() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(123);

        // High duplication rate
        let (rec_tree, _events) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            2.0,
            0.0,
            0.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Per-species with duplication: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );

        // Should have duplications or at least the original 3 leaves
        assert!(
            d > 0 || leaves >= 3,
            "Should have duplications or at least 3 leaves"
        );
    }

    #[test]
    fn test_per_species_with_transfer() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(456);

        // High transfer rate (uniform)
        let (rec_tree, events) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            0.0,
            2.0,
            0.0,
            None,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Per-species with transfer: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );
        println!("Total event count: {}", events.len());
    }

    #[test]
    fn test_per_species_compare_event_counts() {
        // Compare per-gene-copy vs per-species models with same rates
        // Per-species should have fewer events on average when there are many gene copies
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut total_per_gene = 0;
        let mut total_per_species = 0;

        for seed in 0..10 {
            let mut rng1 = StdRng::seed_from_u64(seed);
            let mut rng2 = StdRng::seed_from_u64(seed);

            let (rec1, _) = simulate_dtl(
                &species_tree,
                species_tree.root,
                1.0,
                0.5,
                0.5,
                None,
                None,
                false,
                &mut rng1,
            )
            .unwrap();
            let (rec2, _) = simulate_dtl_per_species(
                &species_tree,
                species_tree.root,
                1.0,
                0.5,
                0.5,
                None,
                None,
                false,
                &mut rng2,
            )
            .unwrap();

            let (_, d1, t1, l1, _) = count_events(&rec1);
            let (_, d2, t2, l2, _) = count_events(&rec2);

            total_per_gene += d1 + t1 + l1;
            total_per_species += d2 + t2 + l2;
        }

        println!("Total DTL events over 10 runs:");
        println!("  Per-gene-copy model: {}", total_per_gene);
        println!("  Per-species model:   {}", total_per_species);

        // Per-species should generally have fewer events because duplications don't increase rate
        // But this isn't guaranteed for all seeds, so we just print for comparison
    }

    #[test]
    fn test_per_species_require_extant() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(999);

        // High loss rate, but require extant genes
        let (rec_tree, _events) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            0.0,
            0.0,
            5.0,
            None,
            None,
            true,
            &mut rng,
        )
        .unwrap();

        let extant = count_extant_genes(&rec_tree);
        assert!(
            extant > 0,
            "Should have at least one extant gene when require_extant=true"
        );
        println!(
            "Per-species with high loss (require_extant=true): {} extant genes",
            extant
        );
    }

    // Tests for replacement transfers

    #[test]
    fn test_replacement_transfer_per_species() {
        let newick = "((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // High duplication + transfer rates with 100% replacement
        let (rec_tree, _events) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            1.0,
            2.0,
            0.0,
            None,
            Some(1.0),
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Replacement per-species: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );

        // With replacement transfers, we should see loss events from replacements
        assert!(t > 0, "Should have transfers");
        assert!(l > 0, "Replacement transfers should cause loss events");
    }

    #[test]
    fn test_replacement_transfer_per_gene() {
        let newick = "((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        let mut rng = StdRng::seed_from_u64(42);

        // High duplication + transfer rates with 100% replacement
        let (rec_tree, _events) = simulate_dtl(
            &species_tree,
            species_tree.root,
            1.0,
            2.0,
            0.0,
            None,
            Some(1.0),
            false,
            &mut rng,
        )
        .unwrap();

        let (s, d, t, l, leaves) = count_events(&rec_tree);
        println!(
            "Replacement per-gene: S={}, D={}, T={}, L={}, Leaves={}",
            s, d, t, l, leaves
        );

        // With replacement transfers, we should see both transfer and loss events
        assert!(t > 0, "Should have transfers");
        assert!(l > 0, "Replacement transfers should cause loss events");
    }

    #[test]
    fn test_replacement_none_backward_compatible() {
        let newick = "((A:1,B:1)AB:1,C:2)root:0;";
        let mut species_tree = parse_tree(newick);
        species_tree.assign_depths();

        // None vs Some(0.0) should produce identical results
        let mut rng1 = StdRng::seed_from_u64(42);
        let mut rng2 = StdRng::seed_from_u64(42);

        let (rec1, _) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            1.0,
            1.0,
            0.5,
            None,
            None,
            false,
            &mut rng1,
        )
        .unwrap();
        let (rec2, _) = simulate_dtl_per_species(
            &species_tree,
            species_tree.root,
            1.0,
            1.0,
            0.5,
            None,
            Some(0.0),
            false,
            &mut rng2,
        )
        .unwrap();

        let events1 = count_events(&rec1);
        let events2 = count_events(&rec2);
        assert_eq!(
            events1, events2,
            "None and Some(0.0) should produce identical results"
        );
    }
}
