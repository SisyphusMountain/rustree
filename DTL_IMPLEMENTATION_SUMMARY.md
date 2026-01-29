# DTL Simulator Implementation Summary

Complete implementation of Duplication-Transfer-Loss (DTL) gene tree simulation with comprehensive benchmarking and XML export capabilities.

## Features Implemented

### 1. Core DTL Simulation

**File**: `src/dtl.rs`

- **Forward-time simulation** along species tree branches
- **Exponential waiting times** for D/T/L events
- **Contemporaneity-based transfers**: Uses species tree topology to find valid transfer recipients
- **Event types**: Speciation, Duplication, Transfer, Loss, Leaf
- **Reconciliation structure**: RecTree with gene-to-species mappings

**Key Algorithm Features:**
- Sum of rates across all active lineages for efficient event sampling
- Uniform selection of transfer recipients among contemporary branches
- Time subdivision tracking for species tree traversal
- Proper handling of gene lineage extinction

### 2. Event Tracking & Logging

**Structure**: `DTLEvent`

Tracks complete event information:
- Time of occurrence (absolute from root)
- Gene node involved
- Event type
- Species branch location
- For transfers: donor and recipient species
- Child gene nodes created

**CSV Export**: `save_events_to_csv()`
- Comma-separated format compatible with data analysis tools
- Full event history for downstream analysis

### 3. XML Export (RecPhyloXML)

**Method**: `RecTree::to_xml()`

Generates standard RecPhyloXML format:
- Complete species tree with branch lengths
- Reconciled gene tree with branch lengths
- Event annotations at each node
- Transfer events with `<branchingOut>` and `<transferBack>` tags
- Compatible with reconciliation visualization tools

### 4. Comprehensive Benchmarking

**File**: `tests/dtl_benchmark.rs`

Six benchmark suites testing different aspects:

1. **Quick Test** - Fast performance check (always runs)
2. **Varying Rates** - Tests different combinations of D/T/L rates
3. **Varying Tree Size** - Scalability with species tree size (5-100 species)
4. **Scalability** - Performance under increasing event rates
5. **Transfer Intensity** - Analysis of HGT patterns
6. **Loss Impact** - Effects of gene loss on simulation and extinction

### 5. Example Code

**File**: `examples/dtl_to_xml.rs`

Demonstrates:
- DTL simulation on custom species trees
- DTL simulation on birth-death generated trees
- XML export
- CSV event logging

## Performance Summary

### Quick Benchmark Results
- **Speed**: ~7,000 simulations/second
- **Average time**: 0.144 ms per simulation (10 species, balanced rates)
- **Average events**: 132 events per simulation
- **Processing rate**: ~850-1000 events/ms

### Performance Characteristics

**Fast scenarios** (<1 ms):
- Small trees (5-20 species)
- Low-to-moderate event rates
- High loss rates (natural limiting factor)

**Moderate scenarios** (1-100 ms):
- Medium trees (20-50 species)
- Balanced event rates
- Typical phylogenetic studies

**Slow scenarios** (>100 ms):
- Large trees (100+ species)
- High duplication + transfer rates
- Low loss rates (gene families explode)

### Key Performance Insights

1. **Loss events speed up simulation** by reducing active lineages
2. **Transfer + Duplication together** causes exponential gene family growth
3. **Processing rate is constant** at ~900 events/ms regardless of event types
4. **Stochastic variation** causes high variance in simulation time
5. **Memory efficient** up to 1.5M events tested

## API Usage

### Basic Simulation

```rust
use rustree::dtl::{simulate_dtl, save_events_to_csv};
use rustree::newick::newick::parse_newick;

// Parse species tree
let newick = "((A:1,B:1)AB:1,C:2)root:0;";
let mut nodes = parse_newick(newick).unwrap();
let mut species_tree = nodes.pop().unwrap().to_flat_tree();
species_tree.assign_depths();

// Simulate gene tree
let (rec_tree, events) = simulate_dtl(
    &species_tree,
    species_tree.root,
    1.0,  // duplication rate
    0.5,  // transfer rate
    0.5,  // loss rate
    &mut rng,
);

// Export results
let xml = rec_tree.to_xml();
std::fs::write("output.xml", &xml).unwrap();
save_events_to_csv(&events, "events.csv").unwrap();
```

### Event Analysis

```rust
use rustree::dtl::{count_events, count_extant_genes};

let (s, d, t, l, leaves) = count_events(&rec_tree);
let extant = count_extant_genes(&rec_tree);

println!("Speciation: {}, Duplication: {}, Transfer: {}, Loss: {}", s, d, t, l);
println!("Total leaves: {}, Extant genes: {}", leaves, extant);
```

## Files Created/Modified

### New Files
- `src/dtl.rs` - Core DTL simulation (660+ lines)
- `tests/dtl_benchmark.rs` - Comprehensive benchmarks (450+ lines)
- `examples/dtl_to_xml.rs` - Example usage
- `DTL_BENCHMARK_RESULTS.md` - Benchmark analysis
- `DTL_IMPLEMENTATION_SUMMARY.md` - This file
- `run_dtl_benchmarks.sh` - Benchmark runner script

### Modified Files
- `src/node.rs` - Added `to_xml()` method for RecTree (140+ lines)
- `src/lib.rs` - Added DTL module export

## Testing

All tests pass:
```bash
cargo test dtl

running 5 tests
test dtl::tests::test_dtl_pure_speciation ... ok
test dtl::tests::test_dtl_with_loss ... ok
test dtl::tests::test_dtl_xml_export ... ok
test dtl::tests::test_dtl_with_duplication ... ok
test dtl::tests::test_dtl_with_transfer ... ok
```

## Running the Code

### Quick Test
```bash
cargo run --example dtl_to_xml
```

Generates:
- `output_rectree.xml` - RecPhyloXML file
- `output_dtl_events.csv` - Event log
- `output_rectree_bd.xml` - Random BD tree example
- `output_dtl_bd_events.csv` - BD tree events

### Run Benchmarks
```bash
# Quick check
cargo test --test dtl_benchmark benchmark_dtl_quick_test -- --nocapture

# All fast benchmarks
./run_dtl_benchmarks.sh

# Individual benchmarks
cargo test --test dtl_benchmark benchmark_dtl_varying_rates -- --ignored --nocapture
cargo test --test dtl_benchmark benchmark_dtl_loss_impact -- --ignored --nocapture
```

## Scientific Applications

This implementation enables:

1. **Gene family evolution studies** - Simulate realistic gene trees within species phylogenies
2. **Horizontal gene transfer research** - Model HGT with contemporaneity constraints
3. **Genome evolution** - Study duplication and loss patterns
4. **Reconciliation analysis** - Generate test data for reconciliation algorithms
5. **Method validation** - Benchmark reconciliation inference methods with known truth

## Future Enhancements (Potential)

- Variable rates along branches (rate shifts)
- Birth-death model for gene tree simulation (analogous to species trees)
- Multiple gene family simulation in parallel
- GPU acceleration for large-scale simulations
- Direct integration with ALeRax/GeneRax formats

## References

The DTL model implementation follows standard reconciliation theory:
- Uses exponential waiting times for event sampling
- Contemporaneity-based transfer constraints
- RecPhyloXML format for standardized output
- Forward-time simulation for efficiency

## Conclusion

This implementation provides a fast, accurate, and well-tested DTL simulator suitable for both research and methods development. The comprehensive benchmarking demonstrates excellent performance characteristics and identifies the key factors affecting simulation speed. The XML export enables integration with existing phylogenetic tools, while the CSV event logging supports detailed analysis of evolutionary patterns.
