# DTL Simulator Benchmark Results

Comprehensive benchmarking of the Duplication-Transfer-Loss (DTL) gene tree simulator.

## Quick Performance Summary

**Quick Benchmark (10 species, balanced rates):**
- Average time per simulation: **0.144 ms**
- Average events per simulation: **132 events**
- Throughput: **~7,000 simulations/second**

## 1. Varying Event Rates

Species tree: 10 species (generated via birth-death process)
Simulations per configuration: 100

| Configuration      | Time(ms) | S   | D    | T    | L   | Leaves | Extant |
|--------------------|----------|-----|------|------|-----|--------|--------|
| Pure speciation    | 0.020    | 11  | 0    | 0    | 0   | 12     | 12     |
| High duplication   | 0.100    | 29  | 26   | 0    | 0   | 57     | 57     |
| High transfer      | 0.150    | 34  | 0    | 30   | 0   | 65     | 65     |
| High loss          | 0.010    | 5   | 0    | 0    | 3   | 2      | 2      |
| Dup+Trans          | 0.770    | 128 | 114  | 115  | 0   | 359    | 359    |
| Dup+Loss           | 0.040    | 10  | 9    | 0    | 9   | 11     | 11     |
| Trans+Loss         | 0.060    | 11  | 0    | 10   | 9   | 12     | 12     |
| All events         | 0.220    | 35  | 31   | 30   | 30  | 67     | 67     |
| Balanced           | 0.070    | 17  | 7    | 7    | 7   | 25     | 25     |
| High D+T, low L    | 9.450    | 940 | 1618 | 1619 | 409 | 3769   | 3769   |

**Key Findings:**
- Pure speciation is fastest (0.020 ms)
- High loss rates reduce simulation time dramatically (fewer lineages to track)
- Combined high duplication and transfer rates can generate very large gene families (3769 leaves)
- Performance scales with the number of active lineages being simulated

## 2. Varying Species Tree Size

Event rates: λ_d=1.0, λ_t=0.5, λ_l=0.5
Simulations per tree size: 50 for small trees, 20 for large

| N Species | Time(ms)  | Events    | S       | D       | T       | L       | Leaves  |
|-----------|-----------|-----------|---------|---------|---------|---------|---------|
| 5         | 0.148     | 143       | 25      | 31      | 14      | 15      | 56      |
| 10        | 29.437    | 30,072    | 5,643   | 6,270   | 3,121   | 3,149   | 11,887  |
| 20        | 3.313     | 3,605     | 683     | 746     | 372     | 376     | 1,427   |
| 50        | 87.932    | 91,698    | 18,528  | 18,191  | 9,129   | 9,103   | 36,745  |
| 100       | 1566.793  | 1,478,687 | 293,470 | 297,161 | 148,710 | 148,678 | 590,665 |

**Key Findings:**
- Performance shows high variability depending on the actual tree topology
- Note: 10-species tree took longer than 20-species (stochastic variation in gene tree size)
- 100-species tree with these rates generates ~1.5 million events
- Time scales roughly linearly with the number of events generated

## 3. Scalability Test

Testing performance with increasing event rates
Species tree: 20 species, 10 replicates per configuration

| Rate Multiplier | Time(ms) | Total Events | Events/ms | Leaves  |
|-----------------|----------|--------------|-----------|---------|
| 0.1x            | 0.105    | 65           | 619.34    | 31      |
| 0.5x            | 0.488    | 481          | 985.77    | 213     |
| 1x              | 3.454    | 2,918        | 844.84    | 1,255   |
| 2x              | 229.106  | 216,243      | 943.86    | 90,548  |

**Key Findings:**
- Processing rate remains fairly constant at ~850-1000 events/ms
- At 2x rates, gene families explode to 90K+ leaves
- Higher multipliers (5x, 10x) cause exponential growth in gene family size

## 4. Transfer Event Analysis

Varying transfer rates to study horizontal gene transfer patterns
Species tree: 15 species, λ_d=1.0, λ_l=0.3, 50 replicates

| λ_transfer | Time(ms) | Transfers | T/(S+D+T) | Leaves | T/Leaf |
|------------|----------|-----------|-----------|--------|--------|
| 0.0        | 0.244    | 0         | 0.000     | 116    | 0.000  |
| 0.1        | 0.275    | 6         | 0.047     | 128    | 0.047  |
| 0.5        | 0.897    | 83        | 0.192     | 382    | 0.217  |
| 1.0        | 3.327    | 545       | 0.323     | 1,524  | 0.358  |
| 2.0        | 46.726   | 10,182    | 0.491     | 19,216 | 0.530  |

**Key Findings:**
- Transfer rate directly correlates with gene family expansion
- At λ_t=2.0, transfers make up ~49% of all diversification events
- High transfer rates lead to extensive gene tree growth
- Transfer efficiency (T/Leaf) increases with transfer rate

## 5. Gene Loss Impact Analysis

Varying loss rates to study gene family extinction
Species tree: 10 species, λ_d=1.0, λ_t=0.5, 100 replicates

| λ_loss | Time(ms) | Losses | Extant | Extinct(%) | L/(D+T+L) |
|--------|----------|--------|--------|------------|-----------|
| 0.0    | 0.863    | 0      | 480    | 0.0        | 0.000     |
| 0.1    | 0.601    | 13     | 324    | 0.0        | 0.063     |
| 0.3    | 0.368    | 24     | 176    | 1.0        | 0.168     |
| 0.5    | 0.209    | 22     | 84     | 2.0        | 0.260     |
| 1.0    | 0.069    | 12     | 17     | 4.0        | 0.422     |
| 2.0    | 0.021    | 4      | 2      | 21.0       | 0.600     |
| 5.0    | 0.011    | 1      | 0      | 48.0       | 0.783     |

**Key Findings:**
- High loss rates dramatically reduce simulation time (fewer active lineages)
- At λ_l=5.0, 48% of gene families go extinct before reaching any extant species
- Loss becomes the dominant event at high rates (78.3% of all DTL events at λ_l=5.0)
- Loss acts as a natural limiter on gene family size
- Average extant genes drops from 480 (no loss) to 0 (extreme loss)

## Performance Characteristics

### Strengths:
1. **Fast for typical scenarios**: Sub-millisecond performance for small trees with moderate rates
2. **Consistent event processing**: ~850-1000 events/ms regardless of event types
3. **Memory efficient**: No memory issues observed up to 1.5M events
4. **Scales well**: Linear scaling with number of events

### Limitations:
1. **Exponential growth with high rates**: Combined high D+T rates can cause gene families to explode
2. **Stochastic variation**: Performance varies significantly due to random tree topology
3. **Large trees**: 100+ species with high rates can take seconds per simulation

### Recommendations:
- For typical phylogenetic studies (10-50 species, moderate rates): Expect <100ms per simulation
- For large-scale simulations (100+ species): Consider reducing rates or using batch processing
- High transfer rates (>2.0) should be used cautiously as they cause rapid gene family expansion
- Loss events naturally limit computation time - consider including moderate loss rates

## Hardware & Configuration

- Test environment: Standard developer machine
- Rust: debug build (unoptimized)
- RNG: StdRng with fixed seeds for reproducibility
- All benchmarks use the same species tree topology per configuration

## Running the Benchmarks

```bash
# Quick test (always runs)
cargo test --test dtl_benchmark benchmark_dtl_quick_test -- --nocapture

# Comprehensive benchmarks (require --ignored flag)
cargo test --test dtl_benchmark benchmark_dtl_varying_rates -- --ignored --nocapture
cargo test --test dtl_benchmark benchmark_dtl_loss_impact -- --ignored --nocapture

# Warning: These can take several minutes
cargo test --test dtl_benchmark benchmark_dtl_varying_species_tree_size -- --ignored --nocapture
cargo test --test dtl_benchmark benchmark_dtl_scalability -- --ignored --nocapture
cargo test --test dtl_benchmark benchmark_dtl_transfer_intensity -- --ignored --nocapture
```
