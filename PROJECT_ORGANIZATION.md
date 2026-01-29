# Project Organization Analysis

## Current Structure

```
rustree/
├── src/
│   ├── lib.rs              (18 lines)   - Module declarations
│   ├── main.rs             (93 lines)   - CLI entry point
│   ├── node.rs            (650 lines)   - Core data structures
│   ├── dtl.rs             (773 lines)   - DTL simulation
│   ├── sampling.rs        (540 lines)   - Tree sampling operations
│   ├── surgery.rs         (338 lines)   - Tree surgery (SPR)
│   ├── bd.rs              (225 lines)   - Birth-death simulation
│   ├── metric_functions.rs(202 lines)   - Tree metrics
│   ├── debug.rs           (131 lines)   - Debug visualization
│   ├── comparison.rs       (46 lines)   - Tree comparison
│   └── newick/
│       ├── mod.rs
│       ├── newick.rs                    - Newick parsing
│       └── newick.pest                  - Newick grammar
└── examples/              (22+ benchmarks and examples)
```

---

## File-by-File Analysis

### 1. **node.rs** (650 lines) - NEEDS SPLITTING

**Contains**:
- `Node` struct (recursive tree representation)
- `FlatNode` struct (flat/array-based tree)
- `FlatTree` struct (collection of FlatNodes)
- `RecTree` struct (reconciled tree for DTL)
- `Event` enum (Speciation, Duplication, Transfer, Loss)
- `TraversalOrder` enum
- `NodeIter` and `FlatTreeIter` (tree iterators)
- Conversion functions: `node_to_flat()`, `flat_to_node()`
- Tree traversal implementations
- XML generation (`to_xml()`)

**Problems**:
- ⚠️ **Too many responsibilities**: Data structures, conversions, iterators, and XML export
- ⚠️ **Large file**: 650 lines makes it hard to navigate
- ⚠️ **Mixed concerns**: XML export doesn't belong with core data structures

**Recommendation**: Split into multiple files

---

### 2. **dtl.rs** (773 lines) - WELL ORGANIZED BUT COULD SPLIT

**Contains**:
- `DTLEvent` struct
- `GeneCopy` struct (internal)
- `simulate_dtl()` - Main simulation
- `simulate_dtl_internal()` - Internal implementation
- `simulate_dtl_batch()` - Batch simulation
- `simulate_dtl_dfs_batch()` - DFS version
- `save_events_to_csv()` - CSV export
- `count_extant_genes()` - Analysis function
- `count_events()` - Analysis function
- `find_time_index()` - Utility

**Problems**:
- ⚠️ **Analysis functions mixed with simulation**: `count_extant_genes()`, `count_events()` should be elsewhere
- ⚠️ **CSV export mixed with simulation**: I/O code should be separate
- ⚠️ **Utility functions**: `find_time_index()` could be in a utils module

**Recommendation**: Separate simulation logic from analysis and I/O

---

### 3. **sampling.rs** (540 lines) - TOO MANY RESPONSIBILITIES

**Contains**:
- `remove_node()` - Tree surgery
- `spr()` - Subtree pruning and regrafting
- `sample_random_leaves()` - Random sampling
- `sample_diversified_leaves()` - Diversified sampling
- `sample_clustered_leaves()` - Clustered sampling
- `find_deepest_nodes()` - Tree analysis
- `sample_all_gene_trees()` - Gene tree sampling
- `species_tree_sample_to_string()` - String formatting
- Many internal helper functions

**Problems**:
- ⚠️ **SPR operations** (`spr()`, `remove_node()`) overlap with `surgery.rs`
- ⚠️ **Three different sampling strategies** could be separate implementations
- ⚠️ **Analysis functions** (`find_deepest_nodes()`) mixed with sampling
- ⚠️ **String formatting** mixed with algorithms

**Recommendation**: Split sampling strategies, move SPR to surgery.rs

---

### 4. **surgery.rs** (338 lines) - CLEAR PURPOSE

**Contains**:
- `is_ancestor()` - Tree relationship check
- `detach()` - Detach subtree
- `attach()` - Attach subtree
- `spr_topology()` - SPR on topology

**Assessment**: ✓ **Well-organized** - Clear purpose, focused functionality

**Minor issue**: Some SPR operations are in `sampling.rs` instead

---

### 5. **bd.rs** (225 lines) - WELL ORGANIZED

**Contains**:
- `TreeEvent` struct
- `simulate_bd_tree()` - Birth-death simulation
- `save_events_to_csv()` - CSV export

**Problems**:
- ⚠️ **CSV export** should be in I/O module

**Assessment**: Mostly good, just I/O separation needed

---

### 6. **metric_functions.rs** (202 lines) - UNCLEAR NAME

**Contains**:
- `give_depth()` - Assign depths to tree

**Problems**:
- ⚠️ **Misleading name**: "metric_functions" suggests multiple metrics
- ⚠️ **Only one function**: Why is this a separate module?
- ⚠️ **Should be in node.rs** or a tree utilities module

**Recommendation**: Merge into tree utilities or node module

---

### 7. **comparison.rs** (46 lines) - SMALL BUT FOCUSED

**Contains**:
- `compare_nodes()` - Compare with branch lengths
- `compare_nodes_topology()` - Compare topology only

**Assessment**: ✓ **Good** - Small, focused, clear purpose

---

### 8. **debug.rs** (131 lines) - CLEAR PURPOSE

**Contains**:
- `diffed_flat_tree_table()` - Debug visualization
- Internal formatting helpers

**Assessment**: ✓ **Good** - Clear debugging purpose

---

### 9. **newick/** - WELL ORGANIZED

**Contains**:
- `parse_newick()` - Parse Newick strings
- Internal parsing helpers

**Assessment**: ✓ **Good** - Separate module for parsing

---

### 10. **examples/** - TOO MANY BENCHMARKS

**Contains**: 22+ benchmark files

**Problems**:
- ⚠️ **Many redundant benchmarks** from optimization session
- ⚠️ **Unclear naming**: `benchmark_1000_batch` vs `benchmark_batch_dtl`
- ⚠️ **No organization**: All in one flat directory

**Recommendation**: Clean up old benchmarks, organize into subdirectories

---

## Major Issues Summary

### 🔴 Critical Issues

1. **node.rs is too large (650 lines)** - Contains data structures, conversions, iterators, and XML export
2. **I/O scattered across modules** - CSV and XML export in multiple files
3. **Sampling and surgery overlap** - SPR operations in both `sampling.rs` and `surgery.rs`
4. **Analysis functions scattered** - `count_events()` in dtl.rs, `find_deepest_nodes()` in sampling.rs

### 🟡 Medium Issues

1. **metric_functions.rs has only one function** - Should be merged elsewhere
2. **Too many benchmark examples** - 22+ files, many redundant
3. **Utility functions scattered** - `find_time_index()` in dtl.rs, helpers everywhere

---

## Recommended New Structure

```
rustree/
├── src/
│   ├── lib.rs                      - Module declarations & re-exports
│   ├── main.rs                     - CLI entry point
│   │
│   ├── core/                       - Core data structures
│   │   ├── mod.rs
│   │   ├── node.rs                 - Node struct only
│   │   ├── flat_tree.rs            - FlatTree & FlatNode
│   │   ├── rec_tree.rs             - RecTree (DTL reconciliation)
│   │   ├── events.rs               - Event enum
│   │   └── iterators.rs            - NodeIter, FlatTreeIter
│   │
│   ├── simulation/                 - Tree simulation
│   │   ├── mod.rs
│   │   ├── birth_death.rs          - Birth-death simulation
│   │   ├── dtl.rs                  - DTL simulation core
│   │   └── dtl_batch.rs            - DTL batch variants (BFS/DFS)
│   │
│   ├── io/                         - Input/Output
│   │   ├── mod.rs
│   │   ├── newick.rs               - Newick parsing (move from newick/)
│   │   ├── csv.rs                  - CSV export (all CSV functions)
│   │   └── xml.rs                  - XML/RecPhyloXML export
│   │
│   ├── operations/                 - Tree operations
│   │   ├── mod.rs
│   │   ├── surgery.rs              - Tree surgery (SPR, detach, attach)
│   │   ├── sampling/               - Sampling strategies
│   │   │   ├── mod.rs
│   │   │   ├── random.rs           - Random sampling
│   │   │   ├── diversified.rs      - Diversified sampling
│   │   │   └── clustered.rs        - Clustered sampling
│   │   └── conversion.rs           - node_to_flat, flat_to_node
│   │
│   ├── analysis/                   - Tree analysis
│   │   ├── mod.rs
│   │   ├── metrics.rs              - Tree metrics (depth, etc.)
│   │   ├── counting.rs             - count_events, count_extant_genes
│   │   ├── comparison.rs           - compare_nodes functions
│   │   └── traversal.rs            - Tree traversal utilities
│   │
│   ├── utils/                      - Utilities
│   │   ├── mod.rs
│   │   └── binary_search.rs        - find_time_index, etc.
│   │
│   └── debug.rs                    - Debug visualization (keep as-is)
│
└── examples/
    ├── basic/                      - Basic usage examples
    │   ├── parse_newick.rs
    │   ├── simulate_bd.rs
    │   └── simulate_dtl.rs
    ├── benchmarks/                 - Performance benchmarks
    │   ├── dtl_batch.rs            - Main DTL benchmark
    │   ├── multicore.rs            - Multi-threading benchmark
    │   └── profiling.rs            - Profiling helper
    └── analysis/                   - Analysis examples
        ├── dtl_growth.rs
        └── event_counting.rs
```

---

## Specific Refactoring Steps

### Phase 1: Split node.rs (HIGHEST PRIORITY)

**Current**: node.rs (650 lines)

**Split into**:
1. `core/node.rs` - `Node` struct
2. `core/flat_tree.rs` - `FlatTree`, `FlatNode`
3. `core/rec_tree.rs` - `RecTree`
4. `core/events.rs` - `Event` enum
5. `core/iterators.rs` - Iterator implementations
6. `io/xml.rs` - XML export (move `to_xml()` methods)
7. `operations/conversion.rs` - `node_to_flat()`, `flat_to_node()`

**Benefit**: Each file <200 lines, clear responsibilities

---

### Phase 2: Consolidate I/O

**Move**:
- `bd.rs::save_events_to_csv()` → `io/csv.rs`
- `dtl.rs::save_events_to_csv()` → `io/csv.rs`
- `node.rs::to_xml()` → `io/xml.rs`
- `newick/` → `io/newick.rs`

**Benefit**: All I/O in one place, easier to maintain

---

### Phase 3: Split sampling.rs

**Current**: sampling.rs (540 lines)

**Split into**:
1. `operations/sampling/random.rs` - Random sampling
2. `operations/sampling/diversified.rs` - Diversified sampling
3. `operations/sampling/clustered.rs` - Clustered sampling
4. `operations/surgery.rs` - Move SPR functions (merge with existing)
5. `analysis/traversal.rs` - `find_deepest_nodes()`, etc.

**Benefit**: Each sampling strategy is independent, easier to understand

---

### Phase 4: Organize analysis functions

**Create** `analysis/` module:
- `analysis/metrics.rs` - `give_depth()` (from metric_functions.rs)
- `analysis/counting.rs` - `count_events()`, `count_extant_genes()` (from dtl.rs)
- `analysis/comparison.rs` - Move comparison.rs here
- `analysis/traversal.rs` - Tree traversal utilities

**Delete**: `metric_functions.rs` (only had one function)

---

### Phase 5: Clean up examples

**Delete redundant benchmarks**:
- Keep: `dtl_batch.rs`, `multicore.rs`, `profiling.rs`
- Archive/delete: All the individual optimization benchmarks

**Organize into subdirectories**:
- `examples/basic/` - Simple usage examples
- `examples/benchmarks/` - Performance tests
- `examples/analysis/` - Analysis examples

---

## Benefits of Reorganization

### ✓ Better Discoverability
- Clear module hierarchy
- File names match functionality
- Each file <300 lines

### ✓ Separation of Concerns
- Data structures separate from operations
- I/O separate from algorithms
- Analysis separate from simulation

### ✓ Easier Maintenance
- Smaller files easier to understand
- Related code grouped together
- Clear dependencies between modules

### ✓ Better Testing
- Each module can be tested independently
- Mock I/O for testing algorithms
- Test sampling strategies separately

### ✓ Future Extensions
- Easy to add new sampling strategies
- Easy to add new export formats
- Easy to add new analysis functions

---

## Migration Strategy

### Option A: Gradual (Recommended)
1. Start with Phase 1 (split node.rs)
2. Test thoroughly
3. Continue with Phases 2-5 over time

### Option B: All at Once
- Create new structure
- Move all code
- Update all imports
- Risk: Many merge conflicts if working with others

**Recommendation**: Use Option A - split one module at a time, test, then continue

---

## Current Assessment

### ✓ What's Good
- Clear separation between simulation types (bd.rs, dtl.rs)
- Debug module is well-focused
- Surgery module has clear purpose
- Comparison module is small and focused

### ⚠️ What Needs Work
- node.rs is too large (650 lines, too many responsibilities)
- sampling.rs is too large (540 lines, multiple responsibilities)
- I/O scattered across modules
- Analysis functions scattered
- Too many benchmark examples (22+)

### 🎯 Priority Order
1. **Split node.rs** (biggest pain point)
2. **Consolidate I/O** (makes testing easier)
3. **Split sampling.rs** (second biggest file)
4. **Organize analysis** (improves discoverability)
5. **Clean examples** (reduces clutter)

---

## Conclusion

The project is **moderately well-organized** but has grown organically. The main issues are:

1. **node.rs is doing too much** - needs to be split into multiple files
2. **I/O is scattered** - CSV and XML export in multiple modules
3. **sampling.rs has too many responsibilities** - needs modularization

**The good news**: The core algorithms are clean and focused. The refactoring is mostly about organizing existing code, not rewriting it.

**Estimated effort**:
- Phase 1 (split node.rs): 2-3 hours
- Phase 2 (consolidate I/O): 1-2 hours
- Phase 3 (split sampling.rs): 2-3 hours
- Phase 4 (organize analysis): 1 hour
- Phase 5 (clean examples): 30 minutes

**Total**: ~8-10 hours for complete reorganization
