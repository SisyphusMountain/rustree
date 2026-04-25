# Robustness, Test Coverage, and Modularity Plan

**Created**: 2026-04-25  
**Scope**: Current `rustree` codebase, including Rust core, Python bindings, R bindings, tests, and CI.  
**Goal**: Increase test coverage, reduce panic paths, make API boundaries more robust, and make module ownership easier to understand.

## 1. Current Baseline

The core Rust library is already in a healthy state mechanically:

- `cargo check` passes.
- `cargo test` passes.
- `cargo clippy --all-targets -- -D warnings` passes.
- `cargo check --features r` passes.
- `cargo test --features r` passes.
- `cargo check --features python` passes.
- `cargo clippy --features python -- -D warnings` passes.

Known local verification gaps:

- `cargo test --features python` fails on macOS at link time because the `python` feature enables PyO3's `extension-module` feature.
- Python pytest tests were not run locally because the active Python 3.13 environment did not have `pytest` installed.
- Some untracked user files are present in the working tree and should not be treated as part of this plan unless intentionally added later.

## 2. Highest Priority Weaknesses

### 2.1 Python feature is not cleanly testable with Cargo

**Files**

- `Cargo.toml`
- `pyproject.toml`
- `.github/workflows/ci.yml`

**Problem**

`pyo3` is declared with `features = ["extension-module"]` under the optional dependency. This is appropriate for building a Python extension with maturin, but it makes plain Cargo test builds fragile on platforms where the Python symbols are not linked into the test artifact. Locally, `cargo test --features python` failed with unresolved Python C API symbols.

CI currently includes `cargo test --features python`, so this configuration is likely either platform-sensitive or relying on Linux behavior that does not match macOS development.

**Plan**

1. Split Python features:
   - `python = ["pyo3", "numpy"]`
   - `python-extension-module = ["python", "pyo3/extension-module"]`
2. Remove `features = ["extension-module"]` from the `pyo3` dependency declaration.
3. Change `pyproject.toml` so maturin uses `python-extension-module`.
4. Change CI:
   - Use `cargo check --features python`, `cargo test --features python`, and `cargo clippy --features python` for Rust-side Python binding tests.
   - Use `maturin build --features python-extension-module` for wheel builds.
5. Add a short comment in `Cargo.toml` explaining why the features are split.

**Acceptance Criteria**

- `cargo test --features python` passes locally on macOS.
- `maturin build --features python-extension-module` still produces an importable wheel.
- CI clearly distinguishes Rust test builds from extension-module packaging builds.

## 3. Panic and Error Handling Weaknesses

### 3.1 `RecTree` constructors panic on invalid mappings

**File**

- `src/node/rectree.rs`

**Problem**

`RecTree::new` and `RecTree::with_dtl_events` use `assert_eq!` and `assert!` to validate mapping lengths and species indices. These constructors sit at a central boundary between simulation, parsing, sampling, Python, and R. A bad mapping should become a recoverable error, not a panic that can crash Python or R.

**Plan**

1. Add checked constructors:
   - `RecTree::try_new(...) -> Result<Self, RustreeError>`
   - `RecTree::try_with_dtl_events(...) -> Result<Self, RustreeError>`
2. Move validation into a shared private helper:
   - gene tree length matches `node_mapping`
   - gene tree length matches `event_mapping`
   - every `Some(species_idx)` is in bounds
3. Keep current `new` constructors temporarily for internal call sites that are proven valid, but document them as unchecked or panicking constructors.
4. Gradually migrate parser and binding paths to checked constructors first:
   - RecPhyloXML parsing
   - Python `PyGeneTree` sampling methods
   - R conversions
   - `GeneForest` constructors
5. Add tests for each invalid constructor case.

**Acceptance Criteria**

- Invalid reconciliation mappings return `RustreeError::Validation` or `RustreeError::Index`.
- Python and R callers receive normal language-level errors.
- No public parser or binding path can panic because of mapping length mismatch.

### 3.2 Metric APIs still contain panic-based preconditions

**File**

- `src/metric_functions.rs`

**Problem**

Several metric functions still assume well-formed input:

- `LcaTable::new` asserts that the tree is non-empty.
- `LcaTable::lca` indexes directly into `first`.
- `make_intervals` can underflow when no depths are assigned.
- `find_contemporaneity` expects all depths to be assigned.
- `distance_between` and its helper paths expect intact parent chains.

The code is mostly fine for internally generated trees, but user-supplied trees and binding inputs should get structured errors.

**Plan**

1. Add checked variants:
   - `LcaTable::try_new(&FlatTree) -> Result<LcaTable, RustreeError>`
   - `LcaTable::try_lca(&self, u, v) -> Result<usize, RustreeError>`
   - `FlatTree::try_make_intervals() -> Result<Vec<f64>, RustreeError>`
   - `FlatTree::try_find_contemporaneity(...) -> Result<Vec<Vec<usize>>, RustreeError>`
2. Keep current fast methods for internal call sites if needed, but route binding-facing methods through checked variants.
3. Change `precompute_lca_depths` to use the checked LCA table.
4. Convert broken parent-chain cases to errors.
5. Add tests for:
   - empty tree
   - invalid root index
   - node index out of bounds
   - disconnected parent chain
   - missing depth
   - single-node tree

**Acceptance Criteria**

- Calling public distance APIs on malformed trees returns an error instead of panicking.
- Python and R pairwise-distance methods surface clear messages.
- Existing performance for valid trees does not materially regress.

### 3.3 DTL helper precomputation can panic

**Files**

- `src/simulation/dtl/utils.rs`
- `src/python/species_tree.rs`

**Problem**

Assortative transfer support calls `precompute_lca_depths().expect(...)` in DTL utility code and Python iterator creation. A malformed species tree or missing depth can panic instead of producing a recoverable simulation error.

**Plan**

1. Change `precompute_lca(...)` to return `Result<Option<Vec<Vec<f64>>>, String>` or `Result<..., RustreeError>`.
2. Update per-gene and per-species DTL iterator constructors to propagate the error.
3. Update Python iterator constructors to map the error to `PyValueError`.
4. Add tests for `transfer_alpha` with missing species depths.

**Acceptance Criteria**

- DTL simulation with invalid species tree metadata fails with a descriptive error.
- No `expect("Failed to precompute LCA depths")` remains in binding-facing paths.

## 4. Test Coverage Weaknesses

### 4.1 Python sampling tests depend heavily on random outcomes and skips

**File**

- `tests/python/test_tree_sampling.py`

**Problem**

Many Python tests skip when simulated trees do not contain enough extant genes. This makes the suite broad but less deterministic. A test suite can pass while skipping large parts of the sampling behavior.

**Plan**

1. Keep a small number of simulation-based smoke tests.
2. Add deterministic fixtures:
   - fixed species Newick strings
   - fixed RecPhyloXML strings
   - hand-built small reconciled trees from XML fixtures
3. For sampling tests, assert exact expected output:
   - leaf names
   - Newick topology
   - event counts
   - node mappings
   - species mappings after pruning
4. Replace permissive assertions like "works or raises" with explicit expected behavior.
5. Add pytest markers:
   - `unit`
   - `integration`
   - `slow`
   - `external`

**Acceptance Criteria**

- Core Python sampling tests do not skip because of random simulation shape.
- Sampling behavior is asserted against deterministic topology and mapping fixtures.
- Random simulation tests are clearly labeled as smoke or integration tests.

### 4.2 Python tests are not validated locally without installing a wheel

**Files**

- `pyproject.toml`
- `.github/workflows/ci.yml`
- `tests/python/*`

**Problem**

The Python tests require a built/importable module and local dependencies. There is no simple documented command that sets up the module and runs the suite from a clean checkout.

**Plan**

1. Add a short `docs/TESTING.md` or update an existing docs page with:
   - Rust core test commands
   - Rust feature test commands
   - Python wheel build and pytest commands
   - R binding build and test commands
2. Consider adding a `justfile`, `Makefile`, or `scripts/test_python.sh` if the project wants command wrappers.
3. Ensure the Python tests do not manually mutate `sys.path` to `target/release` when a wheel install is expected.

**Acceptance Criteria**

- A fresh developer can run the Python tests from documented commands.
- CI and local instructions use the same build mode.

### 4.3 Test suite has no coverage measurement gate

**Files**

- `.github/workflows/ci.yml`
- `pyproject.toml`

**Problem**

There are many tests, but no coverage job or baseline. Without a baseline, it is hard to prioritize coverage work or detect accidental drops.

**Plan**

1. Add a non-blocking Rust coverage job first:
   - Prefer `cargo llvm-cov` if adopted.
   - Generate text and lcov reports.
2. Add Python coverage with `pytest-cov`.
3. Publish coverage as CI artifacts.
4. After two or three stable runs, set minimum thresholds:
   - initially low enough to avoid churn
   - ratchet upward only when coverage is intentionally added

**Acceptance Criteria**

- CI produces Rust and Python coverage reports.
- Coverage is visible without needing local tooling.
- Future PRs can identify uncovered modules before merging.

### 4.4 R tests are split between Rust unit tests and standalone R scripts

**Files**

- `src/r/mod.rs`
- `tests/r/*`
- `R/README.md`

**Problem**

R binding coverage exists, but the workflow is less unified than Rust and Python. Some R tests are standalone scripts that expect a release shared library. This makes them easy to skip unintentionally.

**Plan**

1. Document the exact R test workflow.
2. Add a CI job for `tests/r` if it is not already covered.
3. Standardize temporary output paths in R tests.
4. Add small deterministic tests for:
   - invalid seeds
   - invalid DTL rates
   - replacement transfer parameter
   - RecPhyloXML parse errors
   - sampling with missing names

**Acceptance Criteria**

- R tests can run from a clean checkout with one documented command.
- R binding failures are caught in CI, not only through Rust `#[cfg(feature = "r")]` tests.

## 5. Test Hygiene Weaknesses

### 5.1 Rust tests write files in the repository root

**File**

- `tests/bd_tests.rs`

**Problem**

Some tests write files such as `bd_tree_basic.nwk`, `bd_tree_basic_events.csv`, `bd_tree_pure_birth.nwk`, and `bd_tree_pure_birth_events.csv` into the repository root and then remove them best-effort. This can leave artifacts after failure and can race under parallel test execution.

**Plan**

1. Replace root writes with `tempfile::tempdir()`.
2. Write all generated Newick and CSV files under the temp directory.
3. Remove println-heavy debug output from normal unit tests unless the output is part of failure diagnosis.
4. Add a quick check that `cargo test` leaves no generated files in the repository root.

**Acceptance Criteria**

- `tests/bd_tests.rs` creates no root-level files.
- Re-running `cargo test` after an interrupted test run does not require cleanup.

### 5.2 Some ignored tests require external private data

**Files**

- `tests/test_alerax_real_file.rs`
- `tests/test_real_separate_files.rs`

**Problem**

The ignored tests are appropriate for external datasets, but they do not help default CI coverage. The important parsing behavior should be represented by small committed fixtures.

**Plan**

1. Keep the ignored tests for large external validation.
2. Add minimal committed fixtures that exercise the same parser branches.
3. Ensure default `cargo test` covers:
   - real-ish separate species/gene files
   - missing species mappings
   - malformed file errors
   - multi-family parsing if relevant

**Acceptance Criteria**

- Default CI covers the parser behaviors currently only exercised by external-data tests.
- External-data tests remain documented for manual validation.

## 6. Modularity and API Design Weaknesses

### 6.1 Error handling is still split between `String` and `RustreeError`

**Files**

- `src/error.rs`
- `src/simulation/*`
- `src/metric_functions.rs`
- `src/node/gene_forest.rs`
- `src/external/alerax.rs`

**Problem**

The library has started migrating to `RustreeError`, but many functions still return `Result<_, String>`. This makes error categorization inconsistent and forces bindings to treat all failures as generic value errors or runtime errors.

**Plan**

1. Migrate central internal modules first:
   - metrics
   - DTL simulation
   - BD event generation
   - GeneForest sampling
2. Add variants if needed:
   - `RustreeError::Invariant`
   - `RustreeError::MissingData`
   - `RustreeError::ExternalTool`
3. Keep external tool parsing errors descriptive, but convert them before crossing public boundaries.
4. Update Python and R conversions to map categories to useful host-language errors.

**Acceptance Criteria**

- New public Rust APIs return `RustreeError`, not `String`.
- Binding layers no longer need ad hoc string parsing to choose error types.

### 6.2 Binding modules contain duplicated parsing and CSV logic

**Files**

- `src/python/species_tree.rs`
- `src/python/gene_tree.rs`
- `src/r/analysis.rs`
- `src/r/species.rs`
- `src/bindings_common/mod.rs`

**Problem**

Some common logic has already moved into `bindings_common`, but duplication remains. For example, distance type parsing is centralized for some methods but still duplicated in CSV export paths.

**Plan**

1. Move all distance type parsing through `bindings_common::parse_distance_type`.
2. Move shared CSV serialization helpers out of language bindings.
3. Keep Python/R wrappers thin:
   - validate host-language inputs
   - call Rust core functions
   - convert return values
4. Add parity tests that compare Python, R, and Rust outputs for shared operations.

**Acceptance Criteria**

- No duplicate distance-type parser remains in Python or R modules.
- Binding wrappers are mostly host-language conversion code.

### 6.3 Visualization methods use fixed temp file names

**Files**

- `src/python/species_tree.rs`
- `src/python/gene_tree.rs`
- `src/r/viz.rs`

**Problem**

Visualization helpers write fixed filenames into the system temp directory, such as `rustree_species_temp.nwk` and `rustree_temp.svg`. Parallel calls can collide, and failed calls can leave stale files behind.

**Plan**

1. Use `tempfile::tempdir()` or `NamedTempFile` for all visualization intermediates.
2. Keep paths scoped to the temporary directory lifetime.
3. Add tests that call visualization path-building code in parallel if thirdkind can be mocked or isolated.

**Acceptance Criteria**

- Parallel visualization calls cannot overwrite each other's temporary files.
- Temporary files are cleaned up reliably.

## 7. Behavioral Robustness Weaknesses

### 7.1 Partial invalid-name handling is inconsistent

**Files**

- `src/python/gene_tree.rs`
- `src/python/species_tree.rs`
- `src/sampling.rs`
- `tests/python/test_tree_sampling.py`

**Problem**

Some sampling APIs reject no matches but silently ignore partially invalid name lists. Tests currently document this as acceptable in at least one Python case. Silent partial success can hide data issues.

**Plan**

1. Decide explicit API behavior:
   - strict mode: error if any requested name is missing
   - permissive mode: ignore missing names but return/report ignored names
2. Prefer strict default for binding APIs.
3. Add optional permissive behavior only if there is a real user workflow.
4. Update tests to assert the chosen behavior.

**Acceptance Criteria**

- Sampling APIs have documented behavior for missing and duplicate names.
- Python and R behavior matches Rust behavior.

### 7.2 Tree structural invariants are implicit

**Files**

- `src/node/mod.rs`
- `src/node/conversion.rs`
- `src/metric_functions.rs`
- `src/sampling.rs`

**Problem**

Many functions assume `FlatTree` is well-formed:

- root index is valid
- parent pointers are consistent
- child pointers are in bounds
- every non-root node has one parent
- no cycles exist
- binary child shape is respected

These invariants are not centralized, so each module either trusts the tree or handles only local symptoms.

**Plan**

1. Add `FlatTree::validate() -> Result<(), RustreeError>`.
2. Validate:
   - non-empty nodes
   - root in bounds
   - child indices in bounds
   - parent-child consistency
   - no cycles
   - all nodes reachable from root
   - binary child shape
   - optional depth consistency
3. Use validation in parser tests, constructors, and binding-facing operations.
4. Add malformed-tree tests built manually in Rust.

**Acceptance Criteria**

- Public operations can cheaply validate user-constructed trees before heavy computation.
- Error messages identify the broken invariant and node index.

## 8. Performance and Scale Weaknesses

### 8.1 Pairwise outputs are inherently O(n^2)

**Files**

- `src/metric_functions.rs`
- `src/python/species_tree.rs`
- `src/python/gene_tree.rs`

**Problem**

Pairwise distances and LCA-depth matrices allocate O(n^2) memory. This is expected mathematically, but binding APIs can accidentally create very large pandas DataFrames or matrices.

**Plan**

1. Document O(n^2) behavior in Rust and binding docstrings.
2. Add optional limits or warnings in Python/R for very large trees.
3. Add streaming pairwise-distance iterator or CSV writer for large outputs.
4. Add performance tests for medium-sized trees to detect major regressions.

**Acceptance Criteria**

- Users get clear feedback before generating enormous pairwise outputs.
- Large-output workflows have a streaming option.

### 8.2 External-tool integrations need clearer boundary tests

**Files**

- `src/external/alerax.rs`
- `src/python/alerax.rs`
- `src/python/forest.rs`

**Problem**

ALERax integration touches file generation, command execution, parsing output files, and mapping back into `RecTree`. These boundaries are high-risk and should have focused tests that do not require the external binary.

**Plan**

1. Separate pure parsing from command execution.
2. Add fixtures for:
   - rates files
   - likelihood files
   - event count files
   - transfer files
   - missing or malformed output files
3. Add tests for generated config/families files.
4. Keep end-to-end ALERax tests ignored or external-only.

**Acceptance Criteria**

- Most ALERax parsing and file preparation logic is covered without installing ALERax.
- External binary tests are a small final integration layer.

## 9. Suggested Implementation Order

### Phase 1: Make CI and local test modes reliable

1. Split Python features into testable and extension-module packaging features.
2. Update CI and `pyproject.toml`.
3. Add or update testing documentation.
4. Replace root-writing Rust tests with temp directories.

**Why first**: This makes every later change easier to validate consistently.

### Phase 2: Remove the most dangerous panic paths

1. Add checked `RecTree` constructors.
2. Route parser and binding paths through checked constructors.
3. Add checked metric/LCA variants.
4. Propagate DTL LCA precompute errors.

**Why second**: These are the highest-impact robustness fixes for Python/R users.

### Phase 3: Make tests deterministic and higher signal

1. Add deterministic Python fixtures for sampling and mappings.
2. Add malformed-tree Rust tests.
3. Add committed fixtures replacing behavior hidden behind ignored external-data tests.
4. Add R workflow coverage and deterministic R binding tests.

**Why third**: It prevents future refactors from preserving only happy paths.

### Phase 4: Improve API consistency and modularity

1. Continue migrating `Result<_, String>` to `Result<_, RustreeError>`.
2. Move duplicated binding logic into `bindings_common` or core modules.
3. Add `FlatTree::validate()`.
4. Make visualization temp-file handling robust.

**Why fourth**: These changes are easier and safer after panic paths and tests are improved.

### Phase 5: Add coverage and performance guardrails

1. Add non-blocking Rust and Python coverage jobs.
2. Add performance regression tests for medium-sized trees.
3. Document or guard O(n^2) APIs.
4. Add streaming alternatives where needed.

**Why fifth**: Coverage and performance gates should reflect stable behavior, not churn during structural cleanup.

## 10. Tracking Checklist

- [ ] Split `python` and `python-extension-module` features.
- [ ] Update maturin build feature in `pyproject.toml`.
- [ ] Update CI Python feature commands.
- [ ] Add checked `RecTree` constructors.
- [ ] Route parsers and bindings through checked `RecTree` constructors.
- [ ] Add checked `LcaTable` and metric helpers.
- [ ] Make DTL LCA precompute errors recoverable.
- [ ] Replace repository-root test file writes with temp directories.
- [ ] Add deterministic Python sampling fixtures.
- [ ] Add malformed-tree Rust tests.
- [ ] Add small committed fixtures for external-data parser paths.
- [ ] Document Rust, Python, and R test workflows.
- [ ] Add coverage reporting jobs.
- [ ] Continue `String` to `RustreeError` migration.
- [ ] Add `FlatTree::validate()`.
- [ ] Replace fixed visualization temp filenames with unique temp files.
- [ ] Clarify strict vs permissive behavior for missing sample names.
- [ ] Add O(n^2) output warnings or streaming alternatives.

