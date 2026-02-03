# Python Bindings TODO: Pairwise Distance Functions

This document provides step-by-step instructions for adding pairwise distance functionality to the Python bindings.

## Quick Summary

Add 2 methods to `PySpeciesTree` in `src/python.rs`:
1. `pairwise_distances()` - Returns pandas DataFrame
2. `save_pairwise_distances_csv()` - Saves to CSV file

## Step 1: Add Required Imports

At the top of `src/python.rs`, add these imports:

```rust
use pyo3::types::PyDict;
use crate::metric_functions::{DistanceType, PairwiseDistance};
```

## Step 2: Add Methods to PySpeciesTree

In the `#[pymethods] impl PySpeciesTree` block (around line 27), add these two methods:

### Method 1: pairwise_distances

```rust
/// Compute all pairwise distances between nodes in the tree.
///
/// # Arguments
/// * `distance_type` - Type of distance: "topological" (edge count) or "metric" (branch length sum)
/// * `leaves_only` - If true, only compute distances between leaf nodes
///
/// # Returns
/// A pandas DataFrame with columns: node1, node2, distance
///
/// # Example
/// ```python
/// import rustree
///
/// tree = rustree.simulate_species_tree(10, 1.0, 0.3, seed=42)
/// distances = tree.pairwise_distances("metric", leaves_only=True)
/// print(distances.head())
/// ```
fn pairwise_distances(&self, py: Python, distance_type: &str, leaves_only: bool) -> PyResult<PyObject> {
    let dist_type = match distance_type.to_lowercase().as_str() {
        "topological" => DistanceType::Topological,
        "metric" => DistanceType::Metric,
        _ => return Err(PyValueError::new_err(
            format!("Invalid distance type '{}'. Must be 'topological' or 'metric'", distance_type)
        )),
    };

    let distances = self.tree.pairwise_distances(dist_type, leaves_only);

    // Convert to pandas DataFrame
    let pandas = py.import("pandas")?;
    let dict = PyDict::new(py);

    let node1_list: Vec<String> = distances.iter().map(|d| d.node1.clone()).collect();
    let node2_list: Vec<String> = distances.iter().map(|d| d.node2.clone()).collect();
    let distance_list: Vec<f64> = distances.iter().map(|d| d.distance).collect();

    dict.set_item("node1", node1_list)?;
    dict.set_item("node2", node2_list)?;
    dict.set_item("distance", distance_list)?;

    let df = pandas.call_method1("DataFrame", (dict,))?;
    Ok(df.into())
}
```

### Method 2: save_pairwise_distances_csv

```rust
/// Save pairwise distances to a CSV file.
///
/// # Arguments
/// * `filepath` - Path to the output CSV file
/// * `distance_type` - Type of distance: "topological" or "metric"
/// * `leaves_only` - If true, only include leaf node distances
///
/// # Example
/// ```python
/// import rustree
///
/// tree = rustree.simulate_species_tree(10, 1.0, 0.3, seed=42)
/// tree.save_pairwise_distances_csv("distances.csv", "metric", leaves_only=True)
/// ```
fn save_pairwise_distances_csv(
    &self,
    filepath: &str,
    distance_type: &str,
    leaves_only: bool
) -> PyResult<()> {
    use std::fs::File;
    use std::io::Write;

    let dist_type = match distance_type.to_lowercase().as_str() {
        "topological" => DistanceType::Topological,
        "metric" => DistanceType::Metric,
        _ => return Err(PyValueError::new_err(
            format!("Invalid distance type '{}'. Must be 'topological' or 'metric'", distance_type)
        )),
    };

    let distances = self.tree.pairwise_distances(dist_type, leaves_only);

    let mut file = File::create(filepath)
        .map_err(|e| PyValueError::new_err(format!("Failed to create file '{}': {}", filepath, e)))?;

    // Write header
    writeln!(file, "{}", PairwiseDistance::csv_header())
        .map_err(|e| PyValueError::new_err(format!("Failed to write header: {}", e)))?;

    // Write data rows
    for dist in distances {
        writeln!(file, "{}", dist.to_csv_row())
            .map_err(|e| PyValueError::new_err(format!("Failed to write row: {}", e)))?;
    }

    Ok(())
}
```

## Step 3: Rebuild the Python Module

```bash
cd rustree
cargo build --release
```

The compiled module will be at: `target/release/librustree.so` (or `.dylib` on macOS, `.dll` on Windows)

## Step 4: Test the Implementation

```bash
# Run the test suite
pytest tests/python/test_pairwise_distances.py -v

# Run the demo
python tests/python/demo_pairwise_distances.py
```

## Expected Test Output

```
tests/python/test_pairwise_distances.py::test_pairwise_distances_topological_basic PASSED
tests/python/test_pairwise_distances.py::test_pairwise_distances_topological_self_distance PASSED
tests/python/test_pairwise_distances.py::test_pairwise_distances_topological_symmetry PASSED
...
======================= 35 passed in 2.45s =======================
```

## Verification Checklist

After implementation, verify:

- [ ] Both methods compile without errors
- [ ] Import `rustree` in Python successfully
- [ ] `tree.pairwise_distances()` returns a pandas DataFrame
- [ ] DataFrame has columns: `node1`, `node2`, `distance`
- [ ] Both "topological" and "metric" distance types work
- [ ] Both `leaves_only=True` and `leaves_only=False` work
- [ ] `tree.save_pairwise_distances_csv()` creates a valid CSV file
- [ ] Invalid distance type raises `ValueError`
- [ ] Invalid file path raises appropriate error
- [ ] All 35 tests in `test_pairwise_distances.py` pass
- [ ] Demo script runs without errors

## Example Usage After Implementation

```python
import rustree

# Simulate or parse a tree
tree = rustree.simulate_species_tree(10, 1.0, 0.3, seed=42)

# Get metric distances between leaves
df = tree.pairwise_distances("metric", leaves_only=True)
print(df)
#   node1 node2  distance
# 0    S0    S0      0.00
# 1    S0    S1      2.45
# 2    S0    S2      3.21
# ...

# Get topological distances (edge counts)
df_topo = tree.pairwise_distances("topological", leaves_only=True)

# Include internal nodes
df_all = tree.pairwise_distances("metric", leaves_only=False)

# Save to CSV
tree.save_pairwise_distances_csv("output.csv", "metric", leaves_only=True)
```

## Common Issues and Solutions

### Issue: "pandas not found"
**Solution:** The user needs pandas installed. This is a runtime dependency, not a compilation dependency.
```bash
pip install pandas
```

### Issue: Compilation errors about PyDict
**Solution:** Make sure you added the import: `use pyo3::types::PyDict;`

### Issue: Compilation errors about DistanceType
**Solution:** Make sure you added the import: `use crate::metric_functions::{DistanceType, PairwiseDistance};`

### Issue: Tests fail with "method not found"
**Solution:** Rebuild the module: `cargo build --release`

## File Locations

- **Python bindings:** `rustree/src/python.rs`
- **Rust implementation:** `rustree/src/metric_functions.rs`
- **Test suite:** `rustree/tests/python/test_pairwise_distances.py`
- **Demo script:** `rustree/tests/python/demo_pairwise_distances.py`
- **Documentation:** `rustree/tests/python/README_PAIRWISE_DISTANCES.md`

## Additional Notes

- The implementation already exists in Rust and R bindings
- This is just exposing existing functionality to Python
- No changes needed to `metric_functions.rs`
- The methods follow the same pattern as other `PySpeciesTree` methods
- Error handling uses PyO3's `PyResult` and `PyValueError`
- pandas DataFrame is created dynamically via Python interop
