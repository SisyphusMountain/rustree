# Pairwise Distance Tests

This directory contains comprehensive tests for pairwise distance functionality in rustree.

## Current Status

The pairwise distance functions (`pairwise_distances` and `save_pairwise_distances_csv`) are currently:
- ✅ **Implemented in Rust** (`src/metric_functions.rs`)
- ✅ **Exposed in R bindings** (`src/r.rs`)
- ❌ **Not yet exposed in Python bindings** (`src/python.rs`)

## Required Python Bindings

To make the tests in `test_pairwise_distances.py` work, the following methods need to be added to `PySpeciesTree` in `src/python.rs`:

### 1. `pairwise_distances` method

```rust
/// Compute all pairwise distances between nodes in the tree.
///
/// # Arguments
/// * `distance_type` - Type of distance: "topological" (edge count) or "metric" (branch length sum)
/// * `leaves_only` - If true, only compute distances between leaf nodes
///
/// # Returns
/// A pandas DataFrame with columns: node1, node2, distance
fn pairwise_distances(&self, py: Python, distance_type: &str, leaves_only: bool) -> PyResult<PyObject> {
    use crate::metric_functions::DistanceType;

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

### 2. `save_pairwise_distances_csv` method

```rust
/// Save pairwise distances to a CSV file.
///
/// # Arguments
/// * `filepath` - Path to the output CSV file
/// * `distance_type` - Type of distance: "topological" or "metric"
/// * `leaves_only` - If true, only include leaf node distances
fn save_pairwise_distances_csv(
    &self,
    filepath: &str,
    distance_type: &str,
    leaves_only: bool
) -> PyResult<()> {
    use crate::metric_functions::DistanceType;
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
        .map_err(|e| PyValueError::new_err(format!("Failed to create file: {}", e)))?;

    // Write header
    writeln!(file, "{}", crate::metric_functions::PairwiseDistance::csv_header())
        .map_err(|e| PyValueError::new_err(format!("Failed to write header: {}", e)))?;

    // Write data rows
    for dist in distances {
        writeln!(file, "{}", dist.to_csv_row())
            .map_err(|e| PyValueError::new_err(format!("Failed to write row: {}", e)))?;
    }

    Ok(())
}
```

### Required imports in python.rs

Add these imports at the top of `src/python.rs`:

```rust
use pyo3::types::PyDict;
use crate::metric_functions::{DistanceType, PairwiseDistance};
```

## Test Coverage

The test file `test_pairwise_distances.py` provides comprehensive coverage including:

### Basic Functionality Tests
- ✅ DataFrame structure and column validation
- ✅ Self-distance verification (should be 0)
- ✅ Symmetry checks (distance(A,B) == distance(B,A))
- ✅ Specific distance value verification on known trees

### Distance Type Tests
- ✅ Topological distances (edge count)
  - Correct edge counting between nodes
  - Integer values
  - All non-negative
- ✅ Metric distances (branch length sum)
  - Correct summation of branch lengths
  - Handles varying branch lengths
  - All non-negative
  - Generally >= topological distances

### Parameter Tests
- ✅ `leaves_only=True` (only leaf nodes)
- ✅ `leaves_only=False` (all nodes including internal)
- ✅ Comparison between the two modes

### CSV Export Tests
- ✅ File creation
- ✅ Content matches direct method output
- ✅ Proper CSV header
- ✅ File overwrite behavior
- ✅ Both distance types (topological and metric)

### Edge Cases
- ✅ Single node tree
- ✅ Two node tree
- ✅ Large trees (20+ leaves)
- ✅ Trees with varying branch lengths

### Error Handling
- ✅ Invalid distance type
- ✅ Invalid file path
- ✅ Case sensitivity checks

## Running the Tests

Once the Python bindings are added:

```bash
# Run all tests
pytest rustree/tests/python/test_pairwise_distances.py -v

# Run specific test
pytest rustree/tests/python/test_pairwise_distances.py::test_pairwise_distances_metric_basic -v

# Run with coverage
pytest rustree/tests/python/test_pairwise_distances.py --cov=rustree --cov-report=html
```

## Example Usage

Once implemented, the Python API will work like this:

```python
import rustree

# Create or load a tree
tree = rustree.simulate_species_tree(10, 1.0, 0.3, seed=42)

# Compute metric distances between leaves
distances = tree.pairwise_distances("metric", leaves_only=True)
print(distances.head())

# Output:
#   node1 node2  distance
# 0    S0    S0      0.00
# 1    S0    S1      2.45
# 2    S0    S2      3.21
# ...

# Compute topological distances (edge count)
topo_distances = tree.pairwise_distances("topological", leaves_only=True)

# Include internal nodes
all_distances = tree.pairwise_distances("metric", leaves_only=False)

# Save to CSV
tree.save_pairwise_distances_csv("distances.csv", "metric", leaves_only=True)
```

## Integration with Existing Code

The pairwise distance functionality integrates with existing rustree features:

- Works with both simulated trees (`simulate_species_tree`) and parsed trees (`parse_species_tree`)
- Complements other tree metrics like `tree_height()`, `num_leaves()`, etc.
- CSV export follows similar patterns to gene tree CSV export
- Uses the same `FlatTree` structure as other operations

## Performance Notes

- Complexity: O(n²) for n nodes (unavoidable for pairwise distances)
- For large trees with many nodes, use `leaves_only=True` to reduce computation
- Distance computation uses efficient LCA-based algorithms
- CSV writing is buffered for efficiency
