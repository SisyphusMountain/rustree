# Pairwise Distances Implementation for Python Bindings

## Summary

Successfully implemented the `pairwise_distances` method for `PySpeciesTree` in the Python bindings (`/home/enzo/Documents/git/WP2/rustree/src/python.rs`).

## Implementation Details

### Method Signature

```rust
#[pyo3(signature = (distance_type, leaves_only=true))]
fn pairwise_distances(&self, py: Python, distance_type: &str, leaves_only: bool) -> PyResult<PyObject>
```

### Python API

```python
tree.pairwise_distances(distance_type: str, leaves_only: bool = True) -> pd.DataFrame
```

### Parameters

1. **`distance_type`** (str): Type of distance metric to compute
   - Accepted values: "topological", "topo", "metric", "branch", "patristic"
   - "topological"/"topo": Count of edges between nodes
   - "metric"/"branch"/"patristic": Sum of branch lengths between nodes

2. **`leaves_only`** (bool, default=True): Whether to compute distances only between leaf nodes
   - `True`: Only compute distances between leaf nodes (n² pairs for n leaves)
   - `False`: Compute distances between all nodes (m² pairs for m total nodes)

### Return Value

Returns a pandas DataFrame with three columns:
- **`node1`**: Name of the first node
- **`node2`**: Name of the second node
- **`distance`**: Distance between node1 and node2

The DataFrame includes all pairs (including symmetric pairs and self-distances).

### Error Handling

- Raises `ValueError` if `distance_type` is not one of the accepted values
- Error message: "Invalid distance_type 'X'. Use 'topological' or 'metric'."

## Implementation Pattern

The implementation follows the pattern from the R bindings (`rustree/src/r.rs` lines 865-898):

1. Parse and validate the `distance_type` string
2. Map to `DistanceType` enum (Topological or Metric)
3. Call `crate::metric_functions::pairwise_distances`
4. Convert results to pandas DataFrame

## Bonus Feature

A companion method `save_pairwise_distances_csv` was also implemented (lines 289-317):

```python
tree.save_pairwise_distances_csv(filepath: str, distance_type: str, leaves_only: bool = True)
```

This method saves the pairwise distances directly to a CSV file without creating an intermediate DataFrame.

## Files Modified

1. **`/home/enzo/Documents/git/WP2/rustree/src/python.rs`**
   - Added `pairwise_distances` method to `PySpeciesTree` impl block (lines 218-263)
   - Added `save_pairwise_distances_csv` method (lines 265-317)

2. **`/home/enzo/Documents/git/WP2/rustree/tests/python/test_species_tree.py`**
   - Added 8 comprehensive test cases for `pairwise_distances`:
     - `test_pairwise_distances_metric_leaves`
     - `test_pairwise_distances_topological_leaves`
     - `test_pairwise_distances_all_nodes`
     - `test_pairwise_distances_self_distance_zero`
     - `test_pairwise_distances_symmetric`
     - `test_pairwise_distances_invalid_type`
     - `test_pairwise_distances_aliases`
     - `test_pairwise_distances_default_leaves_only`

3. **`/home/enzo/Documents/git/WP2/rustree/examples/test_pairwise_distances.py`** (created)
   - Comprehensive example script demonstrating all features
   - Shows usage with simulated and parsed trees
   - Demonstrates distance type aliases
   - Shows CSV export functionality

## Usage Examples

### Basic Usage

```python
import rustree

# Simulate a species tree
tree = rustree.simulate_species_tree(5, 1.0, 0.5, seed=42)

# Compute metric distances between leaves
df = tree.pairwise_distances("metric", leaves_only=True)
print(df)
```

### Using Different Distance Types

```python
# Topological distances (edge counts)
df_topo = tree.pairwise_distances("topological")

# Metric distances (branch lengths)
df_metric = tree.pairwise_distances("metric")
```

### Computing Distances for All Nodes

```python
# Include internal nodes
df_all = tree.pairwise_distances("metric", leaves_only=False)
```

### Saving to CSV

```python
# Save directly to CSV file
tree.save_pairwise_distances_csv("distances.csv", "metric", leaves_only=True)
```

### Parsed Tree Example

```python
# Parse a Newick tree
newick = "((A:1.0,B:1.0):1.0,C:2.0):0.0;"
tree = rustree.parse_species_tree(newick)

# Compute distances
df = tree.pairwise_distances("metric")

# Expected distances:
# A-B: 2.0 (1.0 + 1.0 through parent)
# A-C: 4.0 (1.0 + 1.0 + 2.0 through root)
# B-C: 4.0 (1.0 + 1.0 + 2.0 through root)
```

## Testing

The implementation includes comprehensive tests covering:
- Both distance types (topological and metric)
- Both leaves-only and all-nodes modes
- Self-distances (should be zero)
- Symmetry (d(A,B) == d(B,A))
- Invalid distance type error handling
- Distance type aliases
- Default parameter values

Run tests with:
```bash
cd /home/enzo/Documents/git/WP2/rustree
maturin develop --release
python tests/python/test_species_tree.py
```

## Compilation Status

The code compiles successfully with no errors:
```
Compiling rustree v0.1.0 (/home/enzo/Documents/git/WP2/rustree)
Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.60s
```

## Technical Notes

1. The `py: Python` parameter is automatically provided by PyO3 and should not be included in the signature macro
2. String slices are efficiently converted using `.as_str()` to avoid unnecessary allocations
3. The implementation uses case-insensitive matching for distance types via `.to_lowercase()`
4. The function returns a `PyObject` to allow seamless integration with pandas DataFrames
5. Error messages are user-friendly and suggest valid alternatives

## Dependencies

The implementation relies on:
- `crate::metric_functions::DistanceType` - Distance type enum
- `crate::metric_functions::PairwiseDistance` - Distance result struct
- `FlatTree::pairwise_distances()` - Core distance computation
- `pyo3` - Python bindings framework
- `pandas` (Python) - DataFrame creation
