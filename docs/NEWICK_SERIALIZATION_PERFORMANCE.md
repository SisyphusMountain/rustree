# Newick Serialization Performance Notes

This note summarizes the Newick serialization optimization work for large
birth-death trees, using a pure-birth species tree with:

- extant species: `1_000_000`
- birth rate: `lambda = 1.0`
- death rate: `mu = 0.0`
- total nodes: `1_999_999`
- Newick output size: about `33.9 MB`

Timings below are release-mode measurements on the local development machine.
They should be treated as directional rather than portable benchmark claims.

## Implemented Changes

`FlatTree::to_newick()` previously converted the flat tree into a recursive
`Node` tree and then serialized that recursive copy. The optimized path now:

- Serializes directly from `FlatTree` indices.
- Builds a preallocated byte buffer and converts it to `String` at the end.
- Uses a fast fixed-six-decimal branch-length formatter for common finite
  nonnegative branch lengths.
- Falls back to Rust's standard `"{:.6}"` formatter for unusual values and
  half-way rounding cases, preserving byte-for-byte compatibility.
- Uses a three-digit lookup table for faster fractional decimal output.
- Adds `FlatTree::write_newick<W: Write>()`.
- Updates Python and R `save_newick` bindings to use the faster writer path.
- Adds tests comparing optimized flat-tree output to the previous
  `tree.to_node().to_newick()` behavior.

Measured effect:

| Path | Time |
| --- | ---: |
| Previous `FlatTree -> Node -> Newick` construction | `~0.68-0.74 s` |
| Optimized Newick string construction | `~0.226 s` |
| Optimized construction plus file write | `~0.23 s` |

End-to-end generation plus Newick save for the benchmark tree dropped from
roughly `~0.86 s` to roughly `~0.36-0.39 s`.

## Verification

The optimized output was checked against the legacy path with byte-for-byte
comparisons on:

- parsed Newick trees
- simulated pure-birth trees
- simulated birth-death trees
- branch lengths near half-way fixed-decimal rounding cases
- negative and infinite branch-length fallback cases
- single-child error behavior

Verification commands run:

```sh
cargo test
cargo check --features python
cargo check --features r
```

## Profiling Results

Apple `sample` was used on a temporary loop that repeatedly serialized the same
million-species tree. Top-of-stack samples for the optimized serializer were
approximately:

| Area | Share |
| --- | ---: |
| `write_flat_newick_vec_recursive` traversal/output assembly | `~67%` |
| `push_length_fixed6` branch-length formatting | `~21%` |
| `memmove` / `memcpy` | `~11%` |
| UTF-8 validation/conversion | `<1%` |

Component benchmarks showed that flat, contiguous work is cheap in isolation:

| Variant | Time |
| --- | ---: |
| Current library construction | `~0.226 s` |
| Recursive copy without library `Result/String` overhead | `~0.216 s` |
| Iterative traversal | `~0.230 s` |
| Numeric labels from node indices | `~0.188 s` |
| No branch lengths | `~0.165 s` |
| No names | `~0.129 s` |
| Shape only | `~0.068 s` |
| Flat loop formatting all branch lengths | `~0.008 s` |
| Flat loop copying all names | `~0.005 s` |

The main remaining single-thread cost is not simple copying or decimal
formatting alone. It is the combination of tree-order traversal, non-contiguous
node access, and accessing millions of heap-allocated node-name strings.

## Tested But Not Implemented

### Numeric Label Fast Path

For birth-death simulated trees, node names are numeric strings matching node
indices. Emitting labels from the node index instead of reading `node.name`
reduced construction from about `~0.216 s` to `~0.188 s`.

This is promising but needs API care, because user-provided trees and renamed
nodes must still use stored names.

### Newick-Order Operation Vector

A temporary prototype copied the traversal into a contiguous operation vector in
the exact order needed for Newick output, then serialized from that vector.

| Variant | Time |
| --- | ---: |
| Current recursive full construction | `~0.216 s` |
| Build ordered ops only | `~0.075 s` |
| Write from ordered ops only | `~0.048 s` |
| Build ordered ops plus write, stored names | `~0.123 s` |
| Build ordered ops plus write, numeric labels | `~0.100 s` |

This confirms that improving locality can materially help. The tradeoff is
extra temporary memory: about four million operations for a two-million-node
binary tree.

### Parallel Subtree Serialization

A temporary Rayon prototype serialized independent subtrees in parallel and
stitched the resulting chunks in Newick order. Output was byte-identical in the
prototype.

| Split depth | Time |
| --- | ---: |
| 2 | `~0.094 s` |
| 3 | `~0.086 s` |
| 4 | `~0.074 s` |
| 5 | `~0.045 s` |
| 6 | `~0.034 s` |

This is the largest identified opportunity. It should be implemented carefully
with bounded split depth, predictable memory use, and fallback to the
single-thread path for small trees.

## Promising Opportunities Not Yet Tested

- A tree-only birth-death simulator that skips `Vec<TreeEvent>` allocation when
  events are not needed.
- Avoiding per-node `String` allocations for simulated numeric names.
- Storing names in an arena/string table to improve locality.
- Combining a numeric-label fast path with parallel subtree serialization.
- Combining a Newick-order operation vector with parallel chunking.
- Adding a non-Newick binary/compact export format for near-write-speed output.

## Recommended Next Step

The highest-impact next experiment is a production `to_newick_parallel` or
parallel-backed `write_newick` path for large trees:

1. Choose a split depth based on tree size and Rayon thread count.
2. Serialize subtrees independently into byte buffers.
3. Stitch buffers in deterministic Newick order.
4. Keep the current single-thread serializer for small trees and as a fallback.
5. Add byte-equivalence tests against the current optimized serializer.

Expected target for the benchmark tree: `~35-50 ms` Newick construction on the
local machine, before filesystem write time.
