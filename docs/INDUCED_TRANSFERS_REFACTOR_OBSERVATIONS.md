# Induced Transfers Refactor Observations

## Summary

`src/induced_transfers.rs` had two public behaviors sharing the same sampled-tree projection setup:

- `ghost_lengths`
- `induced_transfers` / `induced_transfers_with_algorithm`

The simple projection path and the Damien-style path both needed the same sampled subtree, `NodeMark` state, complete-tree projection, and complete-to-sampled index mapping. The setup was duplicated and used node names to translate kept complete-tree nodes back to sampled-tree nodes.

## Refactor Opportunities

1. Centralize sampled projection setup.
   A single helper should extract the sampled tree, compute marks, compute complete-tree projections, and expose a complete-to-sampled mapping.

2. Prefer index mappings over name lookups.
   `extract_induced_subtree_by_names` already returns an old-to-new mapping. Using that mapping avoids incorrect behavior for unnamed or duplicate internal node names.

3. Validate public event indices.
   Transfer events can come from external bindings. Invalid donor or recipient indices should return a `RustreeError::Index` instead of panicking on vector indexing.

4. Split Damien-style logic into phases.
   The Damien-style implementation mixed event extraction, extended graph construction, extinction masking, sampled projection, path walking, and detectability filtering in one function. Separating those phases makes invariants easier to test and maintain.

5. Make signed transfer-node semantics explicit.
   Damien-style graph nodes used raw signed `i64` values, where negative transfer nodes represent recipient-side split nodes. A small typed wrapper documents and constrains that convention without changing the algorithm.

6. Avoid silent invariant fallbacks.
   Missing transfer metadata should be returned as an error instead of silently assigning `(time=0.0, gene_id=0)`.

7. Remove arbitrary remapping bounds.
   Transfer-node remapping used a hard-coded hop limit. A visited set is clearer and detects cycles directly.

## Expected Benefits

- Fewer duplicate paths for projection-related bugs.
- Safer handling of user-supplied or binding-supplied transfer events.
- More robust sampled-node mapping in trees with duplicate or empty internal names.
- Damien-style code with clearer ownership boundaries and explicit graph invariants.
