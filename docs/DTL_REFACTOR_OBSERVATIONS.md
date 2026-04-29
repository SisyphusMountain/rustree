# DTL Refactor Observations

## Summary

The DTL simulation code already has a useful separation between public entry
points, streaming, mutable simulation state, and the shared Gillespie loop.
The next improvement is to make the boundary code less duplicated and less
panic-prone while preserving the existing simulation behavior.

## Observations

1. `per_gene.rs` and `per_species.rs` duplicate the same setup steps:
   validate rates, clone the species tree, generate species events, compute the
   time subdivision, compute contemporaneity, optionally precompute LCA depths,
   and construct `DtlSimIter`.

2. `DTLConfig` is a public bag of fields. Validation is split between DTL core,
   Python helpers, R helpers, and binding-common functions. This makes it easy
   for Rust and binding paths to drift, especially for `replacement_transfer`.

3. Several DTL paths still return `Result<_, String>` even though the crate has
   a shared `RustreeError` type. Typed errors make Python/R conversion clearer
   and avoid losing error categories.

4. Assortative-transfer setup can panic through `precompute_lca_depths().expect`
   in both core DTL utility code and Python iterator setup. This should be a
   recoverable simulation/setup error.

5. Python's DTL iterator retries indefinitely when `require_extant` is true,
   while Rust's `DtlSimIter` has a retry cap. The two paths should share the
   same failure behavior.

6. `SimulationState::genes_per_species` is an `Option`, but it is always
   initialized as `Some`. Removing the `Option` simplifies state operations and
   removes an impossible internal-error branch.

7. `finalize_simulation` currently relies on invariants and a panicking
   `RecTree::new` constructor. A checked `RecTree` constructor lets DTL
   finalization return recoverable errors at that boundary.

## Implemented Direction

- Add validated `DTLConfig` construction.
- Add a shared DTL preparation helper used by Rust and Python entry points.
- Migrate DTL simulation APIs to `RustreeError`.
- Propagate LCA precompute errors instead of panicking.
- Align Python iterator retry behavior with Rust iterator behavior.
- Simplify `SimulationState` by removing always-present optional state.
- Add checked `RecTree` constructors and use them from DTL finalization.
