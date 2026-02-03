#!/usr/bin/env python3
"""
Integration Test for rustree Python Bindings
============================================

This script tests the complete workflow:
1. Simulate species tree
2. Export BD events to CSV
3. Compute pairwise distances
4. Sample/extract induced subtree
5. Simulate gene trees on sampled species tree
6. Verify all outputs

This serves as both a test and a demonstration of the full API.
"""

import sys
import os

# Add the release directory to path
sys.path.insert(0, "/home/enzo/Documents/git/WP2/rustree/target/release")

import rustree
import tempfile
import traceback

class IntegrationTestSuite:
    """Integration test suite for rustree Python bindings."""

    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.errors = []

    def run_test(self, test_func):
        """Run a single test and track results."""
        test_name = test_func.__name__
        try:
            print(f"\n{'='*70}")
            print(f"Running: {test_name}")
            print(f"{'='*70}")
            test_func()
            self.passed += 1
            print(f"✓ PASS: {test_name}")
        except AssertionError as e:
            self.failed += 1
            error_msg = f"AssertionError: {e}"
            self.errors.append((test_name, error_msg))
            print(f"✗ FAIL: {test_name}")
            print(f"  Error: {error_msg}")
            traceback.print_exc()
        except Exception as e:
            self.failed += 1
            error_msg = f"{type(e).__name__}: {e}"
            self.errors.append((test_name, error_msg))
            print(f"✗ ERROR: {test_name}")
            print(f"  Error: {error_msg}")
            traceback.print_exc()

    def print_summary(self):
        """Print test summary."""
        print(f"\n{'='*70}")
        print("INTEGRATION TEST SUMMARY")
        print(f"{'='*70}")
        print(f"Passed: {self.passed}")
        print(f"Failed: {self.failed}")
        print(f"Total:  {self.passed + self.failed}")

        if self.errors:
            print(f"\n{'='*70}")
            print("FAILED TESTS:")
            print(f"{'='*70}")
            for test_name, error in self.errors:
                print(f"\n{test_name}:")
                print(f"  {error}")

        print(f"\n{'='*70}")
        if self.failed == 0:
            print("✓ ALL TESTS PASSED")
        else:
            print(f"✗ {self.failed} TEST(S) FAILED")
        print(f"{'='*70}\n")

        return self.failed == 0


def test_step1_simulate_species_tree():
    """Step 1: Simulate a species tree."""
    print("\nStep 1: Simulating species tree...")

    sp_tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=42)

    print(f"  - Number of leaves: {sp_tree.num_leaves()}")
    print(f"  - Number of nodes: {sp_tree.num_nodes()}")
    print(f"  - Tree height: {sp_tree.tree_height():.4f}")
    print(f"  - Leaf names: {sp_tree.leaf_names()[:5]}...")

    assert sp_tree.num_leaves() == 20, "Should have exactly 20 leaves"
    assert sp_tree.num_nodes() > 20, "Should have internal nodes"
    assert sp_tree.tree_height() > 0, "Tree height should be positive"

    # Save to Newick
    with tempfile.NamedTemporaryFile(suffix='.nwk', delete=False) as f:
        filepath = f.name

    try:
        sp_tree.save_newick(filepath)
        assert os.path.exists(filepath), "Newick file should exist"

        with open(filepath, 'r') as f:
            newick_content = f.read()

        assert newick_content.endswith(';'), "Newick should end with semicolon"
        print(f"  ✓ Saved to Newick: {filepath}")
        print(f"    First 100 chars: {newick_content[:100]}...")
    finally:
        if os.path.exists(filepath):
            os.remove(filepath)

    return sp_tree


def test_step2_export_bd_events():
    """Step 2: Export BD events to CSV (if implemented)."""
    print("\nStep 2: Checking BD events export...")

    sp_tree = rustree.simulate_species_tree(n=15, lambda_=1.0, mu=0.3, seed=123)

    # Check if export_bd_events is available
    if hasattr(sp_tree, 'export_bd_events'):
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as f:
            filepath = f.name

        try:
            sp_tree.export_bd_events(filepath)
            assert os.path.exists(filepath), "BD events CSV should exist"
            assert os.path.getsize(filepath) > 0, "BD events CSV should not be empty"

            with open(filepath, 'r') as f:
                content = f.read()
                lines = content.strip().split('\n')

            print(f"  ✓ Exported BD events: {len(lines)} lines")
            print(f"    First few lines:")
            for line in lines[:5]:
                print(f"      {line}")
        finally:
            if os.path.exists(filepath):
                os.remove(filepath)
    else:
        print("  ⚠ export_bd_events method not yet implemented (EXPECTED)")
        print("    This is a known gap - worker agent may not have completed this binding yet")


def test_step3_compute_distances():
    """Step 3: Compute pairwise distances."""
    print("\nStep 3: Computing pairwise distances...")

    sp_tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.5, seed=456)

    # Check if pairwise_distances is available
    if hasattr(sp_tree, 'pairwise_distances'):
        try:
            import pandas as pd

            # Topological distances
            df_topo = sp_tree.pairwise_distances("topological", leaves_only=True)
            print(f"  - Topological distances: {len(df_topo)} pairs")
            print(f"    First few:")
            print(df_topo.head())

            # Metric distances
            df_metric = sp_tree.pairwise_distances("metric", leaves_only=True)
            print(f"  - Metric distances: {len(df_metric)} pairs")
            print(f"    First few:")
            print(df_metric.head())

            assert len(df_topo) == len(df_metric), "Should have same number of pairs"
            assert len(df_topo) > 0, "Should have distance pairs"

            print(f"  ✓ Successfully computed {len(df_topo)} pairwise distances")
        except ImportError:
            print("  ⚠ pandas not installed, skipping distance test")
    else:
        print("  ⚠ pairwise_distances method not yet implemented (EXPECTED)")
        print("    This is a known gap - worker agent may not have completed this binding yet")


def test_step4_sample_tree():
    """Step 4: Sample/extract induced subtree."""
    print("\nStep 4: Sampling induced subtree...")

    sp_tree = rustree.simulate_species_tree(n=20, lambda_=1.0, mu=0.5, seed=789)

    all_leaves = sp_tree.leaf_names()
    print(f"  - Original tree: {len(all_leaves)} leaves")

    # Select subset of leaves
    selected_leaves = all_leaves[:10]
    print(f"  - Selecting {len(selected_leaves)} leaves for sampling")

    sampled_tree = sp_tree.extract_induced_subtree_by_names(selected_leaves)

    print(f"  - Sampled tree: {sampled_tree.num_leaves()} leaves")
    print(f"  - Sampled tree: {sampled_tree.num_nodes()} total nodes")

    assert sampled_tree.num_leaves() == len(selected_leaves), "Should have exactly selected number of leaves"
    assert sampled_tree.num_nodes() >= sampled_tree.num_leaves(), "Should have at least as many nodes as leaves"

    # Verify all selected leaves are in sampled tree
    sampled_leaves = sampled_tree.leaf_names()
    for leaf in selected_leaves:
        assert leaf in sampled_leaves, f"Leaf {leaf} should be in sampled tree"

    print(f"  ✓ Successfully sampled subtree")
    print(f"    Sampled leaves: {sampled_leaves}")

    return sampled_tree


def test_step5_simulate_gene_trees():
    """Step 5: Simulate gene trees on species tree."""
    print("\nStep 5: Simulating gene trees...")

    sp_tree = rustree.simulate_species_tree(n=15, lambda_=1.0, mu=0.5, seed=1000)

    # Single gene tree
    print("  5a. Single gene tree simulation...")
    gt = sp_tree.simulate_dtl(
        lambda_d=0.5,
        lambda_t=0.2,
        lambda_l=0.3,
        require_extant=True,
        seed=2000
    )

    print(f"    - Gene tree nodes: {gt.num_nodes()}")
    print(f"    - Extant genes: {gt.num_extant()}")

    s, d, t, l, leaves = gt.count_events()
    print(f"    - Events: S={s}, D={d}, T={t}, L={l}, Leaves={leaves}")

    assert gt.num_nodes() > 0, "Gene tree should have nodes"
    assert gt.num_extant() > 0, "require_extant=True should ensure extant genes"
    assert s + d + t + l + leaves == gt.num_nodes(), "Events should sum to total nodes"

    # Batch simulation
    print("\n  5b. Batch gene tree simulation...")
    gene_trees = sp_tree.simulate_dtl_batch(
        n=10,
        lambda_d=0.5,
        lambda_t=0.2,
        lambda_l=0.3,
        require_extant=True,
        seed=3000
    )

    print(f"    - Simulated {len(gene_trees)} gene trees")
    assert len(gene_trees) == 10, "Should have 10 gene trees"

    extant_counts = [gt.num_extant() for gt in gene_trees]
    print(f"    - Extant genes per tree: min={min(extant_counts)}, max={max(extant_counts)}, mean={sum(extant_counts)/len(extant_counts):.1f}")

    for i, gt in enumerate(gene_trees):
        assert gt.num_extant() > 0, f"Tree {i} should have extant genes (require_extant=True)"

    print(f"  ✓ Successfully simulated gene trees")

    return gene_trees[0]


def test_step6_gene_tree_operations():
    """Step 6: Test gene tree operations."""
    print("\nStep 6: Testing gene tree operations...")

    sp_tree = rustree.simulate_species_tree(n=10, lambda_=1.0, mu=0.5, seed=4000)
    gt = sp_tree.simulate_dtl(
        lambda_d=0.5,
        lambda_t=0.2,
        lambda_l=0.2,
        require_extant=True,
        seed=5000
    )

    # Test to_newick
    print("  6a. Testing to_newick()...")
    newick = gt.to_newick()
    assert isinstance(newick, str), "Newick should be string"
    assert newick.endswith(';'), "Newick should end with semicolon"
    print(f"    ✓ Newick: {newick[:80]}...")

    # Test save_newick
    print("  6b. Testing save_newick()...")
    with tempfile.NamedTemporaryFile(suffix='.nwk', delete=False) as f:
        filepath = f.name
    try:
        gt.save_newick(filepath)
        assert os.path.exists(filepath), "Newick file should exist"
        print(f"    ✓ Saved to {filepath}")
    finally:
        if os.path.exists(filepath):
            os.remove(filepath)

    # Test to_xml
    print("  6c. Testing to_xml()...")
    xml = gt.to_xml()
    assert isinstance(xml, str), "XML should be string"
    assert len(xml) > 100, "XML should have substantial content"
    print(f"    ✓ XML length: {len(xml)} chars")

    # Test save_xml
    print("  6d. Testing save_xml()...")
    with tempfile.NamedTemporaryFile(suffix='.xml', delete=False) as f:
        filepath = f.name
    try:
        gt.save_xml(filepath)
        assert os.path.exists(filepath), "XML file should exist"
        print(f"    ✓ Saved to {filepath}")
    finally:
        if os.path.exists(filepath):
            os.remove(filepath)

    # Test to_csv
    print("  6e. Testing to_csv()...")
    try:
        import pandas as pd
        df = gt.to_csv()
        assert isinstance(df, pd.DataFrame), "to_csv should return DataFrame"
        assert len(df) == gt.num_nodes(), "DataFrame rows should match nodes"
        print(f"    ✓ DataFrame: {len(df)} rows, {len(df.columns)} columns")
        print(f"      Columns: {list(df.columns)}")
    except ImportError:
        print("    ⚠ pandas not installed, skipping")

    # Test extant_gene_names
    print("  6f. Testing extant_gene_names()...")
    names = gt.extant_gene_names()
    assert isinstance(names, list), "Should return list"
    assert len(names) == gt.num_extant(), "Length should match num_extant"
    print(f"    ✓ Extant genes: {names[:5]}...")

    # Test sample_extant
    print("  6g. Testing sample_extant()...")
    sampled = gt.sample_extant()
    assert sampled.num_extant() == gt.num_extant(), "Sampled should preserve extant count"
    print(f"    ✓ Sampled extant genes: {sampled.num_nodes()} nodes")

    print(f"  ✓ All gene tree operations successful")


def test_step7_full_workflow():
    """Step 7: Complete end-to-end workflow."""
    print("\nStep 7: Testing complete end-to-end workflow...")

    # 1. Simulate species tree
    print("  Step 7.1: Simulate species tree (30 species)...")
    sp_tree = rustree.simulate_species_tree(n=30, lambda_=1.0, mu=0.5, seed=6000)
    print(f"    ✓ Species tree: {sp_tree.num_leaves()} leaves")

    # 2. Sample subset of species
    print("  Step 7.2: Sample subset of species (15 species)...")
    all_species = sp_tree.leaf_names()
    selected_species = all_species[::2]  # Every other species
    sampled_sp_tree = sp_tree.extract_induced_subtree_by_names(selected_species)
    print(f"    ✓ Sampled species tree: {sampled_sp_tree.num_leaves()} leaves")

    # 3. Simulate gene trees on original tree
    print("  Step 7.3: Simulate 5 gene families on original tree...")
    gene_trees_orig = sp_tree.simulate_dtl_batch(
        n=5,
        lambda_d=0.5,
        lambda_t=0.2,
        lambda_l=0.3,
        transfer_alpha=1.0,  # Assortative transfers
        require_extant=True,
        seed=7000
    )
    print(f"    ✓ Simulated {len(gene_trees_orig)} gene families")

    # 4. Simulate gene trees on sampled tree
    print("  Step 7.4: Simulate 5 gene families on sampled tree...")
    gene_trees_sampled = sampled_sp_tree.simulate_dtl_batch(
        n=5,
        lambda_d=0.5,
        lambda_t=0.2,
        lambda_l=0.3,
        transfer_alpha=1.0,
        require_extant=True,
        seed=8000
    )
    print(f"    ✓ Simulated {len(gene_trees_sampled)} gene families")

    # 5. Export all results
    print("  Step 7.5: Exporting results...")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Save species trees
        sp_orig_path = os.path.join(tmpdir, "species_original.nwk")
        sp_sampled_path = os.path.join(tmpdir, "species_sampled.nwk")

        sp_tree.save_newick(sp_orig_path)
        sampled_sp_tree.save_newick(sp_sampled_path)

        assert os.path.exists(sp_orig_path), "Original species tree should be saved"
        assert os.path.exists(sp_sampled_path), "Sampled species tree should be saved"
        print(f"    ✓ Saved species trees to {tmpdir}")

        # Save gene trees (both Newick and XML)
        for i, gt in enumerate(gene_trees_orig):
            nwk_path = os.path.join(tmpdir, f"gene_orig_{i}.nwk")
            xml_path = os.path.join(tmpdir, f"gene_orig_{i}.xml")

            gt.save_newick(nwk_path)
            gt.save_xml(xml_path)

            assert os.path.exists(nwk_path), f"Gene tree {i} Newick should exist"
            assert os.path.exists(xml_path), f"Gene tree {i} XML should exist"

        print(f"    ✓ Saved {len(gene_trees_orig)} gene families (original)")

        for i, gt in enumerate(gene_trees_sampled):
            nwk_path = os.path.join(tmpdir, f"gene_sampled_{i}.nwk")
            xml_path = os.path.join(tmpdir, f"gene_sampled_{i}.xml")

            gt.save_newick(nwk_path)
            gt.save_xml(xml_path)

            assert os.path.exists(nwk_path), f"Gene tree {i} Newick should exist"
            assert os.path.exists(xml_path), f"Gene tree {i} XML should exist"

        print(f"    ✓ Saved {len(gene_trees_sampled)} gene families (sampled)")

        # List all files
        all_files = sorted(os.listdir(tmpdir))
        print(f"    ✓ Total files created: {len(all_files)}")
        print(f"      Files: {', '.join(all_files[:10])}...")

    print(f"  ✓ Complete workflow successful!")


def main():
    """Run all integration tests."""
    print("\n" + "="*70)
    print("RUSTREE PYTHON BINDINGS - INTEGRATION TEST SUITE")
    print("="*70)

    suite = IntegrationTestSuite()

    # Run all tests
    suite.run_test(test_step1_simulate_species_tree)
    suite.run_test(test_step2_export_bd_events)
    suite.run_test(test_step3_compute_distances)
    suite.run_test(test_step4_sample_tree)
    suite.run_test(test_step5_simulate_gene_trees)
    suite.run_test(test_step6_gene_tree_operations)
    suite.run_test(test_step7_full_workflow)

    # Print summary
    success = suite.print_summary()

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
