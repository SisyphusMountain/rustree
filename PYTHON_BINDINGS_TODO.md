# Python Bindings: Implementation Status

All planned Python binding features have been implemented in `src/python.rs`.

## Implemented Features

### PySpeciesTree

- `to_newick()` - Convert to Newick format
- `num_nodes()` - Get node count
- `num_leaves()` - Get extant species count
- `tree_height()` - Get total tree height
- `root_index()` - Get root node index
- `leaf_names()` - Get leaf names
- `extract_induced_subtree_by_names(names)` - Extract induced subtree by leaf names
- `save_newick(filepath)` - Save Newick to file
- `save_bd_events_csv(filepath)` - Save birth-death events to CSV
- `get_bd_events()` - Get BD events as a dictionary (pandas-compatible)
- `get_ltt_data()` - Get Lineages Through Time data
- `plot_ltt(...)` - Plot LTT using matplotlib
- `simulate_dtl(...)` - Simulate gene tree (per-gene-copy DTL model)
- `simulate_dtl_batch(...)` - Batch simulate gene trees (per-gene-copy DTL model)
- `simulate_dtl_per_species(...)` - Simulate gene tree (Zombi-style per-species DTL model)
- `simulate_dtl_per_species_batch(...)` - Batch simulate gene trees (per-species DTL model)
- `pairwise_distances(distance_type, leaves_only)` - Compute pairwise distances as pandas DataFrame
- `save_pairwise_distances_csv(filepath, distance_type, leaves_only)` - Save pairwise distances to CSV

### PyGeneTree

- `to_newick()` - Convert to Newick format
- `save_newick(filepath)` - Save Newick to file
- `num_nodes()` - Get node count
- `num_extant()` - Get extant gene count
- `count_events()` - Count events by type (speciations, duplications, transfers, losses, leaves)
- `extant_gene_names()` - Get names of extant genes
- `sample_extant()` - Sample gene tree keeping only extant genes
- `sample_by_names(names)` - Sample gene tree by gene names
- `sample_species_leaves(species_leaf_names)` - Sample by species leaves with LCA-based reconciliation
- `to_xml()` - Export to RecPhyloXML string
- `save_xml(filepath)` - Save RecPhyloXML to file
- `to_svg(filepath, open_browser)` - Generate SVG visualization via thirdkind
- `display()` - Display in Jupyter notebook
- `to_csv(filepath)` - Export gene tree data as pandas DataFrame / CSV

### Module-Level Functions

- `simulate_species_tree(n, lambda_, mu, seed)` - Simulate birth-death species tree
- `parse_species_tree(newick_str)` - Parse Newick string into species tree
- `parse_recphyloxml(filepath)` - Parse RecPhyloXML file into reconciled gene tree
