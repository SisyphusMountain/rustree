# RecPhyloXML Parsing in rustree

## Overview

rustree now supports parsing RecPhyloXML files (e.g., from ALERax) in addition to exporting them. This allows you to read reconciled gene trees with species tree mappings and event information.

## Features

- ✅ Parse RecPhyloXML files containing both species tree and reconciled gene tree
- ✅ Handle all event types: Speciation, Duplication, Transfer, Loss, Leaf
- ✅ Support for missing branch lengths (defaults to 0.0)
- ✅ Handle transfer events with `<branchingOut>` and `<transferBack>` tags
- ✅ Round-trip support: parse → export → parse
- ✅ Python and R bindings

## Quick Start

### Rust

#### Format 1: Complete RecPhyloXML (both trees in one file)

```rust
use rustree::RecTreeOwned;

// Parse from file (contains both species and gene trees)
let rec_tree = RecTreeOwned::from_xml_file("alerax_output.xml")?;

// Access the trees
println!("Species nodes: {}", rec_tree.species_tree.nodes.len());
println!("Gene nodes: {}", rec_tree.gene_tree.nodes.len());

// Query event information
let event = rec_tree.event_for(gene_node_idx);
let species_idx = rec_tree.species_node_for(gene_node_idx);

// Export back to XML
let xml = rec_tree.to_xml();
```

#### Format 2: Separate Files (Newick + XML)

```rust
use rustree::RecTreeOwned;

// Parse from separate files (species tree in Newick, gene tree in XML)
let rec_tree = RecTreeOwned::from_separate_files(
    "species_tree.nwk",
    "gene_tree_rec.xml"
)?;

// Use same API as Format 1
println!("Species nodes: {}", rec_tree.species_tree.nodes.len());
println!("Gene nodes: {}", rec_tree.gene_tree.nodes.len());
```

### Python

```python
import rustree

# Parse RecPhyloXML file
gt = rustree.parse_recphyloxml("alerax_output.xml")

# Access tree information
print(f"Gene tree nodes: {gt.num_nodes()}")
print(f"Extant genes: {gt.num_extant()}")

# Count events
spec, dup, trans, loss, leaf = gt.count_events()
print(f"Duplications: {dup}, Transfers: {trans}")

# Export to XML
xml_str = gt.to_xml()
gt.save_xml("output.xml")

# Export to CSV
gt.to_csv("events.csv")
```

### R

```r
library(rustree)

# Parse RecPhyloXML file
gt <- parse_recphyloxml_r("alerax_output.xml")

# Access tree information
print(paste("Gene tree nodes:", length(gt$gene_tree$nodes)))
print(paste("Species tree nodes:", length(gt$species_tree$nodes)))

# Export to XML
xml_str <- gene_tree_to_xml_r(gt)
save_xml_r(gt, "output.xml")

# Export to CSV
save_csv_r(gt, "events.csv")
```

## Data Structure

### RecTreeOwned

The `RecTreeOwned` structure owns both the species tree and gene tree:

```rust
pub struct RecTreeOwned {
    pub species_tree: FlatTree,
    pub gene_tree: FlatTree,
    pub node_mapping: Vec<usize>,   // Maps gene node idx → species node idx
    pub event_mapping: Vec<Event>,  // Maps gene node idx → Event type
}
```

### Event Types

```rust
pub enum Event {
    Speciation,   // Gene tree follows species tree split
    Duplication,  // Gene duplicated within a species
    Transfer,     // Horizontal gene transfer (donor side)
    Loss,         // Gene loss event
    Leaf,         // Extant gene (leaf node)
}
```

## Example: Parse ALERax Output

```rust
use rustree::{RecTreeOwned, Event};

fn main() -> Result<(), String> {
    // Parse ALERax RecPhyloXML file
    let rec_tree = RecTreeOwned::from_xml_file("alerax_g.xml")?;

    // Count event types
    let duplication_count = rec_tree.event_mapping.iter()
        .filter(|e| **e == Event::Duplication)
        .count();

    println!("Duplications: {}", duplication_count);

    // Find all extant genes
    for (idx, node) in rec_tree.gene_tree.nodes.iter().enumerate() {
        if rec_tree.event_mapping[idx] == Event::Leaf {
            let species_idx = rec_tree.node_mapping[idx];
            let species_name = &rec_tree.species_tree.nodes[species_idx].name;
            println!("Gene {} is in species {}", node.name, species_name);
        }
    }

    Ok(())
}
```

## Running the Example

```bash
# Run the example with an ALERax file
cargo run --example parse_recphyloxml path/to/alerax_g.xml

# Save the exported XML
cargo run --example parse_recphyloxml path/to/alerax_g.xml --save
```

## Format Support

### Supported RecPhyloXML Features

- ✅ **Format 1**: Complete RecPhyloXML with both `<spTree>` and `<recGeneTree>` sections
- ✅ **Format 2**: Gene-tree-only XML (no `<spTree>`) with separate Newick species tree
- ✅ All event types in `<eventsRec>`
- ✅ Branch lengths (optional, defaults to 0.0)
- ✅ Transfer events with `<transferBack>`
- ✅ Node names (including ambiguous names)
- ✅ Alternative event tags:
  - `<leaf>` - Standard leaf/extant gene marker
  - `<P>` - Alternative present-day gene marker
  - `<T>` - Generic internal node marker (treated as speciation)
- ✅ Extra attributes ignored (e.g., `ts` timestamp)

### Format Notes

1. **Missing branch lengths**: If `<branchLength>` tags are missing, values default to 0.0
2. **Ambiguous node names**: Species lookup uses first occurrence in tree traversal
3. **Binary trees only**: Only the first 2 children are processed (warning issued if more)
4. **Transfer events**: `<branchingOut>` marks donor, `<transferBack>` marks recipient

## Testing

Run the test suite:

```bash
# Run all RecPhyloXML tests
cargo test --test recphyloxml_tests

# Test with real ALERax file
cargo test --test test_alerax_real_file -- --nocapture

# Run all tests
cargo test
```

## Error Handling

Parsing errors are returned as `Result<RecTreeOwned, String>`:

```rust
match RecTreeOwned::from_xml_file("file.xml") {
    Ok(rec_tree) => {
        // Success - use rec_tree
    }
    Err(e) => {
        eprintln!("Failed to parse: {}", e);
        // Handle error
    }
}
```

Common errors:
- Missing `<spTree>` or `<recGeneTree>` sections
- Invalid XML format
- Species referenced in gene tree not found in species tree
- I/O errors (file not found, permission denied, etc.)

## Implementation Details

### Parser

- Uses `quick-xml` for efficient XML parsing
- Event-based (streaming) parser for memory efficiency
- Supports large files

### Conversion Process

1. Parse species tree section → build FlatTree + name-to-index HashMap
2. Parse gene tree section → build FlatTree with event information
3. For each gene node, lookup species index using HashMap
4. Create mappings: node_mapping (gene→species) and event_mapping (gene→event)

### Compatibility

The XML parser is designed to work with:
- ALERax output format
- rustree's own XML export format
- Standard RecPhyloXML format

## API Reference

### Main Functions

#### Format 1: Complete RecPhyloXML

```rust
// Parse from string (with both trees)
RecTreeOwned::from_xml(xml_content: &str) -> Result<Self, String>

// Parse from file (with both trees)
RecTreeOwned::from_xml_file(filepath: &str) -> Result<Self, String>
```

#### Format 2: Separate Files

```rust
// Parse from separate files (Newick + XML)
RecTreeOwned::from_separate_files(
    species_newick_path: &str,
    gene_xml_path: &str
) -> Result<Self, String>

// Parse gene-tree-only XML with pre-loaded species tree
RecTreeOwned::from_gene_tree_xml(
    xml_content: &str,
    species_tree: FlatTree
) -> Result<Self, String>

RecTreeOwned::from_gene_tree_xml_file(
    xml_filepath: &str,
    species_tree: FlatTree
) -> Result<Self, String>
```

#### Common Methods

```rust
// Export to XML
rec_tree.to_xml() -> String

// Convert to borrowed RecTree
rec_tree.as_rectree() -> RecTree<'_>

// Query methods
rec_tree.species_node_for(gene_idx: usize) -> usize
rec_tree.event_for(gene_idx: usize) -> &Event
rec_tree.get_full_info(gene_idx: usize) -> (&FlatNode, usize, &Event)
```

### Python Bindings

```python
rustree.parse_recphyloxml(filepath: str) -> PyGeneTree
```

### R Bindings

```r
parse_recphyloxml_r(filepath: str) -> List
```

## See Also

- [Python Tutorial](PYTHON_TUTORIAL.md)
- [R Tutorial](R_TUTORIAL.md)
- [RecPhyloXML Format Specification](http://www.recg.org/)
