use std::process::Command;
use std::fs;
use std::path::{Path, PathBuf};
use pest::Parser;
use rustree::newick::newick::{NewickParser, newick_to_tree, Rule};
use rustree::node::compare_nodes;

#[test]
fn test_extract_extant() {
    // Paths: using an existing species tree file and a temporary output directory.
    let species_tree = "./tests/test_tree_2.nwk";
    let output_dir = "./tests/temp_output";
    
    // Clean and create the output directory.
    let _ = fs::remove_dir_all(output_dir);
    fs::create_dir_all(output_dir).expect("Failed to create temporary output directory");
    
    // Set the number of extant nodes.
    let n_extant = "2";
    
    // Invoke the extract_extant binary with verbose flag.
    let target_executable = "./target/debug/extract_extant";
    let output = Command::new(target_executable)
                        .arg(species_tree)
                        .arg(n_extant)
                        .arg(output_dir)
                        .arg("--verbose")
                        .output()
                        .expect("Failed to execute extract_extant binary");
                        
    assert!(output.status.success(), "extract_extant did not run successfully");
    
    // Verify the expected output file exists and is non-empty.
    let output_file = Path::new(output_dir).join("extant_species_tree.nwk");
    assert!(output_file.exists(), "Output file not found");
    let contents = fs::read_to_string(&output_file).expect("Failed to read output file");
    assert!(!contents.is_empty(), "Output file is empty");

    // Convert the computed Newick string to a Node tree.
    let computed_tree = {
        let mut pairs = NewickParser::parse(Rule::newick, contents.trim())
                          .expect("Failed to parse computed Newick string");
        let mut nodes = newick_to_tree(pairs.next().unwrap());
        nodes.pop().expect("No computed tree found")
    };

    // Read and convert the expected tree from test_tree_2_F.nwk.
    let expected_path = PathBuf::from("./tests/test_tree_2_F.nwk");
    let expected_str = fs::read_to_string(&expected_path)
                         .expect("Failed to read expected Newick file");
    let expected_tree = {
        let mut pairs = NewickParser::parse(Rule::newick, expected_str.trim())
                          .expect("Failed to parse expected Newick string");
        let mut nodes = newick_to_tree(pairs.next().unwrap());
        nodes.pop().expect("No expected tree found")
    };

    // Compare the computed extant tree with the expected tree.
    if !compare_nodes(&computed_tree, &expected_tree) {
        panic!("Computed extant tree does not match expected tree.");
    } else {
        println!("Computed extant tree matches expected tree.");
        println!("Computed tree: {:?}", computed_tree);
        println!("Expected tree: {:?}", expected_tree);
    }

    // Cleanup (optional): remove temporary files if needed.
    let _ = fs::remove_dir_all(output_dir);
}
