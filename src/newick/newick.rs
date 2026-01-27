// rustree/src/newick/newick.rs
// This module contains functions to parse Newick formatted strings into Tree structures
// and to convert Tree structures back into Newick formatted strings.
use crate::node::{Node, FlatTree};

#[derive(Parser)]
#[grammar = "newick/newick.pest"]
pub struct NewickParser; // The parser struct, needs to be public to be used outside
                         // The following two functions are used to convert a Newick string into an arborescent Tree structure.

/// Parse a Newick formatted string into a vector of Node trees.
///
/// # Arguments
/// * `newick_str` - A string slice containing the Newick formatted tree
///
/// # Returns
/// A Result containing a vector of Node trees on success, or an error message on failure.
pub fn parse_newick(newick_str: &str) -> Result<Vec<Node>, String> {
    use pest::Parser;
    let mut pairs = NewickParser::parse(Rule::newick, newick_str)
        .map_err(|e| e.to_string())?;

    Ok(newick_to_tree(pairs.next().ok_or("No tree found in input")?))
}

fn newick_to_tree(pair: pest::iterators::Pair<Rule>) -> Vec<Node> {
    let mut vec_trees: Vec<Node> = Vec::new();
    for inner in pair.into_inner() {
        let tree = handle_pair(inner);
        vec_trees.push(tree.unwrap());
    }
    vec_trees
}

pub fn handle_pair(pair: pest::iterators::Pair<Rule>) -> Option<Node> {
    match pair.as_rule() {
        Rule::newick => {
            for inner in pair.into_inner() {
                if inner.as_rule() == Rule::subtree {
                    return handle_pair(inner);
                }
            }
            None
        }
        Rule::subtree => handle_pair(pair.into_inner().next().unwrap()),
        Rule::leaf => {
            let mut name = String::new();
            let mut length = 0.0;

            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    }
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    }
                    _ => {} // Ignore other rules
                }
            }
            Some(Node {
                name,
                left_child: None,
                right_child: None,
                depth: None,
                length,
            })
        }
        Rule::internal => {
            let mut name = String::new();
            let mut length = 0.0;
            let mut first_subtree = None;
            let mut second_subtree = None;

            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::subtree => {
                        let subtree = handle_pair(inner_pair).unwrap();
                        if first_subtree.is_none() {
                            first_subtree = Some(subtree);
                        } else {
                            second_subtree = Some(subtree);
                        }
                    }
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    }
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    }
                    _ => {}
                }
            }

            Some(Node {
                name,
                left_child: first_subtree.map(Box::new),
                right_child: second_subtree.map(Box::new),
                depth: None,
                length,
            })
        }
        Rule::NAME | Rule::LENGTH => None,
    }
}
impl Node {
    pub fn to_newick(&self) -> String {
        node_to_newick_internal(self)
    }
}

impl FlatTree {
    pub fn to_newick(&self) -> String {
        let node_tree = self.to_node();
        node_to_newick_internal(&node_tree)
    }
}

fn node_to_newick_internal(node: &Node) -> String {
    /* Takes a node and returns the corresponding subtree in Newick format.
        --------------------------------
        INPUT:
            - node: the node to convert to Newick format.
        OUTPUT:
            - the Newick representation of the subtree rooted at node.
        Warning: rounds the lengths to 6 decimal places.
    */
    if let (Some(left_child), Some(right_child)) = (&node.left_child, &node.right_child) {
        // This is an internal node with both left and right children.
        format!(
            "({},{}){}:{:.6}",
            node_to_newick_internal(left_child),
            node_to_newick_internal(right_child),
            node.name,
            node.length
        )
    } else {
        // This is a leaf node.
        format!("{}:{:.6}", node.name, node.length)
    }
}

