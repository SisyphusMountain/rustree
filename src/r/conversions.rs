#![allow(clippy::type_complexity)]
//! R ↔ Rust conversion helpers for tree and event data structures.

use extendr_api::prelude::*;
use std::str::FromStr;

use crate::bd::{BDEvent, TreeEvent};
use crate::dtl::DTLEvent;
use crate::node::{Event, FlatNode, FlatTree};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ApeOrder {
    Cladewise,
    Postorder,
}

impl ApeOrder {
    fn parse(order: &str) -> Result<Self> {
        match order {
            "cladewise" => Ok(Self::Cladewise),
            "postorder" | "pruningwise" => Ok(Self::Postorder),
            other => Err(Error::Other(format!(
                "Unsupported ape edge order '{}'. Expected 'cladewise' or 'postorder'.",
                other
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Cladewise => "cladewise",
            Self::Postorder => "postorder",
        }
    }
}

fn preorder_indices_checked(tree: &FlatTree) -> Result<Vec<usize>> {
    let n = tree.nodes.len();
    if n == 0 {
        return Err(Error::Other(
            "tree must contain at least one node".to_string(),
        ));
    }
    if tree.root >= n {
        return Err(Error::Other(format!(
            "tree root index {} is outside the valid node range 0..{}",
            tree.root,
            n.saturating_sub(1)
        )));
    }

    let mut stack = Vec::with_capacity(n);
    let mut visited = vec![false; n];
    let mut preorder = Vec::with_capacity(n);
    stack.push(tree.root);

    while let Some(idx) = stack.pop() {
        if idx >= n {
            return Err(Error::Other(format!(
                "tree contains child index {} outside the valid node range 0..{}",
                idx,
                n.saturating_sub(1)
            )));
        }
        if visited[idx] {
            return Err(Error::Other(
                "tree contains a cycle or a node with multiple parents".to_string(),
            ));
        }

        visited[idx] = true;
        preorder.push(idx);

        let node = &tree.nodes[idx];
        if let Some(right) = node.right_child {
            stack.push(right);
        }
        if let Some(left) = node.left_child {
            stack.push(left);
        }
    }

    if preorder.len() != n {
        return Err(Error::Other(
            "tree contains nodes that are not reachable from the root".to_string(),
        ));
    }

    Ok(preorder)
}

pub(crate) fn flattree_to_ape_phylo(
    tree: &FlatTree,
    order: &str,
    use_node_labels: bool,
    include_root_edge: bool,
) -> Result<Robj> {
    let order = ApeOrder::parse(order)?;
    let n = tree.nodes.len();
    if n > i32::MAX as usize {
        return Err(Error::Other(format!(
            "ape phylo objects use 32-bit R integer node IDs; tree has {} nodes",
            n
        )));
    }

    let preorder = preorder_indices_checked(tree)?;
    let mut tip_nodes = Vec::new();
    let mut internal_nodes = Vec::new();

    for &idx in &preorder {
        let node = &tree.nodes[idx];
        match (node.left_child, node.right_child) {
            (None, None) => tip_nodes.push(idx),
            (Some(_), Some(_)) => internal_nodes.push(idx),
            _ => {
                return Err(Error::Other(format!(
                    "ape conversion expects binary or leaf nodes; node {} ('{}') has exactly one child",
                    idx, node.name
                )));
            }
        }
    }

    if tip_nodes.is_empty() {
        return Err(Error::Other(
            "tree must contain at least one tip".to_string(),
        ));
    }

    let n_tip = tip_nodes.len();
    let n_node = internal_nodes.len();
    let mut rustree_to_ape = vec![0_i32; n];

    for (pos, &idx) in tip_nodes.iter().enumerate() {
        rustree_to_ape[idx] = (pos + 1) as i32;
    }
    for (pos, &idx) in internal_nodes.iter().enumerate() {
        rustree_to_ape[idx] = (n_tip + pos + 1) as i32;
    }

    let mut derived_parent = vec![None; n];
    for &idx in &preorder {
        let node = &tree.nodes[idx];
        for child in [node.left_child, node.right_child].into_iter().flatten() {
            if child >= n {
                return Err(Error::Other(format!(
                    "tree contains child index {} outside the valid node range 0..{}",
                    child,
                    n.saturating_sub(1)
                )));
            }
            if let Some(existing_parent) = derived_parent[child] {
                return Err(Error::Other(format!(
                    "node {} has multiple parents: {} and {}",
                    child, existing_parent, idx
                )));
            }
            derived_parent[child] = Some(idx);
        }
    }

    for (idx, node) in tree.nodes.iter().enumerate() {
        let expected_parent = if idx == tree.root {
            None
        } else {
            derived_parent[idx]
        };
        if node.parent != expected_parent {
            return Err(Error::Other(format!(
                "parent and child links are inconsistent at node {} ('{}'): parent={:?}, child-derived parent={:?}",
                idx, node.name, node.parent, expected_parent
            )));
        }
    }

    let mut edge_order = preorder.clone();
    if order == ApeOrder::Postorder {
        edge_order.reverse();
    }
    edge_order.retain(|&idx| idx != tree.root);
    let n_edge = edge_order.len();
    if n_edge > i32::MAX as usize {
        return Err(Error::Other(format!(
            "ape phylo edge matrices use 32-bit R dimensions; tree has {} edges",
            n_edge
        )));
    }

    for &idx in &edge_order {
        if derived_parent[idx].is_none() {
            return Err(Error::Other(format!("non-root node {} has no parent", idx)));
        }
    }

    let edge = Robj::from(RMatrix::<i32>::new_matrix(n_edge, 2, |row, col| {
        let idx = edge_order[row];
        match col {
            0 => {
                let parent = derived_parent[idx].expect("edge parent checked before allocation");
                rustree_to_ape[parent]
            }
            _ => rustree_to_ape[idx],
        }
    }));

    let edge_lengths: Vec<f64> = edge_order
        .iter()
        .map(|&idx| tree.nodes[idx].length)
        .collect();
    let tip_labels: Vec<String> = tip_nodes
        .iter()
        .map(|&idx| tree.nodes[idx].name.clone())
        .collect();

    let mut names = vec!["edge", "edge.length", "tip.label", "Nnode"];
    let mut values = vec![
        edge,
        Robj::from(edge_lengths),
        Robj::from(tip_labels),
        Robj::from(n_node as i32),
    ];

    if use_node_labels && n_node > 0 {
        let node_labels: Vec<String> = internal_nodes
            .iter()
            .map(|&idx| tree.nodes[idx].name.clone())
            .collect();
        if node_labels.iter().any(|label| !label.is_empty()) {
            names.push("node.label");
            values.push(Robj::from(node_labels));
        }
    }

    let root_length = tree.nodes[tree.root].length;
    if include_root_edge && root_length.is_finite() && root_length != 0.0 {
        names.push("root.edge");
        values.push(Robj::from(root_length));
    }

    let phylo = List::from_names_and_values(names, values)?;
    let phylo = Robj::from(phylo).set_class(["phylo"])?;
    phylo.set_attrib(Symbol::from_string("order"), order.as_str())
}

pub(crate) fn flattrees_to_ape_multiphylo<'a, I>(
    trees: I,
    order: &str,
    use_node_labels: bool,
    include_root_edge: bool,
) -> Result<Robj>
where
    I: IntoIterator<Item = &'a FlatTree>,
{
    let phylo_values: Vec<Robj> = trees
        .into_iter()
        .map(|tree| flattree_to_ape_phylo(tree, order, use_node_labels, include_root_edge))
        .collect::<Result<Vec<_>>>()?;

    Robj::from(List::from_values(phylo_values)).set_class(["multiPhylo"])
}

// ============================================================================
// FlatTree ↔ R list
// ============================================================================

pub(crate) fn flattree_to_rlist(tree: &FlatTree) -> List {
    let names: Vec<String> = tree.nodes.iter().map(|n| n.name.clone()).collect();
    let parents: Vec<Rint> = tree
        .nodes
        .iter()
        .map(|n| n.parent.map(|p| Rint::from(p as i32)).unwrap_or(Rint::na()))
        .collect();
    let left_children: Vec<Rint> = tree
        .nodes
        .iter()
        .map(|n| {
            n.left_child
                .map(|c| Rint::from(c as i32))
                .unwrap_or(Rint::na())
        })
        .collect();
    let right_children: Vec<Rint> = tree
        .nodes
        .iter()
        .map(|n| {
            n.right_child
                .map(|c| Rint::from(c as i32))
                .unwrap_or(Rint::na())
        })
        .collect();
    let lengths: Vec<f64> = tree.nodes.iter().map(|n| n.length).collect();
    let depths: Vec<Rfloat> = tree
        .nodes
        .iter()
        .map(|n| n.depth.map(Rfloat::from).unwrap_or(Rfloat::na()))
        .collect();
    // bd_event: convert Option<BDEvent> to Rstr for proper NA handling in R
    let bd_events: Vec<Rstr> = tree
        .nodes
        .iter()
        .map(|n| match n.bd_event {
            Some(e) => Rstr::from(e.as_str()),
            None => Rstr::na(),
        })
        .collect();

    list!(
        name = names,
        parent = parents,
        left_child = left_children,
        right_child = right_children,
        length = lengths,
        depth = depths,
        root = tree.root as i32,
        bd_event = bd_events
    )
}

fn derive_missing_depth(idx: usize, nodes: &mut [FlatNode], visiting: &mut [bool]) -> Result<f64> {
    if idx >= nodes.len() {
        return Err(Error::Other(format!(
            "node index {} is outside the valid node range 0..{}",
            idx,
            nodes.len().saturating_sub(1)
        )));
    }
    if let Some(depth) = nodes[idx].depth {
        return Ok(depth);
    }
    if visiting[idx] {
        return Err(Error::Other(
            "tree parent links contain a cycle while deriving node depths".to_string(),
        ));
    }

    visiting[idx] = true;
    let parent = nodes[idx].parent;
    let length = nodes[idx].length;
    let depth = match parent {
        Some(parent) => derive_missing_depth(parent, nodes, visiting)? + length,
        None => 0.0,
    };
    nodes[idx].depth = Some(depth);
    visiting[idx] = false;

    Ok(depth)
}

pub(crate) fn rlist_to_flattree(list: &List) -> Result<FlatTree> {
    let names: Vec<String> = list
        .dollar("name")?
        .as_str_vector()
        .ok_or("Failed to get name column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let parents: Vec<i32> = list
        .dollar("parent")?
        .as_integer_vector()
        .ok_or("Failed to get parent column")?;
    let left_children: Vec<i32> = list
        .dollar("left_child")?
        .as_integer_vector()
        .ok_or("Failed to get left_child column")?;
    let right_children: Vec<i32> = list
        .dollar("right_child")?
        .as_integer_vector()
        .ok_or("Failed to get right_child column")?;
    let lengths: Vec<f64> = list
        .dollar("length")?
        .as_real_vector()
        .ok_or("Failed to get length column")?;
    let depths: Option<Vec<f64>> = match list.dollar("depth") {
        Ok(depths) if depths.is_null() => None,
        Ok(depths) => Some(
            depths
                .as_real_vector()
                .ok_or("Failed to get depth column")?,
        ),
        Err(_) => None,
    };
    let root: i32 = list
        .dollar("root")?
        .as_integer()
        .ok_or("Failed to get root")?;

    // Try to get bd_event (optional, for backward compatibility)
    // Use Rstr to properly handle NA values
    let bd_events_robj = list.dollar("bd_event").ok();

    let mut nodes: Vec<FlatNode> = (0..names.len())
        .map(|i| {
            let bd_event = bd_events_robj.as_ref().and_then(|robj| {
                // Get the i-th element as Rstr to check for NA
                if let Some(str_iter) = robj.as_str_iter() {
                    let strs: Vec<&str> = str_iter.collect();
                    if i < strs.len() {
                        let s = strs[i];
                        // Check if it's NA (empty or "NA" string indicates NA from R)
                        if s.is_empty() || s == "NA" {
                            return None;
                        }
                        return BDEvent::from_str(s).ok();
                    }
                }
                None
            });
            FlatNode {
                name: names[i].clone(),
                parent: if parents[i].is_na() {
                    None
                } else {
                    Some(parents[i] as usize)
                },
                left_child: if left_children[i].is_na() {
                    None
                } else {
                    Some(left_children[i] as usize)
                },
                right_child: if right_children[i].is_na() {
                    None
                } else {
                    Some(right_children[i] as usize)
                },
                length: lengths[i],
                depth: depths
                    .as_ref()
                    .and_then(|values| values.get(i))
                    .and_then(|depth| if depth.is_na() { None } else { Some(*depth) }),
                bd_event,
            }
        })
        .collect();

    if nodes.iter().any(|node| node.depth.is_none()) {
        let mut visiting = vec![false; nodes.len()];
        for idx in 0..nodes.len() {
            derive_missing_depth(idx, &mut nodes, &mut visiting)?;
        }
    }

    Ok(FlatTree {
        nodes,
        root: root as usize,
    })
}

// ============================================================================
// BD events ↔ R list
// ============================================================================

pub(crate) fn bd_events_to_rlist(events: &[TreeEvent]) -> List {
    let times: Vec<f64> = events.iter().map(|e| e.time).collect();
    let node_ids: Vec<i32> = events.iter().map(|e| e.node_id as i32).collect();
    let event_types: Vec<String> = events
        .iter()
        .map(|e| e.event_type.as_str().to_string())
        .collect();
    let child1s: Vec<Rint> = events
        .iter()
        .map(|e| e.child1.map(|c| Rint::from(c as i32)).unwrap_or(Rint::na()))
        .collect();
    let child2s: Vec<Rint> = events
        .iter()
        .map(|e| e.child2.map(|c| Rint::from(c as i32)).unwrap_or(Rint::na()))
        .collect();

    list!(
        time = times,
        node_id = node_ids,
        event_type = event_types,
        child1 = child1s,
        child2 = child2s
    )
}

pub(crate) fn rlist_to_bd_events(list: &List) -> Result<Vec<TreeEvent>> {
    let times: Vec<f64> = list
        .dollar("time")?
        .as_real_vector()
        .ok_or("Failed to get time column")?;
    let node_ids: Vec<i32> = list
        .dollar("node_id")?
        .as_integer_vector()
        .ok_or("Failed to get node_id column")?;
    let event_types: Vec<String> = list
        .dollar("event_type")?
        .as_str_vector()
        .ok_or("Failed to get event_type column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let child1s: Vec<i32> = list
        .dollar("child1")?
        .as_integer_vector()
        .ok_or("Failed to get child1 column")?;
    let child2s: Vec<i32> = list
        .dollar("child2")?
        .as_integer_vector()
        .ok_or("Failed to get child2 column")?;

    let events: Vec<TreeEvent> = (0..times.len())
        .map(|i| {
            let event_type =
                BDEvent::from_str(&event_types[i]).map_err(|e| format!("At index {}: {}", i, e))?;
            Ok(TreeEvent {
                time: times[i],
                node_id: node_ids[i] as usize,
                event_type,
                child1: if child1s[i].is_na() {
                    None
                } else {
                    Some(child1s[i] as usize)
                },
                child2: if child2s[i].is_na() {
                    None
                } else {
                    Some(child2s[i] as usize)
                },
            })
        })
        .collect::<std::result::Result<Vec<TreeEvent>, String>>()
        .map_err(extendr_api::Error::Other)?;

    Ok(events)
}

// ============================================================================
// RecTree (gene tree + reconciliation) ↔ R list
// ============================================================================

pub(crate) fn rectree_to_rlist(
    species_tree: &FlatTree,
    gene_tree: &FlatTree,
    node_mapping: &[Option<usize>],
    event_mapping: &[Event],
) -> Result<List> {
    let names: Vec<String> = gene_tree.nodes.iter().map(|n| n.name.clone()).collect();
    let parents: Vec<Rint> = gene_tree
        .nodes
        .iter()
        .map(|n| n.parent.map(|p| Rint::from(p as i32)).unwrap_or(Rint::na()))
        .collect();
    let left_children: Vec<Rint> = gene_tree
        .nodes
        .iter()
        .map(|n| {
            n.left_child
                .map(|c| Rint::from(c as i32))
                .unwrap_or(Rint::na())
        })
        .collect();
    let right_children: Vec<Rint> = gene_tree
        .nodes
        .iter()
        .map(|n| {
            n.right_child
                .map(|c| Rint::from(c as i32))
                .unwrap_or(Rint::na())
        })
        .collect();
    let lengths: Vec<f64> = gene_tree.nodes.iter().map(|n| n.length).collect();
    let depths: Vec<Rfloat> = gene_tree
        .nodes
        .iter()
        .map(|n| n.depth.map(Rfloat::from).unwrap_or(Rfloat::na()))
        .collect();

    // Validate all node_mapping indices before using them
    for (gene_idx, opt_idx) in node_mapping.iter().enumerate() {
        if let Some(idx) = opt_idx {
            if *idx >= species_tree.nodes.len() {
                return Err(Error::Other(format!(
                    "Invalid species node index {} in node_mapping at gene node {}. Species tree has {} nodes.",
                    idx, gene_idx, species_tree.nodes.len()
                )));
            }
        }
    }

    let species_nodes: Vec<String> = node_mapping
        .iter()
        .map(|opt_idx| match opt_idx {
            Some(idx) => species_tree.nodes[*idx].name.clone(),
            None => "NA".to_string(),
        })
        .collect();

    let events: Vec<String> = event_mapping
        .iter()
        .map(|e| match e {
            Event::Speciation => "Speciation".to_string(),
            Event::Duplication => "Duplication".to_string(),
            Event::Transfer => "Transfer".to_string(),
            Event::Loss => "Loss".to_string(),
            Event::Leaf => "Leaf".to_string(),
        })
        .collect();

    // Also store the species tree for later use
    let species_tree_list = flattree_to_rlist(species_tree);

    Ok(list!(
        name = names,
        parent = parents,
        left_child = left_children,
        right_child = right_children,
        length = lengths,
        depth = depths,
        species_node = species_nodes,
        event = events,
        root = gene_tree.root as i32,
        species_tree = species_tree_list
    ))
}

pub(crate) fn rlist_to_genetree(
    list: &List,
) -> Result<(FlatTree, FlatTree, Vec<Option<usize>>, Vec<Event>)> {
    // Get gene tree data
    let names: Vec<String> = list
        .dollar("name")?
        .as_str_vector()
        .ok_or("Failed to get name column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let parents: Vec<i32> = list
        .dollar("parent")?
        .as_integer_vector()
        .ok_or("Failed to get parent column")?;
    let left_children: Vec<i32> = list
        .dollar("left_child")?
        .as_integer_vector()
        .ok_or("Failed to get left_child column")?;
    let right_children: Vec<i32> = list
        .dollar("right_child")?
        .as_integer_vector()
        .ok_or("Failed to get right_child column")?;
    let lengths: Vec<f64> = list
        .dollar("length")?
        .as_real_vector()
        .ok_or("Failed to get length column")?;
    let depths: Vec<f64> = list
        .dollar("depth")?
        .as_real_vector()
        .ok_or("Failed to get depth column")?;
    let species_node_names: Vec<String> = list
        .dollar("species_node")?
        .as_str_vector()
        .ok_or("Failed to get species_node column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let event_strs: Vec<String> = list
        .dollar("event")?
        .as_str_vector()
        .ok_or("Failed to get event column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let root: i32 = list
        .dollar("root")?
        .as_integer()
        .ok_or("Failed to get root")?;

    // Get species tree
    let species_tree_list: List = list
        .dollar("species_tree")?
        .try_into()
        .map_err(|_| "Failed to get species_tree")?;
    let species_tree = rlist_to_flattree(&species_tree_list)?;

    // Build gene tree
    let gene_nodes: Vec<FlatNode> = (0..names.len())
        .map(|i| FlatNode {
            name: names[i].clone(),
            parent: if parents[i].is_na() {
                None
            } else {
                Some(parents[i] as usize)
            },
            left_child: if left_children[i].is_na() {
                None
            } else {
                Some(left_children[i] as usize)
            },
            right_child: if right_children[i].is_na() {
                None
            } else {
                Some(right_children[i] as usize)
            },
            length: lengths[i],
            depth: if depths[i].is_na() {
                None
            } else {
                Some(depths[i])
            },
            bd_event: None,
        })
        .collect();

    let gene_tree = FlatTree {
        nodes: gene_nodes,
        root: root as usize,
    };

    // Build node mapping (find species node index by name)
    let node_mapping: Vec<Option<usize>> = species_node_names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            if name == "NA" || name.is_empty() {
                Ok(None)
            } else {
                species_tree
                    .nodes
                    .iter()
                    .position(|n| &n.name == name)
                    .map(Some)
                    .ok_or_else(|| {
                        format!(
                        "Species node name '{}' (at gene node index {}) not found in species tree",
                        name, i
                    )
                    })
            }
        })
        .collect::<std::result::Result<Vec<Option<usize>>, String>>()
        .map_err(extendr_api::Error::Other)?;

    // Build event mapping
    let event_mapping: Vec<Event> = event_strs.iter()
        .enumerate()
        .map(|(i, s)| match s.as_str() {
            "Speciation" => Ok(Event::Speciation),
            "Duplication" => Ok(Event::Duplication),
            "Transfer" => Ok(Event::Transfer),
            "Loss" => Ok(Event::Loss),
            "Leaf" => Ok(Event::Leaf),
            unknown => Err(format!(
                "Unknown reconciliation event type '{}' at index {}. Expected one of: Speciation, Duplication, Transfer, Loss, Leaf",
                unknown, i
            )),
        })
        .collect::<std::result::Result<Vec<Event>, String>>()
        .map_err(extendr_api::Error::Other)?;

    Ok((gene_tree, species_tree, node_mapping, event_mapping))
}

// ============================================================================
// DTL events ↔ R list
// ============================================================================

/// Serialize a Vec<DTLEvent> to an R list (data-frame-like columns).
pub(crate) fn dtl_events_to_rlist(events: &[DTLEvent], species_tree: &FlatTree) -> List {
    let n = events.len();
    let mut event_types: Vec<String> = Vec::with_capacity(n);
    let mut times: Vec<f64> = Vec::with_capacity(n);
    let mut gene_ids: Vec<i32> = Vec::with_capacity(n);
    let mut species_names: Vec<String> = Vec::with_capacity(n);
    let mut from_species: Vec<Rstr> = Vec::with_capacity(n);
    let mut to_species: Vec<Rstr> = Vec::with_capacity(n);

    let sp_name = |idx: usize| -> String {
        if idx < species_tree.nodes.len() {
            species_tree.nodes[idx].name.clone()
        } else {
            format!("node_{}", idx)
        }
    };

    for event in events {
        match event {
            DTLEvent::Speciation {
                time,
                gene_id,
                species_id,
                ..
            } => {
                event_types.push("Speciation".to_string());
                times.push(*time);
                gene_ids.push(*gene_id as i32);
                species_names.push(sp_name(*species_id));
                from_species.push(Rstr::na());
                to_species.push(Rstr::na());
            }
            DTLEvent::Duplication {
                time,
                gene_id,
                species_id,
                ..
            } => {
                event_types.push("Duplication".to_string());
                times.push(*time);
                gene_ids.push(*gene_id as i32);
                species_names.push(sp_name(*species_id));
                from_species.push(Rstr::na());
                to_species.push(Rstr::na());
            }
            DTLEvent::Transfer {
                time,
                gene_id,
                species_id,
                from_species: from_sp,
                to_species: to_sp,
                ..
            } => {
                event_types.push("Transfer".to_string());
                times.push(*time);
                gene_ids.push(*gene_id as i32);
                species_names.push(sp_name(*species_id));
                from_species.push(Rstr::from(sp_name(*from_sp)));
                to_species.push(Rstr::from(sp_name(*to_sp)));
            }
            DTLEvent::Loss {
                time,
                gene_id,
                species_id,
            } => {
                event_types.push("Loss".to_string());
                times.push(*time);
                gene_ids.push(*gene_id as i32);
                species_names.push(sp_name(*species_id));
                from_species.push(Rstr::na());
                to_species.push(Rstr::na());
            }
            DTLEvent::Leaf {
                time,
                gene_id,
                species_id,
            } => {
                event_types.push("Leaf".to_string());
                times.push(*time);
                gene_ids.push(*gene_id as i32);
                species_names.push(sp_name(*species_id));
                from_species.push(Rstr::na());
                to_species.push(Rstr::na());
            }
        }
    }

    list!(
        event_type = event_types,
        time = times,
        gene_id = gene_ids,
        species = species_names,
        from_species = from_species,
        to_species = to_species
    )
}

/// Deserialize an R DTL events list back to Vec<DTLEvent>.
pub(crate) fn rlist_to_dtl_events(
    events_list: &List,
    species_tree: &FlatTree,
) -> Result<Vec<DTLEvent>> {
    let event_types: Vec<String> = events_list
        .dollar("event_type")?
        .as_str_vector()
        .ok_or("Failed to get event_type column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let times: Vec<f64> = events_list
        .dollar("time")?
        .as_real_vector()
        .ok_or("Failed to get time column")?;
    let gene_ids: Vec<i32> = events_list
        .dollar("gene_id")?
        .as_integer_vector()
        .ok_or("Failed to get gene_id column")?;
    let species_names: Vec<String> = events_list
        .dollar("species")?
        .as_str_vector()
        .ok_or("Failed to get species column")?
        .iter()
        .map(|s| s.to_string())
        .collect();

    // from_species and to_species may contain NAs
    let from_species_robj = events_list.dollar("from_species")?;
    let to_species_robj = events_list.dollar("to_species")?;
    let from_species_strs: Vec<String> = from_species_robj
        .as_str_vector()
        .ok_or("Failed to get from_species column")?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let to_species_strs: Vec<String> = to_species_robj
        .as_str_vector()
        .ok_or("Failed to get to_species column")?
        .iter()
        .map(|s| s.to_string())
        .collect();

    // Build name→index map
    let name_to_idx: std::collections::HashMap<&str, usize> = species_tree
        .nodes
        .iter()
        .enumerate()
        .map(|(i, n)| (n.name.as_str(), i))
        .collect();

    let resolve = |name: &str| -> Result<usize> {
        name_to_idx
            .get(name)
            .copied()
            .ok_or_else(|| Error::Other(format!("Species '{}' not found in tree", name)))
    };

    let n = event_types.len();
    let mut dtl_events = Vec::with_capacity(n);

    for i in 0..n {
        let sp_idx = resolve(&species_names[i])?;
        let event = match event_types[i].as_str() {
            "Speciation" => DTLEvent::Speciation {
                time: times[i],
                gene_id: gene_ids[i] as usize,
                species_id: sp_idx,
                left_child: 0,
                right_child: 0,
            },
            "Duplication" => DTLEvent::Duplication {
                time: times[i],
                gene_id: gene_ids[i] as usize,
                species_id: sp_idx,
                child1: 0,
                child2: 0,
            },
            "Transfer" => {
                let from_idx = resolve(&from_species_strs[i])?;
                let to_idx = resolve(&to_species_strs[i])?;
                DTLEvent::Transfer {
                    time: times[i],
                    gene_id: gene_ids[i] as usize,
                    species_id: sp_idx,
                    from_species: from_idx,
                    to_species: to_idx,
                    donor_child: 0,
                    recipient_child: 0,
                }
            }
            "Loss" => DTLEvent::Loss {
                time: times[i],
                gene_id: gene_ids[i] as usize,
                species_id: sp_idx,
            },
            "Leaf" => DTLEvent::Leaf {
                time: times[i],
                gene_id: gene_ids[i] as usize,
                species_id: sp_idx,
            },
            other => return Err(Error::Other(format!("Unknown event type: {}", other))),
        };
        dtl_events.push(event);
    }

    Ok(dtl_events)
}
