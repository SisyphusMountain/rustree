// Compute induced transfers: project transfer donors/recipients from a complete
// species tree onto a sampled subtree.
//
// Given a complete species tree, a subset of sampled leaf names, and a list of
// DTL events, we compute for each transfer the "induced donor" and "induced
// recipient" — the branches of the sampled tree that the original donor and
// recipient project onto.

use crate::dtl::DTLEvent;
use crate::error::RustreeError;
use crate::node::FlatTree;
use crate::sampling::{
    extract_induced_subtree_by_names, find_leaf_indices_by_names, mark_nodes_postorder, NodeMark,
};
use std::collections::{HashMap, HashSet};

/// An induced transfer: a transfer event with both complete-tree and
/// sampled-tree donor/recipient species.
#[derive(Clone, Debug)]
pub struct InducedTransfer {
    pub time: f64,
    pub gene_id: usize,
    /// Donor species index in the complete tree.
    pub from_species_complete: usize,
    /// Recipient species index in the complete tree.
    pub to_species_complete: usize,
    /// Induced donor: index in the sampled tree (None if projection fails).
    pub from_species_sampled: Option<usize>,
    /// Induced recipient: index in the sampled tree (None if projection fails).
    pub to_species_sampled: Option<usize>,
}

/// Algorithm mode for induced transfers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InducedTransferAlgorithm {
    /// Current rustree behavior: project every transfer event independently.
    Projection,
    /// Damien-style induced transfer inference on sampled lineages.
    DamienStyle,
}

#[derive(Clone, Copy, Debug)]
struct RawTransferEvent {
    time: f64,
    gene_id: usize,
    from_species: usize,
    to_species: usize,
}

#[derive(Clone, Copy, Debug)]
struct RawInducedTransfer {
    time: f64,
    gene_id: usize,
    donor: i64,
    recipient: i64,
}

/// Returns the sibling of `node_idx` (the other child of its parent).
fn get_sibling(tree: &FlatTree, node_idx: usize) -> Option<usize> {
    let parent_idx = tree.nodes[node_idx].parent?;
    let parent = &tree.nodes[parent_idx];
    if parent.left_child == Some(node_idx) {
        parent.right_child
    } else {
        parent.left_child
    }
}

/// Projects a single node onto the nearest Keep node by walking through the tree.
///
/// - **Keep**: returns itself.
/// - **HasDescendant**: walks down to the child that has kept descendants.
/// - **Discard**: walks to sibling; if sibling is also Discard, walks up to parent.
fn project_node(tree: &FlatTree, marks: &[NodeMark], node_idx: usize) -> Option<usize> {
    match marks[node_idx] {
        NodeMark::Keep => Some(node_idx),
        NodeMark::HasDescendant => {
            let left = tree.nodes[node_idx].left_child;
            let right = tree.nodes[node_idx].right_child;
            let next = match (left, right) {
                (Some(l), _) if marks[l] != NodeMark::Discard => l,
                (_, Some(r)) => r,
                _ => return None,
            };
            project_node(tree, marks, next)
        }
        NodeMark::Discard => {
            let sibling = get_sibling(tree, node_idx);
            match sibling {
                Some(s) if marks[s] != NodeMark::Discard => project_node(tree, marks, s),
                _ => {
                    // Sibling is Discard or doesn't exist — go up to parent
                    let parent = tree.nodes[node_idx].parent?;
                    project_node(tree, marks, parent)
                }
            }
        }
    }
}

/// Computes the projection of every node in the complete tree onto the nearest
/// Keep node (a node that is present in the sampled tree).
///
/// Returns a `Vec<Option<usize>>` indexed by complete-tree node index, where
/// each entry is the index (in the complete tree) of the Keep node it projects
/// onto, or `None` if no projection exists.
pub(crate) fn compute_projection(
    complete_tree: &FlatTree,
    marks: &[NodeMark],
) -> Vec<Option<usize>> {
    (0..complete_tree.nodes.len())
        .map(|idx| project_node(complete_tree, marks, idx))
        .collect()
}

/// Computes the ghost length for each branch of the sampled tree.
///
/// The ghost length of a sampled-tree branch is the sum of branch lengths of all
/// non-sampled nodes (HasDescendant and Discard) in the complete tree that project
/// onto that branch.
///
/// # Arguments
/// * `complete_tree` - The complete species tree
/// * `sampled_leaf_names` - Names of leaves to keep
///
/// # Returns
/// A tuple of `(sampled_tree, ghost_lengths)` where `ghost_lengths` is a `Vec<f64>`
/// indexed by sampled-tree node index.
pub fn ghost_lengths(
    complete_tree: &FlatTree,
    sampled_leaf_names: &[String],
) -> Result<(FlatTree, Vec<f64>), RustreeError> {
    let (sampled_tree, _) = extract_induced_subtree_by_names(complete_tree, sampled_leaf_names)
        .ok_or_else(|| {
            RustreeError::Tree("Failed to extract induced subtree from leaf names".to_string())
        })?;

    let keep_indices = find_leaf_indices_by_names(complete_tree, sampled_leaf_names);
    let mut marks = vec![NodeMark::Discard; complete_tree.nodes.len()];
    mark_nodes_postorder(complete_tree, complete_tree.root, &keep_indices, &mut marks);

    let projection = compute_projection(complete_tree, &marks);

    // Map complete-tree Keep-node names → sampled-tree indices
    let mut keep_name_to_sampled: HashMap<&str, usize> = HashMap::new();
    for (idx, node) in sampled_tree.nodes.iter().enumerate() {
        keep_name_to_sampled.insert(&node.name, idx);
    }

    let mut ghost = vec![0.0; sampled_tree.nodes.len()];

    for (complete_idx, mark) in marks.iter().enumerate() {
        if *mark == NodeMark::Keep {
            continue;
        }
        // Non-Keep node: add its branch length to the projected sampled-tree branch
        if let Some(keep_idx) = projection[complete_idx] {
            let name = &complete_tree.nodes[keep_idx].name;
            if let Some(&sampled_idx) = keep_name_to_sampled.get(name.as_str()) {
                ghost[sampled_idx] += complete_tree.nodes[complete_idx].length;
            }
        }
    }

    Ok((sampled_tree, ghost))
}

/// Computes induced transfers by projecting each transfer's donor and recipient
/// from a complete species tree onto a sampled subtree.
///
/// # Arguments
/// * `complete_tree` - The complete species tree (must have depths assigned)
/// * `sampled_leaf_names` - Names of leaves to keep (must be a subset of complete tree leaves)
/// * `events` - DTL events from a gene tree simulation on the complete tree
///
/// # Returns
/// A vector of `InducedTransfer` — one per Transfer event in `events`.
pub fn induced_transfers(
    complete_tree: &FlatTree,
    sampled_leaf_names: &[String],
    events: &[DTLEvent],
) -> Result<Vec<InducedTransfer>, RustreeError> {
    let (sampled_tree, _) = extract_induced_subtree_by_names(complete_tree, sampled_leaf_names)
        .ok_or_else(|| {
            RustreeError::Tree("Failed to extract induced subtree from leaf names".to_string())
        })?;

    // Step 1: Mark complete tree nodes
    let keep_indices = find_leaf_indices_by_names(complete_tree, sampled_leaf_names);
    let mut marks = vec![NodeMark::Discard; complete_tree.nodes.len()];
    mark_nodes_postorder(complete_tree, complete_tree.root, &keep_indices, &mut marks);

    // Step 2: Compute projection (complete tree index → Keep node index in complete tree)
    let projection = compute_projection(complete_tree, &marks);

    // Step 3: Build mapping from complete tree Keep-node names to sampled tree indices
    let mut keep_name_to_sampled: HashMap<&str, usize> = HashMap::new();
    for (idx, node) in sampled_tree.nodes.iter().enumerate() {
        keep_name_to_sampled.insert(&node.name, idx);
    }

    // Step 4: Compose projection → Keep node name → sampled tree index
    let complete_to_sampled: Vec<Option<usize>> = projection
        .iter()
        .map(|proj| {
            proj.and_then(|keep_idx| {
                let name = &complete_tree.nodes[keep_idx].name;
                keep_name_to_sampled.get(name.as_str()).copied()
            })
        })
        .collect();

    // Step 5: For each transfer event, produce an InducedTransfer
    Ok(events
        .iter()
        .filter_map(|event| {
            if let DTLEvent::Transfer {
                time,
                gene_id,
                from_species,
                to_species,
                ..
            } = event
            {
                Some(InducedTransfer {
                    time: *time,
                    gene_id: *gene_id,
                    from_species_complete: *from_species,
                    to_species_complete: *to_species,
                    from_species_sampled: complete_to_sampled[*from_species],
                    to_species_sampled: complete_to_sampled[*to_species],
                })
            } else {
                None
            }
        })
        .collect())
}

/// Computes induced transfers with explicit algorithm selection.
///
/// `remove_undetectable` currently affects only `DamienStyle` mode.
pub fn induced_transfers_with_algorithm(
    complete_tree: &FlatTree,
    sampled_leaf_names: &[String],
    events: &[DTLEvent],
    algorithm: InducedTransferAlgorithm,
    remove_undetectable: bool,
) -> Result<Vec<InducedTransfer>, RustreeError> {
    match algorithm {
        InducedTransferAlgorithm::Projection => induced_transfers(complete_tree, sampled_leaf_names, events),
        InducedTransferAlgorithm::DamienStyle => induced_transfers_damien_style(
            complete_tree,
            sampled_leaf_names,
            events,
            remove_undetectable,
        ),
    }
}

fn induced_transfers_damien_style(
    complete_tree: &FlatTree,
    sampled_leaf_names: &[String],
    events: &[DTLEvent],
    remove_undetectable: bool,
) -> Result<Vec<InducedTransfer>, RustreeError> {
    let mut transfer_events: Vec<RawTransferEvent> = events
        .iter()
        .filter_map(|e| {
            if let DTLEvent::Transfer {
                time,
                gene_id,
                from_species,
                to_species,
                ..
            } = e
            {
                Some(RawTransferEvent {
                    time: *time,
                    gene_id: *gene_id,
                    from_species: *from_species,
                    to_species: *to_species,
                })
            } else {
                None
            }
        })
        .collect();
    transfer_events.sort_by(|a, b| a.time.total_cmp(&b.time));

    if transfer_events.is_empty() {
        return Ok(Vec::new());
    }

    let sampled_name_set: HashSet<&str> = sampled_leaf_names.iter().map(|s| s.as_str()).collect();

    let mut sampled_leaf_indices = HashSet::new();
    let mut ghost_leaf_indices = HashSet::new();
    let mut non_ghost_leaf_indices = Vec::new();
    let mut leaf_order = Vec::new();
    collect_leaf_order(complete_tree, complete_tree.root, &mut leaf_order);
    for &idx in &leaf_order {
        if sampled_name_set.contains(complete_tree.nodes[idx].name.as_str()) {
            sampled_leaf_indices.insert(idx);
            non_ghost_leaf_indices.push(idx);
        } else {
            ghost_leaf_indices.insert(idx);
        }
    }

    let mut extinct_nodes: HashSet<i64> = ghost_leaf_indices.iter().map(|&i| i as i64).collect();

    // Extended edge matrix (Damien's signed transfer nodes).
    let original_count = complete_tree.nodes.len() as i64;
    let mut edges: Vec<(i64, i64)> = complete_tree
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(child, n)| n.parent.map(|p| (p as i64, child as i64)))
        .collect();

    let mut transfer_node_info: HashMap<i64, (f64, usize)> = HashMap::new();
    let mut curr_transfer_node = original_count;
    for tr in &transfer_events {
        split_edge_with_transfer_nodes(
            &mut edges,
            tr.from_species as i64,
            tr.to_species as i64,
            curr_transfer_node,
        )?;
        transfer_node_info.insert(-curr_transfer_node, (tr.time, tr.gene_id));
        curr_transfer_node += 1;
    }

    let child_to_parent: HashMap<i64, i64> = edges.iter().map(|&(p, c)| (c, p)).collect();

    // Damien-style extinction mask on the extended graph:
    // nodes that are not ancestors of any non-ghost tip are extinct.
    let mut all_extended_nodes: HashSet<i64> = HashSet::new();
    for &(p, c) in &edges {
        all_extended_nodes.insert(p);
        all_extended_nodes.insert(c);
    }
    all_extended_nodes.insert(complete_tree.root as i64);

    let mut non_extinct_extended: HashSet<i64> = HashSet::new();
    for &tip in &non_ghost_leaf_indices {
        let mut cur = tip as i64;
        while let Some(parent) = child_to_parent.get(&cur).copied() {
            non_extinct_extended.insert(parent);
            cur = parent;
        }
    }
    for &node in &all_extended_nodes {
        if !non_extinct_extended.contains(&node)
            && !non_ghost_leaf_indices.contains(&(node as usize))
        {
            extinct_nodes.insert(node);
        }
    }

    let mut already_dealt: HashSet<i64> = HashSet::new();
    already_dealt.insert(complete_tree.root as i64);
    let mut the_path: HashMap<i64, i64> = HashMap::new();
    the_path.insert(complete_tree.root as i64, complete_tree.root as i64);

    let mut tr_lists: Vec<Vec<RawInducedTransfer>> = Vec::new();
    for &tip in &non_ghost_leaf_indices {
        let mut s = tip as i64;
        let mut transfer_mode = false;
        let mut recipient = 0i64;
        let mut extinct_mode = false;
        let mut extinct_path: Vec<i64> = Vec::new();
        let mut local_rows: Vec<RawInducedTransfer> = Vec::new();

        loop {
            let up = if s < 0 {
                let u = -s;
                if !transfer_mode {
                    recipient = s;
                    transfer_mode = true;
                }
                u
            } else if let Some(u) = child_to_parent.get(&s).copied() {
                u
            } else {
                break;
            };

            if !extinct_nodes.contains(&up) {
                if transfer_mode {
                    let (time, gene_id) = transfer_node_info.get(&recipient).copied().unwrap_or((0.0, 0));
                    local_rows.push(RawInducedTransfer {
                        time,
                        gene_id,
                        donor: up,
                        recipient,
                    });
                    transfer_mode = false;
                }
                if extinct_mode {
                    extinct_path.push(up);
                    for &x in &extinct_path {
                        the_path.insert(x.abs(), up);
                    }
                    extinct_mode = false;
                }
            } else if !extinct_mode {
                extinct_path.clear();
                extinct_path.push(up);
                extinct_mode = true;
            } else {
                extinct_path.push(up);
            }

            if already_dealt.contains(&up) {
                if transfer_mode {
                    let donor = the_path.get(&up).copied().unwrap_or(up);
                    let (time, gene_id) = transfer_node_info.get(&recipient).copied().unwrap_or((0.0, 0));
                    local_rows.push(RawInducedTransfer {
                        time,
                        gene_id,
                        donor,
                        recipient,
                    });
                }
                if extinct_mode {
                    let map_to = the_path.get(&up).copied().unwrap_or(up);
                    for &x in &extinct_path {
                        the_path.insert(x.abs(), map_to);
                    }
                }
                break;
            }

            already_dealt.insert(up);
            s = up;
        }

        if !local_rows.is_empty() {
            if remove_undetectable {
                simplify_cross_back_transfers(&mut local_rows);
            }
            if !local_rows.is_empty() {
                tr_lists.push(local_rows);
            }
        }
    }

    if tr_lists.is_empty() {
        return Ok(Vec::new());
    }

    let (sampled_tree, _) = extract_induced_subtree_by_names(complete_tree, sampled_leaf_names)
        .ok_or_else(|| {
            RustreeError::Tree("Failed to extract induced subtree from leaf names".to_string())
        })?;

    let sampled_marks = {
        let keep_indices = find_leaf_indices_by_names(complete_tree, sampled_leaf_names);
        let mut marks = vec![NodeMark::Discard; complete_tree.nodes.len()];
        mark_nodes_postorder(complete_tree, complete_tree.root, &keep_indices, &mut marks);
        marks
    };
    let projection = compute_projection(complete_tree, &sampled_marks);

    let mut keep_name_to_sampled: HashMap<&str, usize> = HashMap::new();
    for (idx, node) in sampled_tree.nodes.iter().enumerate() {
        keep_name_to_sampled.insert(node.name.as_str(), idx);
    }
    let complete_to_sampled: Vec<Option<usize>> = projection
        .iter()
        .map(|proj| {
            proj.and_then(|keep_idx| {
                keep_name_to_sampled
                    .get(complete_tree.nodes[keep_idx].name.as_str())
                    .copied()
            })
        })
        .collect();

    let transfer_desc_map = closest_existing_desc_map(&edges, original_count);

    let mut all_rows: Vec<RawInducedTransfer> = tr_lists.into_iter().flatten().collect();
    all_rows.sort_by(|a, b| a.time.total_cmp(&b.time));

    let mut out = Vec::with_capacity(all_rows.len());
    for row in all_rows {
        let from_complete = map_damien_node_to_complete_keep(row.donor, original_count, &transfer_desc_map, &projection)?;
        let to_complete = map_damien_node_to_complete_keep(row.recipient, original_count, &transfer_desc_map, &projection)?;

        let from_sampled = complete_to_sampled[from_complete];
        let to_sampled = complete_to_sampled[to_complete];

        if remove_undetectable {
            if from_sampled.is_none() || to_sampled.is_none() {
                continue;
            }
            let fs = from_sampled.unwrap();
            let ts = to_sampled.unwrap();
            if fs == ts {
                continue;
            }
            if are_sister_in_sampled_tree(&sampled_tree, fs, ts)
                || are_direct_edge_in_sampled_tree(&sampled_tree, fs, ts)
            {
                continue;
            }
        }

        out.push(InducedTransfer {
            time: row.time,
            gene_id: row.gene_id,
            from_species_complete: from_complete,
            to_species_complete: to_complete,
            from_species_sampled: from_sampled,
            to_species_sampled: to_sampled,
        });
    }

    Ok(out)
}

fn collect_leaf_order(tree: &FlatTree, idx: usize, out: &mut Vec<usize>) {
    let n = &tree.nodes[idx];
    if n.left_child.is_none() && n.right_child.is_none() {
        out.push(idx);
        return;
    }
    if let Some(l) = n.left_child {
        collect_leaf_order(tree, l, out);
    }
    if let Some(r) = n.right_child {
        collect_leaf_order(tree, r, out);
    }
}

fn split_edge_with_transfer_nodes(
    edges: &mut Vec<(i64, i64)>,
    from: i64,
    to: i64,
    transfer_node: i64,
) -> Result<(), RustreeError> {
    // Recipient split: p->to => p->(-k), (-k)->to
    let rec_pos = edges
        .iter()
        .position(|&(_, c)| c == to)
        .ok_or_else(|| RustreeError::Tree(format!("Cannot find recipient edge for node {}", to)))?;
    let (rec_parent, _) = edges[rec_pos];
    edges.swap_remove(rec_pos);
    edges.push((rec_parent, -transfer_node));
    edges.push((-transfer_node, to));

    // Donor split: p->from => p->k, k->from
    let don_pos = edges
        .iter()
        .position(|&(_, c)| c == from)
        .ok_or_else(|| RustreeError::Tree(format!("Cannot find donor edge for node {}", from)))?;
    let (don_parent, _) = edges[don_pos];
    edges.swap_remove(don_pos);
    edges.push((don_parent, transfer_node));
    edges.push((transfer_node, from));

    Ok(())
}

fn closest_existing_desc_map(edges: &[(i64, i64)], transfer_offset: i64) -> HashMap<i64, i64> {
    let sub: Vec<(i64, i64)> = edges
        .iter()
        .copied()
        .filter(|(p, _)| *p >= transfer_offset || *p <= -transfer_offset)
        .collect();
    let next: HashMap<i64, i64> = sub.iter().map(|(p, c)| (*p, *c)).collect();

    let mut map = HashMap::new();
    for (start, first_child) in sub {
        let mut end = first_child;
        while let Some(n) = next.get(&end).copied() {
            end = n;
        }
        map.insert(start, end);
    }
    map
}

fn map_damien_node_to_complete_keep(
    node: i64,
    transfer_offset: i64,
    transfer_desc_map: &HashMap<i64, i64>,
    projection: &[Option<usize>],
) -> Result<usize, RustreeError> {
    let mut cur = node;
    for _ in 0..8 {
        if cur.abs() >= transfer_offset {
            if let Some(next) = transfer_desc_map.get(&cur).copied() {
                cur = next;
                continue;
            }
            return Err(RustreeError::Tree(format!(
                "No descendant mapping for transfer node {}",
                cur
            )));
        }
        if cur < 0 {
            return Err(RustreeError::Tree(format!(
                "Unexpected negative non-transfer node {}",
                cur
            )));
        }
        let idx = cur as usize;
        return projection[idx].ok_or_else(|| {
            RustreeError::Tree(format!("Node {} did not project onto sampled tree", idx))
        });
    }
    Err(RustreeError::Tree(
        "Too many transfer-node remapping hops".to_string(),
    ))
}

fn simplify_cross_back_transfers(rows: &mut Vec<RawInducedTransfer>) {
    if rows.len() < 2 {
        return;
    }
    let mut keep = vec![true; rows.len()];
    for i in 0..rows.len() - 1 {
        if rows[i].donor == rows[i + 1].recipient && rows[i].recipient == rows[i + 1].donor {
            keep[i] = false;
        }
    }
    let mut new_rows = Vec::with_capacity(rows.len());
    for (i, r) in rows.iter().enumerate() {
        if keep[i] {
            new_rows.push(*r);
        }
    }
    *rows = new_rows;
}

fn are_sister_in_sampled_tree(tree: &FlatTree, from_idx: usize, to_idx: usize) -> bool {
    tree.nodes[from_idx].parent.is_some() && tree.nodes[from_idx].parent == tree.nodes[to_idx].parent
}

fn are_direct_edge_in_sampled_tree(tree: &FlatTree, from_idx: usize, to_idx: usize) -> bool {
    tree.nodes[from_idx].parent == Some(to_idx) || tree.nodes[to_idx].parent == Some(from_idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::newick::parse_newick;
    use crate::sampling::extract_induced_subtree_by_names;
    use std::collections::HashSet;
    use std::fs;

    fn make_tree(newick: &str) -> FlatTree {
        let mut nodes = parse_newick(newick).unwrap();
        let mut tree = nodes.pop().unwrap().to_flat_tree();
        tree.assign_depths();
        tree
    }

    #[test]
    fn test_projection_balanced_tree() {
        // Tree: ((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;
        // Sample: {A, C}
        // Expected projections:
        //   A → A (Keep), B → A (sibling), C → C (Keep), D → C (sibling)
        //   AB → A (HasDescendant, left child A is Keep)
        //   CD → C (HasDescendant, left child C is Keep)
        //   root → root (Keep, both subtrees have kept descendants)
        let tree = make_tree("((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "C".into()];

        let keep_indices = find_leaf_indices_by_names(&tree, &sampled_names);
        let mut marks = vec![NodeMark::Discard; tree.nodes.len()];
        mark_nodes_postorder(&tree, tree.root, &keep_indices, &mut marks);

        let projection = compute_projection(&tree, &marks);

        // Get node indices by name
        let idx = |name: &str| tree.find_node_index(name).unwrap();

        // Keep nodes project to themselves
        assert_eq!(projection[idx("A")], Some(idx("A")));
        assert_eq!(projection[idx("C")], Some(idx("C")));
        assert_eq!(projection[idx("root")], Some(idx("root")));

        // B projects to its sibling A
        assert_eq!(projection[idx("B")], Some(idx("A")));
        // D projects to its sibling C
        assert_eq!(projection[idx("D")], Some(idx("C")));

        // AB (HasDescendant) projects down to A
        assert_eq!(projection[idx("AB")], Some(idx("A")));
        // CD (HasDescendant) projects down to C
        assert_eq!(projection[idx("CD")], Some(idx("C")));
    }

    #[test]
    fn test_projection_one_subtree_sampled() {
        // Tree: ((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;
        // Sample: {A, B} — only left subtree
        // root is HasDescendant, CD subtree is all Discard
        let tree = make_tree("((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "B".into()];

        let keep_indices = find_leaf_indices_by_names(&tree, &sampled_names);
        let mut marks = vec![NodeMark::Discard; tree.nodes.len()];
        mark_nodes_postorder(&tree, tree.root, &keep_indices, &mut marks);

        let projection = compute_projection(&tree, &marks);

        let idx = |name: &str| tree.find_node_index(name).unwrap();

        // AB is Keep (MRCA of A, B)
        assert_eq!(projection[idx("AB")], Some(idx("AB")));
        assert_eq!(projection[idx("A")], Some(idx("A")));
        assert_eq!(projection[idx("B")], Some(idx("B")));

        // root is HasDescendant → projects down to AB
        assert_eq!(projection[idx("root")], Some(idx("AB")));

        // CD is Discard → sibling AB has kept descendants → projects to AB
        assert_eq!(projection[idx("CD")], Some(idx("AB")));

        // C is Discard → sibling D is Discard → go up to CD (Discard) →
        // sibling of CD is AB → projects to AB
        assert_eq!(projection[idx("C")], Some(idx("AB")));
        assert_eq!(projection[idx("D")], Some(idx("AB")));
    }

    #[test]
    fn test_induced_transfers_basic() {
        let complete_tree = make_tree("((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "C".into()];

        let idx = |name: &str| complete_tree.find_node_index(name).unwrap();

        // A transfer from B (complete) to D (complete)
        let events = vec![DTLEvent::Transfer {
            time: 0.5,
            gene_id: 0,
            species_id: idx("B"),
            from_species: idx("B"),
            to_species: idx("D"),
            donor_child: 1,
            recipient_child: 2,
        }];

        let induced = induced_transfers(&complete_tree, &sampled_names, &events).unwrap();
        assert_eq!(induced.len(), 1);

        let t = &induced[0];
        assert_eq!(t.from_species_complete, idx("B"));
        assert_eq!(t.to_species_complete, idx("D"));

        // B projects to A in the complete tree, then A maps to sampled tree
        let sampled_tree = extract_induced_subtree_by_names(&complete_tree, &sampled_names)
            .unwrap()
            .0;
        let sampled_a = sampled_tree.find_node_index("A").unwrap();
        let sampled_c = sampled_tree.find_node_index("C").unwrap();
        assert_eq!(t.from_species_sampled, Some(sampled_a));
        assert_eq!(t.to_species_sampled, Some(sampled_c));
    }

    #[test]
    fn test_ghost_lengths_balanced() {
        // Tree: ((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;
        // Sample: {A, C}
        //
        // Non-Keep nodes and their projections:
        //   B  (length 1) → projects to A
        //   D  (length 1) → projects to C
        //   AB (length 1) → projects to A
        //   CD (length 1) → projects to C
        //
        // Ghost lengths: A = 1 + 1 = 2, C = 1 + 1 = 2, root = 0
        let complete = make_tree("((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "C".into()];

        let (sampled, ghost) = ghost_lengths(&complete, &sampled_names).unwrap();

        let sampled_a = sampled.find_node_index("A").unwrap();
        let sampled_c = sampled.find_node_index("C").unwrap();

        assert!(
            (ghost[sampled_a] - 2.0).abs() < 1e-10,
            "Ghost length of A should be 2.0, got {}",
            ghost[sampled_a]
        );
        assert!(
            (ghost[sampled_c] - 2.0).abs() < 1e-10,
            "Ghost length of C should be 2.0, got {}",
            ghost[sampled_c]
        );
    }

    #[test]
    fn test_ghost_lengths_one_subtree() {
        // Tree: ((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;
        // Sample: {A, B} — only left subtree
        //
        // Non-Keep nodes: root (length 0), CD (length 1), C (length 1), D (length 1)
        // All project to AB.
        // Ghost length of AB = 0 + 1 + 1 + 1 = 3
        let complete = make_tree("((A:1,B:1)AB:1,(C:1,D:1)CD:1)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "B".into()];

        let (sampled, ghost) = ghost_lengths(&complete, &sampled_names).unwrap();

        let sampled_ab = sampled.find_node_index("AB").unwrap();
        let sampled_a = sampled.find_node_index("A").unwrap();
        let sampled_b = sampled.find_node_index("B").unwrap();

        assert!(
            (ghost[sampled_ab] - 3.0).abs() < 1e-10,
            "Ghost length of AB should be 3.0, got {}",
            ghost[sampled_ab]
        );
        assert!(
            (ghost[sampled_a]).abs() < 1e-10,
            "Ghost length of A should be 0.0"
        );
        assert!(
            (ghost[sampled_b]).abs() < 1e-10,
            "Ghost length of B should be 0.0"
        );
    }

    #[test]
    fn test_induced_transfers_non_transfer_events_filtered() {
        let complete_tree = make_tree("((A:1,B:1)AB:1,C:2)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "C".into()];

        let events = vec![
            DTLEvent::Duplication {
                time: 0.3,
                gene_id: 0,
                species_id: 0,
                child1: 1,
                child2: 2,
            },
            DTLEvent::Loss {
                time: 0.4,
                gene_id: 1,
                species_id: 0,
            },
            DTLEvent::Leaf {
                time: 1.0,
                gene_id: 2,
                species_id: 0,
            },
        ];

        let induced = induced_transfers(&complete_tree, &sampled_names, &events).unwrap();
        assert_eq!(
            induced.len(),
            0,
            "Non-transfer events should be filtered out"
        );
    }

    #[test]
    fn test_damien_style_matches_expected_fixture() {
        let base = "testdata/induced_tr";
        let complete_newick = fs::read_to_string(format!("{}/complete_tree.nwk", base)).unwrap();
        let sampled_newick = fs::read_to_string(format!("{}/sampled_tree.nwk", base)).unwrap();
        let events_tsv = fs::read_to_string(format!("{}/complete_events.tsv", base)).unwrap();
        let expected_tsv = fs::read_to_string(format!("{}/expected_induced.tsv", base)).unwrap();

        let mut complete_nodes = parse_newick(complete_newick.trim()).unwrap();
        let complete = complete_nodes.pop().unwrap().to_flat_tree();

        let mut sampled_nodes = parse_newick(sampled_newick.trim()).unwrap();
        let sampled = sampled_nodes.pop().unwrap().to_flat_tree();
        let sampled_names_set: HashSet<String> = sampled
            .nodes
            .iter()
            .filter(|n| n.left_child.is_none() && n.right_child.is_none())
            .map(|n| n.name.split('_').next().unwrap_or(n.name.as_str()).to_string())
            .collect();
        let mut sampled_names: Vec<String> = sampled_names_set.into_iter().collect();
        sampled_names.sort();

        let mut events: Vec<DTLEvent> = Vec::new();
        for (i, line) in events_tsv.lines().enumerate() {
            if i == 0 || line.trim().is_empty() {
                continue;
            }
            let mut cols = line.split('\t');
            let time_s = cols.next().unwrap();
            let event_type = cols.next().unwrap();
            let nodes_s = cols.next().unwrap();
            if event_type != "T" {
                continue;
            }
            let p: Vec<&str> = nodes_s.split(';').collect();
            let from_species_name = p[0];
            let gene_id: usize = p[1].parse().unwrap();
            let to_species_name = p[4];
            let time: f64 = time_s.parse().unwrap();

            let from_species = complete.find_node_index(from_species_name).unwrap();
            let to_species = complete.find_node_index(to_species_name).unwrap();

            events.push(DTLEvent::Transfer {
                time,
                gene_id,
                species_id: from_species,
                from_species,
                to_species,
                donor_child: 0,
                recipient_child: 0,
            });
        }

        let induced = induced_transfers_with_algorithm(
            &complete,
            &sampled_names,
            &events,
            InducedTransferAlgorithm::DamienStyle,
            false,
        )
        .unwrap();

        let observed: Vec<(String, String)> = induced
            .iter()
            .map(|t| {
                (
                    complete.nodes[t.from_species_complete].name.clone(),
                    complete.nodes[t.to_species_complete].name.clone(),
                )
            })
            .collect();

        let expected: Vec<(String, String)> = expected_tsv
            .lines()
            .skip(1)
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let mut c = l.split('\t');
                (
                    c.next().unwrap().to_string(),
                    c.next().unwrap().to_string(),
                )
            })
            .collect();

        assert_eq!(observed, expected, "Damien-style induced transfers should match expected_induced.tsv exactly");

        let reduced = induced_transfers_with_algorithm(
            &complete,
            &sampled_names,
            &events,
            InducedTransferAlgorithm::DamienStyle,
            true,
        )
        .unwrap();
        assert_eq!(
            reduced.len(),
            15,
            "With remove_undetectable=true this fixture should yield 15 transfers"
        );
    }
}
