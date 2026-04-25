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
    extract_induced_subtree, find_leaf_indices_by_names, mark_nodes_postorder, NodeMark,
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

struct SampleProjection {
    sampled_tree: FlatTree,
    marks: Vec<NodeMark>,
    projection: Vec<Option<usize>>,
    complete_to_sampled: Vec<Option<usize>>,
}

impl SampleProjection {
    fn new(complete_tree: &FlatTree, sampled_leaf_names: &[String]) -> Result<Self, RustreeError> {
        let keep_indices = find_leaf_indices_by_names(complete_tree, sampled_leaf_names);
        let (sampled_tree, old_to_new) = extract_induced_subtree(complete_tree, &keep_indices)
            .ok_or_else(|| {
                RustreeError::Tree("Failed to extract induced subtree from leaf names".to_string())
            })?;

        let mut marks = vec![NodeMark::Discard; complete_tree.nodes.len()];
        mark_nodes_postorder(complete_tree, complete_tree.root, &keep_indices, &mut marks);

        let projection = compute_projection(complete_tree, &marks);
        let complete_to_sampled = projection
            .iter()
            .map(|proj| proj.and_then(|keep_idx| old_to_new[keep_idx]))
            .collect();

        Ok(Self {
            sampled_tree,
            marks,
            projection,
            complete_to_sampled,
        })
    }

    fn sampled_index(&self, complete_idx: usize) -> Result<Option<usize>, RustreeError> {
        self.complete_to_sampled
            .get(complete_idx)
            .copied()
            .ok_or_else(|| {
                RustreeError::Index(format!(
                    "species index {} is out of bounds for complete tree projection",
                    complete_idx
                ))
            })
    }
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
    donor: ExtendedNodeId,
    recipient: ExtendedNodeId,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct ExtendedNodeId(i64);

impl ExtendedNodeId {
    fn original(idx: usize) -> Self {
        Self(idx as i64)
    }

    fn donor_transfer(idx: i64) -> Self {
        Self(idx)
    }

    fn recipient_transfer(idx: i64) -> Self {
        Self(-idx)
    }

    fn is_recipient_side(self) -> bool {
        self.0 < 0
    }

    fn counterpart(self) -> Self {
        Self(-self.0)
    }

    fn canonical_path_key(self) -> Self {
        Self(self.0.abs())
    }

    fn is_transfer(self, transfer_offset: i64) -> bool {
        self.0.abs() >= transfer_offset
    }

    fn as_original_index(self) -> Result<usize, RustreeError> {
        if self.0 < 0 {
            return Err(RustreeError::Tree(format!(
                "Unexpected negative non-transfer node {}",
                self
            )));
        }
        Ok(self.0 as usize)
    }
}

impl std::fmt::Display for ExtendedNodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
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
    let projection = SampleProjection::new(complete_tree, sampled_leaf_names)?;
    let mut ghost = vec![0.0; projection.sampled_tree.nodes.len()];

    for (complete_idx, mark) in projection.marks.iter().enumerate() {
        if *mark == NodeMark::Keep {
            continue;
        }
        if let Some(sampled_idx) = projection.complete_to_sampled[complete_idx] {
            ghost[sampled_idx] += complete_tree.nodes[complete_idx].length;
        }
    }

    Ok((projection.sampled_tree, ghost))
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
    let projection = SampleProjection::new(complete_tree, sampled_leaf_names)?;
    let transfer_events = transfer_events_from_dtl(events, complete_tree.nodes.len())?;
    let mut out = Vec::with_capacity(transfer_events.len());

    for transfer in transfer_events {
        out.push(InducedTransfer {
            time: transfer.time,
            gene_id: transfer.gene_id,
            from_species_complete: transfer.from_species,
            to_species_complete: transfer.to_species,
            from_species_sampled: projection.sampled_index(transfer.from_species)?,
            to_species_sampled: projection.sampled_index(transfer.to_species)?,
        });
    }

    Ok(out)
}

fn transfer_events_from_dtl(
    events: &[DTLEvent],
    node_count: usize,
) -> Result<Vec<RawTransferEvent>, RustreeError> {
    let mut transfers = Vec::new();
    for event in events {
        if let DTLEvent::Transfer {
            time,
            gene_id,
            from_species,
            to_species,
            ..
        } = event
        {
            validate_species_index(*from_species, node_count, "donor species")?;
            validate_species_index(*to_species, node_count, "recipient species")?;
            transfers.push(RawTransferEvent {
                time: *time,
                gene_id: *gene_id,
                from_species: *from_species,
                to_species: *to_species,
            });
        }
    }
    Ok(transfers)
}

fn validate_species_index(idx: usize, node_count: usize, label: &str) -> Result<(), RustreeError> {
    if idx >= node_count {
        return Err(RustreeError::Index(format!(
            "{} index {} is out of bounds for complete tree with {} nodes",
            label, idx, node_count
        )));
    }
    Ok(())
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
        InducedTransferAlgorithm::Projection => {
            induced_transfers(complete_tree, sampled_leaf_names, events)
        }
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
    let mut transfer_events = transfer_events_from_dtl(events, complete_tree.nodes.len())?;
    transfer_events.sort_by(|a, b| a.time.total_cmp(&b.time));

    if transfer_events.is_empty() {
        return Ok(Vec::new());
    }

    let projection = SampleProjection::new(complete_tree, sampled_leaf_names)?;
    let leaf_partition = LeafPartition::new(complete_tree, sampled_leaf_names);
    let graph = ExtendedTransferGraph::new(complete_tree, &transfer_events)?;
    let extinct_nodes = compute_extinct_nodes(complete_tree, &graph, &leaf_partition);
    let tr_lists = collect_damien_rows(
        complete_tree.root,
        &leaf_partition.non_ghost_leaf_indices,
        &graph,
        &extinct_nodes,
        remove_undetectable,
    )?;

    if tr_lists.is_empty() {
        return Ok(Vec::new());
    }

    let mut all_rows: Vec<RawInducedTransfer> = tr_lists.into_iter().flatten().collect();
    induced_transfers_from_damien_rows(&mut all_rows, &graph, &projection, remove_undetectable)
}

struct LeafPartition {
    ghost_leaf_indices: HashSet<usize>,
    non_ghost_leaf_indices: Vec<usize>,
}

impl LeafPartition {
    fn new(tree: &FlatTree, sampled_leaf_names: &[String]) -> Self {
        let sampled_name_set: HashSet<&str> =
            sampled_leaf_names.iter().map(|s| s.as_str()).collect();
        let mut ghost_leaf_indices = HashSet::new();
        let mut non_ghost_leaf_indices = Vec::new();
        let mut leaf_order = Vec::new();

        collect_leaf_order(tree, tree.root, &mut leaf_order);
        for idx in leaf_order {
            if sampled_name_set.contains(tree.nodes[idx].name.as_str()) {
                non_ghost_leaf_indices.push(idx);
            } else {
                ghost_leaf_indices.insert(idx);
            }
        }

        Self {
            ghost_leaf_indices,
            non_ghost_leaf_indices,
        }
    }
}

struct ExtendedTransferGraph {
    edges: Vec<(ExtendedNodeId, ExtendedNodeId)>,
    child_to_parent: HashMap<ExtendedNodeId, ExtendedNodeId>,
    transfer_node_info: HashMap<ExtendedNodeId, (f64, usize)>,
    original_count: i64,
}

impl ExtendedTransferGraph {
    fn new(tree: &FlatTree, transfer_events: &[RawTransferEvent]) -> Result<Self, RustreeError> {
        let original_count = tree.nodes.len() as i64;
        let mut edges: Vec<(ExtendedNodeId, ExtendedNodeId)> = tree
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(child, n)| {
                n.parent
                    .map(|p| (ExtendedNodeId::original(p), ExtendedNodeId::original(child)))
            })
            .collect();

        let mut transfer_node_info = HashMap::new();
        for (curr_transfer_node, tr) in (original_count..).zip(transfer_events.iter()) {
            split_edge_with_transfer_nodes(
                &mut edges,
                tr.from_species,
                tr.to_species,
                curr_transfer_node,
            )?;
            transfer_node_info.insert(
                ExtendedNodeId::recipient_transfer(curr_transfer_node),
                (tr.time, tr.gene_id),
            );
        }

        let child_to_parent = edges.iter().map(|&(p, c)| (c, p)).collect();

        Ok(Self {
            edges,
            child_to_parent,
            transfer_node_info,
            original_count,
        })
    }

    fn all_nodes(&self, root: usize) -> HashSet<ExtendedNodeId> {
        let mut nodes = HashSet::new();
        for &(parent, child) in &self.edges {
            nodes.insert(parent);
            nodes.insert(child);
        }
        nodes.insert(ExtendedNodeId::original(root));
        nodes
    }

    fn transfer_metadata(&self, recipient: ExtendedNodeId) -> Result<(f64, usize), RustreeError> {
        self.transfer_node_info
            .get(&recipient)
            .copied()
            .ok_or_else(|| {
                RustreeError::Tree(format!(
                    "Missing transfer metadata for recipient-side transfer node {}",
                    recipient
                ))
            })
    }
}

fn compute_extinct_nodes(
    complete_tree: &FlatTree,
    graph: &ExtendedTransferGraph,
    leaf_partition: &LeafPartition,
) -> HashSet<ExtendedNodeId> {
    let all_extended_nodes = graph.all_nodes(complete_tree.root);
    let non_ghost_tips: HashSet<ExtendedNodeId> = leaf_partition
        .non_ghost_leaf_indices
        .iter()
        .copied()
        .map(ExtendedNodeId::original)
        .collect();
    let mut extinct_nodes: HashSet<ExtendedNodeId> = leaf_partition
        .ghost_leaf_indices
        .iter()
        .copied()
        .map(ExtendedNodeId::original)
        .collect();

    // Damien-style extinction mask on the extended graph:
    // nodes that are not ancestors of any non-ghost tip are extinct.
    let mut non_extinct_extended: HashSet<ExtendedNodeId> = HashSet::new();
    for &tip in &leaf_partition.non_ghost_leaf_indices {
        let mut cur = ExtendedNodeId::original(tip);
        while let Some(parent) = graph.child_to_parent.get(&cur).copied() {
            non_extinct_extended.insert(parent);
            cur = parent;
        }
    }

    for &node in &all_extended_nodes {
        if !non_extinct_extended.contains(&node) && !non_ghost_tips.contains(&node) {
            extinct_nodes.insert(node);
        }
    }

    extinct_nodes
}

fn collect_damien_rows(
    root: usize,
    non_ghost_leaf_indices: &[usize],
    graph: &ExtendedTransferGraph,
    extinct_nodes: &HashSet<ExtendedNodeId>,
    remove_undetectable: bool,
) -> Result<Vec<Vec<RawInducedTransfer>>, RustreeError> {
    let root = ExtendedNodeId::original(root);
    let mut already_dealt: HashSet<ExtendedNodeId> = HashSet::new();
    already_dealt.insert(root);
    let mut the_path: HashMap<ExtendedNodeId, ExtendedNodeId> = HashMap::new();
    the_path.insert(root, root);
    let mut tr_lists = Vec::new();

    for &tip in non_ghost_leaf_indices {
        let mut local_rows = collect_damien_rows_for_tip(
            ExtendedNodeId::original(tip),
            graph,
            extinct_nodes,
            &mut already_dealt,
            &mut the_path,
        )?;

        if !local_rows.is_empty() {
            if remove_undetectable {
                simplify_cross_back_transfers(&mut local_rows);
            }
            if !local_rows.is_empty() {
                tr_lists.push(local_rows);
            }
        }
    }

    Ok(tr_lists)
}

fn collect_damien_rows_for_tip(
    tip: ExtendedNodeId,
    graph: &ExtendedTransferGraph,
    extinct_nodes: &HashSet<ExtendedNodeId>,
    already_dealt: &mut HashSet<ExtendedNodeId>,
    the_path: &mut HashMap<ExtendedNodeId, ExtendedNodeId>,
) -> Result<Vec<RawInducedTransfer>, RustreeError> {
    let mut s = tip;
    let mut transfer_mode = false;
    let mut recipient = None;
    let mut extinct_mode = false;
    let mut extinct_path: Vec<ExtendedNodeId> = Vec::new();
    let mut local_rows = Vec::new();

    loop {
        let up = if s.is_recipient_side() {
            let u = s.counterpart();
            if !transfer_mode {
                recipient = Some(s);
                transfer_mode = true;
            }
            u
        } else if let Some(u) = graph.child_to_parent.get(&s).copied() {
            u
        } else {
            break;
        };

        if !extinct_nodes.contains(&up) {
            if transfer_mode {
                push_damien_transfer_row(&mut local_rows, graph, up, recipient)?;
                transfer_mode = false;
                recipient = None;
            }
            if extinct_mode {
                extinct_path.push(up);
                for &node in &extinct_path {
                    the_path.insert(node.canonical_path_key(), up);
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
                push_damien_transfer_row(&mut local_rows, graph, donor, recipient)?;
            }
            if extinct_mode {
                let map_to = the_path.get(&up).copied().unwrap_or(up);
                for &node in &extinct_path {
                    the_path.insert(node.canonical_path_key(), map_to);
                }
            }
            break;
        }

        already_dealt.insert(up);
        s = up;
    }

    Ok(local_rows)
}

fn push_damien_transfer_row(
    rows: &mut Vec<RawInducedTransfer>,
    graph: &ExtendedTransferGraph,
    donor: ExtendedNodeId,
    recipient: Option<ExtendedNodeId>,
) -> Result<(), RustreeError> {
    let recipient = recipient.ok_or_else(|| {
        RustreeError::Tree("Transfer mode had no recipient-side transfer node".to_string())
    })?;
    let (time, gene_id) = graph.transfer_metadata(recipient)?;
    rows.push(RawInducedTransfer {
        time,
        gene_id,
        donor,
        recipient,
    });
    Ok(())
}

fn induced_transfers_from_damien_rows(
    all_rows: &mut [RawInducedTransfer],
    graph: &ExtendedTransferGraph,
    projection: &SampleProjection,
    remove_undetectable: bool,
) -> Result<Vec<InducedTransfer>, RustreeError> {
    let transfer_desc_map = closest_existing_desc_map(&graph.edges, graph.original_count);
    all_rows.sort_by(|a, b| a.time.total_cmp(&b.time));
    let mut out = Vec::with_capacity(all_rows.len());

    for row in all_rows.iter().copied() {
        let from_complete = map_damien_node_to_complete_keep(
            row.donor,
            graph.original_count,
            &transfer_desc_map,
            &projection.projection,
        )?;
        let to_complete = map_damien_node_to_complete_keep(
            row.recipient,
            graph.original_count,
            &transfer_desc_map,
            &projection.projection,
        )?;

        let from_sampled = projection.sampled_index(from_complete)?;
        let to_sampled = projection.sampled_index(to_complete)?;

        if remove_undetectable {
            let (Some(fs), Some(ts)) = (from_sampled, to_sampled) else {
                continue;
            };
            if fs == ts {
                continue;
            }
            if are_sister_in_sampled_tree(&projection.sampled_tree, fs, ts)
                || are_direct_edge_in_sampled_tree(&projection.sampled_tree, fs, ts)
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
    edges: &mut Vec<(ExtendedNodeId, ExtendedNodeId)>,
    from: usize,
    to: usize,
    transfer_node: i64,
) -> Result<(), RustreeError> {
    let from = ExtendedNodeId::original(from);
    let to = ExtendedNodeId::original(to);
    let recipient_node = ExtendedNodeId::recipient_transfer(transfer_node);
    let donor_node = ExtendedNodeId::donor_transfer(transfer_node);

    // Recipient split: p->to => p->(-k), (-k)->to
    let rec_pos = edges
        .iter()
        .position(|&(_, c)| c == to)
        .ok_or_else(|| RustreeError::Tree(format!("Cannot find recipient edge for node {}", to)))?;
    let (rec_parent, _) = edges[rec_pos];
    edges.swap_remove(rec_pos);
    edges.push((rec_parent, recipient_node));
    edges.push((recipient_node, to));

    // Donor split: p->from => p->k, k->from
    let don_pos = edges
        .iter()
        .position(|&(_, c)| c == from)
        .ok_or_else(|| RustreeError::Tree(format!("Cannot find donor edge for node {}", from)))?;
    let (don_parent, _) = edges[don_pos];
    edges.swap_remove(don_pos);
    edges.push((don_parent, donor_node));
    edges.push((donor_node, from));

    Ok(())
}

fn closest_existing_desc_map(
    edges: &[(ExtendedNodeId, ExtendedNodeId)],
    transfer_offset: i64,
) -> HashMap<ExtendedNodeId, ExtendedNodeId> {
    let sub: Vec<(ExtendedNodeId, ExtendedNodeId)> = edges
        .iter()
        .copied()
        .filter(|(p, _)| p.is_transfer(transfer_offset))
        .collect();
    let next: HashMap<ExtendedNodeId, ExtendedNodeId> = sub.iter().map(|(p, c)| (*p, *c)).collect();

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
    node: ExtendedNodeId,
    transfer_offset: i64,
    transfer_desc_map: &HashMap<ExtendedNodeId, ExtendedNodeId>,
    projection: &[Option<usize>],
) -> Result<usize, RustreeError> {
    let mut cur = node;
    let mut seen = HashSet::new();

    loop {
        if !seen.insert(cur) {
            return Err(RustreeError::Tree(format!(
                "Cycle while remapping transfer node {}",
                cur
            )));
        }

        if cur.is_transfer(transfer_offset) {
            if let Some(next) = transfer_desc_map.get(&cur).copied() {
                cur = next;
                continue;
            }
            return Err(RustreeError::Tree(format!(
                "No descendant mapping for transfer node {}",
                cur
            )));
        }
        let idx = cur.as_original_index()?;
        if idx >= projection.len() {
            return Err(RustreeError::Index(format!(
                "mapped complete-tree node {} is out of bounds for projection with {} nodes",
                idx,
                projection.len()
            )));
        }
        return projection[idx].ok_or_else(|| {
            RustreeError::Tree(format!("Node {} did not project onto sampled tree", idx))
        });
    }
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
    tree.nodes[from_idx].parent.is_some()
        && tree.nodes[from_idx].parent == tree.nodes[to_idx].parent
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
    use std::path::Path;

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
    fn test_induced_transfers_duplicate_internal_names_use_index_mapping() {
        let complete_tree = make_tree("((A:1,B:1)X:1,(C:1,D:1)X:1)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "B".into(), "C".into(), "D".into()];

        let complete_a = complete_tree.find_node_index("A").unwrap();
        let complete_c = complete_tree.find_node_index("C").unwrap();
        let left_x = complete_tree.nodes[complete_a].parent.unwrap();
        let right_x = complete_tree.nodes[complete_c].parent.unwrap();
        assert_ne!(left_x, right_x);
        assert_eq!(complete_tree.nodes[left_x].name, "X");
        assert_eq!(complete_tree.nodes[right_x].name, "X");

        let events = vec![DTLEvent::Transfer {
            time: 0.5,
            gene_id: 0,
            species_id: left_x,
            from_species: left_x,
            to_species: right_x,
            donor_child: 1,
            recipient_child: 2,
        }];

        let induced = induced_transfers(&complete_tree, &sampled_names, &events).unwrap();
        assert_eq!(induced.len(), 1);

        let sampled_tree = extract_induced_subtree_by_names(&complete_tree, &sampled_names)
            .unwrap()
            .0;
        let sampled_a = sampled_tree.find_node_index("A").unwrap();
        let sampled_c = sampled_tree.find_node_index("C").unwrap();
        let sampled_left_x = sampled_tree.nodes[sampled_a].parent.unwrap();
        let sampled_right_x = sampled_tree.nodes[sampled_c].parent.unwrap();

        assert_ne!(sampled_left_x, sampled_right_x);
        assert_eq!(sampled_tree.nodes[sampled_left_x].name, "X");
        assert_eq!(sampled_tree.nodes[sampled_right_x].name, "X");
        assert_eq!(induced[0].from_species_sampled, Some(sampled_left_x));
        assert_eq!(induced[0].to_species_sampled, Some(sampled_right_x));
    }

    #[test]
    fn test_induced_transfers_invalid_species_index_returns_error() {
        let complete_tree = make_tree("((A:1,B:1)AB:1,C:2)root:0;");
        let sampled_names: Vec<String> = vec!["A".into(), "C".into()];
        let invalid_idx = complete_tree.nodes.len();
        let to_species = complete_tree.find_node_index("C").unwrap();
        let events = vec![DTLEvent::Transfer {
            time: 0.5,
            gene_id: 0,
            species_id: invalid_idx,
            from_species: invalid_idx,
            to_species,
            donor_child: 1,
            recipient_child: 2,
        }];

        let err = induced_transfers(&complete_tree, &sampled_names, &events).unwrap_err();
        assert!(matches!(
            err,
            RustreeError::Index(msg) if msg.contains("donor species index")
        ));
    }

    #[test]
    fn test_damien_style_matches_expected_fixture() {
        let base = "testdata/induced_tr";
        assert!(
            Path::new(base).exists(),
            "missing fixture directory: {base}"
        );

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
            .map(|n| {
                n.name
                    .split('_')
                    .next()
                    .unwrap_or(n.name.as_str())
                    .to_string()
            })
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
                (c.next().unwrap().to_string(), c.next().unwrap().to_string())
            })
            .collect();

        assert_eq!(
            observed, expected,
            "Damien-style induced transfers should match expected_induced.tsv exactly"
        );

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
            0,
            "With remove_undetectable=true this compact fixture should yield no detectable transfers"
        );
    }
}
