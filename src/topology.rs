use crate::error::RustreeError;
use crate::node::FlatTree;
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::hash::{Hash, Hasher};

type LeafOrder = (Vec<String>, Vec<Option<usize>>, Vec<bool>);
type RootedTraversal = (Vec<Option<usize>>, Vec<usize>);

/// Exact canonical key for an unrooted topology with labeled leaves.
///
/// Leaf labels are stored in sorted order, and each non-trivial split is encoded
/// as the normalized smaller side of the bipartition using a bitset over that
/// leaf order.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct UnrootedTopologyKey {
    leaf_names: Vec<String>,
    splits: Vec<Vec<u64>>,
}

impl UnrootedTopologyKey {
    /// Sorted leaf names used as the canonical bit ordering.
    pub fn leaf_names(&self) -> &[String] {
        &self.leaf_names
    }

    /// Canonical normalized non-trivial splits as bitsets.
    pub fn splits(&self) -> &[Vec<u64>] {
        &self.splits
    }

    /// Number of leaves in the topology.
    pub fn leaf_count(&self) -> usize {
        self.leaf_names.len()
    }

    /// Stable 64-bit hash of this exact key.
    pub fn hash64(&self) -> u64 {
        stable_hash64(self)
    }
}

/// Canonical recursive shape for an unlabeled rooted tree.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum UnlabeledShape {
    Leaf,
    Node(Vec<UnlabeledShape>),
}

/// Exact canonical key for an unlabeled unrooted topology.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct UnrootedShapeKey {
    shape: UnlabeledShape,
}

impl UnrootedShapeKey {
    /// Canonical recursive shape.
    pub fn shape(&self) -> &UnlabeledShape {
        &self.shape
    }

    /// Stable 64-bit hash of this exact key.
    pub fn hash64(&self) -> u64 {
        stable_hash64(self)
    }
}

impl FlatTree {
    /// Build an exact canonical key for the labeled unrooted topology.
    pub fn unrooted_topology_key(&self) -> Result<UnrootedTopologyKey, RustreeError> {
        unrooted_topology_key(self)
    }

    /// Build a stable 64-bit hash for the labeled unrooted topology.
    pub fn unrooted_topology_hash(&self) -> Result<u64, RustreeError> {
        self.unrooted_topology_key().map(|key| key.hash64())
    }

    /// Build an exact canonical key for the unlabeled unrooted topology.
    pub fn unrooted_shape_key(&self) -> Result<UnrootedShapeKey, RustreeError> {
        unrooted_shape_key(self)
    }

    /// Build a stable 64-bit hash for the unlabeled unrooted topology.
    pub fn unrooted_shape_hash(&self) -> Result<u64, RustreeError> {
        self.unrooted_shape_key().map(|key| key.hash64())
    }
}

/// Build an exact canonical key for the labeled unrooted topology.
pub fn unrooted_topology_key(tree: &FlatTree) -> Result<UnrootedTopologyKey, RustreeError> {
    validate_nonempty(tree)?;

    let (leaf_names, leaf_ids, is_leaf) = canonical_leaf_order(tree)?;
    let (adjacency, active) = contracted_adjacency(tree, &is_leaf)?;
    let root = first_active(&active).ok_or_else(|| {
        RustreeError::Tree("Failed to find an active node after topology contraction".to_string())
    })?;
    let (parent, order) = rooted_order(&adjacency, &active, root)?;

    let leaf_count = leaf_names.len();
    let word_count = bit_word_count(leaf_count);
    let mut descendant_sets = vec![vec![0_u64; word_count]; tree.nodes.len()];

    for &node_idx in order.iter().rev() {
        let mut bits = vec![0_u64; word_count];
        if let Some(leaf_id) = leaf_ids[node_idx] {
            set_bit(&mut bits, leaf_id);
        }
        for &neighbor in &adjacency[node_idx] {
            if !active[neighbor] || parent[node_idx] == Some(neighbor) {
                continue;
            }
            or_assign(&mut bits, &descendant_sets[neighbor]);
        }
        descendant_sets[node_idx] = bits;
    }

    let mut splits = Vec::new();
    for &node_idx in &order {
        if node_idx == root {
            continue;
        }
        if let Some(split) = normalized_split(&descendant_sets[node_idx], leaf_count) {
            splits.push(split);
        }
    }
    splits.sort();
    splits.dedup();

    Ok(UnrootedTopologyKey { leaf_names, splits })
}

/// Build a stable 64-bit hash for the labeled unrooted topology.
pub fn unrooted_topology_hash(tree: &FlatTree) -> Result<u64, RustreeError> {
    unrooted_topology_key(tree).map(|key| key.hash64())
}

/// Build an exact canonical key for the unlabeled unrooted topology.
pub fn unrooted_shape_key(tree: &FlatTree) -> Result<UnrootedShapeKey, RustreeError> {
    validate_nonempty(tree)?;

    let contractible = vec![true; tree.nodes.len()];
    let (adjacency, active) = contract_adjacency(build_adjacency(tree), contractible);

    let centers = tree_centers(&adjacency, &active)?;
    let shape = match centers.as_slice() {
        [center] => rooted_shape(&adjacency, *center, None),
        [left, right] => {
            let mut children = vec![
                rooted_shape(&adjacency, *left, Some(*right)),
                rooted_shape(&adjacency, *right, Some(*left)),
            ];
            children.sort();
            UnlabeledShape::Node(children)
        }
        _ => {
            return Err(RustreeError::Tree(format!(
                "Expected 1 or 2 centers in contracted tree, found {}",
                centers.len()
            )));
        }
    };

    Ok(UnrootedShapeKey { shape })
}

/// Build a stable 64-bit hash for the unlabeled unrooted topology.
pub fn unrooted_shape_hash(tree: &FlatTree) -> Result<u64, RustreeError> {
    unrooted_shape_key(tree).map(|key| key.hash64())
}

fn validate_nonempty(tree: &FlatTree) -> Result<(), RustreeError> {
    if tree.nodes.is_empty() {
        return Err(RustreeError::Tree("Cannot hash an empty tree".to_string()));
    }
    if tree.root >= tree.nodes.len() {
        return Err(RustreeError::Tree(format!(
            "Tree root index {} is out of bounds for {} nodes",
            tree.root,
            tree.nodes.len()
        )));
    }
    Ok(())
}

fn canonical_leaf_order(tree: &FlatTree) -> Result<LeafOrder, RustreeError> {
    let mut leaves = Vec::new();
    let mut is_leaf = vec![false; tree.nodes.len()];

    for (idx, node) in tree.nodes.iter().enumerate() {
        if node.left_child.is_none() && node.right_child.is_none() {
            is_leaf[idx] = true;
            leaves.push((node.name.clone(), idx));
        }
    }

    leaves.sort_by(|left, right| left.0.cmp(&right.0));

    for pair in leaves.windows(2) {
        if pair[0].0 == pair[1].0 {
            return Err(RustreeError::Validation(format!(
                "Duplicate leaf label '{}' is not supported for labeled topology hashing",
                pair[0].0
            )));
        }
    }

    let mut leaf_names = Vec::with_capacity(leaves.len());
    let mut leaf_ids = vec![None; tree.nodes.len()];
    for (leaf_id, (name, node_idx)) in leaves.into_iter().enumerate() {
        leaf_names.push(name);
        leaf_ids[node_idx] = Some(leaf_id);
    }

    Ok((leaf_names, leaf_ids, is_leaf))
}

fn contracted_adjacency(
    tree: &FlatTree,
    is_leaf: &[bool],
) -> Result<(Vec<Vec<usize>>, Vec<bool>), RustreeError> {
    let contractible: Vec<bool> = is_leaf.iter().map(|is_leaf| !is_leaf).collect();
    let (adjacency, active) = contract_adjacency(build_adjacency(tree), contractible);
    if first_active(&active).is_none() {
        return Err(RustreeError::Tree(
            "Topology contraction removed every node from the tree".to_string(),
        ));
    }
    Ok((adjacency, active))
}

fn build_adjacency(tree: &FlatTree) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); tree.nodes.len()];
    for (idx, node) in tree.nodes.iter().enumerate() {
        if let Some(left) = node.left_child {
            adjacency[idx].push(left);
            adjacency[left].push(idx);
        }
        if let Some(right) = node.right_child {
            adjacency[idx].push(right);
            adjacency[right].push(idx);
        }
    }
    adjacency
}

fn contract_adjacency(
    mut adjacency: Vec<Vec<usize>>,
    contractible: Vec<bool>,
) -> (Vec<Vec<usize>>, Vec<bool>) {
    let mut active = vec![true; adjacency.len()];
    let mut queue = VecDeque::new();

    for node_idx in 0..adjacency.len() {
        if contractible[node_idx] && adjacency[node_idx].len() == 2 {
            queue.push_back(node_idx);
        }
    }

    while let Some(node_idx) = queue.pop_front() {
        if !active[node_idx] || !contractible[node_idx] || adjacency[node_idx].len() != 2 {
            continue;
        }

        let left = adjacency[node_idx][0];
        let right = adjacency[node_idx][1];
        if !active[left] || !active[right] || left == right {
            continue;
        }

        remove_edge(&mut adjacency, left, node_idx);
        remove_edge(&mut adjacency, right, node_idx);
        adjacency[node_idx].clear();
        active[node_idx] = false;
        add_edge(&mut adjacency, left, right);

        if contractible[left] && adjacency[left].len() == 2 {
            queue.push_back(left);
        }
        if contractible[right] && adjacency[right].len() == 2 {
            queue.push_back(right);
        }
    }

    (adjacency, active)
}

fn remove_edge(adjacency: &mut [Vec<usize>], from: usize, to: usize) {
    if let Some(pos) = adjacency[from].iter().position(|&neighbor| neighbor == to) {
        adjacency[from].swap_remove(pos);
    }
}

fn add_edge(adjacency: &mut [Vec<usize>], left: usize, right: usize) {
    if !adjacency[left].contains(&right) {
        adjacency[left].push(right);
    }
    if !adjacency[right].contains(&left) {
        adjacency[right].push(left);
    }
}

fn first_active(active: &[bool]) -> Option<usize> {
    active.iter().position(|is_active| *is_active)
}

fn rooted_order(
    adjacency: &[Vec<usize>],
    active: &[bool],
    root: usize,
) -> Result<RootedTraversal, RustreeError> {
    let active_count = active.iter().filter(|is_active| **is_active).count();
    let mut parent = vec![None; adjacency.len()];
    let mut order = Vec::with_capacity(active_count);
    let mut stack = vec![root];
    parent[root] = Some(root);

    while let Some(node_idx) = stack.pop() {
        order.push(node_idx);
        for &neighbor in &adjacency[node_idx] {
            if !active[neighbor] || Some(neighbor) == parent[node_idx] {
                continue;
            }
            if parent[neighbor].is_some() {
                continue;
            }
            parent[neighbor] = Some(node_idx);
            stack.push(neighbor);
        }
    }

    parent[root] = None;

    if order.len() != active_count {
        return Err(RustreeError::Tree(format!(
            "Expected a connected tree after contraction, but visited {} of {} active nodes",
            order.len(),
            active_count
        )));
    }

    Ok((parent, order))
}

fn bit_word_count(bit_count: usize) -> usize {
    (bit_count.saturating_add(63)) / 64
}

fn set_bit(words: &mut [u64], bit: usize) {
    if words.is_empty() {
        return;
    }
    words[bit / 64] |= 1_u64 << (bit % 64);
}

fn or_assign(dst: &mut [u64], src: &[u64]) {
    for (dst_word, src_word) in dst.iter_mut().zip(src.iter()) {
        *dst_word |= *src_word;
    }
}

fn popcount(words: &[u64]) -> usize {
    words.iter().map(|word| word.count_ones() as usize).sum()
}

fn complement_bitset(words: &[u64], leaf_count: usize) -> Vec<u64> {
    let mut complement: Vec<u64> = words.iter().map(|word| !word).collect();
    if let Some(last_word) = complement.last_mut() {
        let used_bits = leaf_count % 64;
        if used_bits != 0 {
            *last_word &= (1_u64 << used_bits) - 1;
        }
    }
    complement
}

fn normalized_split(words: &[u64], leaf_count: usize) -> Option<Vec<u64>> {
    let size = popcount(words);
    if size <= 1 || size + 1 >= leaf_count {
        return None;
    }

    let complement = complement_bitset(words, leaf_count);
    let complement_size = leaf_count - size;

    let normalized = match size.cmp(&complement_size) {
        Ordering::Less => words.to_vec(),
        Ordering::Greater => complement,
        Ordering::Equal => {
            if subset_lex_cmp(words, &complement) == Ordering::Less {
                words.to_vec()
            } else {
                complement
            }
        }
    };

    Some(normalized)
}

fn subset_lex_cmp(left: &[u64], right: &[u64]) -> Ordering {
    for (left_word, right_word) in left.iter().zip(right.iter()) {
        let diff = left_word ^ right_word;
        if diff == 0 {
            continue;
        }
        let bit_idx = diff.trailing_zeros() as usize;
        let mask = 1_u64 << bit_idx;
        return if (left_word & mask) != 0 {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    Ordering::Equal
}

fn tree_centers(adjacency: &[Vec<usize>], active: &[bool]) -> Result<Vec<usize>, RustreeError> {
    let active_count = active.iter().filter(|is_active| **is_active).count();
    if active_count == 0 {
        return Err(RustreeError::Tree(
            "Cannot compute the center of an empty contracted tree".to_string(),
        ));
    }
    if active_count <= 2 {
        return Ok(active
            .iter()
            .enumerate()
            .filter_map(|(idx, is_active)| is_active.then_some(idx))
            .collect());
    }

    let mut degree: Vec<usize> = adjacency.iter().map(Vec::len).collect();
    let mut peeled = active
        .iter()
        .map(|is_active| !*is_active)
        .collect::<Vec<_>>();
    let mut leaves = VecDeque::new();

    for (idx, is_active) in active.iter().enumerate() {
        if *is_active && degree[idx] <= 1 {
            leaves.push_back(idx);
        }
    }

    let mut remaining = active_count;
    while remaining > 2 {
        let layer_len = leaves.len();
        if layer_len == 0 {
            return Err(RustreeError::Tree(
                "Failed to peel leaves while searching for tree centers".to_string(),
            ));
        }
        remaining -= layer_len;
        for _ in 0..layer_len {
            let leaf_idx = leaves
                .pop_front()
                .expect("queue length checked before pop_front");
            if peeled[leaf_idx] {
                continue;
            }
            peeled[leaf_idx] = true;
            for &neighbor in &adjacency[leaf_idx] {
                if peeled[neighbor] {
                    continue;
                }
                degree[neighbor] = degree[neighbor].saturating_sub(1);
                if degree[neighbor] == 1 {
                    leaves.push_back(neighbor);
                }
            }
        }
    }

    Ok(active
        .iter()
        .enumerate()
        .filter_map(|(idx, is_active)| (*is_active && !peeled[idx]).then_some(idx))
        .collect())
}

fn rooted_shape(
    adjacency: &[Vec<usize>],
    node_idx: usize,
    parent: Option<usize>,
) -> UnlabeledShape {
    let mut children = Vec::new();
    for &neighbor in &adjacency[node_idx] {
        if Some(neighbor) == parent {
            continue;
        }
        children.push(rooted_shape(adjacency, neighbor, Some(node_idx)));
    }

    if children.is_empty() {
        UnlabeledShape::Leaf
    } else {
        children.sort();
        UnlabeledShape::Node(children)
    }
}

fn stable_hash64<T: Hash>(value: &T) -> u64 {
    let mut hasher = Fnv1a64::new();
    value.hash(&mut hasher);
    hasher.finish()
}

struct Fnv1a64 {
    state: u64,
}

impl Fnv1a64 {
    const OFFSET_BASIS: u64 = 0xcbf29ce484222325;
    const PRIME: u64 = 0x100000001b3;

    fn new() -> Self {
        Self {
            state: Self::OFFSET_BASIS,
        }
    }

    fn write_bytes(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.state ^= u64::from(*byte);
            self.state = self.state.wrapping_mul(Self::PRIME);
        }
    }
}

impl Hasher for Fnv1a64 {
    fn finish(&self) -> u64 {
        self.state
    }

    fn write(&mut self, bytes: &[u8]) {
        self.write_bytes(bytes);
    }

    fn write_u8(&mut self, value: u8) {
        self.write_bytes(&[value]);
    }

    fn write_u16(&mut self, value: u16) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_u32(&mut self, value: u32) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_u64(&mut self, value: u64) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_u128(&mut self, value: u128) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_usize(&mut self, value: usize) {
        self.write_bytes(&(value as u64).to_le_bytes());
    }

    fn write_i8(&mut self, value: i8) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_i16(&mut self, value: i16) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_i32(&mut self, value: i32) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_i64(&mut self, value: i64) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_i128(&mut self, value: i128) {
        self.write_bytes(&value.to_le_bytes());
    }

    fn write_isize(&mut self, value: isize) {
        self.write_bytes(&(value as i64).to_le_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::{unrooted_shape_key, unrooted_topology_key};
    use crate::newick::parse_newick;

    fn flat_tree(newick: &str) -> crate::node::FlatTree {
        let mut nodes = parse_newick(newick).expect("failed to parse test tree");
        nodes.pop().expect("missing root").to_flat_tree()
    }

    #[test]
    fn labeled_key_is_root_invariant() {
        let tree1 = flat_tree("((A:1,B:1):1,(C:1,D:1):1):0;");
        let tree2 = flat_tree("(A:1,(B:1,(C:1,D:1):1):1):0;");

        let key1 = unrooted_topology_key(&tree1).unwrap();
        let key2 = unrooted_topology_key(&tree2).unwrap();

        assert_eq!(key1, key2);
        assert_eq!(key1.hash64(), key2.hash64());
    }

    #[test]
    fn labeled_key_tracks_leaf_labels() {
        let tree1 = flat_tree("((A:1,B:1):1,(C:1,D:1):1):0;");
        let tree2 = flat_tree("((A:1,C:1):1,(B:1,D:1):1):0;");

        let key1 = unrooted_topology_key(&tree1).unwrap();
        let key2 = unrooted_topology_key(&tree2).unwrap();

        assert_ne!(key1, key2);
        assert_ne!(key1.hash64(), key2.hash64());
    }

    #[test]
    fn labeled_key_rejects_duplicate_leaf_labels() {
        let tree = flat_tree("((A:1,A:1):1,(C:1,D:1):1):0;");
        let err = unrooted_topology_key(&tree).unwrap_err();
        assert!(err.to_string().contains("Duplicate leaf label"));
    }

    #[test]
    fn unlabeled_key_ignores_rooting_and_leaf_names() {
        let tree1 = flat_tree("((A:1,B:1):1,(C:1,D:1):1):0;");
        let tree2 = flat_tree("(w:1,(x:1,(y:1,z:1):1):1):0;");

        let key1 = unrooted_shape_key(&tree1).unwrap();
        let key2 = unrooted_shape_key(&tree2).unwrap();

        assert_eq!(key1, key2);
        assert_eq!(key1.hash64(), key2.hash64());
    }

    #[test]
    fn unlabeled_key_distinguishes_different_shapes() {
        let symmetric = flat_tree("((A:1,B:1):1,((C:1,D:1):1,(E:1,F:1):1):1):0;");
        let caterpillar = flat_tree("((((A:1,B:1):1,C:1):1,D:1):1,(E:1,F:1):1):0;");

        let key1 = unrooted_shape_key(&symmetric).unwrap();
        let key2 = unrooted_shape_key(&caterpillar).unwrap();

        assert_ne!(key1, key2);
        assert_ne!(key1.hash64(), key2.hash64());
    }
}
