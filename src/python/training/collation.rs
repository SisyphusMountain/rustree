//! Batch-level operations: collation, GCN normalization, and inference batch construction.
//!
//! Contains CollatedTask struct, gcn_norm_internal, compute_varlen_metadata,
//! collate_task_tensors, and the build_otf_batch / compute_gcn_norm / build_inference_batch pyfunctions.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

use crate::bd::simulate_bd_tree_bwd;
use crate::dtl::simulate_dtl_per_species;
use crate::node::{remap_gene_tree_indices, Event, RecTree};
use crate::sampling::extract_extant_subtree;
use crate::sampling::extract_induced_subtree;

use crate::python::validate_dtl_rates;

use super::extraction::{extract_base_sample_internal, parse_base_sample_dict, RustBaseSample};
use super::tensors::{build_task_tensors_internal, TaskTensors};

/// GCN normalization: add self-loops, compute D^{-1/2} A D^{-1/2} edge weights.
///
/// Matches PyG's `gcn_norm(edge_index, None, num_nodes, improved=False,
/// add_self_loops=True, flow="source_to_target")`.
///
/// Returns (new_src, new_dst, edge_weights) with self-loops included.
pub(super) fn gcn_norm_internal(
    src: &[i32],
    dst: &[i32],
    num_nodes: usize,
) -> (Vec<i32>, Vec<i32>, Vec<f32>) {
    let n_orig = src.len();
    let n_total = n_orig + num_nodes; // original edges + self-loops

    let mut new_src = Vec::with_capacity(n_total);
    let mut new_dst = Vec::with_capacity(n_total);
    new_src.extend_from_slice(src);
    new_dst.extend_from_slice(dst);
    for i in 0..num_nodes {
        new_src.push(i as i32);
        new_dst.push(i as i32);
    }

    // Degree = number of edges targeting each node (flow="source_to_target")
    let mut deg = vec![0u32; num_nodes];
    for &d in &new_dst {
        deg[d as usize] += 1;
    }

    // D^{-1/2}
    let inv_sqrt: Vec<f32> = deg
        .iter()
        .map(|&d| {
            if d == 0 {
                0.0f32
            } else {
                (d as f32).powf(-0.5)
            }
        })
        .collect();

    // Edge weights: inv_sqrt[src] * inv_sqrt[dst]
    let weights: Vec<f32> = new_src
        .iter()
        .zip(new_dst.iter())
        .map(|(&s, &d)| inv_sqrt[s as usize] * inv_sqrt[d as usize])
        .collect();

    (new_src, new_dst, weights)
}

/// Compute varlen attention metadata from sizes array.
/// Returns (prefix_sums_with_leading_zero, max_size).
fn compute_varlen_metadata(sizes: &[i32]) -> (Vec<i32>, i32) {
    let mut cu = Vec::with_capacity(sizes.len() + 1);
    cu.push(0);
    let mut cumsum: i32 = 0;
    let mut max_val: i32 = 0;
    for &s in sizes {
        cumsum += s;
        cu.push(cumsum);
        if s > max_val {
            max_val = s;
        }
    }
    (cu, max_val)
}

#[derive(Clone, Copy, Default)]
struct TreeMetadataRequest {
    topological: bool,
    branch_length: bool,
}

impl TreeMetadataRequest {
    fn from_metrics(metrics: Option<Vec<String>>) -> Result<Self, String> {
        let mut request = Self::default();
        for metric in metrics.unwrap_or_default() {
            match metric.as_str() {
                "topological" => request.topological = true,
                "branch_length" => request.branch_length = true,
                other => {
                    return Err(format!(
                        "tree metadata metric must be 'topological' or 'branch_length', got '{}'",
                        other
                    ));
                }
            }
        }
        Ok(request)
    }

    fn any(self) -> bool {
        self.topological || self.branch_length
    }
}

struct SegmentTreeMetadata {
    token_depth: Vec<f32>,
    block_starts: Vec<i64>,
    block_lengths: Vec<i64>,
    block_lca_depth: Vec<f32>,
    local_lca_depth_padded: Vec<f32>,
    max_local_block: usize,
    num_blocks: usize,
}

struct PackedTreeMetadata {
    token_depth: Vec<f32>,
    cu_tokens: Vec<i64>,
    block_starts: Vec<i64>,
    block_lengths: Vec<i64>,
    block_ptr: Vec<i64>,
    block_segment: Vec<i64>,
    block_local_index: Vec<i64>,
    block_lca_depth_padded: Vec<f32>,
    block_lca_shape: (usize, usize, usize),
    local_lca_depth_padded: Vec<f32>,
    local_lca_shape: (usize, usize, usize),
}

#[derive(Default)]
struct TreeMetadataBundle {
    topological: Option<PackedTreeMetadata>,
    branch_length: Option<PackedTreeMetadata>,
}

fn next_power_of_two_or_zero(value: usize) -> usize {
    if value == 0 {
        0
    } else {
        value.next_power_of_two()
    }
}

fn build_segment_tree_metadata_fast(
    parents: &[i32],
    edge_lengths: Option<&[f32]>,
    max_block_size: usize,
) -> Result<SegmentTreeMetadata, String> {
    if parents.is_empty() {
        return Err("tree metadata segment cannot be empty".to_string());
    }
    if max_block_size == 0 {
        return Err("tree metadata max_block_size must be positive".to_string());
    }
    if let Some(lengths) = edge_lengths {
        if lengths.len() != parents.len() {
            return Err("edge_lengths must have the same length as parents".to_string());
        }
    }

    let n = parents.len();
    let mut root: Option<usize> = None;
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut topo_depth = vec![0usize; n];
    let mut token_depth = vec![0.0f32; n];

    for node in 0..n {
        let p = parents[node];
        if p < 0 {
            if root.replace(node).is_some() {
                return Err("parents must describe exactly one rooted tree".to_string());
            }
            token_depth[node] = edge_lengths.map(|lengths| lengths[node]).unwrap_or(0.0);
            continue;
        }
        let parent = p as usize;
        if parent >= n {
            return Err(format!("parent index {} out of bounds for {} nodes", p, n));
        }
        if parent >= node {
            return Err(
                "fast tree metadata builder expects preorder parents before children".to_string(),
            );
        }
        children[parent].push(node);
        topo_depth[node] = topo_depth[parent] + 1;
        token_depth[node] = if let Some(lengths) = edge_lengths {
            token_depth[parent] + lengths[node]
        } else {
            topo_depth[node] as f32
        };
    }

    let root = root.ok_or_else(|| "parents must describe exactly one rooted tree".to_string())?;
    let mut subtree_size = vec![1usize; n];
    for node in (0..n).rev() {
        let p = parents[node];
        if p >= 0 {
            subtree_size[p as usize] += subtree_size[node];
        }
    }

    let max_log = (usize::BITS - n.leading_zeros()) as usize + 1;
    let mut up = vec![vec![usize::MAX; n]; max_log.max(1)];
    for node in 0..n {
        if parents[node] >= 0 {
            up[0][node] = parents[node] as usize;
        }
    }
    for level in 1..up.len() {
        for node in 0..n {
            let prev = up[level - 1][node];
            if prev != usize::MAX {
                up[level][node] = up[level - 1][prev];
            }
        }
    }

    fn lca(mut a: usize, mut b: usize, topo_depth: &[usize], up: &[Vec<usize>]) -> usize {
        if topo_depth[a] < topo_depth[b] {
            std::mem::swap(&mut a, &mut b);
        }
        let diff = topo_depth[a] - topo_depth[b];
        for level in 0..up.len() {
            if ((diff >> level) & 1) != 0 {
                a = up[level][a];
            }
        }
        if a == b {
            return a;
        }
        for level in (0..up.len()).rev() {
            let ua = up[level][a];
            let ub = up[level][b];
            if ua != ub {
                a = ua;
                b = ub;
            }
        }
        up[0][a]
    }

    let mut block_roots: Vec<usize> = Vec::new();
    let mut block_starts: Vec<i64> = Vec::new();
    let mut block_lengths: Vec<i64> = Vec::new();

    fn partition(
        node: usize,
        children: &[Vec<usize>],
        subtree_size: &[usize],
        max_block_size: usize,
        block_roots: &mut Vec<usize>,
        block_starts: &mut Vec<i64>,
        block_lengths: &mut Vec<i64>,
    ) {
        if subtree_size[node] <= max_block_size {
            block_roots.push(node);
            block_starts.push(node as i64);
            block_lengths.push(subtree_size[node] as i64);
            return;
        }
        block_roots.push(node);
        block_starts.push(node as i64);
        block_lengths.push(1);
        for &child in &children[node] {
            partition(
                child,
                children,
                subtree_size,
                max_block_size,
                block_roots,
                block_starts,
                block_lengths,
            );
        }
    }
    partition(
        root,
        &children,
        &subtree_size,
        max_block_size,
        &mut block_roots,
        &mut block_starts,
        &mut block_lengths,
    );

    let num_blocks = block_roots.len();
    let mut block_lca_depth = vec![0.0f32; num_blocks * num_blocks];
    for qb in 0..num_blocks {
        for kb in 0..num_blocks {
            let ancestor = lca(block_roots[qb], block_roots[kb], &topo_depth, &up);
            block_lca_depth[qb * num_blocks + kb] = token_depth[ancestor];
        }
    }

    let max_block_len = block_lengths
        .iter()
        .map(|&length| length as usize)
        .max()
        .unwrap_or(0);
    let max_local_block = next_power_of_two_or_zero(max_block_len);
    let mut local_lca_depth_padded = vec![0.0f32; num_blocks * max_local_block * max_local_block];
    if max_local_block > 0 {
        for block in 0..num_blocks {
            let start = block_starts[block] as usize;
            let len = block_lengths[block] as usize;
            for qi in 0..len {
                let q_node = start + qi;
                for ki in 0..len {
                    let k_node = start + ki;
                    let ancestor = lca(q_node, k_node, &topo_depth, &up);
                    let offset = (block * max_local_block + qi) * max_local_block + ki;
                    local_lca_depth_padded[offset] = token_depth[ancestor];
                }
            }
        }
    }

    Ok(SegmentTreeMetadata {
        token_depth,
        block_starts,
        block_lengths,
        block_lca_depth,
        local_lca_depth_padded,
        max_local_block,
        num_blocks,
    })
}

fn build_packed_tree_metadata_fast(
    parents: &[i32],
    cu_tokens: &[i32],
    edge_lengths: Option<&[f32]>,
    max_block_size: usize,
) -> Result<PackedTreeMetadata, String> {
    if cu_tokens.len() < 2 {
        return Err("cu_tokens must contain at least one segment".to_string());
    }
    if cu_tokens.first().copied() != Some(0)
        || cu_tokens.last().copied() != Some(parents.len() as i32)
    {
        return Err("cu_tokens must start at 0 and end at len(parents)".to_string());
    }
    if let Some(lengths) = edge_lengths {
        if lengths.len() != parents.len() {
            return Err("edge_lengths must have the same length as parents".to_string());
        }
    }

    let mut segments = Vec::with_capacity(cu_tokens.len() - 1);
    for window in cu_tokens.windows(2) {
        let start = window[0] as usize;
        let end = window[1] as usize;
        let local_parents: Vec<i32> = parents[start..end]
            .iter()
            .map(|&p| {
                if p >= start as i32 && p < end as i32 {
                    p - start as i32
                } else {
                    -1
                }
            })
            .collect();
        let local_lengths = edge_lengths.map(|lengths| &lengths[start..end]);
        segments.push(build_segment_tree_metadata_fast(
            &local_parents,
            local_lengths,
            max_block_size,
        )?);
    }

    let num_segments = segments.len();
    let total_tokens = parents.len();
    let total_blocks: usize = segments.iter().map(|meta| meta.num_blocks).sum();
    let max_blocks = segments
        .iter()
        .map(|meta| meta.num_blocks)
        .max()
        .unwrap_or(0);
    let max_local_block = segments
        .iter()
        .map(|meta| meta.max_local_block)
        .max()
        .unwrap_or(0);

    let mut token_depth = Vec::with_capacity(total_tokens);
    let mut block_starts = Vec::with_capacity(total_blocks);
    let mut block_lengths = Vec::with_capacity(total_blocks);
    let mut block_ptr = Vec::with_capacity(num_segments + 1);
    let mut block_segment = Vec::with_capacity(total_blocks);
    let mut block_local_index = Vec::with_capacity(total_blocks);
    let mut block_lca_depth_padded = vec![0.0f32; num_segments * max_blocks * max_blocks];
    let mut local_lca_depth_padded = vec![0.0f32; total_blocks * max_local_block * max_local_block];

    let mut token_offset = 0i64;
    let mut block_offset = 0usize;
    block_ptr.push(0);
    for (segment_idx, meta) in segments.iter().enumerate() {
        token_depth.extend_from_slice(&meta.token_depth);
        for block in 0..meta.num_blocks {
            block_starts.push(meta.block_starts[block] + token_offset);
            block_lengths.push(meta.block_lengths[block]);
            block_segment.push(segment_idx as i64);
            block_local_index.push(block as i64);
        }
        for qb in 0..meta.num_blocks {
            for kb in 0..meta.num_blocks {
                let dst = (segment_idx * max_blocks + qb) * max_blocks + kb;
                let src = qb * meta.num_blocks + kb;
                block_lca_depth_padded[dst] = meta.block_lca_depth[src];
            }
        }
        if max_local_block > 0 {
            for block in 0..meta.num_blocks {
                for qi in 0..meta.max_local_block {
                    for ki in 0..meta.max_local_block {
                        let dst_block = block_offset + block;
                        let dst = (dst_block * max_local_block + qi) * max_local_block + ki;
                        let src = (block * meta.max_local_block + qi) * meta.max_local_block + ki;
                        local_lca_depth_padded[dst] = meta.local_lca_depth_padded[src];
                    }
                }
            }
        }
        token_offset += meta.token_depth.len() as i64;
        block_offset += meta.num_blocks;
        block_ptr.push(block_offset as i64);
    }

    Ok(PackedTreeMetadata {
        token_depth,
        cu_tokens: cu_tokens.iter().map(|&v| v as i64).collect(),
        block_starts,
        block_lengths,
        block_ptr,
        block_segment,
        block_local_index,
        block_lca_depth_padded,
        block_lca_shape: (num_segments, max_blocks, max_blocks),
        local_lca_depth_padded,
        local_lca_shape: (total_blocks, max_local_block, max_local_block),
    })
}

fn repeat_segment_tree_metadata(
    meta: &SegmentTreeMetadata,
    segment_len: usize,
    num_segments: usize,
) -> PackedTreeMetadata {
    let num_blocks = meta.num_blocks;
    let max_local_block = meta.max_local_block;
    let total_blocks = num_segments * num_blocks;
    let mut token_depth = Vec::with_capacity(num_segments * segment_len);
    let mut cu_tokens = Vec::with_capacity(num_segments + 1);
    let mut block_starts = Vec::with_capacity(total_blocks);
    let mut block_lengths = Vec::with_capacity(total_blocks);
    let mut block_ptr = Vec::with_capacity(num_segments + 1);
    let mut block_segment = Vec::with_capacity(total_blocks);
    let mut block_local_index = Vec::with_capacity(total_blocks);
    let mut block_lca_depth_padded = vec![0.0f32; num_segments * num_blocks * num_blocks];
    let mut local_lca_depth_padded = vec![0.0f32; total_blocks * max_local_block * max_local_block];

    cu_tokens.push(0);
    block_ptr.push(0);
    for segment in 0..num_segments {
        let token_offset = (segment * segment_len) as i64;
        token_depth.extend_from_slice(&meta.token_depth);
        cu_tokens.push(((segment + 1) * segment_len) as i64);
        for block in 0..num_blocks {
            block_starts.push(meta.block_starts[block] + token_offset);
            block_lengths.push(meta.block_lengths[block]);
            block_segment.push(segment as i64);
            block_local_index.push(block as i64);
        }
        block_ptr.push(((segment + 1) * num_blocks) as i64);

        let block_lca_dst = segment * num_blocks * num_blocks;
        block_lca_depth_padded[block_lca_dst..block_lca_dst + meta.block_lca_depth.len()]
            .copy_from_slice(&meta.block_lca_depth);

        if max_local_block > 0 {
            let local_len = meta.local_lca_depth_padded.len();
            let local_dst = segment * local_len;
            local_lca_depth_padded[local_dst..local_dst + local_len]
                .copy_from_slice(&meta.local_lca_depth_padded);
        }
    }

    PackedTreeMetadata {
        token_depth,
        cu_tokens,
        block_starts,
        block_lengths,
        block_ptr,
        block_segment,
        block_local_index,
        block_lca_depth_padded,
        block_lca_shape: (num_segments, num_blocks, num_blocks),
        local_lca_depth_padded,
        local_lca_shape: (total_blocks, max_local_block, max_local_block),
    }
}

fn build_repeated_tree_metadata_fast(
    single_parents: &[i32],
    single_edge_lengths: Option<&[f32]>,
    num_segments: usize,
    max_block_size: usize,
) -> Result<PackedTreeMetadata, String> {
    let meta =
        build_segment_tree_metadata_fast(single_parents, single_edge_lengths, max_block_size)?;
    Ok(repeat_segment_tree_metadata(
        &meta,
        single_parents.len(),
        num_segments,
    ))
}

fn build_tree_metadata_bundle(
    parents: &[i32],
    cu_tokens: &[i32],
    edge_lengths: &[f32],
    request: TreeMetadataRequest,
) -> Result<TreeMetadataBundle, String> {
    if !request.any() {
        return Ok(TreeMetadataBundle::default());
    }
    Ok(TreeMetadataBundle {
        topological: if request.topological {
            Some(build_packed_tree_metadata_fast(
                parents, cu_tokens, None, 64,
            )?)
        } else {
            None
        },
        branch_length: if request.branch_length {
            Some(build_packed_tree_metadata_fast(
                parents,
                cu_tokens,
                Some(edge_lengths),
                64,
            )?)
        } else {
            None
        },
    })
}

fn build_repeated_tree_metadata_bundle(
    single_parents: &[i32],
    single_edge_lengths: &[f32],
    num_segments: usize,
    request: TreeMetadataRequest,
) -> Result<TreeMetadataBundle, String> {
    if !request.any() {
        return Ok(TreeMetadataBundle::default());
    }
    Ok(TreeMetadataBundle {
        topological: if request.topological {
            Some(build_repeated_tree_metadata_fast(
                single_parents,
                None,
                num_segments,
                64,
            )?)
        } else {
            None
        },
        branch_length: if request.branch_length {
            Some(build_repeated_tree_metadata_fast(
                single_parents,
                Some(single_edge_lengths),
                num_segments,
                64,
            )?)
        } else {
            None
        },
    })
}

#[cfg(feature = "python")]
fn write_tree_metadata_arrays(
    py: Python<'_>,
    result: &Bound<'_, pyo3::types::PyDict>,
    prefix: &str,
    metric: &str,
    metadata: &PackedTreeMetadata,
) -> PyResult<()> {
    use numpy::ndarray::Array3;
    use numpy::PyArray1;
    use numpy::ToPyArray;

    let key = |field: &str| format!("{}_tree_metadata_{}_{}", prefix, metric, field);
    result.set_item(
        key("token_depth"),
        PyArray1::from_slice(py, &metadata.token_depth),
    )?;
    result.set_item(
        key("cu_tokens"),
        PyArray1::from_slice(py, &metadata.cu_tokens),
    )?;
    result.set_item(
        key("block_starts"),
        PyArray1::from_slice(py, &metadata.block_starts),
    )?;
    result.set_item(
        key("block_lengths"),
        PyArray1::from_slice(py, &metadata.block_lengths),
    )?;
    result.set_item(
        key("block_ptr"),
        PyArray1::from_slice(py, &metadata.block_ptr),
    )?;
    result.set_item(
        key("block_segment"),
        PyArray1::from_slice(py, &metadata.block_segment),
    )?;
    result.set_item(
        key("block_local_index"),
        PyArray1::from_slice(py, &metadata.block_local_index),
    )?;
    let block_lca = Array3::from_shape_vec(
        metadata.block_lca_shape,
        metadata.block_lca_depth_padded.clone(),
    )
    .map_err(|e| PyValueError::new_err(format!("block_lca_depth_padded: {}", e)))?;
    result.set_item(key("block_lca_depth_padded"), block_lca.to_pyarray(py))?;
    let local_lca = Array3::from_shape_vec(
        metadata.local_lca_shape,
        metadata.local_lca_depth_padded.clone(),
    )
    .map_err(|e| PyValueError::new_err(format!("local_lca_depth_padded: {}", e)))?;
    result.set_item(key("local_lca_depth_padded"), local_lca.to_pyarray(py))?;
    Ok(())
}

#[cfg(feature = "python")]
fn write_tree_metadata_bundle(
    py: Python<'_>,
    result: &Bound<'_, pyo3::types::PyDict>,
    prefix: &str,
    bundle: &TreeMetadataBundle,
) -> PyResult<()> {
    if let Some(ref metadata) = bundle.topological {
        write_tree_metadata_arrays(py, result, prefix, "topological", metadata)?;
    }
    if let Some(ref metadata) = bundle.branch_length {
        write_tree_metadata_arrays(py, result, prefix, "branch_length", metadata)?;
    }
    Ok(())
}

/// Collated tensors for one task across a full batch.
struct CollatedTask {
    gene_x: Vec<i32>,
    gene_y: Vec<i32>,
    event: Vec<i32>,
    event_true: Vec<i32>,
    frontier_mask: Vec<u8>,
    is_leaf: Vec<u8>,
    mask_label_node: Vec<u8>,
    // GCN-normalized gene edges (includes self-loops)
    gene_ei_src: Vec<i32>,
    gene_ei_dst: Vec<i32>,
    gene_ew: Vec<f32>,
    g_dir_src: Vec<i32>,
    g_dir_dst: Vec<i32>,
    gene_tree_parent: Vec<i32>,
    gene_tree_branch_length: Vec<f32>,
    true_root_index: Vec<i32>,
    g_sizes: Vec<i32>,
    g_ptr: Vec<i32>,
    g_batch: Vec<i32>,
    root_edges_ptr: Vec<i32>,
    // Varlen attention metadata
    cu_g: Vec<i32>,
    max_g: i32,
}

/// Collate per-sample task tensors into a single batch.
///
/// Also applies GCN normalization to gene edges and computes varlen attention metadata.
fn collate_task_tensors(tensors: &[TaskTensors]) -> CollatedTask {
    let mut gene_x = Vec::new();
    let mut gene_y = Vec::new();
    let mut event = Vec::new();
    let mut event_true = Vec::new();
    let mut frontier_mask = Vec::new();
    let mut is_leaf = Vec::new();
    let mut mask_label_node = Vec::new();
    let mut g_edge_src_raw = Vec::new();
    let mut g_edge_dst_raw = Vec::new();
    let mut g_dir_src = Vec::new();
    let mut g_dir_dst = Vec::new();
    let mut gene_tree_parent = Vec::new();
    let mut gene_tree_branch_length = Vec::new();
    let mut true_root_index = Vec::new();
    let mut g_sizes: Vec<i32> = Vec::new();
    let mut g_ptr: Vec<i32> = vec![0];
    let mut g_batch: Vec<i32> = Vec::new();
    let mut root_edges_ptr: Vec<i32> = vec![0];

    let mut cum_g: i32 = 0;
    let mut cum_re: i32 = 0;

    for (s_idx, t) in tensors.iter().enumerate() {
        let ng = t.x_gene.len() as i32;
        let ndir = t.g_dir_src.len() as i32;

        gene_x.extend_from_slice(&t.x_gene);
        gene_y.extend(t.g_true_sp.iter().map(|&v| v - 1));
        event.extend_from_slice(&t.event_input);
        event_true.extend_from_slice(&t.event_true);
        frontier_mask.extend_from_slice(&t.frontier_mask);
        is_leaf.extend_from_slice(&t.is_leaf);
        mask_label_node.extend_from_slice(&t.mask_label_node);

        for &s in &t.g_edge_src {
            g_edge_src_raw.push(s + cum_g);
        }
        for &d in &t.g_edge_dst {
            g_edge_dst_raw.push(d + cum_g);
        }
        for &s in &t.g_dir_src {
            g_dir_src.push(s + cum_g);
        }
        for &d in &t.g_dir_dst {
            g_dir_dst.push(d + cum_g);
        }
        for &p in &t.tree_parent {
            gene_tree_parent.push(if p >= 0 { p + cum_g } else { -1 });
        }
        gene_tree_branch_length.extend_from_slice(&t.tree_branch_length);

        true_root_index.push(t.root_edge_target);
        g_sizes.push(ng);
        for _ in 0..ng {
            g_batch.push(s_idx as i32);
        }
        cum_g += ng;
        g_ptr.push(cum_g);
        cum_re += ndir;
        root_edges_ptr.push(cum_re);
    }

    // Apply GCN normalization to gene edges
    let total_gene_nodes = cum_g as usize;
    let (gene_ei_src, gene_ei_dst, gene_ew) =
        gcn_norm_internal(&g_edge_src_raw, &g_edge_dst_raw, total_gene_nodes);

    // Compute varlen attention metadata
    let (cu_g, max_g) = compute_varlen_metadata(&g_sizes);

    CollatedTask {
        gene_x,
        gene_y,
        event,
        event_true,
        frontier_mask,
        is_leaf,
        mask_label_node,
        gene_ei_src,
        gene_ei_dst,
        gene_ew,
        g_dir_src,
        g_dir_dst,
        gene_tree_parent,
        gene_tree_branch_length,
        true_root_index,
        g_sizes,
        g_ptr,
        g_batch,
        root_edges_ptr,
        cu_g,
        max_g,
    }
}

/// Build a complete batched dict of numpy arrays for on-the-fly training.
///
/// Simulates one species tree and n_gene_trees gene trees, builds all three
/// task tensors (map / evt / root), collates them, and returns flat numpy arrays
/// ready to be passed to batch_to_heterodata() in Python.
///
/// Parameters
/// ----------
/// n_sp : int          Target number of extant species.
/// lambda_birth : float  Speciation rate (must be > mu_death).
/// mu_death : float    Extinction rate (must be >= 0).
/// sp_seed : int       Seed for species tree simulation.
/// n_gene_trees : int  Number of gene trees per batch.
/// gt_seeds : list\[int\]  Per-gene-tree simulation seeds (length n_gene_trees).
/// lambda_d/t/l : float  DTL rates for gene tree simulation.
/// enable_event : bool   Include evt_* keys in returned dict.
/// enable_root : bool    Include root_* keys in returned dict.
/// sample_order : str    "random" or "bottom_up".
/// map/evt/root_coloring_seeds : list\[int\]  Per-sample RNG seeds for tensor construction.
/// min_gene_leaves : int   Minimum extant gene leaves required.
/// max_gene_nodes : int    Maximum total gene nodes (0 = no limit).
/// max_retries : int       Retry attempts per gene tree slot on simulation failure.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (n_sp, lambda_birth, mu_death, sp_seed, n_gene_trees, gt_seeds, lambda_d, lambda_t, lambda_l, enable_event, enable_root, sample_order, map_coloring_seeds, evt_coloring_seeds, root_coloring_seeds, min_gene_leaves, max_gene_nodes, max_retries, species_tree_metadata_metrics=None, gene_tree_metadata_metrics=None))]
#[allow(clippy::too_many_arguments)]
pub fn build_otf_batch(
    py: Python,
    n_sp: usize,
    lambda_birth: f64,
    mu_death: f64,
    sp_seed: u64,
    n_gene_trees: usize,
    gt_seeds: Vec<u64>,
    lambda_d: f64,
    lambda_t: f64,
    lambda_l: f64,
    enable_event: bool,
    enable_root: bool,
    sample_order: &str,
    map_coloring_seeds: Vec<u64>,
    evt_coloring_seeds: Vec<u64>,
    root_coloring_seeds: Vec<u64>,
    min_gene_leaves: usize,
    max_gene_nodes: usize,
    max_retries: usize,
    species_tree_metadata_metrics: Option<Vec<String>>,
    gene_tree_metadata_metrics: Option<Vec<String>>,
) -> PyResult<PyObject> {
    use crate::node::TraversalOrder;
    use numpy::PyArray1;
    use pyo3::types::PyDict;

    // Validate
    if gt_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "gt_seeds length {} != n_gene_trees {}",
            gt_seeds.len(),
            n_gene_trees
        )));
    }
    if map_coloring_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "map_coloring_seeds length {} != n_gene_trees {}",
            map_coloring_seeds.len(),
            n_gene_trees
        )));
    }
    if enable_event && evt_coloring_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "evt_coloring_seeds length {} != n_gene_trees {}",
            evt_coloring_seeds.len(),
            n_gene_trees
        )));
    }
    if enable_root && root_coloring_seeds.len() != n_gene_trees {
        return Err(PyValueError::new_err(format!(
            "root_coloring_seeds length {} != n_gene_trees {}",
            root_coloring_seeds.len(),
            n_gene_trees
        )));
    }
    validate_dtl_rates(lambda_d, lambda_t, lambda_l)?;
    let species_tree_metadata_request =
        TreeMetadataRequest::from_metrics(species_tree_metadata_metrics)
            .map_err(PyValueError::new_err)?;
    let gene_tree_metadata_request = TreeMetadataRequest::from_metrics(gene_tree_metadata_metrics)
        .map_err(PyValueError::new_err)?;

    // ---- GIL-free section ------------------------------------------------
    // All simulation, tensor building, and collation is pure Rust with no
    // Python object interaction.  Releasing the GIL here lets other Python
    // threads (prefetcher producers, GPU transfer thread) run concurrently.
    let sample_order_owned = sample_order.to_owned();
    let result_or_err = py.allow_threads(move || -> Result<_, String> {
        let sample_order: &str = &sample_order_owned;

        // 1. Simulate species tree
        let mut sp_rng = StdRng::seed_from_u64(sp_seed);
        let (mut sp_tree_raw, _) = simulate_bd_tree_bwd(n_sp, lambda_birth, mu_death, &mut sp_rng)?;
        sp_tree_raw.assign_depths();
        let sp_tree_raw_arc = Arc::new(sp_tree_raw);
        let (extant_sp, _) = extract_extant_subtree(&sp_tree_raw_arc).ok_or_else(|| {
            "build_otf_batch: failed to extract extant species subtree".to_string()
        })?;
        let extant_sp_arc = Arc::new(extant_sp);
        let n_sp_nodes = extant_sp_arc.nodes.len();

        // 2. Simulate and filter gene trees (one per slot, individual seeds)
        //    Parallelised with rayon: each slot is independent (own seed, own RNG).
        let species_identity: Vec<Option<usize>> = (0..n_sp_nodes).map(Some).collect();

        let valid_gene_trees: Vec<RecTree> = (0..n_gene_trees)
            .into_par_iter()
            .map(|i| {
                for attempt in 0..=(max_retries as u64) {
                    let gt_seed = gt_seeds[i].wrapping_add(attempt);
                    let mut gt_rng = StdRng::seed_from_u64(gt_seed);

                    let (mut rec_tree, events) = match simulate_dtl_per_species(
                        &extant_sp_arc,
                        extant_sp_arc.root,
                        lambda_d,
                        lambda_t,
                        lambda_l,
                        None,
                        None,
                        false,
                        &mut gt_rng,
                    ) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };
                    rec_tree.species_tree = Arc::clone(&extant_sp_arc);
                    rec_tree.dtl_events = Some(events);

                    // Extract extant gene leaves
                    let extant_gene_indices: std::collections::HashSet<usize> = rec_tree
                        .gene_tree
                        .nodes
                        .iter()
                        .enumerate()
                        .filter(|(idx, n)| {
                            n.left_child.is_none()
                                && n.right_child.is_none()
                                && rec_tree.event_mapping[*idx] == Event::Leaf
                        })
                        .map(|(idx, _)| idx)
                        .collect();

                    if extant_gene_indices.is_empty() {
                        continue;
                    }
                    if extant_gene_indices.len() < min_gene_leaves {
                        continue;
                    }

                    let (extant_gene_tree, gene_old_to_new) =
                        match extract_induced_subtree(&rec_tree.gene_tree, &extant_gene_indices) {
                            Some(r) => r,
                            None => continue,
                        };

                    if max_gene_nodes > 0 && extant_gene_tree.nodes.len() > max_gene_nodes {
                        continue;
                    }

                    let (new_node_mapping, new_event_mapping) = match remap_gene_tree_indices(
                        &extant_gene_tree,
                        &gene_old_to_new,
                        &rec_tree.node_mapping,
                        &rec_tree.event_mapping,
                        &species_identity,
                    ) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };

                    return Ok(RecTree::new(
                        Arc::clone(&extant_sp_arc),
                        extant_gene_tree,
                        new_node_mapping,
                        new_event_mapping,
                    ));
                }
                Err(format!(
                    "build_otf_batch: gene tree slot {} failed after {} retries",
                    i, max_retries
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // 3. Species tree tensors (computed once, replicated B times during collation)
        let sp_preorder: Vec<usize> = extant_sp_arc
            .iter_indices(TraversalOrder::PreOrder)
            .collect();
        let species_names: Vec<String> = sp_preorder
            .iter()
            .map(|&i| extant_sp_arc.nodes[i].name.clone())
            .collect();

        let mut sp_name_to_id: HashMap<String, i32> = HashMap::new();
        let mut sp_name_to_idx_map: HashMap<String, usize> = HashMap::new();
        for (i, name) in species_names.iter().enumerate() {
            sp_name_to_id.insert(name.clone(), (i + 1) as i32);
            sp_name_to_idx_map.insert(name.clone(), i);
        }

        let x_sp_single: Vec<i32> = species_names
            .iter()
            .map(|name| sp_name_to_id[name])
            .collect();

        let mut sp_tree_parent_s: Vec<i32> = Vec::with_capacity(n_sp_nodes);
        let mut sp_tree_branch_length_s: Vec<f32> = Vec::with_capacity(n_sp_nodes);
        for &idx in &sp_preorder {
            let node = &extant_sp_arc.nodes[idx];
            sp_tree_parent_s.push(
                node.parent
                    .map(|p| sp_name_to_idx_map[&extant_sp_arc.nodes[p].name] as i32)
                    .unwrap_or(-1),
            );
            sp_tree_branch_length_s.push(if node.parent.is_some() {
                node.length as f32
            } else {
                0.0
            });
        }

        let mut sp_child_src_s: Vec<i32> = Vec::new();
        let mut sp_child_dst_s: Vec<i32> = Vec::new();
        let mut sp_parent_src_s: Vec<i32> = Vec::new();
        let mut sp_parent_dst_s: Vec<i32> = Vec::new();
        for &idx in &sp_preorder {
            let node = &extant_sp_arc.nodes[idx];
            let pi = sp_name_to_idx_map[&node.name] as i32;
            for child in [node.left_child, node.right_child].into_iter().flatten() {
                let ci = sp_name_to_idx_map[&extant_sp_arc.nodes[child].name] as i32;
                sp_child_src_s.push(pi);
                sp_child_dst_s.push(ci);
                sp_parent_src_s.push(ci);
                sp_parent_dst_s.push(pi);
            }
        }

        // 4. Extract base samples from gene trees (parallel)
        let base_samples: Vec<RustBaseSample> = valid_gene_trees
            .par_iter()
            .map(extract_base_sample_internal)
            .collect::<Result<Vec<_>, _>>()?;

        // 5. Build task tensors (parallel -- each sample has its own seed/RNG)
        let map_tensors: Vec<TaskTensors> = base_samples
            .par_iter()
            .enumerate()
            .map(|(i, s)| {
                let mut rng = StdRng::seed_from_u64(map_coloring_seeds[i]);
                build_task_tensors_internal(s, &sp_name_to_id, &mut rng, false, false, sample_order)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let evt_tensors: Vec<TaskTensors> = if enable_event {
            base_samples
                .par_iter()
                .enumerate()
                .map(|(i, s)| {
                    let mut rng = StdRng::seed_from_u64(evt_coloring_seeds[i]);
                    build_task_tensors_internal(
                        s,
                        &sp_name_to_id,
                        &mut rng,
                        true,
                        false,
                        sample_order,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            Vec::new()
        };

        let root_tensors: Vec<TaskTensors> = if enable_root {
            base_samples
                .par_iter()
                .enumerate()
                .map(|(i, s)| {
                    let mut rng = StdRng::seed_from_u64(root_coloring_seeds[i]);
                    build_task_tensors_internal(
                        s,
                        &sp_name_to_id,
                        &mut rng,
                        false,
                        true,
                        sample_order,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            Vec::new()
        };

        // 6. Compute species GCN norm once on single-tree edges, then replicate B times
        let (sp_child_norm_src_s, sp_child_norm_dst_s, sp_child_ew_s) =
            gcn_norm_internal(&sp_child_src_s, &sp_child_dst_s, n_sp_nodes);
        let (sp_parent_norm_src_s, sp_parent_norm_dst_s, sp_parent_ew_s) =
            gcn_norm_internal(&sp_parent_src_s, &sp_parent_dst_s, n_sp_nodes);

        let b = n_gene_trees;
        let n_sp_child_e = sp_child_norm_src_s.len(); // includes self-loops
        let n_sp_parent_e = sp_parent_norm_src_s.len();
        let mut sp_x: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_child_src: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_dst: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_ew: Vec<f32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_parent_src: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_dst: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_ew: Vec<f32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_tree_parent: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_tree_branch_length: Vec<f32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_batch: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_sizes: Vec<i32> = Vec::with_capacity(b);
        let mut sp_ptr: Vec<i32> = Vec::with_capacity(b + 1);
        sp_ptr.push(0);

        for s_idx in 0..b {
            let offset = (s_idx * n_sp_nodes) as i32;
            sp_x.extend_from_slice(&x_sp_single);
            for (&s, &d) in sp_child_norm_src_s.iter().zip(sp_child_norm_dst_s.iter()) {
                sp_child_src.push(s + offset);
                sp_child_dst.push(d + offset);
            }
            sp_child_ew.extend_from_slice(&sp_child_ew_s);
            for (&s, &d) in sp_parent_norm_src_s.iter().zip(sp_parent_norm_dst_s.iter()) {
                sp_parent_src.push(s + offset);
                sp_parent_dst.push(d + offset);
            }
            sp_parent_ew.extend_from_slice(&sp_parent_ew_s);
            for &p in &sp_tree_parent_s {
                sp_tree_parent.push(if p >= 0 { p + offset } else { -1 });
            }
            sp_tree_branch_length.extend_from_slice(&sp_tree_branch_length_s);
            for _ in 0..n_sp_nodes {
                sp_batch.push(s_idx as i32);
            }
            sp_sizes.push(n_sp_nodes as i32);
            sp_ptr.push((s_idx + 1) as i32 * n_sp_nodes as i32);
        }

        // Compute species varlen attention metadata
        let (cu_sp, max_sp) = compute_varlen_metadata(&sp_sizes);

        // 7. Collate task tensors (includes gene GCN norm + varlen metadata)
        let map_c = collate_task_tensors(&map_tensors);
        let evt_c = if enable_event {
            Some(collate_task_tensors(&evt_tensors))
        } else {
            None
        };
        let root_c = if enable_root {
            Some(collate_task_tensors(&root_tensors))
        } else {
            None
        };

        Ok((
            sp_x,
            sp_child_src,
            sp_child_dst,
            sp_child_ew,
            sp_parent_src,
            sp_parent_dst,
            sp_parent_ew,
            sp_tree_parent,
            sp_tree_branch_length,
            sp_ptr,
            sp_sizes,
            sp_batch,
            cu_sp,
            max_sp,
            map_c,
            evt_c,
            root_c,
        ))
    });
    // ---- End GIL-free section --------------------------------------------

    let (
        sp_x,
        sp_child_src,
        sp_child_dst,
        sp_child_ew,
        sp_parent_src,
        sp_parent_dst,
        sp_parent_ew,
        sp_tree_parent,
        sp_tree_branch_length,
        sp_ptr,
        sp_sizes,
        sp_batch,
        cu_sp,
        max_sp,
        map_c,
        evt_c,
        root_c,
    ) = result_or_err.map_err(PyValueError::new_err)?;

    let (species_tree_metadata, map_tree_metadata, evt_tree_metadata, root_tree_metadata) = py
        .allow_threads(|| -> Result<_, String> {
            let species_segment_len = sp_sizes.first().copied().unwrap_or(0) as usize;
            let species_segments = sp_sizes.len();
            let species_parent_single = &sp_tree_parent[..species_segment_len];
            let species_length_single = &sp_tree_branch_length[..species_segment_len];
            Ok((
                build_repeated_tree_metadata_bundle(
                    species_parent_single,
                    species_length_single,
                    species_segments,
                    species_tree_metadata_request,
                )?,
                build_tree_metadata_bundle(
                    &map_c.gene_tree_parent,
                    &map_c.cu_g,
                    &map_c.gene_tree_branch_length,
                    gene_tree_metadata_request,
                )?,
                if let Some(ref c) = evt_c {
                    Some(build_tree_metadata_bundle(
                        &c.gene_tree_parent,
                        &c.cu_g,
                        &c.gene_tree_branch_length,
                        gene_tree_metadata_request,
                    )?)
                } else {
                    None
                },
                if let Some(ref c) = root_c {
                    Some(build_tree_metadata_bundle(
                        &c.gene_tree_parent,
                        &c.cu_g,
                        &c.gene_tree_branch_length,
                        gene_tree_metadata_request,
                    )?)
                } else {
                    None
                },
            ))
        })
        .map_err(PyValueError::new_err)?;

    // 8. Build result dict of numpy arrays
    let result = PyDict::new(py);

    let edge_2d = |src: &[i32], dst: &[i32], name: &str| -> PyResult<PyObject> {
        use numpy::ndarray::Array2;
        use numpy::ToPyArray;
        let n = src.len();
        let mut flat = Vec::with_capacity(2 * n);
        flat.extend_from_slice(src);
        flat.extend_from_slice(dst);
        let arr = Array2::from_shape_vec([2, n], flat)
            .map_err(|e| PyValueError::new_err(format!("{}: {}", name, e)))?;
        Ok(arr.to_pyarray(py).into_any().unbind())
    };

    // Species (GCN-normalized edges with self-loops)
    result.set_item("sp_x", PyArray1::from_slice(py, &sp_x))?;
    result.set_item(
        "sp_child_ei",
        edge_2d(&sp_child_src, &sp_child_dst, "sp_child_ei")?,
    )?;
    result.set_item("sp_child_ew", PyArray1::from_slice(py, &sp_child_ew))?;
    result.set_item(
        "sp_parent_ei",
        edge_2d(&sp_parent_src, &sp_parent_dst, "sp_parent_ei")?,
    )?;
    result.set_item("sp_parent_ew", PyArray1::from_slice(py, &sp_parent_ew))?;
    result.set_item("sp_tree_parent", PyArray1::from_slice(py, &sp_tree_parent))?;
    result.set_item(
        "sp_tree_branch_length",
        PyArray1::from_slice(py, &sp_tree_branch_length),
    )?;
    result.set_item("sp_ptr", PyArray1::from_slice(py, &sp_ptr))?;
    result.set_item("sp_sizes", PyArray1::from_slice(py, &sp_sizes))?;
    result.set_item("sp_batch", PyArray1::from_slice(py, &sp_batch))?;
    result.set_item("cu_sp", PyArray1::from_slice(py, &cu_sp))?;
    result.set_item("max_sp", max_sp)?;
    write_tree_metadata_bundle(py, &result, "species", &species_tree_metadata)?;

    // Helper: write one collated task into result dict under a given prefix
    let write_task =
        |c: &CollatedTask, prefix: &str, tree_metadata: &TreeMetadataBundle| -> PyResult<()> {
            result.set_item(
                format!("{}_gene_x", prefix),
                PyArray1::from_slice(py, &c.gene_x),
            )?;
            result.set_item(
                format!("{}_gene_y", prefix),
                PyArray1::from_slice(py, &c.gene_y),
            )?;
            result.set_item(
                format!("{}_event", prefix),
                PyArray1::from_slice(py, &c.event),
            )?;
            result.set_item(
                format!("{}_event_true", prefix),
                PyArray1::from_slice(py, &c.event_true),
            )?;
            result.set_item(
                format!("{}_is_leaf", prefix),
                PyArray1::from_slice(py, &c.is_leaf),
            )?;
            result.set_item(
                format!("{}_frontier_mask", prefix),
                PyArray1::from_slice(py, &c.frontier_mask),
            )?;
            result.set_item(
                format!("{}_mask_label_node", prefix),
                PyArray1::from_slice(py, &c.mask_label_node),
            )?;
            result.set_item(
                format!("{}_gene_ei", prefix),
                edge_2d(
                    &c.gene_ei_src,
                    &c.gene_ei_dst,
                    &format!("{}_gene_ei", prefix),
                )?,
            )?;
            result.set_item(
                format!("{}_gene_ew", prefix),
                PyArray1::from_slice(py, &c.gene_ew),
            )?;
            result.set_item(
                format!("{}_g_dir_edge", prefix),
                edge_2d(
                    &c.g_dir_src,
                    &c.g_dir_dst,
                    &format!("{}_g_dir_edge", prefix),
                )?,
            )?;
            result.set_item(
                format!("{}_gene_tree_parent", prefix),
                PyArray1::from_slice(py, &c.gene_tree_parent),
            )?;
            result.set_item(
                format!("{}_gene_tree_branch_length", prefix),
                PyArray1::from_slice(py, &c.gene_tree_branch_length),
            )?;
            write_tree_metadata_bundle(py, &result, &format!("{}_gene", prefix), tree_metadata)?;
            result.set_item(
                format!("{}_true_root_index", prefix),
                PyArray1::from_slice(py, &c.true_root_index),
            )?;
            result.set_item(
                format!("{}_g_ptr", prefix),
                PyArray1::from_slice(py, &c.g_ptr),
            )?;
            result.set_item(
                format!("{}_g_sizes", prefix),
                PyArray1::from_slice(py, &c.g_sizes),
            )?;
            result.set_item(
                format!("{}_g_batch", prefix),
                PyArray1::from_slice(py, &c.g_batch),
            )?;
            result.set_item(
                format!("{}_root_edges_ptr", prefix),
                PyArray1::from_slice(py, &c.root_edges_ptr),
            )?;
            result.set_item(
                format!("{}_cu_g", prefix),
                PyArray1::from_slice(py, &c.cu_g),
            )?;
            result.set_item(format!("{}_max_g", prefix), c.max_g)?;
            Ok(())
        };

    write_task(&map_c, "map", &map_tree_metadata)?;
    if let (Some(ref c), Some(ref metadata)) = (&evt_c, &evt_tree_metadata) {
        write_task(c, "evt", metadata)?;
    }
    if let (Some(ref c), Some(ref metadata)) = (&root_c, &root_tree_metadata) {
        write_task(c, "root", metadata)?;
    }

    Ok(result.into())
}

/// Compute GCN normalization on edge indices.
///
/// Adds self-loops and computes D^{-1/2} A D^{-1/2} edge weights.
/// Matches PyG's `gcn_norm(edge_index, None, num_nodes, improved=False,
/// add_self_loops=True, flow="source_to_target")`.
///
/// Parameters
/// ----------
/// edge_index : numpy array of shape (2, E), int32
/// num_nodes : int
///
/// Returns (edge_index_with_self_loops, edge_weights) as numpy arrays.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (edge_index, num_nodes))]
pub fn compute_gcn_norm(
    py: Python,
    edge_index: &Bound<'_, numpy::PyArray2<i32>>,
    num_nodes: usize,
) -> PyResult<(PyObject, PyObject)> {
    use numpy::ndarray::Array2;
    use numpy::ToPyArray;
    use numpy::{PyArray1, PyArrayMethods, PyUntypedArrayMethods};

    let shape = edge_index.shape();
    if shape[0] != 2 {
        return Err(PyValueError::new_err(format!(
            "compute_gcn_norm: edge_index must have shape (2, E), got ({}, {})",
            shape[0], shape[1]
        )));
    }

    let ei = unsafe { edge_index.as_array() };
    let src: Vec<i32> = ei.row(0).to_vec();
    let dst: Vec<i32> = ei.row(1).to_vec();

    for (&s, &d) in src.iter().zip(dst.iter()) {
        if s < 0 || d < 0 || s as usize >= num_nodes || d as usize >= num_nodes {
            return Err(PyValueError::new_err(format!(
                "compute_gcn_norm: edge indices must be in [0, {}), got ({}, {})",
                num_nodes, s, d
            )));
        }
    }

    let (new_src, new_dst, weights) = gcn_norm_internal(&src, &dst, num_nodes);

    let n = new_src.len();
    let mut flat = Vec::with_capacity(2 * n);
    flat.extend_from_slice(&new_src);
    flat.extend_from_slice(&new_dst);
    let arr = Array2::from_shape_vec([2, n], flat)
        .map_err(|e| PyValueError::new_err(format!("compute_gcn_norm: {}", e)))?;

    let ei_out = arr.to_pyarray(py).into_any().unbind();
    let ew_out = PyArray1::from_slice(py, &weights).into_any().unbind();
    Ok((ei_out, ew_out))
}

/// Build a batched dict of numpy arrays for inference from a base_sample dict.
///
/// Replaces the Python get_sample() + collate() path. Builds task tensors,
/// collates n_copies, applies GCN norm and computes varlen metadata.
/// Returns a flat numpy dict in the same format as build_otf_batch().
///
/// Only mapping task tensors are returned (inference doesn't need separate evt/root batches).
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (base_sample, n_copies, seed, sample_order))]
pub fn build_inference_batch(
    py: Python,
    base_sample: &Bound<'_, pyo3::types::PyDict>,
    n_copies: usize,
    seed: u64,
    sample_order: &str,
) -> PyResult<PyObject> {
    use numpy::PyArray1;
    use pyo3::types::PyDict;

    // Parse the base_sample dict into RustBaseSample
    let base = parse_base_sample_dict(base_sample)?;

    // Build species name -> 1-based ID mapping
    let mut sp_name_to_id: HashMap<String, i32> = HashMap::new();
    for (i, name) in base.species_names.iter().enumerate() {
        sp_name_to_id.insert(name.clone(), (i + 1) as i32);
    }

    let n_sp_nodes = base.species_names.len();
    let sample_order_owned = sample_order.to_owned();

    let result_or_err = py.allow_threads(move || -> Result<_, String> {
        let sample_order: &str = &sample_order_owned;

        // Build species tensors (single copy)
        let x_sp_single: Vec<i32> = base
            .species_names
            .iter()
            .map(|name| sp_name_to_id[name])
            .collect();

        let mut sp_child_src_s: Vec<i32> = Vec::new();
        let mut sp_child_dst_s: Vec<i32> = Vec::new();
        let mut sp_parent_src_s: Vec<i32> = Vec::new();
        let mut sp_parent_dst_s: Vec<i32> = Vec::new();
        let mut sp_tree_parent_s: Vec<i32> = vec![-1; n_sp_nodes];
        let mut sp_tree_branch_length_s: Vec<f32> = vec![0.0; n_sp_nodes];

        // Build species name to index mapping
        let mut sp_name_to_idx: HashMap<String, usize> = HashMap::new();
        for (i, name) in base.species_names.iter().enumerate() {
            sp_name_to_idx.insert(name.clone(), i);
        }

        for (parent_name, children) in &base.sp_children {
            if let Some(&pi) = sp_name_to_idx.get(parent_name) {
                for child_name in children {
                    if let Some(&ci) = sp_name_to_idx.get(child_name) {
                        sp_child_src_s.push(pi as i32);
                        sp_child_dst_s.push(ci as i32);
                        sp_parent_src_s.push(ci as i32);
                        sp_parent_dst_s.push(pi as i32);
                        sp_tree_parent_s[ci] = pi as i32;
                        sp_tree_branch_length_s[ci] =
                            *base.sp_branch_length.get(child_name).unwrap_or(&1.0);
                    }
                }
            }
        }
        for (name, parent_name) in &base.sp_parent {
            if let (Some(&ci), Some(&pi)) =
                (sp_name_to_idx.get(name), sp_name_to_idx.get(parent_name))
            {
                sp_tree_parent_s[ci] = pi as i32;
                sp_tree_branch_length_s[ci] = *base.sp_branch_length.get(name).unwrap_or(&1.0);
            }
        }

        // GCN norm on species edges (single copy)
        let (sp_child_norm_src_s, sp_child_norm_dst_s, sp_child_ew_s) =
            gcn_norm_internal(&sp_child_src_s, &sp_child_dst_s, n_sp_nodes);
        let (sp_parent_norm_src_s, sp_parent_norm_dst_s, sp_parent_ew_s) =
            gcn_norm_internal(&sp_parent_src_s, &sp_parent_dst_s, n_sp_nodes);

        // Build mapping task tensors for each copy
        let mut map_tensors: Vec<TaskTensors> = Vec::with_capacity(n_copies);
        for copy_idx in 0..n_copies {
            let copy_seed = seed.wrapping_add(copy_idx as u64 * 1000003);
            let mut rng = StdRng::seed_from_u64(copy_seed);
            // Inference must always operate on the full reconstructed tree
            // (unrooted, all internal nodes present). Using the training-time
            // partial-coloring path here creates variable-size copies, which
            // breaks batched inference reshapes on the Python side.
            let t = build_task_tensors_internal(
                &base,
                &sp_name_to_id,
                &mut rng,
                false,
                true,
                sample_order,
            )?;
            map_tensors.push(t);
        }

        // Collate
        let map_c = collate_task_tensors(&map_tensors);

        // Replicate species B times
        let b = n_copies;
        let n_sp_child_e = sp_child_norm_src_s.len();
        let n_sp_parent_e = sp_parent_norm_src_s.len();
        let mut sp_x: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_child_src: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_dst: Vec<i32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_child_ew: Vec<f32> = Vec::with_capacity(b * n_sp_child_e);
        let mut sp_parent_src: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_dst: Vec<i32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_parent_ew: Vec<f32> = Vec::with_capacity(b * n_sp_parent_e);
        let mut sp_tree_parent: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_tree_branch_length: Vec<f32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_batch: Vec<i32> = Vec::with_capacity(b * n_sp_nodes);
        let mut sp_sizes: Vec<i32> = Vec::with_capacity(b);
        let mut sp_ptr: Vec<i32> = Vec::with_capacity(b + 1);
        sp_ptr.push(0);

        for s_idx in 0..b {
            let offset = (s_idx * n_sp_nodes) as i32;
            sp_x.extend_from_slice(&x_sp_single);
            for (&s, &d) in sp_child_norm_src_s.iter().zip(sp_child_norm_dst_s.iter()) {
                sp_child_src.push(s + offset);
                sp_child_dst.push(d + offset);
            }
            sp_child_ew.extend_from_slice(&sp_child_ew_s);
            for (&s, &d) in sp_parent_norm_src_s.iter().zip(sp_parent_norm_dst_s.iter()) {
                sp_parent_src.push(s + offset);
                sp_parent_dst.push(d + offset);
            }
            sp_parent_ew.extend_from_slice(&sp_parent_ew_s);
            for &p in &sp_tree_parent_s {
                sp_tree_parent.push(if p >= 0 { p + offset } else { -1 });
            }
            sp_tree_branch_length.extend_from_slice(&sp_tree_branch_length_s);
            for _ in 0..n_sp_nodes {
                sp_batch.push(s_idx as i32);
            }
            sp_sizes.push(n_sp_nodes as i32);
            sp_ptr.push((s_idx + 1) as i32 * n_sp_nodes as i32);
        }

        let (cu_sp, max_sp) = compute_varlen_metadata(&sp_sizes);

        Ok((
            sp_x,
            sp_child_src,
            sp_child_dst,
            sp_child_ew,
            sp_parent_src,
            sp_parent_dst,
            sp_parent_ew,
            sp_tree_parent,
            sp_tree_branch_length,
            sp_ptr,
            sp_sizes,
            sp_batch,
            cu_sp,
            max_sp,
            map_c,
        ))
    });

    let (
        sp_x,
        sp_child_src,
        sp_child_dst,
        sp_child_ew,
        sp_parent_src,
        sp_parent_dst,
        sp_parent_ew,
        sp_tree_parent,
        sp_tree_branch_length,
        sp_ptr,
        sp_sizes,
        sp_batch,
        cu_sp,
        max_sp,
        map_c,
    ) = result_or_err.map_err(PyValueError::new_err)?;

    let result = PyDict::new(py);

    let edge_2d = |src: &[i32], dst: &[i32], name: &str| -> PyResult<PyObject> {
        use numpy::ndarray::Array2;
        use numpy::ToPyArray;
        let n = src.len();
        let mut flat = Vec::with_capacity(2 * n);
        flat.extend_from_slice(src);
        flat.extend_from_slice(dst);
        let arr = Array2::from_shape_vec([2, n], flat)
            .map_err(|e| PyValueError::new_err(format!("{}: {}", name, e)))?;
        Ok(arr.to_pyarray(py).into_any().unbind())
    };

    // Species
    result.set_item("sp_x", PyArray1::from_slice(py, &sp_x))?;
    result.set_item(
        "sp_child_ei",
        edge_2d(&sp_child_src, &sp_child_dst, "sp_child_ei")?,
    )?;
    result.set_item("sp_child_ew", PyArray1::from_slice(py, &sp_child_ew))?;
    result.set_item(
        "sp_parent_ei",
        edge_2d(&sp_parent_src, &sp_parent_dst, "sp_parent_ei")?,
    )?;
    result.set_item("sp_parent_ew", PyArray1::from_slice(py, &sp_parent_ew))?;
    result.set_item("sp_tree_parent", PyArray1::from_slice(py, &sp_tree_parent))?;
    result.set_item(
        "sp_tree_branch_length",
        PyArray1::from_slice(py, &sp_tree_branch_length),
    )?;
    result.set_item("sp_ptr", PyArray1::from_slice(py, &sp_ptr))?;
    result.set_item("sp_sizes", PyArray1::from_slice(py, &sp_sizes))?;
    result.set_item("sp_batch", PyArray1::from_slice(py, &sp_batch))?;
    result.set_item("cu_sp", PyArray1::from_slice(py, &cu_sp))?;
    result.set_item("max_sp", max_sp)?;

    // Mapping task (only task needed for inference)
    let c = &map_c;
    result.set_item("map_gene_x", PyArray1::from_slice(py, &c.gene_x))?;
    result.set_item("map_gene_y", PyArray1::from_slice(py, &c.gene_y))?;
    result.set_item("map_event", PyArray1::from_slice(py, &c.event))?;
    result.set_item("map_event_true", PyArray1::from_slice(py, &c.event_true))?;
    result.set_item("map_is_leaf", PyArray1::from_slice(py, &c.is_leaf))?;
    result.set_item(
        "map_frontier_mask",
        PyArray1::from_slice(py, &c.frontier_mask),
    )?;
    result.set_item(
        "map_mask_label_node",
        PyArray1::from_slice(py, &c.mask_label_node),
    )?;
    result.set_item(
        "map_gene_ei",
        edge_2d(&c.gene_ei_src, &c.gene_ei_dst, "map_gene_ei")?,
    )?;
    result.set_item("map_gene_ew", PyArray1::from_slice(py, &c.gene_ew))?;
    result.set_item(
        "map_g_dir_edge",
        edge_2d(&c.g_dir_src, &c.g_dir_dst, "map_g_dir_edge")?,
    )?;
    result.set_item(
        "map_gene_tree_parent",
        PyArray1::from_slice(py, &c.gene_tree_parent),
    )?;
    result.set_item(
        "map_gene_tree_branch_length",
        PyArray1::from_slice(py, &c.gene_tree_branch_length),
    )?;
    result.set_item(
        "map_true_root_index",
        PyArray1::from_slice(py, &c.true_root_index),
    )?;
    result.set_item("map_g_ptr", PyArray1::from_slice(py, &c.g_ptr))?;
    result.set_item("map_g_sizes", PyArray1::from_slice(py, &c.g_sizes))?;
    result.set_item("map_g_batch", PyArray1::from_slice(py, &c.g_batch))?;
    result.set_item(
        "map_root_edges_ptr",
        PyArray1::from_slice(py, &c.root_edges_ptr),
    )?;
    result.set_item("map_cu_g", PyArray1::from_slice(py, &c.cu_g))?;
    result.set_item("map_max_g", c.max_g)?;

    Ok(result.into())
}
