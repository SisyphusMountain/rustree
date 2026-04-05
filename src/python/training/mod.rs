//! ML/training tensor construction for Python bindings.
//!
//! Contains training sample creation, tensor building, GCN normalization,
//! batch collation, and inference batch construction.

mod collation;
mod extraction;
mod tensors;

pub use collation::{build_inference_batch, build_otf_batch, compute_gcn_norm};
pub use extraction::{
    create_training_sample, create_training_sample_from_sim, from_reconciliation,
};
pub use tensors::build_training_tensors;
