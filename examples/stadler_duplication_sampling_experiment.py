#!/usr/bin/env python3
"""Empirical check of Stadler's sampling-rate transform for DTL duplications.

The transform used here is:

    lambda' = rho * lambda
    mu'     = mu - lambda * (1 - rho)

For two birth-death-sampling triplets with the same (lambda', mu'), Stadler's
result predicts the same sampled-tree distribution. With only gene duplication
events and no transfers/losses, the copy-count distribution at sampled leaves
should therefore be the same up to Monte Carlo error.

This script uses rustree's Python bindings. It simulates complete species trees
conditioned on a fixed number of extant leaves, samples a fixed number of those
leaves, extracts the induced sampled tree, and runs per-gene DTL simulations on
the sampled tree by default. The fixed sampled leaf count keeps both experiments
directly comparable.
"""

from __future__ import annotations

import argparse
import math
import os
import statistics
import sys
from collections import Counter
from dataclasses import dataclass
from typing import Iterable


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET_RELEASE = os.path.join(REPO_ROOT, "target", "release")
if TARGET_RELEASE not in sys.path:
    sys.path.insert(0, TARGET_RELEASE)

import rustree  # noqa: E402


@dataclass(frozen=True)
class Triplet:
    name: str
    lambda_: float
    mu: float
    rho: float

    @property
    def lambda_prime(self) -> float:
        return self.rho * self.lambda_

    @property
    def mu_prime(self) -> float:
        return self.mu - self.lambda_ * (1.0 - self.rho)

    def validate(self) -> None:
        if not (0.0 < self.rho <= 1.0):
            raise ValueError(f"{self.name}: rho must be in (0, 1], got {self.rho}")
        if self.lambda_ <= 0.0:
            raise ValueError(f"{self.name}: lambda must be positive, got {self.lambda_}")
        if self.mu < 0.0:
            raise ValueError(f"{self.name}: mu must be non-negative, got {self.mu}")
        if self.lambda_ <= self.mu:
            raise ValueError(
                f"{self.name}: rustree's BD simulator requires lambda > mu, "
                f"got lambda={self.lambda_}, mu={self.mu}"
            )
        if self.mu_prime < -1e-12:
            raise ValueError(
                f"{self.name}: transformed mu' is negative ({self.mu_prime}); "
                "choose a different apparent rate or rho"
            )


@dataclass
class ExperimentResult:
    triplet: Triplet
    n_complete: int
    counts: list[int]
    per_tree_means: list[float]
    sampled_tree_heights: list[float]


def make_triplets(apparent_lambda: float, apparent_mu: float, second_rho: float) -> tuple[Triplet, Triplet]:
    """Make rho=1 and rho=second_rho triplets with identical apparent rates."""
    first = Triplet("A_rho_1", apparent_lambda, apparent_mu, 1.0)
    second_lambda = apparent_lambda / second_rho
    second_mu = apparent_mu + second_lambda * (1.0 - second_rho)
    second = Triplet("B_incomplete", second_lambda, second_mu, second_rho)
    first.validate()
    second.validate()
    return first, second


def complete_leaf_count(n_sampled: int, rho: float) -> int:
    """Fixed-count approximation to uniform sampling fraction rho."""
    return max(n_sampled, int(round(n_sampled / rho)))


def simulate_one_triplet(
    triplet: Triplet,
    *,
    n_sampled: int,
    species_reps: int,
    gene_reps: int,
    lambda_d: float,
    dtl_on: str,
    base_seed: int,
) -> ExperimentResult:
    n_complete = complete_leaf_count(n_sampled, triplet.rho)
    counts: list[int] = []
    per_tree_means: list[float] = []
    sampled_tree_heights: list[float] = []

    for species_rep in range(species_reps):
        species_seed = base_seed + 10_000 * species_rep + 101
        sample_seed = base_seed + 10_000 * species_rep + 202

        complete_tree = rustree.simulate_species_tree(
            n=n_complete,
            lambda_=triplet.lambda_,
            mu=triplet.mu,
            seed=species_seed,
        )
        sampled_species = complete_tree.sample_leaf_names(n=n_sampled, seed=sample_seed)
        sampled_tree = complete_tree.extract_induced_subtree_by_names(sampled_species)
        sampled_tree_heights.append(sampled_tree.tree_height())

        tree_counts: list[int] = []
        for gene_rep in range(gene_reps):
            gene_seed = base_seed + 10_000 * species_rep + gene_rep + 1_000_003
            if dtl_on == "sampled":
                gene_tree = sampled_tree.simulate_dtl(
                    lambda_d=lambda_d,
                    lambda_t=0.0,
                    lambda_l=0.0,
                    seed=gene_seed,
                )
            elif dtl_on == "complete":
                complete_gene_tree = complete_tree.simulate_dtl(
                    lambda_d=lambda_d,
                    lambda_t=0.0,
                    lambda_l=0.0,
                    seed=gene_seed,
                )
                gene_tree = complete_gene_tree.sample_by_species_names(sampled_species)
            else:
                raise ValueError(f"unknown dtl_on mode: {dtl_on}")
            n_copies = gene_tree.num_extant()
            counts.append(n_copies)
            tree_counts.append(n_copies)

        per_tree_means.append(statistics.fmean(tree_counts))

    return ExperimentResult(
        triplet=triplet,
        n_complete=n_complete,
        counts=counts,
        per_tree_means=per_tree_means,
        sampled_tree_heights=sampled_tree_heights,
    )


def mean(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.fmean(values) if values else float("nan")


def std(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.pstdev(values) if len(values) > 1 else float("nan")


def quantiles(values: list[float]) -> tuple[float, float, float]:
    ordered = sorted(values)
    if not ordered:
        return (float("nan"), float("nan"), float("nan"))

    def q(prob: float) -> float:
        pos = prob * (len(ordered) - 1)
        lo = math.floor(pos)
        hi = math.ceil(pos)
        if lo == hi:
            return ordered[lo]
        weight = pos - lo
        return ordered[lo] * (1.0 - weight) + ordered[hi] * weight

    return q(0.05), q(0.50), q(0.95)


def ks_distance(xs: list[int], ys: list[int]) -> float:
    """Two-sample Kolmogorov-Smirnov distance, without a p-value dependency."""
    x_counts = Counter(xs)
    y_counts = Counter(ys)
    x_total = len(xs)
    y_total = len(ys)
    x_cum = 0
    y_cum = 0
    distance = 0.0
    for value in sorted(set(x_counts) | set(y_counts)):
        x_cum += x_counts[value]
        y_cum += y_counts[value]
        distance = max(distance, abs(x_cum / x_total - y_cum / y_total))
    return distance


def histogram_table(a: list[int], b: list[int], *, max_rows: int = 45) -> str:
    a_counts = Counter(a)
    b_counts = Counter(b)
    values = sorted(set(a_counts) | set(b_counts))

    if len(values) > max_rows:
        # Use integer bins when exact support is too wide.
        lo = min(values)
        hi = max(values)
        width = max(1, math.ceil((hi - lo + 1) / max_rows))
        bins = list(range(lo, hi + 1, width))
        rows = []
        for start in bins:
            end = min(start + width - 1, hi)
            a_bin = sum(count for value, count in a_counts.items() if start <= value <= end)
            b_bin = sum(count for value, count in b_counts.items() if start <= value <= end)
            label = f"{start}-{end}" if start != end else str(start)
            rows.append((label, a_bin, b_bin))
    else:
        rows = [(str(value), a_counts[value], b_counts[value]) for value in values]

    a_total = len(a)
    b_total = len(b)
    out = ["copies  A_count  A_freq    B_count  B_freq"]
    for label, a_count, b_count in rows:
        out.append(
            f"{label:>6}  {a_count:>7}  {a_count / a_total:>7.4f}  "
            f"{b_count:>7}  {b_count / b_total:>7.4f}"
        )
    return "\n".join(out)


def print_summary(result: ExperimentResult, n_sampled: int) -> None:
    q05, q50, q95 = quantiles([float(x) for x in result.counts])
    mean_total = mean(result.counts)
    tree_q05, tree_q50, tree_q95 = quantiles(result.per_tree_means)
    height_q05, height_q50, height_q95 = quantiles(result.sampled_tree_heights)

    print(f"\n{result.triplet.name}")
    print(
        f"  lambda={result.triplet.lambda_:.6g}, mu={result.triplet.mu:.6g}, "
        f"rho={result.triplet.rho:.6g}"
    )
    print(
        f"  lambda'={result.triplet.lambda_prime:.6g}, "
        f"mu'={result.triplet.mu_prime:.6g}"
    )
    print(
        f"  complete extant leaves={result.n_complete}, "
        f"sampled leaves={n_sampled}, simulations={len(result.counts)}"
    )
    print(
        f"  copy count: mean={mean_total:.4f}, sd={std(result.counts):.4f}, "
        f"q05={q05:.1f}, median={q50:.1f}, q95={q95:.1f}, "
        f"mean per sampled leaf={mean_total / n_sampled:.4f}"
    )
    print(
        f"  per-tree expected copies: mean={mean(result.per_tree_means):.4f}, "
        f"sd={std(result.per_tree_means):.4f}, "
        f"q05={tree_q05:.3f}, median={tree_q50:.3f}, q95={tree_q95:.3f}"
    )
    print(
        f"  sampled-tree height: mean={mean(result.sampled_tree_heights):.4f}, "
        f"q05={height_q05:.4f}, median={height_q50:.4f}, q95={height_q95:.4f}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare DTL duplication copy-count distributions for two "
        "birth-death-sampling triplets with the same Stadler transformed rates."
    )
    parser.add_argument("--n-sampled", type=int, default=20)
    parser.add_argument("--species-reps", type=int, default=120)
    parser.add_argument("--gene-reps", type=int, default=120)
    parser.add_argument("--lambda-d", type=float, default=0.12)
    parser.add_argument("--apparent-lambda", type=float, default=1.0)
    parser.add_argument("--apparent-mu", type=float, default=0.2)
    parser.add_argument("--second-rho", type=float, default=0.5)
    parser.add_argument(
        "--dtl-on",
        choices=("sampled", "complete"),
        default="sampled",
        help=(
            "Run DTL on the induced sampled tree, or on the complete tree and "
            "then prune the gene tree to sampled species."
        ),
    )
    parser.add_argument("--seed", type=int, default=12345)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_sampled < 2:
        raise ValueError("--n-sampled must be at least 2")
    if args.species_reps < 1 or args.gene_reps < 1:
        raise ValueError("--species-reps and --gene-reps must be positive")
    if args.lambda_d <= 0.0:
        raise ValueError("--lambda-d must be positive")

    triplet_a, triplet_b = make_triplets(
        args.apparent_lambda,
        args.apparent_mu,
        args.second_rho,
    )

    result_a = simulate_one_triplet(
        triplet_a,
        n_sampled=args.n_sampled,
        species_reps=args.species_reps,
        gene_reps=args.gene_reps,
        lambda_d=args.lambda_d,
        dtl_on=args.dtl_on,
        base_seed=args.seed,
    )
    result_b = simulate_one_triplet(
        triplet_b,
        n_sampled=args.n_sampled,
        species_reps=args.species_reps,
        gene_reps=args.gene_reps,
        lambda_d=args.lambda_d,
        dtl_on=args.dtl_on,
        base_seed=args.seed + 50_000_000,
    )

    print("Stadler sampled-tree DTL duplication experiment")
    print(
        f"n_sampled={args.n_sampled}, species_reps={args.species_reps}, "
        f"gene_reps={args.gene_reps}, lambda_d={args.lambda_d}, "
        f"dtl_on={args.dtl_on}"
    )
    print_summary(result_a, args.n_sampled)
    print_summary(result_b, args.n_sampled)

    mean_diff = mean(result_b.counts) - mean(result_a.counts)
    per_leaf_diff = mean_diff / args.n_sampled
    pooled_mean = (mean(result_a.counts) + mean(result_b.counts)) / 2.0
    rel_diff = mean_diff / pooled_mean if pooled_mean else float("nan")
    ks = ks_distance(result_a.counts, result_b.counts)

    print("\nCopy-count distribution comparison")
    print(f"  mean(B) - mean(A) = {mean_diff:.4f} total copies")
    print(f"  per sampled leaf difference = {per_leaf_diff:.6f}")
    print(f"  relative mean difference = {100.0 * rel_diff:.3f}%")
    print(f"  empirical KS distance = {ks:.4f}")

    print("\nHistogram of total copy counts at sampled leaves")
    print(histogram_table(result_a.counts, result_b.counts))

    print("\nInterpretation")
    if abs(rel_diff) < 0.03 and ks < 0.05:
        print(
            "  The two empirical distributions are not meaningfully different "
            "at this Monte Carlo resolution."
        )
    else:
        print(
            "  The two empirical distributions show a visible difference. "
            "With this script, note that rho<1 is represented by fixed-count "
            "leaf sampling from a complete tree, which is a practical "
            "approximation to the sampling-fraction setup."
        )


if __name__ == "__main__":
    main()
