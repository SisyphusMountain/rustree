"""Compute induced transfers from a TSV event trace.

Reads the testdata/induced_tr/ files and verifies the projected induced
transfers against expected_induced.tsv.

The projection itself runs entirely inside rustree via
``PySpeciesTree.compute_induced_transfers``; this script only handles the
TSV parsing and the comparison with the expected output.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import rustree

DATA_DIR = Path(__file__).resolve().parent.parent / "testdata" / "induced_tr"
COMPLETE_TREE_FILE = DATA_DIR / "complete_tree.nwk"
SAMPLED_TREE_FILE = DATA_DIR / "sampled_tree.nwk"
EVENTS_FILE = DATA_DIR / "complete_events.tsv"
EXPECTED_FILE = DATA_DIR / "expected_induced.tsv"


@dataclass(frozen=True)
class TransferEvent:
    time: float
    donor_gene: int
    donor_species: str
    recipient_species: str


def parse_transfers(path: Path) -> list[TransferEvent]:
    """Parse only the transfer (``T``) rows from a Zombi-style events TSV.

    Format:
        TIME\tEVENT\tNODES
    where, for a transfer, NODES is
        donor_species;donor_gene;donor_child_species;donor_child_gene;recipient_species;recipient_gene
    """
    out: list[TransferEvent] = []
    with path.open() as f:
        header = next(f).rstrip("\n").split("\t")
        assert header == ["TIME", "EVENT", "NODES"], f"unexpected header: {header}"
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            time_s, event, nodes_s = line.split("\t")
            if event != "T":
                continue
            p = nodes_s.split(";")
            out.append(
                TransferEvent(
                    time=float(time_s),
                    donor_gene=int(p[1]),
                    donor_species=p[0],
                    recipient_species=p[4],
                )
            )
    return out


def sampled_species_from_gene_tree(newick: str) -> list[str]:
    """Strip ``_<gene_id>`` suffix from each gene-tree leaf to get species names."""
    leaf_labels = re.findall(r"[,(]([A-Za-z0-9_]+):", newick)
    return sorted({lbl.split("_", 1)[0] for lbl in leaf_labels})


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Parse the complete species tree
    # ------------------------------------------------------------------
    complete_tree = rustree.parse_species_tree(COMPLETE_TREE_FILE.read_text().strip())
    print(
        f"complete tree: {complete_tree.num_leaves()} leaves, "
        f"{complete_tree.num_nodes()} nodes"
    )

    # ------------------------------------------------------------------
    # 2. Parse the sampled gene tree, extract sampled species names
    # ------------------------------------------------------------------
    sampled_newick = SAMPLED_TREE_FILE.read_text().strip()
    sampled_gene_tree = rustree.parse_species_tree(sampled_newick)
    sampled_species = sampled_species_from_gene_tree(sampled_newick)
    print(
        f"sampled gene tree: {sampled_gene_tree.num_leaves()} leaves, "
        f"{sampled_gene_tree.num_nodes()} nodes"
    )
    print(f"sampled species (from gene-tree leaves): {len(sampled_species)}")

    # ------------------------------------------------------------------
    # 3. Parse complete_events.tsv into the list of transfer events
    #    (the list of events is enough; we don't need a full reconciled tree)
    # ------------------------------------------------------------------
    transfers = parse_transfers(EVENTS_FILE)
    print(f"transfer events parsed: {len(transfers)}")

    # ------------------------------------------------------------------
    # 4. Compute induced transfers via rustree
    # ------------------------------------------------------------------
    transfer_tuples = [
        (t.time, t.donor_gene, t.donor_species, t.recipient_species) for t in transfers
    ]
    induced_df: pd.DataFrame = complete_tree.compute_induced_transfers(
        sampled_species, transfer_tuples
    )
    print(f"\ninduced transfers DataFrame: {induced_df.shape}")
    print(induced_df.head())

    # ------------------------------------------------------------------
    # 5. Verify against expected_induced.tsv
    # ------------------------------------------------------------------
    expected_df = pd.read_csv(EXPECTED_FILE, sep="\t").dropna(how="all")
    expected_pairs: list[tuple[str, str]] = list(
        map(tuple, expected_df[["FROM", "TO"]].astype(str).values)
    )

    computed_pairs: list[tuple[str, str]] = list(
        zip(induced_df["from_species_sampled"], induced_df["to_species_sampled"])
    )

    computed_counter = Counter(computed_pairs)
    expected_counter = Counter(expected_pairs)

    missing = expected_counter - computed_counter
    extra = computed_counter - expected_counter

    print("\n--- comparison with expected_induced.tsv ---")
    print(f"computed pairs:               {len(computed_pairs)}")
    print(f"unique computed pairs:        {len(computed_counter)}")
    print(f"expected pairs:               {len(expected_pairs)}")
    print(f"exact in-order match?         {computed_pairs == expected_pairs}")
    print(f"multiset match?               {computed_counter == expected_counter}")
    print(f"expected ⊆ computed?          {not missing}")

    if missing:
        print(f"\nexpected pairs NOT produced by the projection ({sum(missing.values())}):")
        for pair, count in missing.items():
            print(f"  {pair} (x{count})")

    # Per-expected-pair check: which complete-tree donors/recipients
    # produced each expected projected pair, if any.
    rows = []
    for exp_from, exp_to in expected_pairs:
        m = induced_df[
            (induced_df["from_species_sampled"] == exp_from)
            & (induced_df["to_species_sampled"] == exp_to)
        ]
        rows.append(
            {
                "from": exp_from,
                "to": exp_to,
                "found": len(m) > 0,
                "num_source_transfers": len(m),
                "source_from_complete": list(m["from_species_complete"].unique())[:5],
                "source_to_complete": list(m["to_species_complete"].unique())[:5],
            }
        )
    print("\nper-expected-pair check:")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
