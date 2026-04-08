"""Convert a Zombi-style events TSV + species Newick into a RecPhyloXML.

Reads ``testdata/induced_tr/complete_tree.nwk`` and ``complete_events.tsv``,
builds the reconciled gene tree from the event trace, and writes the
RecPhyloXML to disk via :py:func:`rustree.from_reconciliation` /
:py:meth:`PyGeneTree.save_xml`.

Event TSV format (Zombi):
    TIME    EVENT    NODES
where EVENT is one of:
    O                 origin (single species, gene id is the root)
    S / D / T         binary events: parent_species;parent_gene;c1_sp;c1_g;c2_sp;c2_g
    L                 loss leaf:     species;gene
    F                 final / extant leaf at present time: species;gene
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import rustree

DATA_DIR = Path(__file__).resolve().parent.parent / "testdata" / "induced_tr"
COMPLETE_TREE_FILE = DATA_DIR / "complete_tree.nwk"
EVENTS_FILE = DATA_DIR / "complete_events.tsv"
DEFAULT_OUTPUT = DATA_DIR / "complete_reconciled.xml"

# Event codes accepted by rustree.from_reconciliation
EVT_SPECIATION = 0
EVT_DUPLICATION = 1
EVT_TRANSFER = 2
EVT_LEAF = 3
EVT_LOSS = 4


@dataclass
class GeneNode:
    gene_id: str
    species: str
    event_code: int
    time: float
    children: list[str]
    parent: str | None = None
    parent_time: float = 0.0  # set in a second pass after we know the parent


def parse_events(path: Path) -> tuple[dict[str, GeneNode], str]:
    """Build a gene-id-keyed map of GeneNode plus the root gene id.

    The root gene id is taken from the parent of the first ``S/D/T`` event
    (Zombi's ``O`` row only carries the species, not the gene id).
    """
    nodes: dict[str, GeneNode] = {}
    root_gene: str | None = None

    with path.open() as f:
        header = next(f).rstrip("\n").split("\t")
        assert header == ["TIME", "EVENT", "NODES"], f"unexpected header: {header}"
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            time_s, event, nodes_s = line.split("\t")
            time = float(time_s)
            parts = nodes_s.split(";")

            if event == "O":
                # Origin: just records the species; the gene id is the parent
                # of the first speciation/duplication/transfer event.
                continue

            if event in ("S", "D", "T"):
                p_sp, p_g = parts[0], parts[1]
                c1_sp, c1_g = parts[2], parts[3]
                c2_sp, c2_g = parts[4], parts[5]

                code = {
                    "S": EVT_SPECIATION,
                    "D": EVT_DUPLICATION,
                    "T": EVT_TRANSFER,
                }[event]
                nodes[p_g] = GeneNode(
                    gene_id=p_g,
                    species=p_sp,
                    event_code=code,
                    time=time,
                    children=[c1_g, c2_g],
                )
                if root_gene is None:
                    root_gene = p_g
            elif event in ("L", "F"):
                sp, g = parts[0], parts[1]
                code = EVT_LOSS if event == "L" else EVT_LEAF
                nodes[g] = GeneNode(
                    gene_id=g,
                    species=sp,
                    event_code=code,
                    time=time,
                    children=[],
                )
            else:
                raise ValueError(f"unsupported event type {event!r}")

    if root_gene is None:
        raise ValueError("no internal events found; cannot determine root gene")

    # Wire up parent / parent_time so we can compute branch lengths
    for n in nodes.values():
        for c in n.children:
            child = nodes.get(c)
            if child is None:
                raise ValueError(
                    f"event references gene {c!r} which never appears as a parent or leaf"
                )
            child.parent = n.gene_id
            child.parent_time = n.time

    return nodes, root_gene


def build_newick(nodes: dict[str, GeneNode], root: str) -> str:
    """Emit Newick with branch lengths = (node.time - parent.time)."""
    sb: list[str] = []

    def emit(gid: str) -> None:
        n = nodes[gid]
        if n.children:
            sb.append("(")
            for i, c in enumerate(n.children):
                if i:
                    sb.append(",")
                emit(c)
            sb.append(")")
        sb.append(gid)
        sb.append(f":{max(0.0, n.time - n.parent_time)}")

    emit(root)
    sb.append(";")
    return "".join(sb)


def main(argv: list[str]) -> int:
    output = Path(argv[1]) if len(argv) > 1 else DEFAULT_OUTPUT

    sp_newick = COMPLETE_TREE_FILE.read_text().strip()

    nodes, root_gene = parse_events(EVENTS_FILE)
    print(f"parsed {len(nodes)} gene nodes; root gene id = {root_gene}")

    counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    for n in nodes.values():
        counts[n.event_code] += 1
    print(
        f"  speciations={counts[EVT_SPECIATION]}  duplications={counts[EVT_DUPLICATION]}  "
        f"transfers={counts[EVT_TRANSFER]}  leaves={counts[EVT_LEAF]}  losses={counts[EVT_LOSS]}"
    )

    g_newick = build_newick(nodes, root_gene)
    node_species = {gid: n.species for gid, n in nodes.items()}
    node_events = {gid: n.event_code for gid, n in nodes.items()}

    gene_tree = rustree.from_reconciliation(
        sp_newick=sp_newick,
        g_newick=g_newick,
        node_species=node_species,
        node_events=node_events,
    )
    print(
        f"reconciled gene tree: {gene_tree.num_nodes()} nodes, "
        f"{gene_tree.num_extant()} extant"
    )

    gene_tree.save_xml(str(output))
    print(f"\nwrote RecPhyloXML to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
