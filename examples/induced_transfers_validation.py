from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import rustree


DATA_DIR = Path('/home/enzo/Documents/git/rustree/testdata/induced_tr')
COMPLETE_TREE_FILE = DATA_DIR / 'complete_tree.nwk'
SAMPLED_TREE_FILE = DATA_DIR / 'sampled_tree.nwk'
EVENTS_FILE = DATA_DIR / 'complete_events.tsv'
EXPECTED_FILE = DATA_DIR / 'expected_induced.tsv'


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
    parent_time: float = 0.0


@dataclass(frozen=True)
class TransferEvent:
    time: float
    donor_gene: int
    donor_species: str
    recipient_species: str


def parse_events(path: Path) -> tuple[dict[str, GeneNode], str, list[TransferEvent]]:
    nodes: dict[str, GeneNode] = {}
    root_gene: str | None = None
    transfers: list[TransferEvent] = []

    with path.open() as f:
        header = next(f).rstrip('\n').split('\t')
        assert header == ['TIME', 'EVENT', 'NODES'], f'unexpected header: {header}'

        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue

            time_s, event, nodes_s = line.split('\t')
            time = float(time_s)
            p = nodes_s.split(';')

            if event == 'O':
                continue

            if event in ('S', 'D', 'T'):
                parent_sp, parent_g = p[0], p[1]
                _c1_sp, c1_g = p[2], p[3]
                c2_sp, c2_g = p[4], p[5]

                code = {'S': EVT_SPECIATION, 'D': EVT_DUPLICATION, 'T': EVT_TRANSFER}[event]
                nodes[parent_g] = GeneNode(
                    gene_id=parent_g,
                    species=parent_sp,
                    event_code=code,
                    time=time,
                    children=[c1_g, c2_g],
                )

                if root_gene is None:
                    root_gene = parent_g

                if event == 'T':
                    transfers.append(
                        TransferEvent(
                            time=time,
                            donor_gene=int(parent_g),
                            donor_species=parent_sp,
                            recipient_species=c2_sp,
                        )
                    )

            elif event in ('L', 'F'):
                sp, g = p[0], p[1]
                code = EVT_LOSS if event == 'L' else EVT_LEAF
                nodes[g] = GeneNode(
                    gene_id=g,
                    species=sp,
                    event_code=code,
                    time=time,
                    children=[],
                )
            else:
                raise ValueError(f'unsupported event type: {event!r}')

    if root_gene is None:
        raise ValueError('No internal event found, cannot determine root gene id')

    for n in nodes.values():
        for c in n.children:
            if c not in nodes:
                raise ValueError(f'child gene {c} referenced but never defined')
            nodes[c].parent = n.gene_id
            nodes[c].parent_time = n.time

    return nodes, root_gene, transfers


def build_gene_newick(nodes: dict[str, GeneNode], root: str) -> str:
    out: list[str] = []

    def emit(gid: str) -> None:
        n = nodes[gid]
        if n.children:
            out.append('(')
            for i, child in enumerate(n.children):
                if i:
                    out.append(',')
                emit(child)
            out.append(')')
        out.append(gid)
        out.append(f':{max(0.0, n.time - n.parent_time)}')

    emit(root)
    out.append(';')
    return ''.join(out)


def main() -> None:
    assert COMPLETE_TREE_FILE.exists(), COMPLETE_TREE_FILE
    assert SAMPLED_TREE_FILE.exists(), SAMPLED_TREE_FILE
    assert EVENTS_FILE.exists(), EVENTS_FILE
    assert EXPECTED_FILE.exists(), EXPECTED_FILE

    print('Data directory:', DATA_DIR)

    complete_newick = COMPLETE_TREE_FILE.read_text().strip()
    sampled_newick = SAMPLED_TREE_FILE.read_text().strip()

    complete_tree = rustree.parse_species_tree(complete_newick)
    sampled_gene_tree = rustree.parse_species_tree(sampled_newick)

    print('complete tree:', complete_tree.num_leaves(), 'leaves,', complete_tree.num_nodes(), 'nodes')
    print('sampled tree :', sampled_gene_tree.num_leaves(), 'leaves,', sampled_gene_tree.num_nodes(), 'nodes')

    leaf_labels = re.findall(r'[,(]([A-Za-z0-9_]+):', sampled_newick)
    sampled_species = sorted({lbl.split('_', 1)[0] for lbl in leaf_labels})

    print('sampled species names:', len(sampled_species))
    print('first 10 sampled species:', sampled_species[:10])

    nodes, root_gene, transfers = parse_events(EVENTS_FILE)
    gene_newick = build_gene_newick(nodes, root_gene)

    node_species = {gid: n.species for gid, n in nodes.items()}
    node_events = {gid: n.event_code for gid, n in nodes.items()}

    reconciled_gene_tree = rustree.from_reconciliation(
        sp_newick=complete_newick,
        g_newick=gene_newick,
        node_species=node_species,
        node_events=node_events,
    )

    print('parsed events:', len(nodes))
    print('transfer events:', len(transfers))
    print('reconciled gene tree nodes:', reconciled_gene_tree.num_nodes())
    print('reconciled gene tree extant leaves:', reconciled_gene_tree.num_extant())

    transfer_tuples = [
        (t.time, t.donor_gene, t.donor_species, t.recipient_species)
        for t in transfers
    ]

    induced_df: pd.DataFrame = complete_tree.compute_induced_transfers(
        sampled_species, transfer_tuples
    )

    print('induced transfer rows:', induced_df.shape[0])
    print(induced_df.head().to_string(index=False))

    expected_df = pd.read_csv(EXPECTED_FILE, sep='\t').dropna(how='all')
    expected_pairs = list(map(tuple, expected_df[['FROM', 'TO']].astype(str).values))
    computed_pairs = list(zip(induced_df['from_species_sampled'], induced_df['to_species_sampled']))

    expected_counter = Counter(expected_pairs)
    computed_counter = Counter(computed_pairs)

    missing = expected_counter - computed_counter
    extra = computed_counter - expected_counter

    print('expected rows :', len(expected_pairs))
    print('computed rows :', len(computed_pairs))
    print('exact order match   :', computed_pairs == expected_pairs)
    print('multiset exact match:', computed_counter == expected_counter)
    print('missing count:', sum(missing.values()))
    print('extra count  :', sum(extra.values()))

    if missing:
        print('\nMissing pairs:')
        for pair, count in missing.items():
            print(f'  {pair} x{count}')

    if extra:
        print('\nExtra pairs (first 30 shown):')
        for i, (pair, count) in enumerate(extra.items()):
            if i >= 30:
                print(f'  ... and {len(extra)-30} more distinct extra pairs')
                break
            print(f'  {pair} x{count}')

    assert computed_counter == expected_counter, 'Induced transfers do not match expected_induced.tsv'
    print('\n✅ Validation passed: induced transfers match expected_induced.tsv (multiset equality).')


if __name__ == '__main__':
    main()
