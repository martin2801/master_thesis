#!/usr/bin/env python3
"""
match_alignment_to_tree.py

prepare_relax_trees.py prunes the tree per sweep down to that sweep's
genomes + 'no_sweep' genomes (dropping every other sweep's genomes). The
alignment it's paired with needs the same taxon set. This subsets a FASTA
alignment down to exactly the tip names present in a given tree file.

Tip names are extracted from the Newick with a regex rather than a full
parser, because prepare_relax_trees.py writes HyPhy's non-standard
"{Foreground}" branch tag after branch lengths, which generic Newick
parsers (including ete3's reader) don't understand. In Newick grammar, a
leaf name can only appear immediately after '(' or ',' -- internal node
labels appear after ')' -- so this reliably captures only leaf tips
regardless of the trailing HyPhy tag.

Usage:
    python3 match_alignment_to_tree.py \
        --alignment concatenated_core.fasta \
        --tree relax_inputs_concat/sweep_7/concatenated_core.nwk \
        --out subset.fasta
"""

import argparse
import re
from pathlib import Path

from Bio import SeqIO

LEAF_NAME_RE = re.compile(r'[(,]([A-Za-z0-9_.\-]+):')


def get_tree_tip_names(nwk_path: Path) -> set:
    text = nwk_path.read_text()
    return set(LEAF_NAME_RE.findall(text))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--tree", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    tip_names = get_tree_tip_names(args.tree)
    print(f"{len(tip_names)} tips found in {args.tree}")

    records = list(SeqIO.parse(args.alignment, "fasta"))
    kept = [r for r in records if r.id in tip_names]
    missing_from_aln = tip_names - {r.id for r in records}

    with open(args.out, "w") as f:
        for r in kept:
            f.write(f">{r.id}\n{str(r.seq)}\n")

    print(f"Alignment had {len(records)} sequences; kept {len(kept)} matching the tree")
    if missing_from_aln:
        print(f"WARNING: {len(missing_from_aln)} tree tips have no matching sequence "
              f"in the alignment: {sorted(missing_from_aln)[:10]}{' ...' if len(missing_from_aln) > 10 else ''}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()