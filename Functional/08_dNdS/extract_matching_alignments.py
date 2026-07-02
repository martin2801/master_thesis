#!/usr/bin/env python3
"""
For every (sweep, gene) pair marked 'ok' in the manifest, extract just the
sequences whose names match that gene's PRUNED TREE tip set, from the full
(all-genomes) per-gene alignment. Writes the subset alignment alongside the
tree file, so relax_inputs/<sweep>/<gene>.fasta and .nwk have matching tip
sets by construction.

This uses the tree file's own tip names as the source of truth (rather than
recomputing sweep membership independently), guaranteeing tree and alignment
always agree - if the tree changes for any reason, rerunning this will stay
in sync automatically.

Usage:
    python3 extract_matching_alignments.py \
        --manifest relax_inputs/_manifest.tsv \
        --gene-dir per_gene_alignments_nostop \
        --tree-dir relax_inputs \
        --outdir relax_inputs
"""

import argparse
import csv
import os
import re


def get_tree_tip_names(newick_path):
    """
    Extract tip (leaf) names from a Newick string without a full parser -
    tip names are any token immediately preceded by '(' or ',' and followed
    by ':' (branch length), which is safe here since internal nodes in our
    trees are unnamed (ete3's writer only emits names for leaves).
    """
    with open(newick_path) as fh:
        newick = fh.read()
    # tip tokens: preceded by '(' or ',', followed by ':'
    tips = re.findall(r'[(,]([A-Za-z0-9_]+):', newick)
    return set(tips)


def read_fasta(path):
    """Return dict header -> sequence."""
    seqs = {}
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            seqs[header] = "".join(seq_lines)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--gene-dir", required=True, help="Full per-gene alignments (all genomes, stop-stripped)")
    ap.add_argument("--tree-dir", required=True, help="Dir containing <sweep>/<gene>.nwk")
    ap.add_argument("--outdir", required=True, help="Where to write <sweep>/<gene>.fasta (usually same as tree-dir)")
    args = ap.parse_args()

    n_ok = 0
    n_mismatch = 0
    mismatch_log = []

    with open(args.manifest) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [r for r in reader if r["status"] == "ok"]

    print(f"Processing {len(rows)} (sweep, gene) pairs...")

    # cache full gene alignments so we don't re-read the same file 7 times
    gene_cache = {}

    for row in rows:
        sweep = row["sweep"]
        gene_id = row["gene_id"]
        tree_path = row["tree_path"]

        if gene_id not in gene_cache:
            full_fasta = os.path.join(args.gene_dir, f"{gene_id}.fasta")
            gene_cache[gene_id] = read_fasta(full_fasta)
        full_seqs = gene_cache[gene_id]

        tip_names = get_tree_tip_names(tree_path)

        out_dir = os.path.join(args.outdir, sweep)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{gene_id}.fasta")

        missing = tip_names - set(full_seqs.keys())
        if missing:
            n_mismatch += 1
            mismatch_log.append((sweep, gene_id, len(missing)))
            continue  # don't write a broken alignment

        with open(out_path, "w") as out_fh:
            for tip in sorted(tip_names):
                out_fh.write(f">{tip}\n{full_seqs[tip]}\n")

        n_ok += 1

    print(f"\n{n_ok} alignment subsets written")
    print(f"{n_mismatch} pairs skipped due to tip names missing from the source alignment "
          f"(should be 0 - investigate if not)")
    if mismatch_log:
        for sweep, gene_id, n_missing in mismatch_log[:20]:
            print(f"  {sweep}\t{gene_id}\t{n_missing} tips not found in source alignment")


if __name__ == "__main__":
    main()