#!/usr/bin/env python3
"""
For each core gene alignment and each sweep, prune the full phylogeny down
to the genomes present in that gene's alignment AND belonging to either the
target sweep or 'no_sweep', then label the sweep-clade branches for HyPhy
RELAX (Newick branch-label syntax: tag test/foreground branches with
'{Foreground}' after the branch length).

Labeling rule: a branch is Foreground if EVERY leaf descending from it is a
member of the target sweep. This means:
  - all sweep-clade tip branches are labeled Foreground
  - internal branches within the pure sweep clade are labeled Foreground
  - the branch leading INTO the sweep clade (stem branch) is also Foreground,
    since HyPhy's convention for testing "did selection change on this
    lineage" includes the branch where the clade actually originated
  - everything else (no_sweep tips + backbone) is left unlabeled = Reference

Genes with too few surviving sweep tips (default < 3) are skipped and logged,
since RELAX needs enough foreground branches to estimate anything meaningful.

Usage:
    python3 prepare_relax_trees.py \
        --tree full_tree.treefile \
        --labels genome_sweep_labels.txt \
        --gene-dir per_gene_alignments \
        --outdir relax_inputs \
        --min-test-tips 3
"""

import argparse
import csv
import os
import glob
from ete3 import Tree


def load_sweep_labels(labels_path):
    """
    Returns dict: tip_name (prefixed with 'out_') -> sweep label
    (e.g. 'sweep_1', ..., 'sweep_7', 'no_sweep').
    """
    mapping = {}
    with open(labels_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_id = row["genome"].strip()
            sweep = row["sweep"].strip()
            tip_name = f"out_{genome_id}"
            mapping[tip_name] = sweep
    return mapping


def get_gene_tips(fasta_path):
    """Return set of sequence IDs (tip names) present in a gene fasta."""
    tips = set()
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                tips.add(line[1:].strip())
    return tips


def label_and_prune(full_tree, keep_tips, sweep_tips):
    """
    Prune a copy of full_tree to keep_tips, then label branches whose
    entire descendant leaf set is a subset of sweep_tips as Foreground.
    Returns (newick_string_with_labels, n_test_tips_used) or (None, 0)
    if the sweep clade collapses entirely.
    """
    t = full_tree.copy()
    t.prune(list(keep_tips), preserve_branch_length=True)

    surviving_sweep_tips = set(leaf.name for leaf in t.get_leaves()) & sweep_tips
    if len(surviving_sweep_tips) == 0:
        return None, 0

    # Label: any node (leaf or internal) whose full leaf-descendant set is
    # a subset of sweep_tips is Foreground.
    for node in t.traverse("postorder"):
        descendant_leaves = set(leaf.name for leaf in node.get_leaves())
        if descendant_leaves and descendant_leaves.issubset(sweep_tips):
            node.add_feature("is_foreground", True)
        else:
            node.add_feature("is_foreground", False)

    # Build Newick manually with HyPhy-style {Foreground} tags, since
    # ete3's built-in writer doesn't support this custom tag syntax.
    def render(node):
        if node.is_leaf():
            label = node.name
        else:
            children_str = ",".join(render(child) for child in node.children)
            label = f"({children_str})"
            if node.name:
                label += node.name  # support value, if present

        bl = f":{node.dist:.10f}"
        tag = "{Foreground}" if getattr(node, "is_foreground", False) else ""
        return f"{label}{bl}{tag}"

    newick = render(t) + ";"
    return newick, len(surviving_sweep_tips)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True, help="Full IQ-TREE Newick (473 tips)")
    ap.add_argument("--labels", required=True, help="genome_sweep_labels.txt (genome<TAB>sweep)")
    ap.add_argument("--gene-dir", required=True, help="Directory of per-gene FASTA alignments")
    ap.add_argument("--outdir", required=True, help="Output directory for labeled trees")
    ap.add_argument("--min-test-tips", type=int, default=3,
                     help="Skip gene if fewer than this many sweep tips survive pruning")
    ap.add_argument("--drop-reference", action="store_true", default=True,
                     help="Exclude the 'Reference' tip from all trees (default: True)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading full tree...")
    full_tree = Tree(args.tree, format=1)
    all_tip_names = set(leaf.name for leaf in full_tree.get_leaves())
    print(f"  {len(all_tip_names)} tips in full tree")

    if args.drop_reference and "Reference" in all_tip_names:
        full_tree.prune([t for t in all_tip_names if t != "Reference"], preserve_branch_length=True)
        all_tip_names.discard("Reference")
        print(f"  Dropped 'Reference' tip -> {len(all_tip_names)} tips remain")

    print("Loading sweep labels...")
    label_map = load_sweep_labels(args.labels)  # tip_name -> sweep
    sweeps = sorted(set(v for v in label_map.values() if v != "no_sweep"))
    print(f"  {len(sweeps)} sweeps found: {sweeps}")

    no_sweep_tips = set(t for t, s in label_map.items() if s == "no_sweep")
    print(f"  {len(no_sweep_tips)} no_sweep genomes")

    gene_files = sorted(glob.glob(os.path.join(args.gene_dir, "*.fasta")))
    print(f"Found {len(gene_files)} gene alignments\n")

    manifest_rows = []

    for sweep in sweeps:
        sweep_tips = set(t for t, s in label_map.items() if s == sweep)
        sweep_outdir = os.path.join(args.outdir, sweep)
        os.makedirs(sweep_outdir, exist_ok=True)

        n_written = 0
        n_skipped = 0

        for gene_fasta in gene_files:
            gene_id = os.path.splitext(os.path.basename(gene_fasta))[0]
            gene_tips = get_gene_tips(gene_fasta)
            gene_tips.discard("Reference")

            # keep = genes tips that belong to this sweep OR no_sweep
            keep_tips = gene_tips & (sweep_tips | no_sweep_tips)

            if len(keep_tips) < 4:  # need at least a few tips total to do anything
                n_skipped += 1
                manifest_rows.append({
                    "sweep": sweep, "gene_id": gene_id, "status": "skipped_too_few_tips",
                    "n_test_tips": 0, "n_background_tips": 0, "tree_path": ""
                })
                continue

            newick, n_test = label_and_prune(full_tree, keep_tips, sweep_tips)

            if newick is None or n_test < args.min_test_tips:
                n_skipped += 1
                manifest_rows.append({
                    "sweep": sweep, "gene_id": gene_id, "status": "skipped_too_few_sweep_tips",
                    "n_test_tips": n_test, "n_background_tips": len(keep_tips) - n_test,
                    "tree_path": ""
                })
                continue

            tree_path = os.path.join(sweep_outdir, f"{gene_id}.nwk")
            with open(tree_path, "w") as fh:
                fh.write(newick + "\n")

            n_written += 1
            manifest_rows.append({
                "sweep": sweep, "gene_id": gene_id, "status": "ok",
                "n_test_tips": n_test, "n_background_tips": len(keep_tips) - n_test,
                "tree_path": tree_path
            })

        print(f"{sweep}: {n_written} trees written, {n_skipped} genes skipped "
              f"(sweep has {len(sweep_tips)} total genomes)")

    manifest_path = os.path.join(args.outdir, "_manifest.tsv")
    with open(manifest_path, "w", newline="") as fh:
        fieldnames = ["sweep", "gene_id", "status", "n_test_tips", "n_background_tips", "tree_path"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)

    print(f"\nManifest written: {manifest_path}")


if __name__ == "__main__":
    main()