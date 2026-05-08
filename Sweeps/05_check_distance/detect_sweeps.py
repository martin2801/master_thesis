#!/usr/bin/env python3
"""
detect_sweeps.py
Identifies selective sweeps from a phylogenetic tree and pairwise SNP distance list.

For each internal node (top-down, non-overlapping):
  1. Get all tips in the clade
  2. Get the immediate sister clade from tree topology
  3. Pick the single closest sister tip (lowest mean SNP dist to clade members)
  4. Compute avg pairwise SNPs within clade
  5. Compute avg SNPs from clade members to that sister tip
  6. Ratio = sister_mean / within_mean
  7. If ratio >= threshold: call sweep, skip all child nodes

Usage:
    python3 detect_sweeps.py --tree full_tree.treefile --snps snp_list.tab [--threshold 5] [--min-tips 3] [--out sweeps.txt]
"""

import argparse
import itertools
from Bio import Phylo

# --- Arguments ---
parser = argparse.ArgumentParser()
parser.add_argument("--tree",      required=True,  help="Newick tree file")
parser.add_argument("--snps",      required=True,  help="Tab-separated pairwise SNP list: sample1, sample2, distance")
parser.add_argument("--threshold", type=float, default=5.0, help="Ratio threshold to call a sweep (default: 5)")
parser.add_argument("--min-tips",  type=int,   default=3,   help="Minimum clade size to test (default: 3)")
parser.add_argument("--out",       default="sweeps.txt",    help="Output file (default: sweeps.txt)")
args = parser.parse_args()

# --- Load SNP distances ---
print("Loading SNP distances...")
snp = {}
with open(args.snps) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        a, b, d = parts[0], parts[1], int(parts[2])
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        snp[key] = d

def get_snp(a, b):
    if a == b:
        return None
    return snp.get((min(a, b), max(a, b)), None)

# --- Load and midpoint-root tree ---
print("Loading tree...")
tree = Phylo.read(args.tree, "newick")
tree.root_at_midpoint()

# --- Build parent lookup ---
parent = {}
for clade in tree.find_clades(order="level"):
    for child in clade.clades:
        parent[id(child)] = clade

def get_sister(clade):
    """Return the sister clade (other child of the parent node)."""
    p = parent.get(id(clade))
    if p is None:
        return None
    for child in p.clades:
        if child is not clade:
            return child
    return None

def get_closest_sister_tip(clade_tips, sister_clade):
    """
    From all tips in the sister clade, return the single tip
    with the lowest mean SNP distance to the clade members.
    """
    sister_tips_names = [c.name for c in sister_clade.get_terminals()]
    best_tip = None
    best_mean = float("inf")
    for tip in sister_tips_names:
        dists = [get_snp(t, tip) for t in clade_tips]
        dists = [d for d in dists if d is not None]
        if not dists:
            continue
        mean_dist = sum(dists) / len(dists)
        if mean_dist < best_mean:
            best_mean = mean_dist
            best_tip = tip
    return best_tip, best_mean

# --- Traverse tree top-down, skipping children of called sweeps ---
print(f"Scanning tree (min clade size: {args.min_tips}, threshold: {args.threshold}x)...")

sweeps = []
swept_clades = []  # store clade objects that were called as sweeps

def is_inside_sweep(clade):
    """Check if this clade is a descendant of any already-called sweep."""
    for swept in swept_clades:
        if clade in swept.find_clades():
            return True
    return False

for clade in tree.find_clades(order="level"):
    if clade.is_terminal():
        continue
    if is_inside_sweep(clade):
        continue

    tips = [c.name for c in clade.get_terminals()]
    if len(tips) < args.min_tips:
        continue

    sister = get_sister(clade)
    if sister is None:
        continue

    # pick closest single tip from sister clade
    sister_tip, sister_mean = get_closest_sister_tip(tips, sister)
    if sister_tip is None:
        continue

    # within-clade pairwise SNPs (exclude self-comparisons)
    within_vals = [get_snp(a, b) for a, b in itertools.combinations(tips, 2) if a != b]
    within_vals = [v for v in within_vals if v is not None]
    if not within_vals:
        continue

    within_mean = sum(within_vals) / len(within_vals)

    # handle fully clonal clade (within_mean == 0)
    if within_mean == 0:
        ratio = float("inf")
    else:
        ratio = sister_mean / within_mean

    if ratio >= args.threshold:
        sweeps.append({
            "clade_tips":   tips,
            "sister_tip":   sister_tip,
            "sister_mean":  sister_mean,
            "within_mean":  within_mean,
            "within_pairs": len(within_vals),
            "ratio":        ratio,
        })
        swept_clades.append(clade)

# sort by ratio descending for output
sweeps.sort(key=lambda x: x["ratio"], reverse=True)

# --- Print summary ---
print(f"\nSweeps detected: {len(sweeps)}\n")
for i, s in enumerate(sweeps):
    ratio_str = f"{s['ratio']:.2f}x" if s['ratio'] != float("inf") else "inf (fully clonal)"
    print(f"--- Sweep {i+1} ---")
    print(f"  Clade size:         {len(s['clade_tips'])} genomes")
    print(f"  Sister tip:         {s['sister_tip']}")
    print(f"  Within-group mean:  {s['within_mean']:.2f} SNPs ({s['within_pairs']} pairs)")
    print(f"  Avg SNPs to sister: {s['sister_mean']:.2f} SNPs")
    print(f"  Ratio:              {ratio_str}")
    print(f"  Genomes:            {', '.join(s['clade_tips'])}")
    print()

# --- Write output table ---
with open(args.out, "w") as out:
    out.write("sweep_id\tgenome\tsister_tip\twithin_mean\tsister_mean\tratio\tclade_size\n")
    for i, s in enumerate(sweeps):
        ratio_str = f"{s['ratio']:.2f}" if s['ratio'] != float("inf") else "inf"
        for genome in s["clade_tips"]:
            out.write(
                f"{i+1}\t{genome}\t{s['sister_tip']}\t"
                f"{s['within_mean']:.2f}\t{s['sister_mean']:.2f}\t"
                f"{ratio_str}\t{len(s['clade_tips'])}\n"
            )

print(f"Results written to: {args.out}")