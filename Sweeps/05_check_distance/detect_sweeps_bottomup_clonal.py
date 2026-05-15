#!/usr/bin/env python3

"""
detect_sweeps_bottomup_clonal.py
Identifies selective sweeps de novo from a phylogenetic tree and a pairwise
SNP distance list, using a bottom-up traversal strategy with clonality-first
conflict resolution.

Background:
    A selective sweep occurs when positive selection drives a single haplotype
    to high frequency, reducing genetic diversity within the affected lineage.
    This creates a detectable pattern: genomes within a swept clade are nearly
    identical to each other, but clearly distinct from all outside genomes.
    This script detects that pattern using the 5x rule (modified from Birky
    et al. 2005): a clade is called a sweep if the mean SNP distance to the
    nearest outside genome is at least N times greater than the mean pairwise
    SNP distance within the clade (default N=5).

Difference from detect_sweeps_bottomup.py:
    This script is identical to detect_sweeps_bottomup.py in every respect
    except one: the priority rule used when two overlapping candidate clades
    compete for the same genomes in the greedy assignment step.

    detect_sweeps_bottomup.py    → highest ratio wins
    detect_sweeps_bottomup_clonal.py → lowest within_mean wins (this script)

    Prioritizing lowest within_mean rewards internal clonality directly. When
    a tight sub-clade (e.g. 3 genomes, within_mean=1.6 SNPs) is nested inside
    a larger clade (e.g. 6 genomes, within_mean=14 SNPs), the tighter sub-clade
    wins even if the larger clade has a higher ratio. This tends to produce
    smaller, more precisely bounded sweeps at the cost of potentially
    over-splitting larger biologically meaningful sweep groups.

    Use this version when you want to prioritize finding the most clonal,
    tightest possible sweep boundaries. Use detect_sweeps_bottomup.py when
    you want the strongest overall signal (ratio) to define the boundary.

Algorithm:
    1. Load all pairwise SNP distances into a dictionary for O(1) lookup
    2. Load the phylogenetic tree and root it at the midpoint
    3. Traverse all internal nodes bottom-up (leaves first, root last)
    4. For each node:
         a. Get all tips in the clade
         b. Find the immediate sister clade from tree topology
         c. Pick the single closest tip in the sister clade by mean SNP distance
         d. Compute mean pairwise SNP distance within the clade
         e. Compute mean SNP distance from clade members to the sister tip
         f. Compute ratio = sister_mean / within_mean
         g. If ratio >= threshold, save as a candidate
    5. Sort candidates by within_mean ascending (most clonal first)
    6. Greedily assign sweeps, skipping candidates whose genomes overlap
       with an already-called sweep (mutual exclusivity)
    7. Re-sort final sweeps by ratio descending for output
    8. Write results to a tab-separated output file

Usage:
    python3 detect_sweeps_bottomup_clonal.py \\
        --tree  full_tree.treefile \\
        --snps  snp_list.tab \\
        --threshold 5 \\
        --min-tips  3 \\
        --out   sweeps_clonal.txt

Arguments:
    --tree        Newick tree file, e.g. IQ-TREE .treefile output
    --snps        Tab-separated pairwise SNP distance list (from snp-dists -m):
                  sample1 <tab> sample2 <tab> distance
    --threshold   Ratio threshold to call a sweep (default: 5.0)
                  Based on Birky et al. 2005: a 4x ratio corresponds to 95%
                  confidence of independent evolution under neutral coalescent
    --min-tips    Minimum number of genomes a clade must contain to be tested
                  (default: 3, filtering out pairs which may be spurious)
    --out         Output file path (default: sweeps.txt)

Output:
    Tab-separated file with one row per genome per sweep:
    sweep_id | genome | sister_tip | within_mean | sister_mean | ratio | clade_size

    sweep_id    : Integer sweep identifier, ranked by ratio descending
    genome      : Genome name
    sister_tip  : Name of the closest sister genome outside the clade
    within_mean : Mean pairwise SNP distance within the sweep clade
    sister_mean : Mean SNP distance from clade members to the sister tip
    ratio       : sister_mean / within_mean
    clade_size  : Number of genomes in the sweep clade

Notes:
    - Self-comparisons (genome vs itself) are excluded from all calculations
    - Pairs missing from the SNP list are treated as None and excluded from means
    - Fully clonal clades (within_mean == 0) receive a ratio of inf, representing
      the strongest possible sweep signal, and will always win conflicts
    - The tree is midpoint-rooted before traversal; the input treefile is not modified
    - Sweeps are mutually exclusive: each genome appears in at most one sweep
    - The greedy assignment prioritizes lowest within_mean, so the most clonal
      clade wins when nested clades compete for the same genomes
    - Final output is sorted by ratio descending regardless of assignment order

Dependencies:
    - Python 3.6+
    - biopython (pip install biopython)

Example:
    snp-dists -m alignment.fasta > snp_list.tab
    python3 detect_sweeps_bottomup_clonal.py \\
        --tree  full_tree.treefile \\
        --snps  snp_list.tab \\
        --threshold 5 \\
        --min-tips  3 \\
        --out   sweeps_clonal_5x.txt
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
# All pairwise SNP distances are loaded into a dictionary for O(1) lookup.
# The key is always stored as (smaller_name, larger_name) using alphabetical
# ordering — this ensures the same pair is stored under the same key regardless
# of which orientation (A,B) or (B,A) it appears in the file, so lookups never
# need to check both orientations.
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
    """
    Retrieve the SNP distance between two genomes from the lookup dictionary.

    Uses the same (min, max) key format as the storage step, so orientation
    of the pair does not matter. Returns None if the pair is not found in the
    SNP list, which can happen if a genome was absent from the snp-dists run.

    Args:
        a (str): Name of the first genome
        b (str): Name of the second genome

    Returns:
        int or None: SNP distance, or None if not found or if a == b
    """
    if a == b:
        return None
    return snp.get((min(a, b), max(a, b)), None)

def safe_mean(vals):
    """
    Compute the mean of a list while ignoring None values.

    Used wherever distances may be missing due to gaps in the SNP list.
    Returns None if the list is empty or contains only None values.

    Args:
        vals (list): List of numeric values or None

    Returns:
        float or None: Mean of non-None values, or None if no valid values
    """
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None

# --- Load and midpoint-root tree ---
# IQ-TREE outputs unrooted trees. Midpoint rooting places the root at the
# point that minimizes the maximum distance from the root to any tip.
# This is required because parent-child relationships (used for sister
# detection) are only defined on a rooted tree.
print("Loading tree...")
tree = Phylo.read(args.tree, "newick")
tree.root_at_midpoint()

# --- Build parent lookup ---
# BioPython has no built-in get_parent() method, so we construct a dictionary
# mapping id(child) -> parent clade by traversing the tree top-down.
# id() is Python's unique memory address for each object, used as the key
# since clade objects cannot be used directly as dictionary keys.
parent = {}
for clade in tree.find_clades(order="level"):
    for child in clade.clades:
        parent[id(child)] = clade

def get_sister(clade):
    """
    Return the sister clade of a given node (the other child of its parent).

    The sister is defined purely by tree topology — it is whatever clade
    shares the same parent node. This is used to find the reference group
    for comparison rather than using SNP distance to define the sister.
    Returns None if the clade is the root (no parent exists).

    Args:
        clade: BioPython clade object

    Returns:
        BioPython clade object or None
    """
    p = parent.get(id(clade))
    if p is None:
        return None     # root node has no sister
    for child in p.clades:
        if child is not clade:
            return child
    return None

def get_closest_sister_tip(clade_tips, sister_clade):
    """
    From all tips in the sister clade, return the single tip with the lowest
    mean SNP distance to the clade members.

    The sister clade may contain one genome or hundreds. Rather than averaging
    across all of them (which would dilute the signal), we pick the single
    closest tip as the representative sister. This mirrors the manual approach
    of choosing one known sister genome for comparison, but automated via SNP
    distances.

    For each candidate sister tip, the mean is computed across all clade
    members — not just the closest one — so that the result reflects how close
    that tip is to the whole clade, not just one member.

    Args:
        clade_tips (list of str): Genome names in the sweep candidate clade
        sister_clade: BioPython clade object representing the sister

    Returns:
        tuple: (best_tip_name, best_mean_distance)
               Returns (None, inf) if no valid distances are found
    """
    sister_tip_names = [c.name for c in sister_clade.get_terminals()]
    best_tip, best_mean = None, float("inf")  # inf ensures any real distance replaces it
    for tip in sister_tip_names:
        dists = [get_snp(t, tip) for t in clade_tips]  # distances from each clade member to this tip
        dists = [d for d in dists if d is not None]     # remove missing distances
        if not dists:
            continue        # skip this tip entirely if no distances are available
        m = sum(dists) / len(dists)
        if m < best_mean:
            best_mean, best_tip = m, tip  # update if this is the closest tip so far
    return best_tip, best_mean

# --- Test ALL internal nodes bottom-up ---
# find_clades(order="level") traverses top-down (root first).
# Reversing gives bottom-up order so small, tight clades are tested before
# the larger clades that contain them. This ensures the tightest sweep
# boundary is identified before broader encompassing clades compete for
# the same genomes in the greedy assignment step.
print(f"Scanning tree bottom-up (min clade size: {args.min_tips}, threshold: {args.threshold}x)...")

candidates = []  # all nodes that pass the threshold

# bottom-up = reverse of level-order
all_clades = list(tree.find_clades(order="level"))
for clade in reversed(all_clades):
    if clade.is_terminal():
        continue  # tips have no children, nothing to test

    tips = [c.name for c in clade.get_terminals()]
    if len(tips) < args.min_tips:
        continue  # clade too small

    sister = get_sister(clade)
    if sister is None:
        continue  # root node, no sister exists

    sister_tip, sister_mean = get_closest_sister_tip(tips, sister)
    if sister_tip is None:
        continue  # no valid SNP distances to any sister tip

    # --- Sister merging step ---
    # If the closest sister tip has a mean SNP distance of 0 to the clade,
    # it is identical to the clade members. This is a tree resolution artefact:
    # with zero genetic distance there is no phylogenetic signal to place the
    # genome inside or outside the clade, so its placement is arbitrary.
    # We resolve this by merging the sister tip into the clade and climbing
    # one level up the tree to find the next distinct sister. We keep climbing
    # as long as the new sister tip is also at distance 0, stopping when we
    # find a genuinely distinct outside genome or reach the root.
    merged_tips = list(tips)
    current_clade = clade
    while sister_mean == 0:
        merged_tips.append(sister_tip)                     # absorb the identical sister into the clade
        current_clade = parent.get(id(current_clade))     # climb one level up
        if current_clade is None:
            break                                          # reached the root, no further sister exists
        new_sister = get_sister(current_clade)
        if new_sister is None:
            break                                          # no sister at this level either
        sister_tip, sister_mean = get_closest_sister_tip(merged_tips, new_sister)
        if sister_tip is None:
            break                                          # no valid distances to any tip in the new sister
    tips = merged_tips                                     # use the expanded tip list going forward

    # itertools.combinations generates all unique pairs without repetition
    # e.g. for [A, B, C]: (A,B), (A,C), (B,C) — no (A,A) or (B,A) duplicates
    within_vals = [get_snp(a, b) for a, b in itertools.combinations(tips, 2) if a != b]
    within_vals = [v for v in within_vals if v is not None]
    if not within_vals:
        continue

    within_mean = sum(within_vals) / len(within_vals)

    # A fully clonal clade (all genomes identical) has within_mean == 0.
    # This is the strongest possible sweep signal, so ratio is set to inf
    # rather than crashing with a division by zero error. Fully clonal clades
    # will always sort first under the clonality-first priority rule.
    if within_mean == 0:
        ratio = float("inf")
    else:
        ratio = sister_mean / within_mean

    if ratio >= args.threshold:
        candidates.append({
            "clade_tips":   tips,
            "tip_set":      set(tips), # set for fast intersection testing in greedy step
            "sister_tip":   sister_tip,
            "sister_mean":  sister_mean,
            "within_mean":  within_mean,
            "within_pairs": len(within_vals),
            "ratio":        ratio,
        })

# --- Greedy assignment: lowest within_mean first, enforce mutual exclusivity ---
# This is the only difference from detect_sweeps_bottomup.py.
# Instead of sorting by ratio descending (strongest signal wins), we sort by
# within_mean ascending (most clonal clade wins). When a tight sub-clade and
# a larger parent clade both pass the threshold and overlap, the sub-clade
# wins if it has lower internal diversity — even if the parent has a higher ratio.
# Note: no reverse=True needed since ascending is Python's default sort order.

# candidates.sort(key=lambda x: x["ratio"], reverse=True)
candidates.sort(key=lambda x: x["within_mean"])

assigned_tips = set()
sweeps = []

for c in candidates:
    if c["tip_set"] & assigned_tips:  # overlap with already-assigned genomes
        continue
    sweeps.append(c)
    assigned_tips |= c["tip_set"]     # mark all genomes in this sweep as assigned

# Re-sort by ratio descending for final output so sweep IDs reflect signal strength
# regardless of which priority rule was used during assignment
sweeps.sort(key=lambda x: x["ratio"], reverse=True)

# re-sort by ratio descending for final output (already sorted, but explicit)
sweeps.sort(key=lambda x: x["ratio"], reverse=True)

# --- Print summary ---
print(f"\nCandidates passing threshold: {len(candidates)}")
print(f"Sweeps after mutual exclusivity: {len(sweeps)}\n")

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

# --- Write output table in long format ---
# One row per genome per sweep, with sweep-level statistics repeated on every
# row. This makes it easy to load into R or join against metadata tables
# since every genome has its sweep ID and statistics directly attached.
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