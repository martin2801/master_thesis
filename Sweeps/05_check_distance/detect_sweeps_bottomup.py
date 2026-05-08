#!/usr/bin/env python3
"""
detect_sweeps_bottomup.py
Identifies selective sweeps de novo from a phylogenetic tree and a pairwise
SNP distance list, using a bottom-up traversal strategy.

Background:
    A selective sweep occurs when positive selection drives a single haplotype
    to high frequency, reducing genetic diversity within the affected lineage.
    This creates a detectable pattern: genomes within a swept clade are nearly
    identical to each other, but clearly distinct from all outside genomes.
    This script detects that pattern using the 5x rule (modified from Birky
    et al. 2005): a clade is called a sweep if the mean SNP distance to the
    nearest outside genome is at least N times greater than the mean pairwise
    SNP distance within the clade (default N=5).

Strategy:
    Unlike a top-down approach (which would stop at the first node passing the
    threshold), this script tests ALL internal nodes bottom-up (smallest clades
    first). This ensures the tightest, most specific sweep boundary is found
    before larger encompassing clades are considered. After all candidates are
    identified, a greedy mutual exclusivity step assigns each genome to at most
    one sweep, with the highest-ratio candidate winning any conflict.

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
    5. Sort candidates by ratio descending (strongest signal first)
    6. Greedily assign sweeps, skipping candidates whose genomes overlap
       with an already-called sweep (mutual exclusivity)
    7. Write results to a tab-separated output file

Usage:
    python3 detect_sweeps_bottomup.py \\
        --tree  full_tree.treefile \\
        --snps  snp_list.tab \\
        --threshold 5 \\
        --min-tips  3 \\
        --out   sweeps.txt

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
      the strongest possible sweep signal
    - The tree is midpoint-rooted before traversal; the input treefile is not modified
    - Sweeps are mutually exclusive: each genome appears in at most one sweep
    - The greedy assignment prioritizes highest ratio, so the strongest signal
      wins when nested clades compete for the same genomes

Dependencies:
    - Python 3.6+
    - biopython (pip install biopython)

Example:
    snp-dists -m alignment.fasta > snp_list.tab
    python3 detect_sweeps_bottomup.py \\
        --tree  full_tree.treefile \\
        --snps  snp_list.tab \\
        --threshold 5 \\
        --min-tips  3 \\
        --out   sweeps_5x.txt
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
        # Expecting lines of the form: sample1, sample2, distance
        # Skip malformed lines
        if len(parts) < 3:
            continue
        a, b, d = parts[0], parts[1], int(parts[2])
        # Avoid self-comparisons
        if a == b:
            continue
        # Store distances in a way that (a, b) and (b, a) are the same key
        # No matter what orientation the pair appears int he file, it will be stored under the same key
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
    return snp.get((min(a, b), max(a, b)), None) # uses the same key format as stored to get the distance, None if not found

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
# Root at midpoint, the point to minimize the max dist from root to any tip
tree.root_at_midpoint()

# --- Build parent lookup ---
# BioPython has no built-in get_parent() method, so we construct a dictionary
# mapping id(child) -> parent clade by traversing the tree top-down.
# id() is Python's unique memory address for each object, used as the key
# since clade objects cannot be used directly as dictionary keys.
parent = {}
for clade in tree.find_clades(order="level"):
    for child in clade.clades:
        parent[id(child)] = clade # id(): Py unique memory adress for the object

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
        return None
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
    sister_tip_names = [c.name for c in sister_clade.get_terminals()] # get all the tips in the sister clade
    best_tip, best_mean = None, float("inf")          # initialize best mean to infinity so any real mean will be lower
    for tip in sister_tip_names: 
        dists = [get_snp(t, tip) for t in clade_tips] # get the SNP distances from each genome in the sweep clade to that tip -> list of dist
        dists = [d for d in dists if d is not None]   # filter out any missing distances (None)
        if not dists:
            continue                                  # if no valid distances, skip this tip
        m = sum(dists) / len(dists)                   # compute mean distance from all clade members to this sister tip  
        if m < best_mean:
            best_mean, best_tip = m, tip              # if this is the lowest mean so far, update best_mean and best_tip
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
# we want to test smaller clades first, so we reverse the level-order traversal which processes nodes from root to tips
all_clades = list(tree.find_clades(order="level"))
for clade in reversed(all_clades):                  # here the order gets reversed
    if clade.is_terminal():
        continue                                    # skip tips, we only want to test internal nodes

    tips = [c.name for c in clade.get_terminals()]  # skip small clades that don't meet the minimum tip requirement
    if len(tips) < args.min_tips:
        continue                                   

    sister = get_sister(clade)                      # find the sister clade to compare against. If no sister (e.g. root node), skip
    if sister is None:
        continue

    sister_tip, sister_mean = get_closest_sister_tip(tips, sister)
    if sister_tip is None:
        continue                                    # skip clades with no sister (i.e.root)

    # Compute within-clade mean SNP distance
    # itertools.combinations generates all unique pairs without repetition
    # e.g. for [A, B, C]: (A,B), (A,C), (B,C) — no (A,A) or (B,A) duplicates
    within_vals = [get_snp(a, b) for a, b in itertools.combinations(tips, 2) if a != b] # get all pairwise combinations of tips in the clade, compute their SNP distances, and filter out any missing distances (None)
    within_vals = [v for v in within_vals if v is not None]
    if not within_vals:
        continue

    within_mean = sum(within_vals) / len(within_vals)

    # A fully clonal clade (all genomes identical) has within_mean == 0.
    # This is the strongest possible sweep signal, so ratio is set to inf
    # rather than crashing with a division by zero error.
    if within_mean == 0:            # 0 if clade is fully clonal
        ratio = float("inf")
    else:
        ratio = sister_mean / within_mean

    if ratio >= args.threshold:     # if this clade passes the threshold, save it as a candidate sweep
        candidates.append({
            "clade_tips":   tips,
            "tip_set":      set(tips),
            "sister_tip":   sister_tip,
            "sister_mean":  sister_mean,
            "within_mean":  within_mean,
            "within_pairs": len(within_vals),
            "ratio":        ratio,
        })

# --- Greedy assignment: highest ratio first, enforce mutual exclusivity ---
# All candidates that passed the threshold are sorted by ratio descending.
# We then iterate through them and assign each as a sweep only if none of
# its genomes have already been assigned to a previously called sweep.
# The & operator checks for set intersection — if any tips overlap, skip.
# The |= operator adds all tips from the new sweep to the assigned set.
# This ensures each genome appears in at most one sweep, and the strongest
# signal (highest ratio) wins any conflict between overlapping candidates.
candidates.sort(key=lambda x: x["ratio"], reverse=True)

assigned_tips = set()
sweeps = []

for c in candidates:
    # skip if any tip in this clade is already assigned to a sweep
    if c["tip_set"] & assigned_tips:
        continue
    sweeps.append(c)
    assigned_tips |= c["tip_set"]

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