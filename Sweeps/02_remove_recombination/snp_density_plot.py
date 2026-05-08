#!/usr/bin/env python3
"""
snp_density_plot.py
 
Visualize SNP density across core genome alignments before and after
recombination removal. Recombinant regions from BED files are shown
only in a separate bottom track for reference.
 
Usage:
    python3 snp_density_plot.py -r core.full.aln -c masked.aln \
        -g gubbins.bed -f cfml.bed -o snp_density.png

Dependencies:
    biopython, matplotlib, numpy
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import AlignIO

# =============================================================================
# FUNCTIONS
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot SNP density before and after recombination removal."
    )
    parser.add_argument(
        "-r", "--raw",
        required=True,
        help="Path to raw input alignment (FASTA format)"
    )
    parser.add_argument(
        "-c", "--clean",
        required=True,
        help="Path to cleaned/masked alignment (FASTA format)"
    )
    parser.add_argument(
        "-g", "--gubbins-bed",
        required=True,
        help="Path to Gubbins BED file of recombinant regions"
    )
    parser.add_argument(
        "-f", "--cfml-bed",
        required=True,
        help="Path to CFML BED file of recombinant regions"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path for output plot (e.g. snp_density.png)"
    )
    parser.add_argument(
        "-w", "--window-size",
        type=int,
        default=500,
        help="Sliding window size in bp (default: 500)"
    )
    parser.add_argument(
        "-d", "--dpi",
        type=int,
        default=200,
        help="Output image DPI (default: 200)"
    )
    return parser.parse_args()


def load_bed(path):
    """Parse a BED file and return a list of (start, end) tuples."""
    regions = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            regions.append((int(parts[1]), int(parts[2])))
    return regions


def snp_density(aln, window=500):
    """
    Count variable sites (SNPs) per sliding window across the alignment.
    Uses a fully vectorized numpy approach for speed on large alignments.
    Returns arrays of window midpoint positions and SNP counts.
    """
    # Convert alignment to a 2D numpy byte array (sequences x positions)
    # View as uint8 for fast comparisons — avoids per-column Python loops
    print("    Converting alignment to numpy matrix...")
    mat = np.frombuffer(
        b"".join(str(rec.seq).upper().encode() for rec in aln),
        dtype=np.uint8
    ).reshape(len(aln), -1)

    n = mat.shape[1]

    # Boolean mask: True where base is A, C, G, or T (exclude N, -, etc.)
    acgt = np.array([ord(b) for b in "ACGT"], dtype=np.uint8)
    is_acgt = np.isin(mat, acgt)  # shape: (n_seqs, n_positions)

    positions, counts = [], []

    for start in range(0, n, window):
        end = min(start + window, n)
        block      = mat[:, start:end]       # (n_seqs, window_size)
        block_acgt = is_acgt[:, start:end]   # valid base mask

        # For each column, find the first valid base and compare all others to it
        # A column is variable if any sequence differs from the first valid base
        col_snps = 0
        for i in range(block.shape[1]):
            col      = block[:, i]
            col_mask = block_acgt[:, i]
            valid    = col[col_mask]
            if len(valid) > 1 and not np.all(valid == valid[0]):
                col_snps += 1

        counts.append(col_snps)
        positions.append((start + end) / 2)

    return np.array(positions), np.array(counts)


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    # --- Load alignments ---
    print("Loading raw alignment...")
    raw = AlignIO.read(args.raw, "fasta")
    print("Loading clean alignment...")
    clean = AlignIO.read(args.clean, "fasta")

    aln_len = raw.get_alignment_length()
    print(f"Alignment length: {aln_len:,} bp")
    print(f"Number of sequences: {len(raw)}")

    # --- Load recombinant regions ---
    print("Loading BED files...")
    gubbins_regions = load_bed(args.gubbins_bed)
    cfml_regions    = load_bed(args.cfml_bed)
    print(f"  Gubbins regions: {len(gubbins_regions)}")
    print(f"  CFML regions:    {len(cfml_regions)}")

    # --- Compute SNP density ---
    print(f"Computing SNP density (window = {args.window_size} bp)...")
    print("  Processing raw alignment (this may take a few minutes)...")
    pos_raw,   dens_raw   = snp_density(raw,   window=args.window_size)
    print("  Processing clean alignment...")
    pos_clean, dens_clean = snp_density(clean, window=args.window_size)

    # --- Plot ---
    print("Plotting...")
    fig, axes = plt.subplots(
        3, 1,
        figsize=(16, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 3, 1]}
    )

    # Panel 1: raw alignment
    ax = axes[0]
    ax.fill_between(pos_raw, dens_raw, alpha=0.6, color="#2196F3", linewidth=0.5)
    ax.set_ylabel("SNPs per window")
    ax.set_title(f"Before recombination removal  (raw snippy alignment, window = {args.window_size} bp)")
    ax.set_xlim(0, aln_len)

    # Panel 2: clean alignment
    ax = axes[1]
    ax.fill_between(pos_clean, dens_clean, alpha=0.6, color="#4CAF50", linewidth=0.5)
    ax.set_ylabel("SNPs per window")
    ax.set_title(f"After recombination removal  (masked alignment, window = {args.window_size} bp)")
    ax.set_xlim(0, aln_len)

    # Panel 3: recombinant region track
    ax = axes[2]
    for start, end in gubbins_regions:
        ax.barh(0.5, end - start, left=start, height=0.4, color="red",    alpha=0.7)
    for start, end in cfml_regions:
        ax.barh(0.5, end - start, left=start, height=0.4, color="orange", alpha=0.7)
    ax.set_yticks([])
    ax.set_xlabel("Alignment position (bp)")
    ax.set_title("Recombinant regions")
    ax.set_xlim(0, aln_len)

    patches = [
        mpatches.Patch(color="red",    alpha=0.7, label="Gubbins regions"),
        mpatches.Patch(color="orange", alpha=0.7, label="CFML regions"),
    ]
    axes[2].legend(handles=patches, loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()