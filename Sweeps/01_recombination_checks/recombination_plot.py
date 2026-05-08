#!/usr/bin/env python3
"""
recombination_plot.py
=====================
Visualise and compare recombinant genomic regions predicted by Gubbins and
ClonalFrameML for a Salmonella Infantis population.

Usage
-----
    python3 recombination_plot.py -g gubbins_fixed.bed -f cfml_fixed.bed \
        [-o recombination_comparison.png] [-l 4900000] [--dpi 300]

Dependencies
------------
    pandas, matplotlib
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot and compare recombinant regions from Gubbins and ClonalFrameML."
    )
    parser.add_argument(
        "-g", "--gubbins",
        required=True,
        metavar="FILE",
        help="Gubbins recombinant regions BED file"
    )
    parser.add_argument(
        "-f", "--cfml",
        required=True,
        metavar="FILE",
        help="ClonalFrameML recombinant regions BED file"
    )
    parser.add_argument(
        "-o", "--output",
        default="recombination_comparison.png",
        metavar="FILE",
        help="Output image file (default: recombination_comparison.png)"
    )
    parser.add_argument(
        "-l", "--genome-length",
        type=int,
        default=4_900_000,
        metavar="INT",
        help="Approximate chromosome length in bp for x-axis limits "
             "(default: 4900000 for GCF_000506925.1 S. Infantis SI119944)"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        metavar="INT",
        help="Output image resolution in DPI (default: 300)"
    )
    return parser.parse_args()


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    # --- Load BED files ---
    gubbins = pd.read_csv(args.gubbins, sep="\t", header=None,
                          names=["seq", "start", "end", "info"])
    cfml    = pd.read_csv(args.cfml,    sep="\t", header=None,
                          names=["seq", "start", "end", "info"])

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(16, 3))

    for _, r in gubbins.iterrows():
        ax.barh(1, r.end - r.start, left=r.start, height=0.5,
                color="steelblue", alpha=0.5, linewidth=0)

    for _, r in cfml.iterrows():
        ax.barh(0, r.end - r.start, left=r.start, height=0.5,
                color="coral", alpha=0.5, linewidth=0)

    # --- Axes formatting ---
    ax.set_xlim(0, args.genome_length)
    ax.set_ylim(-0.5, 1.8)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["ClonalFrameML", "Gubbins"], fontsize=11)
    ax.set_xlabel("Genome position (bp)", fontsize=11)
    ax.set_title(
        "Recombinant regions — Gubbins vs ClonalFrameML\n"
        "(Salmonella Infantis, Snippy alignment)",
        fontsize=12
    )
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.1f} Mbp"))

    # --- Legend ---
    patches = [mpatches.Patch(color="steelblue", alpha=0.5, label="Gubbins"),
               mpatches.Patch(color="coral",     alpha=0.5, label="ClonalFrameML")]
    ax.legend(handles=patches, loc="upper right", fontsize=10)

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()