#!/usr/bin/env python3
"""
mask_recombination.py

Mask recombinant regions across ALL sequences in a multiple sequence alignment
by replacing bases at BED-specified coordinates with Ns.

Unlike bedtools maskfasta, this script operates on every sequence in the
alignment simultaneously, which is the correct behaviour for multi-FASTA
core genome alignments.

Usage:
    python3 mask_recombination.py -i core.full.aln -b gubbins.bed -o masked.aln

Dependencies:
    biopython
"""

import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import MutableSeq
import subprocess

# =============================================================================
# FUNCTIONS
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Mask recombinant regions in a core genome alignment using a BED file."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input alignment (FASTA format)"
    )
    parser.add_argument(
        "-b", "--bed",
        required=True,
        help="Path to Gubbins BED file of recombinant regions"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path for output masked alignment (FASTA format)"
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

def count_ns_in_file(path):
    """Count total Ns in a FASTA file using grep/tr (fast, avoids re-parsing)."""
    result = subprocess.run(
        f"grep -v '>' {path} | tr -cd 'Nn' | wc -c",
        shell=True, capture_output=True, text=True
    )
    return int(result.stdout.strip())
 
 
def merged_bed_coverage(path):
    """Return total bp covered by BED file after merging overlapping regions."""
    result = subprocess.run(
        f"bedtools merge -i {path} | awk '{{sum += $3 - $2}} END{{print sum}}'",
        shell=True, capture_output=True, text=True
    )
    return int(result.stdout.strip())

# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    # --- Load BED regions ---
    print("Loading BED file...")
    regions = load_bed(args.bed)
    total_region_bp = sum(e - s for s, e in regions)
    print(f"  Recombinant regions: {len(regions)}")
    print(f"  Total region bp (raw, overlaps included): {total_region_bp:,}")

    # --- Load alignment ---
    print("Loading alignment...")
    aln = AlignIO.read(args.input, "fasta")
    n_seqs  = len(aln)
    aln_len = aln.get_alignment_length()
    print(f"  Sequences:        {n_seqs}")
    print(f"  Alignment length: {aln_len:,} bp")

    # --- Mask every sequence at every recombinant region ---
    print("Masking recombinant regions across all sequences...")
    records = []
    for i, rec in enumerate(aln):
        if i % 10 == 0:
            print(f"  Processing sequence {i + 1}/{n_seqs}: {rec.id}")
        seq = MutableSeq(str(rec.seq).upper())
        for start, end in regions:
            for pos in range(start, min(end, aln_len)):
                seq[pos] = 'N'
        rec.seq = seq
        records.append(rec)

    # --- Write output ---
    print("Writing masked alignment...")
    SeqIO.write(records, args.output, "fasta")
    print(f"Saved: {args.output}")

    # --- Sanity check ---
    print("\n=== Alignment Info ===")
    print(f"Alignment length:    {aln_len:,} bp")
    print(f"Number of sequences: {n_seqs}")
 
    print("\n=== N Count Comparison ===")
    raw_n   = count_ns_in_file(args.input)
    clean_n = count_ns_in_file(args.output)
    new_ns  = clean_n - raw_n
    print(f"Ns in raw alignment:   {raw_n:,}")
    print(f"Ns in clean alignment: {clean_n:,}")
    print(f"New Ns added:          {new_ns:,}")
 
    print("\n=== BED File Verification ===")
    merged_bases = merged_bed_coverage(args.bed)
    expected_new_ns = merged_bases * n_seqs
    print(f"Merged BED covered bases: {merged_bases:,} bp  (overlaps removed)")
    print(f"Expected new Ns (merged x {n_seqs} sequences): {expected_new_ns:,}")
 
    print("\n=== Sanity Check ===")
    diff = abs(new_ns - expected_new_ns)
    pct_diff = 100 * diff / expected_new_ns if expected_new_ns > 0 else 0
    if pct_diff < 5:
        print(f"PASS: New Ns ({new_ns:,}) are within 5% of expected ({expected_new_ns:,}).")
        print(f"      Difference of {diff:,} bp is explained by positions already N in the raw alignment.")
    else:
        print(f"WARN: New Ns ({new_ns:,}) differ from expected ({expected_new_ns:,}) by {pct_diff:.1f}%.")
        print(f"      Difference: {diff:,} bp — investigate further.")

if __name__ == "__main__":
    main()