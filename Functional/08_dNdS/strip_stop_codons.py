#!/usr/bin/env python3
"""
Strip the terminal codon (last 3 nt) from every sequence in every per-gene
FASTA alignment. GFF3 CDS coordinates include the stop codon, but HyPhy
(and most codon-model tools) require it removed before fitting - the
substitution model only covers the 61 sense codons.

Since all genomes were sliced using the same reference-derived coordinates,
the stop codon sits at the same trailing position for every sequence in a
given gene file, so a uniform 3-nt trim is safe and consistent across the
whole file.

Writes trimmed output to a NEW directory (does not overwrite originals),
so you can compare / re-run cleanly if needed.

Usage:
    python3 strip_stop_codons.py \
        --in-dir per_gene_alignments \
        --out-dir per_gene_alignments_nostop
"""

import argparse
import glob
import os

STOP_CODONS = {"TAA", "TAG", "TGA"}


def process_fasta(in_path, out_path):
    """
    Trim the last 3 nt from every sequence. Logs (but does not block on)
    cases where the trimmed codon wasn't actually a canonical stop -
    e.g. genomes with masked N's at the gene's very end - since those are
    still safe to trim positionally, just worth knowing about.
    """
    records = []  # (header, seq)
    with open(in_path) as fh:
        header = None
        seq_lines = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "".join(seq_lines)))

    non_stop_trimmed = []
    with open(out_path, "w") as out_fh:
        for header, seq in records:
            last_codon = seq[-3:].upper()
            trimmed_seq = seq[:-3]
            out_fh.write(f">{header}\n{trimmed_seq}\n")
            if last_codon not in STOP_CODONS and "N" not in last_codon and "-" not in last_codon:
                non_stop_trimmed.append((header, last_codon))

    return non_stop_trimmed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    fasta_files = sorted(glob.glob(os.path.join(args.in_dir, "*.fasta")))
    print(f"Processing {len(fasta_files)} gene alignments...")

    total_anomalies = 0
    anomaly_log_path = os.path.join(args.out_dir, "_non_stop_trims.tsv")
    with open(anomaly_log_path, "w") as log_fh:
        log_fh.write("gene_id\tgenome_id\ttrimmed_codon\n")
        for fasta_path in fasta_files:
            gene_id = os.path.splitext(os.path.basename(fasta_path))[0]
            out_path = os.path.join(args.out_dir, os.path.basename(fasta_path))
            anomalies = process_fasta(fasta_path, out_path)
            for genome_id, codon in anomalies:
                log_fh.write(f"{gene_id}\t{genome_id}\t{codon}\n")
                total_anomalies += 1

    print(f"Done. Trimmed alignments written to {args.out_dir}/")
    print(f"{total_anomalies} genome/gene combos had a non-stop, non-N/gap codon trimmed "
          f"(logged in {anomaly_log_path}) - worth a quick look if this number is large.")


if __name__ == "__main__":
    main()