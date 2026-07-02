#!/usr/bin/env python3
"""
Slice a masked whole-genome core alignment (snippy-core, Gubbins-masked)
into per-gene codon alignments, using CDS coordinates from a Prokka GFF
that are LOCAL to each contig, converted into GLOBAL alignment coordinates
via a .fai index (which reflects the same contig order snippy/snippy-core
used to build the concatenated 'Reference' alignment record).

Assumes:
  - alignment is one concatenated sequence per genome (snippy-core core.full.aln,
    Gubbins-masked), collinear with the reference, contigs in .fai order,
    no columns dropped (only masked in place with N/-).
  - GFF coordinates are 1-based inclusive, local to each contig.

Usage:
    python3 slice_genes.py \
        --aln masked.aln \
        --gff ref.gff \
        --fai reference_genome.fna.fai \
        --outdir per_gene_alignments
"""

import argparse
import csv
import os
from Bio import SeqIO
from Bio.Seq import Seq


def load_fai_offsets(fai_path):
    """
    Read .fai (NAME LENGTH OFFSET LINEBASES LINEWIDTH), return:
      offsets: dict contig_name -> 0-based cumulative offset in the
               concatenated alignment (i.e. global_pos = offset + local_pos)
      lengths: dict contig_name -> contig length
    Order in the file = order contigs were concatenated.
    """
    offsets = {}
    lengths = {}
    cumulative = 0
    with open(fai_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            name, length = fields[0], int(fields[1])
            offsets[name] = cumulative
            lengths[name] = length
            cumulative += length
    return offsets, lengths


def parse_gff(gff_path, feature_type="CDS"):
    """Parse GFF3 CDS features. Returns list of dicts with LOCAL (per-contig) coords."""
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            contig, source, ftype, start, end, score, strand, frame, attrs = fields
            if ftype != feature_type:
                continue
            attr_dict = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v
            gene_id = attr_dict.get("locus_tag") or attr_dict.get("ID") or f"{contig}:{start}-{end}"
            genes.append({
                "gene_id": gene_id,
                "contig": contig,
                "start_local": int(start),  # 1-based inclusive, local to contig
                "end_local": int(end),
                "strand": strand,
            })
    return genes


def to_global_coords(gene, offsets, lengths):
    """
    Convert local (contig-relative, 1-based inclusive) coords to global
    0-based half-open slice coordinates into the concatenated alignment.
    Returns (global_start_0based, global_end_exclusive) or None if contig
    not found (e.g. dropped by Prokka's mincontiglen filter) or out of range.
    """
    contig = gene["contig"]
    if contig not in offsets:
        return None
    contig_len = lengths[contig]
    if gene["end_local"] > contig_len or gene["start_local"] < 1:
        return None
    global_start = offsets[contig] + (gene["start_local"] - 1)  # 0-based
    global_end = offsets[contig] + gene["end_local"]             # 0-based exclusive
    return global_start, global_end


def load_alignment(aln_path):
    records = {}
    for rec in SeqIO.parse(aln_path, "fasta"):
        records[rec.id] = rec
    return records


def slice_gene(records, global_start, global_end, strand):
    out = {}
    for genome_id, rec in records.items():
        seq = rec.seq[global_start:global_end]
        if strand == "-":
            seq = seq.reverse_complement()
        out[genome_id] = str(seq)
    return out


def check_gene(seq_dict):
    issues = []
    lengths = {len(s) for s in seq_dict.values()}
    if len(lengths) != 1:
        issues.append(f"unequal lengths across genomes: {lengths}")
        return False, {"issues": issues, "length": None, "n_genomes": len(seq_dict),
                        "mean_missing_frac": None, "max_missing_frac": None,
                        "genomes_with_premature_stop": None}

    length = lengths.pop()
    frame_ok = (length % 3 == 0)
    if not frame_ok:
        issues.append(f"length {length} not multiple of 3")

    n_frac_per_genome = {}
    stop_hits = {}
    for genome_id, seq in seq_dict.items():
        n_count = seq.upper().count("N") + seq.count("-")
        n_frac_per_genome[genome_id] = n_count / max(len(seq), 1)
        if frame_ok and len(seq) >= 3:
            try:
                protein = str(Seq(seq).translate(table=11, to_stop=False))
                premature = protein[:-1].count("*") if protein else 0
                if premature > 0:
                    stop_hits[genome_id] = premature
            except Exception as e:
                issues.append(f"{genome_id}: translation error {e}")

    mean_n = sum(n_frac_per_genome.values()) / len(n_frac_per_genome) if n_frac_per_genome else 1.0
    max_n = max(n_frac_per_genome.values()) if n_frac_per_genome else 1.0

    report = {
        "length": length,
        "n_genomes": len(seq_dict),
        "mean_missing_frac": round(mean_n, 4),
        "max_missing_frac": round(max_n, 4),
        "genomes_with_premature_stop": len(stop_hits),
        "issues": issues,
    }
    passed = frame_ok and (len(stop_hits) == 0) and (mean_n < 0.5)
    return passed, report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True, help="Masked core genome alignment (FASTA, one concatenated seq/genome)")
    ap.add_argument("--gff", required=True, help="Prokka GFF3 with CDS features (local per-contig coords)")
    ap.add_argument("--fai", required=True, help=".fai index of the reference FASTA (defines contig order/offsets)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading .fai offsets...")
    offsets, lengths = load_fai_offsets(args.fai)
    print(f"  {len(offsets)} contigs, total length {sum(lengths.values())}")

    print("Loading alignment...")
    records = load_alignment(args.aln)
    aln_len = len(next(iter(records.values())).seq)
    print(f"  {len(records)} genomes, alignment length {aln_len}")
    if aln_len != sum(lengths.values()):
        print(f"  WARNING: alignment length ({aln_len}) != sum of contig lengths "
              f"({sum(lengths.values())}). Coordinates may be wrong — verify before proceeding.")

    print("Parsing GFF...")
    genes = parse_gff(args.gff)
    print(f"  {len(genes)} CDS features found")

    qc_rows = []
    n_pass = 0
    n_skipped_contig = 0
    for gene in genes:
        coords = to_global_coords(gene, offsets, lengths)
        if coords is None:
            n_skipped_contig += 1
            qc_rows.append({
                "gene_id": gene["gene_id"], "contig": gene["contig"],
                "start_local": gene["start_local"], "end_local": gene["end_local"],
                "strand": gene["strand"], "length": "", "n_genomes": "",
                "mean_missing_frac": "", "max_missing_frac": "",
                "genomes_with_premature_stop": "", "passed": False,
                "issues": "contig not found in .fai (likely dropped by Prokka mincontiglen) or out of range",
            })
            continue

        global_start, global_end = coords
        seq_dict = slice_gene(records, global_start, global_end, gene["strand"])
        passed, report = check_gene(seq_dict)

        row = {
            "gene_id": gene["gene_id"], "contig": gene["contig"],
            "start_local": gene["start_local"], "end_local": gene["end_local"],
            "strand": gene["strand"],
            "length": report["length"], "n_genomes": report["n_genomes"],
            "mean_missing_frac": report["mean_missing_frac"],
            "max_missing_frac": report["max_missing_frac"],
            "genomes_with_premature_stop": report["genomes_with_premature_stop"],
            "passed": passed,
            "issues": ";".join(report["issues"]),
        }
        qc_rows.append(row)

        if passed:
            n_pass += 1
            out_path = os.path.join(args.outdir, f"{gene['gene_id']}.fasta")
            with open(out_path, "w") as fh:
                for genome_id, seq in seq_dict.items():
                    fh.write(f">{genome_id}\n{seq}\n")

    qc_path = os.path.join(args.outdir, "_qc_summary.tsv")
    fieldnames = ["gene_id", "contig", "start_local", "end_local", "strand",
                  "length", "n_genomes", "mean_missing_frac", "max_missing_frac",
                  "genomes_with_premature_stop", "passed", "issues"]
    with open(qc_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in qc_rows:
            writer.writerow(row)

    print(f"\nDone. {n_pass}/{len(genes)} genes passed QC and were written to {args.outdir}/")
    print(f"  {n_skipped_contig} genes skipped (contig missing from .fai / out of range)")
    print(f"QC summary: {qc_path}")


if __name__ == "__main__":
    main()