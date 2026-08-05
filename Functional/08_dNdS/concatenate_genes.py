#!/usr/bin/env python3
"""
concatenate_genes.py

Concatenate a directory of per-gene, in-frame, stop-codon-free alignments
(the output of strip_stop_codons.py in the dN/dS pipeline) into a single
supermatrix, so RELAX can be run once genome-wide instead of once per gene.

Per-gene RELAX runs on this dataset are underpowered: branch lengths within
a single ~100-400bp gene are close to zero substitutions on most branches
(confirmed directly from a sample JSON: 205/206 Test-clade branches had an
estimated codon-model branch length of exactly 0). Concatenating genes pools
enough substitutions across the genome to actually estimate a rate.

Only taxa present in every included gene are kept (strict intersection) by
default, so the resulting alignment is a clean rectangular matrix. Genes
whose aligned length isn't a multiple of 3 (a frame-shift artifact) are
skipped and reported, since that would corrupt the reading frame of every
gene concatenated after it.

Usage:
    python3 concatenate_genes.py \
        --gene-dir per_gene_alignments_nostop \
        --out-fasta concatenated_core.fasta \
        --out-manifest concatenated_core.manifest.tsv \
        --pattern "*.fasta"

Requires Biopython (same environment as the rest of the pipeline).
"""

import argparse
from pathlib import Path

from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-dir", required=True, type=Path,
                     help="Directory of per-gene alignments (e.g. per_gene_alignments_nostop)")
    ap.add_argument("--out-fasta", required=True, type=Path, help="Output concatenated FASTA")
    ap.add_argument("--out-manifest", required=True, type=Path,
                     help="Output TSV: gene, start, end, length, n_taxa_in_gene (1-based, inclusive coords in the concatenated alignment)")
    ap.add_argument("--pattern", default="*.fasta", help="Glob pattern for per-gene alignment files")
    ap.add_argument("--min-presence-frac", type=float, default=1.0,
                     help="Minimum fraction of included genes a taxon must appear in to be kept "
                          "(default 1.0 = strict intersection, no missing data)")
    args = ap.parse_args()

    gene_files = sorted(args.gene_dir.glob(args.pattern))
    if not gene_files:
        raise SystemExit(f"No files matching {args.pattern} in {args.gene_dir}")

    per_gene_seqs = {}   # gene_name -> {taxon: seq_str}
    per_gene_len = {}    # gene_name -> alignment length
    taxon_gene_count = {}  # taxon -> number of included genes it appears in
    skipped = []

    for gf in gene_files:
        gene = gf.stem
        records = list(SeqIO.parse(gf, "fasta"))
        if not records:
            skipped.append((gene, "empty file"))
            continue
        lengths = {len(r.seq) for r in records}
        if len(lengths) != 1:
            skipped.append((gene, f"unequal sequence lengths within gene: {lengths}"))
            continue
        length = lengths.pop()
        if length % 3 != 0:
            skipped.append((gene, f"length {length} not a multiple of 3 (frame-shift risk)"))
            continue

        seqs = {r.id: str(r.seq) for r in records}
        per_gene_seqs[gene] = seqs
        per_gene_len[gene] = length
        for taxon in seqs:
            taxon_gene_count[taxon] = taxon_gene_count.get(taxon, 0) + 1

    n_included_genes = len(per_gene_seqs)
    if n_included_genes == 0:
        raise SystemExit("No usable gene alignments after filtering -- check --pattern and skip reasons above")

    min_count = args.min_presence_frac * n_included_genes
    keep_taxa = sorted(t for t, c in taxon_gene_count.items() if c >= min_count)
    dropped_taxa = sorted(t for t in taxon_gene_count if t not in keep_taxa)

    gene_order = sorted(per_gene_seqs.keys())  # fixed, reproducible order

    concat = {taxon: [] for taxon in keep_taxa}
    manifest_rows = []
    pos = 1
    for gene in gene_order:
        seqs = per_gene_seqs[gene]
        length = per_gene_len[gene]
        n_taxa_in_gene = sum(1 for t in keep_taxa if t in seqs)
        for taxon in keep_taxa:
            seq = seqs.get(taxon)
            if seq is None:
                # shouldn't happen at min_presence_frac=1.0; for looser thresholds, gap-fill
                seq = "-" * length
            concat[taxon].append(seq)
        manifest_rows.append((gene, pos, pos + length - 1, length, n_taxa_in_gene))
        pos += length

    with open(args.out_fasta, "w") as f:
        for taxon in keep_taxa:
            f.write(f">{taxon}\n")
            f.write("".join(concat[taxon]) + "\n")

    with open(args.out_manifest, "w") as f:
        f.write("gene\tstart\tend\tlength\tn_taxa_in_gene\n")
        for row in manifest_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    total_len = pos - 1
    print(f"Genes found:            {len(gene_files)}")
    print(f"Genes included:         {n_included_genes}")
    print(f"Genes skipped:          {len(skipped)}")
    for gene, reason in skipped[:20]:
        print(f"  - {gene}: {reason}")
    if len(skipped) > 20:
        print(f"  ... and {len(skipped) - 20} more")
    print(f"Taxa kept:              {len(keep_taxa)}")
    print(f"Taxa dropped:           {len(dropped_taxa)} (present in <{args.min_presence_frac:.0%} of included genes)")
    print(f"Concatenated length:    {total_len} bp ({total_len // 3} codons)")
    print(f"Wrote alignment to:     {args.out_fasta}")
    print(f"Wrote manifest to:      {args.out_manifest}")


if __name__ == "__main__":
    main()