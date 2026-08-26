#!/usr/bin/env python3
"""
annotate_fst_genes.py

Overlap vcftools Fst output (per-site or windowed) with gene coordinates from a
Bakta or Prokka GFF3 annotation file, to find which genes contain high-Fst SNPs.

Usage:
    python annotate_fst_genes.py \
        --fst fst_results.weir.fst \
        --gff annotation.gff3 \
        --out-snps fst_annotated_snps.tsv \
        --out-genes fst_gene_summary.tsv

Optional:
    --feature-types CDS         Comma-separated GFF3 feature types to use (default: CDS)
    --fst-threshold 0.95        Threshold used to count "high-Fst SNPs" per gene in the summary
    --floor-negative             Clamp genuine negative Fst values to 0 before summarizing
    --min-variants N             (windowed Fst only) drop windows with fewer than N variants

Notes:
    - Works with either per-site (POS, WEIR_AND_COCKERHAM_FST) or windowed
      (BIN_START, BIN_END, WEIGHTED_FST, N_VARIANTS) vcftools output; auto-detected.
    - vcftools' "-nan" (undefined Fst) is parsed correctly and excluded automatically;
      never hand-edit the raw Fst file (e.g. with sed) before running this.
    - CHROM/seqid names in the Fst file and the GFF3 must match exactly. The script
      will warn you if they don't overlap at all.
    - SNPs/windows that don't fall inside any annotated feature are labeled "intergenic".
"""

import argparse
import re
import sys
import bisect
import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Fst loading (same robust logic as plot_fst.py: handles per-site/windowed,
# "-nan" undefined values, and optional flooring of genuine negative Fst)
# ---------------------------------------------------------------------------

def load_fst(path, min_variants=1, floor_negative=False):
    df = pd.read_csv(path, sep="\t")
    cols = set(df.columns)

    if {"BIN_START", "BIN_END", "WEIGHTED_FST", "MEAN_FST", "N_VARIANTS"}.issubset(cols):
        mode = "windowed"
        df["POS_START"] = df["BIN_START"]
        df["POS_END"] = df["BIN_END"]
        fst_col = "WEIGHTED_FST"
    elif {"POS", "WEIR_AND_COCKERHAM_FST"}.issubset(cols):
        mode = "persite"
        df["POS_START"] = df["POS"]
        df["POS_END"] = df["POS"]
        fst_col = "WEIR_AND_COCKERHAM_FST"
    else:
        sys.exit(
            f"Could not recognize Fst file format. Columns found: {list(df.columns)}\n"
            "Expected either windowed vcftools output (BIN_START, BIN_END, WEIGHTED_FST, ...) "
            "or per-site output (POS, WEIR_AND_COCKERHAM_FST)."
        )

    df["FST"] = pd.to_numeric(df[fst_col], errors="coerce")

    n_total = len(df)
    n_undefined = df["FST"].isna().sum()
    if n_undefined:
        print(f"Note: {n_undefined} of {n_total} sites/windows had undefined Fst "
              f"('-nan'/'nan') and will be excluded.")
    df = df.dropna(subset=["FST"])

    if mode == "windowed" and "N_VARIANTS" in df.columns:
        before = len(df)
        df = df[df["N_VARIANTS"] >= min_variants]
        dropped = before - len(df)
        if dropped:
            print(f"Dropped {dropped} windows with fewer than {min_variants} variants.")

    if floor_negative:
        n_neg = (df["FST"] < 0).sum()
        if n_neg:
            print(f"Flooring {n_neg} genuine negative Fst values to 0.")
            df.loc[df["FST"] < 0, "FST"] = 0.0

    df = df.sort_values(["CHROM", "POS_START"]).reset_index(drop=True)
    return df, mode


# ---------------------------------------------------------------------------
# GFF3 loading
# ---------------------------------------------------------------------------

ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")

def parse_attributes(attr_str):
    return dict(ATTR_RE.findall(attr_str))


def load_gff(path, feature_types):
    rows = []
    contig_lengths = {}   # seqid -> length, from 'region' features, in file order
    contig_order = []     # order contigs first appear in the file
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, frame, attrs = parts[:9]

            if seqid not in contig_lengths:
                contig_order.append(seqid)
                contig_lengths[seqid] = 0
            if ftype == "region":
                # 'region' features give the authoritative full contig length
                contig_lengths[seqid] = max(contig_lengths[seqid], int(end))

            if ftype not in feature_types:
                continue
            attrd = parse_attributes(attrs)
            gene_name = attrd.get("gene") or attrd.get("Name") or attrd.get("ID") or "unknown"
            rows.append({
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "feature_type": ftype,
                "gene": gene_name,
                "locus_tag": attrd.get("locus_tag", attrd.get("ID", "")),
                "product": attrd.get("product", ""),
            })
    if not rows:
        sys.exit(f"No features of type {feature_types} found in {path}. "
                  "Check --feature-types matches column 3 of your GFF3.")
    gff = pd.DataFrame(rows).sort_values(["seqid", "start"]).reset_index(drop=True)

    # Fall back to max feature end per contig if no explicit 'region' feature was seen
    for seqid in contig_order:
        if contig_lengths.get(seqid, 0) == 0:
            sub = gff[gff["seqid"] == seqid]
            if len(sub):
                contig_lengths[seqid] = int(sub["end"].max())

    contig_info = {"order": contig_order, "lengths": contig_lengths}
    return gff, contig_info


def build_concat_offsets(contig_info):
    """Cumulative start offset (0-based) for each contig, in file-appearance order."""
    offsets = {}
    running = 0
    for seqid in contig_info["order"]:
        offsets[seqid] = running
        running += contig_info["lengths"].get(seqid, 0)
    return offsets, running


def concatenate_gff_coords(gff_df, contig_info):
    """Return a copy of gff_df with 'start'/'end' shifted into concatenated
    genome coordinates (assuming contigs are laid end-to-end, in the order
    they first appear in the GFF3, with no gaps/padding between them)."""
    offsets, total_len = build_concat_offsets(contig_info)
    out = gff_df.copy()
    out["concat_start"] = out.apply(lambda r: r["start"] + offsets.get(r["seqid"], 0), axis=1)
    out["concat_end"] = out.apply(lambda r: r["end"] + offsets.get(r["seqid"], 0), axis=1)
    return out, total_len


# ---------------------------------------------------------------------------
# Overlap join: for each Fst row (a point, or a window), find any gene(s)
# whose [start, end] interval overlaps it.
# ---------------------------------------------------------------------------

def annotate_overlaps_concat(fst_df, gff_concat_df, total_len):
    """Overlap join in concatenated-genome coordinate space: every Fst row is
    matched against gene 'concat_start'/'concat_end' intervals directly,
    ignoring CHROM/seqid (since the Fst file uses one pseudo-chromosome that
    represents all contigs laid end-to-end)."""
    gff_sorted = gff_concat_df.sort_values("concat_start").reset_index(drop=True)
    starts = gff_sorted["concat_start"].to_numpy()

    max_fst_pos = fst_df["POS_END"].max()
    if max_fst_pos > total_len * 1.05:
        print(f"WARNING: max Fst position ({max_fst_pos}) exceeds the concatenated "
              f"GFF3 genome length ({total_len}) by more than 5%. The concatenation "
              "assumption (contig order/no gaps) may not hold - results may be unreliable.")

    annotated_rows = []
    for _, row in fst_df.iterrows():
        pos_start, pos_end = row["POS_START"], row["POS_END"]
        match_idx = []
        if len(gff_sorted):
            i = bisect.bisect_right(starts, pos_end) - 1
            for j in range(max(0, i - 2), min(len(gff_sorted), i + 3)):
                if gff_sorted["concat_start"].iat[j] <= pos_end and gff_sorted["concat_end"].iat[j] >= pos_start:
                    match_idx.append(j)

        if match_idx:
            for j in match_idx:
                g = gff_sorted.iloc[j]
                annotated_rows.append({
                    **row.to_dict(),
                    "gene": g["gene"],
                    "locus_tag": g["locus_tag"],
                    "product": g["product"],
                    "contig": g["seqid"],
                    "gene_start": g["start"],
                    "gene_end": g["end"],
                    "strand": g["strand"],
                })
        else:
            annotated_rows.append({
                **row.to_dict(),
                "gene": "intergenic",
                "locus_tag": "",
                "product": "",
                "contig": "",
                "gene_start": np.nan,
                "gene_end": np.nan,
                "strand": "",
            })

    return pd.DataFrame(annotated_rows)


def annotate_overlaps(fst_df, gff_df):
    fst_chroms = set(fst_df["CHROM"].unique())
    gff_chroms = set(gff_df["seqid"].unique())
    if not (fst_chroms & gff_chroms):
        print("WARNING: no contig/chromosome names in common between the Fst file and the "
              f"GFF3.\n  Fst CHROM examples: {list(fst_chroms)[:5]}\n"
              f"  GFF3 seqid examples: {list(gff_chroms)[:5]}\n"
              "  Check that both were generated against the same reference and that names "
              "match exactly (e.g. 'contig_1' vs '1').")

    annotated_rows = []
    for chrom, fst_sub in fst_df.groupby("CHROM"):
        gff_sub = gff_df[gff_df["seqid"] == chrom].sort_values("start").reset_index(drop=True)
        starts = gff_sub["start"].to_numpy()
        ends = gff_sub["end"].to_numpy()

        for _, row in fst_sub.iterrows():
            pos_start, pos_end = row["POS_START"], row["POS_END"]
            match_idx = []
            if len(gff_sub):
                # candidate gene(s) whose start <= pos_end (i.e. could overlap);
                # check a small neighborhood around the insertion point for safety
                # in case of overlapping/nested annotations
                i = bisect.bisect_right(starts, pos_end) - 1
                for j in range(max(0, i - 2), min(len(gff_sub), i + 3)):
                    if starts[j] <= pos_end and ends[j] >= pos_start:
                        match_idx.append(j)

            if match_idx:
                for j in match_idx:
                    g = gff_sub.iloc[j]
                    annotated_rows.append({
                        **row.to_dict(),
                        "gene": g["gene"],
                        "locus_tag": g["locus_tag"],
                        "product": g["product"],
                        "gene_start": g["start"],
                        "gene_end": g["end"],
                        "strand": g["strand"],
                    })
            else:
                annotated_rows.append({
                    **row.to_dict(),
                    "gene": "intergenic",
                    "locus_tag": "",
                    "product": "",
                    "gene_start": np.nan,
                    "gene_end": np.nan,
                    "strand": "",
                })

    return pd.DataFrame(annotated_rows)


def summarize_by_gene(annotated_df, fst_threshold):
    genic = annotated_df[annotated_df["gene"] != "intergenic"].copy()
    if genic.empty:
        print("No SNPs/windows overlapped any annotated gene.")
        return pd.DataFrame()

    summary = (
        genic.groupby(["gene", "locus_tag", "product", "gene_start", "gene_end", "strand"])
        .agg(
            n_variants=("FST", "size"),
            mean_fst=("FST", "mean"),
            max_fst=("FST", "max"),
            n_high_fst=("FST", lambda x: (x >= fst_threshold).sum()),
        )
        .reset_index()
        .sort_values(["n_high_fst", "max_fst"], ascending=False)
    )
    return summary


def main():
    ap = argparse.ArgumentParser(description="Overlap Fst results with gene annotations (GFF3).")
    ap.add_argument("--fst", required=True, help="vcftools Fst output file (per-site or windowed)")
    ap.add_argument("--gff", required=True, help="Bakta or Prokka GFF3 annotation file")
    ap.add_argument("--out-snps", default="fst_annotated_snps.tsv",
                     help="Output: every Fst site/window with its overlapping gene (if any)")
    ap.add_argument("--out-genes", default="fst_gene_summary.tsv",
                     help="Output: one row per gene, summarizing overlapping Fst values")
    ap.add_argument("--feature-types", default="CDS",
                     help="Comma-separated GFF3 column-3 feature types to include (default: CDS)")
    ap.add_argument("--fst-threshold", type=float, default=0.95,
                     help="Threshold for counting 'high-Fst' variants per gene (default: 0.95)")
    ap.add_argument("--floor-negative", action="store_true",
                     help="Clamp genuine negative Fst values to 0 before summarizing")
    ap.add_argument("--min-variants", type=int, default=1,
                     help="(windowed Fst only) drop windows with fewer than N variants")
    args = ap.parse_args()

    feature_types = set(t.strip() for t in args.feature_types.split(","))

    print(f"Loading Fst results from {args.fst} ...")
    fst_df, mode = load_fst(args.fst, min_variants=args.min_variants,
                             floor_negative=args.floor_negative)
    print(f"  Mode: {mode}, {len(fst_df)} rows after filtering.")

    print(f"Loading gene annotations from {args.gff} (feature types: {feature_types}) ...")
    gff_df, contig_info = load_gff(args.gff, feature_types)
    print(f"  {len(gff_df)} features loaded across {len(contig_info['order'])} contigs.")

    fst_chroms = set(fst_df["CHROM"].unique())
    gff_chroms = set(contig_info["order"])
    use_concat_mode = fst_chroms.isdisjoint(gff_chroms)

    if use_concat_mode:
        print("Fst CHROM names don't match any GFF3 contig names "
              f"(Fst: {sorted(fst_chroms)}; GFF3 has {len(gff_chroms)} contigs) - "
              "switching to CONCATENATED-COORDINATE mode.")
        print("Assuming contigs were laid end-to-end in the order they appear in the "
              "GFF3, with no gaps, to build the alignment/VCF coordinate system.")
        gff_concat, total_len = concatenate_gff_coords(gff_df, contig_info)
        print(f"  Concatenated GFF3 genome length: {total_len} bp")
        print("Overlapping Fst positions with gene intervals (concatenated coordinates) ...")
        annotated = annotate_overlaps_concat(fst_df, gff_concat, total_len)
    else:
        print("Overlapping Fst positions with gene intervals (matching by contig name) ...")
        annotated = annotate_overlaps(fst_df, gff_df)

    keep_cols = ["CHROM", "POS_START", "POS_END", "FST", "gene", "locus_tag",
                 "product", "contig", "gene_start", "gene_end", "strand"]
    keep_cols = [c for c in keep_cols if c in annotated.columns]
    annotated[keep_cols].to_csv(args.out_snps, sep="\t", index=False)
    print(f"Saved per-site/window annotation to {args.out_snps}")

    n_genic = (annotated["gene"] != "intergenic").sum()
    n_intergenic = (annotated["gene"] == "intergenic").sum()
    print(f"  {n_genic} sites/windows fell inside an annotated {args.feature_types} feature; "
          f"{n_intergenic} were intergenic.")

    summary = summarize_by_gene(annotated, args.fst_threshold)
    if not summary.empty:
        summary.to_csv(args.out_genes, sep="\t", index=False)
        print(f"Saved per-gene summary to {args.out_genes}")
        print(f"\nTop genes by number of high-Fst (>= {args.fst_threshold}) variants:")
        print(summary.head(10).to_string(index=False))


if __name__ == "__main__":
    main()