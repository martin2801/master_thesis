#!/usr/bin/env python3
"""
fisher_gene_test.py

Naive (lineage-UNCORRECTED) gene presence/absence association test,
for direct comparison against a lineage-corrected pyseer run.

For each gene in a Panaroo Rtab matrix, this builds the 2x2 table:

                pESI+   pESI-
    gene present  a       b
    gene absent   c       d

and runs Fisher's exact test, then applies Benjamini-Hochberg FDR
correction across all genes tested (equivalent in spirit to pyseer's
count of unique patterns, though not identical).

Usage:
    python3 fisher_gene_test.py \
        --rtab      gene_presence_absence.Rtab \
        --phenotype pESI_phenotype.tsv \
        --out       fisher_results.tsv \
        [--pyseer   pyseer_results.txt]   # optional: merge in for comparison

Input formats:
    --rtab       Panaroo/Roary style Rtab: first column = gene name,
                 remaining columns = one per genome, header row = genome IDs,
                 values 0/1 (absent/present).
    --phenotype  Two columns, tab-separated, header "samples\tpESI":
                 genome_id \t 0-or-1

Output:
    One row per gene: gene, a, b, c, d, present_in_pESI+, present_in_pESI-,
    odds_ratio, pvalue, qvalue (BH-FDR). Sorted by pvalue ascending.
    If --pyseer is given, its lrt-pvalue and filter-pvalue columns are
    merged in alongside (by gene name) so you can compare naive vs.
    lineage-corrected results directly in one table.
"""

import argparse
import csv
from scipy.stats import fisher_exact

def benjamini_hochberg(pvals):
    """Return BH-FDR adjusted q-values, same order as input pvals."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = order[n - rank]
        bh = pvals[i] * n / (n - rank + 1)
        prev = min(prev, bh)
        q[i] = prev
    return q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--out", default="fisher_results.tsv")
    ap.add_argument("--pyseer", default=None,
                     help="optional pyseer results.txt to merge in for comparison")
    args = ap.parse_args()

    # --- load phenotype ---
    pheno = {}
    with open(args.phenotype) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            sample, status = parts[0], parts[1]
            pheno[sample] = int(status)

    n_pos = sum(1 for v in pheno.values() if v == 1)
    n_neg = sum(1 for v in pheno.values() if v == 0)
    print(f"Loaded phenotype for {len(pheno)} genomes: {n_pos} pESI+, {n_neg} pESI-")

    # --- load Rtab ---
    with open(args.rtab) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        genome_cols = header[1:]

        # map genome column -> phenotype, warn about mismatches
        missing = [g for g in genome_cols if g not in pheno]
        if missing:
            print(f"WARNING: {len(missing)} genomes in Rtab have no phenotype entry "
                  f"and will be skipped, e.g. {missing[:5]}")

        col_status = [pheno.get(g) for g in genome_cols]

        results = []
        for row in reader:
            gene = row[0]
            vals = row[1:]
            a = b = c = d = 0  # a=present+pos b=present+neg c=absent+pos d=absent+neg
            for v, status in zip(vals, col_status):
                if status is None:
                    continue
                present = v not in ("0", "", "0.0")
                if status == 1:
                    if present: a += 1
                    else:       c += 1
                else:
                    if present: b += 1
                    else:       d += 1

            # skip invariant genes (no info): present/absent in everyone
            if (a + b) == 0 or (c + d) == 0:
                continue

            table = [[a, b], [c, d]]
            odds_ratio, pvalue = fisher_exact(table)
            pct_pos = 100 * a / (a + c) if (a + c) else float("nan")
            pct_neg = 100 * b / (b + d) if (b + d) else float("nan")
            results.append({
                "gene": gene, "a_present_pos": a, "b_present_neg": b,
                "c_absent_pos": c, "d_absent_neg": d,
                "pct_present_in_pESI_pos": pct_pos,
                "pct_present_in_pESI_neg": pct_neg,
                "odds_ratio": odds_ratio, "pvalue": pvalue,
            })

    print(f"Tested {len(results)} variable genes (invariant genes skipped)")

    pvals = [r["pvalue"] for r in results]
    qvals = benjamini_hochberg(pvals)
    for r, q in zip(results, qvals):
        r["qvalue"] = q

    results.sort(key=lambda r: r["pvalue"])

    # --- optional merge with pyseer results ---
    pyseer_map = {}
    if args.pyseer:
        with open(args.pyseer) as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                pyseer_map[row["variant"]] = row

    fieldnames = ["gene", "a_present_pos", "b_present_neg", "c_absent_pos", "d_absent_neg",
                  "pct_present_in_pESI_pos", "pct_present_in_pESI_neg",
                  "odds_ratio", "pvalue", "qvalue"]
    if args.pyseer:
        fieldnames += ["pyseer_lrt_pvalue", "pyseer_filter_pvalue"]

    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in results:
            row = dict(r)
            if args.pyseer:
                pr = pyseer_map.get(r["gene"])
                row["pyseer_lrt_pvalue"] = pr["lrt-pvalue"] if pr else "NA"
                row["pyseer_filter_pvalue"] = pr["filter-pvalue"] if pr else "NA"
            w.writerow(row)

    n_sig_raw = sum(1 for r in results if r["pvalue"] < 0.05)
    n_sig_fdr = sum(1 for r in results if r["qvalue"] < 0.05)
    print(f"\nNaive Fisher's exact (NO lineage correction):")
    print(f"  genes with raw p < 0.05:        {n_sig_raw} / {len(results)}")
    print(f"  genes with BH-FDR q < 0.05:     {n_sig_fdr} / {len(results)}")
    print(f"\nResults written to: {args.out}")

if __name__ == "__main__":
    main()