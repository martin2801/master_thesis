#!/usr/bin/env python3
"""
merge_annotations.py

Merges Panaroo's own gene annotations (from gene_presence_absence.csv)
into your Fisher/pyseer results table, keyed on gene cluster ID
(the "group_XXXX" / gene-symbol IDs shared across all Panaroo outputs).

This is step 1 before falling back to fresh BLAST/eggNOG-mapper re-annotation:
Panaroo's gene_presence_absence.csv already carries an "Annotation" column
(pulled from your original Prokka/Bakta run) for every gene cluster, which
often already has a real product name even when the "Gene" column itself
is a generic group_XXXX ID.

Usage:
    python3 merge_annotations.py \
        --results     fisher_results.tsv \
        --panaroo-csv gene_presence_absence.csv \
        --out         fisher_results_annotated.tsv
"""

import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True, help="fisher_results.tsv (or pyseer results.txt)")
ap.add_argument("--panaroo-csv", required=True, help="Panaroo's gene_presence_absence.csv")
ap.add_argument("--out", default="results_annotated.tsv")
ap.add_argument("--gene-col", default="gene",
                 help="name of the gene-ID column in --results (default: 'gene'; "
                      "use 'variant' for a raw pyseer results.txt)")
args = ap.parse_args()

# --- load Panaroo annotations, keyed by gene cluster ID ---
annot = {}
with open(args.panaroo_csv, newline="") as f:
    reader = csv.DictReader(f)
    # Panaroo's csv usually has "Gene" (the cluster ID) and "Annotation"
    # (and sometimes "Non-unique Gene name"); be tolerant of naming.
    cols = reader.fieldnames
    gene_col = "Gene" if "Gene" in cols else cols[0]
    annot_col = "Annotation" if "Annotation" in cols else None
    nonuniq_col = "Non-unique Gene name" if "Non-unique Gene name" in cols else None

    for row in reader:
        gid = row[gene_col]
        annot[gid] = {
            "annotation": row.get(annot_col, "") if annot_col else "",
            "non_unique_gene_name": row.get(nonuniq_col, "") if nonuniq_col else "",
        }

print(f"Loaded annotations for {len(annot)} gene clusters from {args.panaroo_csv}")

# --- merge into results table ---
with open(args.results, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames + ["panaroo_annotation", "panaroo_non_unique_gene_name"]
    rows = list(reader)

n_found = 0
for row in rows:
    gid = row[args.gene_col]
    a = annot.get(gid)
    if a:
        row["panaroo_annotation"] = a["annotation"]
        row["panaroo_non_unique_gene_name"] = a["non_unique_gene_name"]
        if a["annotation"] and "hypothetical" not in a["annotation"].lower():
            n_found += 1
    else:
        row["panaroo_annotation"] = "NOT_FOUND_IN_PANAROO_CSV"
        row["panaroo_non_unique_gene_name"] = ""

with open(args.out, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"{n_found} / {len(rows)} genes got a non-'hypothetical' annotation from Panaroo's own csv")
print(f"Written to: {args.out}")
print("Anything still 'hypothetical protein' / blank / NOT_FOUND needs re-annotation "
      "(see extract_seqs_for_reannotation.py)")