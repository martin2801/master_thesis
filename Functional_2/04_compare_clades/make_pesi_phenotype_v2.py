#!/usr/bin/env python3
"""
Build the pESI phenotype for pyseer / Fisher from 07_plasmid_id genome_summary.tsv.

- Only genomes with a Bakta chromosome GFF3 are included.
- pESI_pos -> 1, pESI_neg -> 0
- pESI_partial -> excluded (default), or forced to 1 / 0 with --partial pos|neg
Also writes a plain list of pESI+ genome IDs (for tree plotting) and a short report.
"""
import argparse
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", required=True, help="07_plasmid_id .../consensus/genome_summary.tsv")
    ap.add_argument("--bakta-dir", required=True, help="bakta_chromosome dir with <g>/<g>.gff3")
    ap.add_argument("--out", required=True, help="phenotype TSV (samples, pESI)")
    ap.add_argument("--pos-list", required=True, help="output list of pESI+ IDs")
    ap.add_argument("--partial", choices=["exclude", "pos", "neg"], default="exclude")
    a = ap.parse_args()

    status = {}
    with open(a.summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            status[row["genome"]] = row

    bakta = sorted(d for d in os.listdir(a.bakta_dir)
                   if os.path.isfile(os.path.join(a.bakta_dir, d, f"{d}.gff3")))
    if not bakta:
        sys.exit(f"No Bakta GFF3s found in {a.bakta_dir}")

    pheno, excluded, no_summary = {}, [], []
    for g in bakta:
        row = status.get(g)
        if row is None:
            no_summary.append(g)
            continue
        s = row["status"]
        if s == "pESI_pos":
            pheno[g] = 1
        elif s == "pESI_neg":
            pheno[g] = 0
        elif a.partial == "pos":
            pheno[g] = 1
        elif a.partial == "neg":
            pheno[g] = 0
        else:
            excluded.append(g)

    with open(a.out, "w") as out:
        out.write("samples\tpESI\n")
        for g in sorted(pheno):
            out.write(f"{g}\t{pheno[g]}\n")
    with open(a.pos_list, "w") as out:
        for g in sorted(pheno):
            if pheno[g] == 1:
                out.write(g + "\n")

    n_pos = sum(pheno.values())
    print(f"Bakta genomes: {len(bakta)}")
    print(f"Phenotype written: {n_pos} pESI+ / {len(pheno) - n_pos} pESI- -> {a.out}")
    print(f"pESI_partial ({a.partial}): {len(excluded) if a.partial == 'exclude' else 'recoded'}")
    if excluded:
        print("  excluded: " + ", ".join(excluded))
    if no_summary:
        print(f"WARNING: {len(no_summary)} Bakta genomes not in summary: " + ", ".join(no_summary))

    # changes vs. old abricate call, if the summary has that column
    changed = [(g, r["abricate_old"], r["status"]) for g, r in status.items()
               if r.get("agree") == "False"]
    if changed:
        print(f"Status changed vs. abricate for {len(changed)} genomes:")
        for g, old, new in sorted(changed):
            print(f"  {g}\tabricate={old}\tnew={new}")


if __name__ == "__main__":
    main()
