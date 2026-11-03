#!/usr/bin/env python3
"""
make_metadata.py

Build the metadata TSV for amrfinder_summary.py (and pyseer / Fisher tests)
from the genome folder plus the ID lists:

    sample      pESI    pESI_detail   group    pESI_binary
    394361      pESI+   partial       AGES     1
    10133415    pESI+   full          Sweep7   1
    11448455    pESI+   partial       other    1
    12345678    pESI-   absent        other    0

  sample       file name of each genome in --genome-dir minus .fna
               (= the sample name your AMRFinderPlus loop used)
  pESI         pESI+ if the genome is in the full OR partial pESI list,
               otherwise pESI-
  pESI_detail  full / partial / absent
  group        AGES / Sweep7 / other
  pESI_binary  1 / 0  (same information as `pESI`, for pyseer / Fisher)

Every genome in --genome-dir gets a row; anything not in a pESI list is
pESI-. The script stops if a genome is in both pESI lists or in both
AGES and Sweep7, and it warns about IDs in the lists that have no
matching .fna file (typo / renamed / excluded genome).

Usage (defaults = your server paths; override any of them):
    python3 make_metadata.py \
        --genome-dir /home/senekowitsch/Salmonella/AGES/Sweep/00_Data/infantis/all_genomes \
        --pesi-full     pESI_genomes_infantis_minimap2.txt \
        --pesi-partial  pESI_partial_genomes_infantis_minimap2.txt \
        --ages          AGES_genomes.txt \
        --sweep7        sweep7_genomes.txt \
        --out           metadata.tsv
"""

import argparse
import glob
import os
import re
import sys

EXTENSIONS = (".fna", ".fasta", ".fa")


def read_ids(path):
    """One ID per line; ignores blank lines, CRLF line endings, whitespace,
    and anything after the first tab/space (in case a column was added)."""
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.append(line.split()[0])
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        print(f"  NOTE: {path}: {len(dup)} duplicated IDs (kept once), "
              f"e.g. {sorted(dup)[:3]}", file=sys.stderr)
    return set(ids)


def genome_samples(genome_dir):
    files = []
    for ext in EXTENSIONS:
        files += glob.glob(os.path.join(genome_dir, f"*{ext}"))
    if not files:
        sys.exit(f"No {'/'.join(EXTENSIONS)} files found in {genome_dir}")
    samples = sorted({re.sub(r"\.(fna|fasta|fa)$", "", os.path.basename(f))
                      for f in files})
    return samples


def match_ids(ids, samples, label):
    """Map list IDs onto sample names. Exact match first; if the file names
    carry a prefix/suffix (e.g. 'GCA_x_10133415' or '10133415_contigs'),
    fall back to the ID appearing as a separate token in the name."""
    sample_set = set(samples)
    matched, missing = {}, []
    token_index = {}
    for s in samples:
        for tok in re.split(r"[^A-Za-z0-9]+", s):
            token_index.setdefault(tok, []).append(s)
    for i in ids:
        if i in sample_set:
            matched[i] = i
        elif len(token_index.get(i, [])) == 1:
            matched[i] = token_index[i][0]
        else:
            missing.append(i)
    if missing:
        print(f"  WARNING: {len(missing)}/{len(ids)} IDs in the {label} list "
              f"have no genome file, e.g. {sorted(missing)[:5]}",
              file=sys.stderr)
    fuzzy = sum(1 for k, v in matched.items() if k != v)
    if fuzzy:
        print(f"  NOTE: {fuzzy} {label} IDs matched a file name by token "
              f"(not exact) - check a few by eye", file=sys.stderr)
    return set(matched.values())

def main():
    base = "/home/senekowitsch/Salmonella/AGES/Sweep/00_Data/infantis"
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--genome-dir", default=f"{base}/all_genomes")
    ap.add_argument("--pesi-full", default="pESI_genomes_infantis_minimap2.txt")
    ap.add_argument("--pesi-partial",
                    default="pESI_partial_genomes_infantis_minimap2.txt")
    ap.add_argument("--ages", default=None,
                    help="optional AGES ID list; omit to skip the AGES group")
    ap.add_argument("--sweep7", default="sweep7_genomes.txt")
    ap.add_argument("--out", default="metadata.tsv")
    args = ap.parse_args()

    samples = genome_samples(args.genome_dir)
    print(f"{len(samples)} genomes in {args.genome_dir}", file=sys.stderr)

    full = match_ids(read_ids(args.pesi_full), samples, "pESI full")
    partial = match_ids(read_ids(args.pesi_partial), samples, "pESI partial")
    if args.ages:
        ages = match_ids(read_ids(args.ages), samples, "AGES")
    else:
        ages = set()
        print("  NOTE: no --ages list given; group = Sweep7 / other only",
              file=sys.stderr)
    sweep7 = match_ids(read_ids(args.sweep7), samples, "Sweep7")

    # consistency checks
    checks = [(full, partial, "pESI full", "pESI partial")]
    if ages:
        checks.append((ages, sweep7, "AGES", "Sweep7"))
    for a, b, na, nb in checks:
        both = a & b
        if both:
            sys.exit(f"ERROR: {len(both)} genomes are in both the {na} and "
                     f"{nb} lists, e.g. {sorted(both)[:5]} - fix the lists")

    rows = []
    for s in samples:
        if s in full:
            detail = "full"
        elif s in partial:
            detail = "partial"
        else:
            detail = "absent"
        pesi = "pESI-" if detail == "absent" else "pESI+"
        group = "AGES" if s in ages else "Sweep7" if s in sweep7 else "other"
        rows.append((s, pesi, detail, group, "0" if pesi == "pESI-" else "1"))

    with open(args.out, "w") as out:
        out.write("sample\tpESI\tpESI_detail\tgroup\tpESI_binary\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    # summary
    from collections import Counter
    c = Counter((r[3], r[2]) for r in rows)
    groups = (["AGES"] if args.ages else []) + ["Sweep7", "other"]
    print(f"\nwrote {args.out} ({len(rows)} genomes)\n", file=sys.stderr)
    print(f"{'group':<8}{'full':>7}{'partial':>9}{'absent':>8}{'total':>7}",
          file=sys.stderr)
    for g in groups:
        n = [c[(g, d)] for d in ("full", "partial", "absent")]
        print(f"{g:<8}{n[0]:>7}{n[1]:>9}{n[2]:>8}{sum(n):>7}",
              file=sys.stderr)
    tot = [sum(c[(g, d)] for g in groups) for d in ("full", "partial",
                                                     "absent")]
    print(f"{'total':<8}{tot[0]:>7}{tot[1]:>9}{tot[2]:>8}{sum(tot):>7}",
          file=sys.stderr)
    print(f"\npESI+ = {tot[0] + tot[1]}, pESI- = {tot[2]}", file=sys.stderr)

if __name__ == "__main__":
    main()
