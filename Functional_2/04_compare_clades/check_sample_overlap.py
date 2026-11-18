#!/usr/bin/env python3
"""Check sample overlap between a pyseer phenotype file, a Panaroo Rtab
and a pyseer similarity matrix."""
import argparse


def load_ids_from_phenotype(path):
    ids = set()
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if parts and parts[0]:
                ids.add(parts[0])
    return ids


def load_ids_from_rtab(path):
    with open(path) as f:
        header = f.readline().strip('\n').split('\t')
    return set(header[1:])  # first column is gene name, rest are sample IDs


def load_ids_from_similarity(path):
    with open(path) as f:
        header = f.readline().strip('\n').split('\t')
    return set(header[1:]) if not header[0].strip() else set(header)


def report(label, ids, show):
    print(f"{label}: {len(ids)}")
    if show and ids:
        for i in sorted(ids)[:show]:
            print(f"    {i}")
        if len(ids) > show:
            print(f"    ... ({len(ids) - show} more)")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phenotype", required=True,
                   help="pyseer phenotype TSV (with header, sample IDs in column 1)")
    p.add_argument("--rtab", required=True,
                   help="Panaroo gene_presence_absence.Rtab")
    p.add_argument("--similarity", required=True,
                   help="pyseer similarity matrix (phylogeny_distance.py --lmm output)")
    p.add_argument("--show", type=int, default=0, metavar="N",
                   help="print up to N example IDs for each mismatch (default: 0)")
    args = p.parse_args()

    pheno_ids = load_ids_from_phenotype(args.phenotype)
    rtab_ids = load_ids_from_rtab(args.rtab)
    sim_ids = load_ids_from_similarity(args.similarity)

    print(f"phenotype: {len(pheno_ids)} samples")
    print(f"gene presence/absence: {len(rtab_ids)} samples")
    print(f"similarity matrix: {len(sim_ids)} samples")

    common = pheno_ids & rtab_ids & sim_ids
    print(f"\nCommon to all three: {len(common)}\n")

    report("In phenotype but not in Rtab", pheno_ids - rtab_ids, args.show)
    report("In phenotype but not in similarity", pheno_ids - sim_ids, args.show)
    report("In Rtab but not in phenotype", rtab_ids - pheno_ids, args.show)
    report("In Rtab but not in similarity", rtab_ids - sim_ids, args.show)
    report("In similarity but not in Rtab", sim_ids - rtab_ids, args.show)


if __name__ == "__main__":
    main()