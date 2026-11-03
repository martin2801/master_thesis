#!/usr/bin/env python3
"""
Consensus contig labelling: pESI / other plasmid / chromosome.

Evidence per contig:
  - minimap2 PAF vs closed pESI reference (query coverage, identity)
  - MOB-suite contig_report.txt (molecule_type == plasmid)
  - Platon  <name>.plasmid.fasta
  - geNomad <name>_summary/<name>_plasmid_summary.tsv

Labels:
  pESI              contig >= min_len and >= min_qcov of its length aligns to pESI
  pESI_short        same, but contig < min_len (often IS26/IS elements -> ambiguous,
                    can be chromosomal; kept OUT of chromosome by default, see flag)
  plasmid_other     >= 2 of 3 tools (MOB, Platon, geNomad) call plasmid
  uncertain         exactly 1 tool calls plasmid
  chromosome        no evidence

Genome status (from pESI reference coverage):
  pESI_pos >= --pos-cov, pESI_partial >= --partial-cov, else pESI_neg

Outputs (in --out):
  contig_labels.tsv              one row per contig
  genome_summary.tsv             one row per genome
  chromosomes/<g>.fasta          chromosome contigs (always written)
  pESI/<g>.fasta                 pESI contigs (only for genomes that have any)
  plasmid_other/<g>.fasta        other plasmid contigs (only for genomes that have any)

Standard library only.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path, records):
    with open(path, "w") as out:
        for cid, seq in records:
            out.write(f">{cid}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


def fasta_ids(path):
    if not os.path.isfile(path):
        return None
    return {n for n, _ in read_fasta(path)}


def merged_length(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def parse_paf(path, min_ident, min_aln):
    """Return {contig: [(qs,qe),...]}, {contig: [(ts,te),...]}, ref_len."""
    q_iv, t_iv, ref_len = defaultdict(list), defaultdict(list), None
    if not os.path.isfile(path):
        return None, None, None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qs, qe = f[0], int(f[2]), int(f[3])
            tlen, ts, te = int(f[6]), int(f[7]), int(f[8])
            matches, alnlen = int(f[9]), int(f[10])
            ref_len = tlen  # single-sequence reference assumed
            if alnlen < min_aln or matches / alnlen < min_ident:
                continue
            q_iv[qname].append((qs, qe))
            t_iv[qname].append((ts, te))
    return q_iv, t_iv, ref_len


def mob_plasmid_ids(mob_dir, name):
    cands = glob.glob(os.path.join(mob_dir, f"{name}_mob_recon_output", "contig_report.txt"))
    if not cands:
        return None
    ids = set()
    with open(cands[0]) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("molecule_type", "").lower() == "plasmid":
                ids.add(row["contig_id"].split()[0])
    return ids


def genomad_plasmid_ids(genomad_dir, name):
    cands = glob.glob(os.path.join(genomad_dir, name, "*_summary", "*_plasmid_summary.tsv"))
    if not cands:
        return None
    ids = set()
    with open(cands[0]) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ids.add(row["seq_name"].split()[0])
    return ids


def platon_plasmid_ids(platon_dir, name):
    return fasta_ids(os.path.join(platon_dir, name, f"{name}.plasmid.fasta"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--ext", default=".fna")
    ap.add_argument("--paf-dir", required=True)
    ap.add_argument("--platon-dir", required=True)
    ap.add_argument("--genomad-dir", required=True)
    ap.add_argument("--mob-dir", required=True)
    ap.add_argument("--abricate-list", help="old pESI+ list, for comparison")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-ident", type=float, default=0.90)
    ap.add_argument("--min-aln", type=int, default=500)
    ap.add_argument("--min-qcov", type=float, default=0.5, help="fraction of contig aligned to pESI")
    ap.add_argument("--min-len", type=int, default=5000, help="shorter pESI-matching contigs -> pESI_short")
    ap.add_argument("--pos-cov", type=float, default=0.80, help="pESI ref coverage for pESI_pos")
    ap.add_argument("--partial-cov", type=float, default=0.20, help="pESI ref coverage for pESI_partial")
    ap.add_argument("--keep-pesi-short", action="store_true", help="keep pESI_short contigs in chromosome output")
    ap.add_argument("--drop-uncertain", action="store_true", help="remove 'uncertain' contigs from chromosome output")
    ap.add_argument("--pesi-include-short", action="store_true",
                    help="also write pESI_short contigs to pESI/<g>.fasta (default: only 'pESI')")
    ap.add_argument("--pesi-status", default="pESI_pos,pESI_partial",
                    help="comma-separated genome statuses for which pESI/<g>.fasta is written "
                         "(default: pESI_pos,pESI_partial)")
    a = ap.parse_args()

    for sub in ("chromosomes", "pESI", "plasmid_other"):
        os.makedirs(os.path.join(a.out, sub), exist_ok=True)

    old_pos = None
    if a.abricate_list:
        with open(a.abricate_list) as fh:
            old_pos = {l.strip() for l in fh if l.strip()}

    keep_labels = {"chromosome", "uncertain"}
    if a.drop_uncertain:
        keep_labels.discard("uncertain")
    if a.keep_pesi_short:
        keep_labels.add("pESI_short")

    pesi_labels = {"pESI"}
    if a.pesi_include_short:
        pesi_labels.add("pESI_short")
    pesi_status = {s.strip() for s in a.pesi_status.split(",") if s.strip()}

    genomes = sorted(glob.glob(os.path.join(a.genomes, f"*{a.ext}")))
    if not genomes:
        sys.exit(f"No *{a.ext} files in {a.genomes}")

    lab_fh = open(os.path.join(a.out, "contig_labels.tsv"), "w")
    lab = csv.writer(lab_fh, delimiter="\t")
    lab.writerow(["genome", "contig", "length", "pesi_qcov", "mob", "platon", "genomad", "label"])

    sum_fh = open(os.path.join(a.out, "genome_summary.tsv"), "w")
    summ = csv.writer(sum_fh, delimiter="\t")
    summ.writerow(["genome", "n_contigs", "total_len", "pesi_ref_coverage", "status", "abricate_old",
                   "agree", "n_pESI", "n_pESI_short", "n_plasmid_other", "n_uncertain",
                   "chromosome_len_kept", "pesi_len_written", "plasmid_other_len",
                   "missing_tools"])

    n_pesi_files, n_other_files = 0, 0
    for fna in genomes:
        name = os.path.basename(fna)[: -len(a.ext)]
        q_iv, t_iv, ref_len = parse_paf(os.path.join(a.paf_dir, f"{name}.paf"), a.min_ident, a.min_aln)
        mob = mob_plasmid_ids(a.mob_dir, name)
        pla = platon_plasmid_ids(a.platon_dir, name)
        gen = genomad_plasmid_ids(a.genomad_dir, name)
        missing = [t for t, v in (("minimap2", q_iv), ("mob", mob), ("platon", pla), ("genomad", gen)) if v is None]
        q_iv, t_iv = q_iv or {}, t_iv or {}
        mob, pla, gen = mob or set(), pla or set(), gen or set()

        counts = defaultdict(int)
        ref_ivs, total_len = [], 0
        kept, pesi, other = [], [], []
        for cid, seq in read_fasta(fna):
            L = len(seq)
            total_len += L
            qcov = merged_length(q_iv.get(cid, [])) / L if L else 0.0
            votes = (cid in mob) + (cid in pla) + (cid in gen)
            if qcov >= a.min_qcov:
                label = "pESI" if L >= a.min_len else "pESI_short"
                ref_ivs.extend(t_iv.get(cid, []))
            elif votes >= 2:
                label = "plasmid_other"
            elif votes == 1:
                label = "uncertain"
            else:
                label = "chromosome"
            counts[label] += 1
            lab.writerow([name, cid, L, f"{qcov:.3f}", int(cid in mob), int(cid in pla), int(cid in gen), label])
            if label in keep_labels:
                kept.append((cid, seq))
            if label in pesi_labels:
                pesi.append((cid, seq))
            if label == "plasmid_other":
                other.append((cid, seq))

        cov = merged_length(ref_ivs) / ref_len if ref_len else 0.0
        status = "pESI_pos" if cov >= a.pos_cov else "pESI_partial" if cov >= a.partial_cov else "pESI_neg"
        old = "NA" if old_pos is None else ("pos" if name in old_pos else "neg")
        agree = "NA" if old_pos is None else str((old == "pos") == (status == "pESI_pos"))

        # chromosome: always written
        write_fasta(os.path.join(a.out, "chromosomes", f"{name}.fasta"), kept)

        # pESI: only if the genome has the right status and any pESI contigs
        pesi_out = os.path.join(a.out, "pESI", f"{name}.fasta")
        if pesi and status in pesi_status:
            write_fasta(pesi_out, pesi)
            n_pesi_files += 1
            pesi_len = sum(len(s) for _, s in pesi)
        else:
            if os.path.exists(pesi_out):  # stale file from an earlier run
                os.remove(pesi_out)
            pesi_len = 0

        # other plasmids: only if any
        other_out = os.path.join(a.out, "plasmid_other", f"{name}.fasta")
        if other:
            write_fasta(other_out, other)
            n_other_files += 1
        elif os.path.exists(other_out):
            os.remove(other_out)

        summ.writerow([name, sum(counts.values()), total_len, f"{cov:.3f}", status, old, agree,
                       counts["pESI"], counts["pESI_short"], counts["plasmid_other"], counts["uncertain"],
                       sum(len(s) for _, s in kept), pesi_len, sum(len(s) for _, s in other),
                       ",".join(missing) or "-"])

    lab_fh.close()
    sum_fh.close()
    print(f"Done: {len(genomes)} genomes -> {a.out}")
    print(f"  pESI FASTA written:          {n_pesi_files}")
    print(f"  plasmid_other FASTA written: {n_other_files}")


if __name__ == "__main__":
    main()