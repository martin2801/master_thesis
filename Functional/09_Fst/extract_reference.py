#!/usr/bin/env python3
"""
extract_reference.py

Extract a single sample's sequence from a whole-genome alignment FASTA,
preserving its coordinate system so downstream tools (e.g. Bakta) can be
run on it and produce gene coordinates that line up 1-to-1 with the
positions in a VCF/Fst file derived from the same alignment (e.g. via
snp-sites).

CRITICAL: only alignment gap characters ('-') are stripped. Any 'N' bases
are left in place, since they represent real (masked/unknown) positions in
the reference's own coordinate system - removing them would shift every
downstream coordinate and break correspondence with your VCF positions.

Usage:
    python extract_reference.py --alignment alignment.fasta --sample Reference \
        --out reference_extracted.fasta --header 1

Optional:
    --sample NAME     FASTA header (or header prefix) of the sequence to extract
                       (default: "Reference")
    --header NAME     Header to write in the output FASTA (default: same as --sample).
                       Set this to match your VCF's CHROM value (e.g. "1") so the
                       resulting Bakta GFF3 contig name matches directly - no
                       concatenation-mode guessing needed downstream.
"""

import argparse
import sys


def read_fasta(path):
    """Yield (header, sequence) tuples from a (possibly multi-line) FASTA file."""
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    ap = argparse.ArgumentParser(description="Extract one sample's sequence from an alignment FASTA.")
    ap.add_argument("--alignment", required=True, help="Whole-genome alignment FASTA file")
    ap.add_argument("--sample", default="Reference",
                     help="Header (or header prefix) of the sequence to extract (default: 'Reference')")
    ap.add_argument("--out", required=True, help="Output FASTA path")
    ap.add_argument("--header", default=None,
                     help="Header to write in the output FASTA. Set this to match your VCF's "
                          "CHROM value (e.g. '1') so downstream contig names line up directly. "
                          "Defaults to --sample.")
    ap.add_argument("--wrap", type=int, default=70,
                     help="Line-wrap width for the output FASTA sequence (default: 70)")
    args = ap.parse_args()

    match = None
    all_headers = []
    for header, seq in read_fasta(args.alignment):
        all_headers.append(header)
        if header == args.sample or header.split()[0] == args.sample:
            match = (header, seq)
            break  # take the first match

    if match is None:
        sys.exit(
            f"Could not find a sequence with header '{args.sample}' in {args.alignment}.\n"
            f"First few headers found in the file: {all_headers[:10]}"
        )

    header, seq = match
    n_gap = seq.count("-")
    n_n = seq.upper().count("N")
    cleaned = seq.replace("-", "")  # strip alignment gaps ONLY

    print(f"Found sequence '{header}', length {len(seq)} (alignment columns).")
    print(f"  Gap characters ('-') removed: {n_gap}")
    print(f"  'N' bases retained (masked/unknown, NOT removed): {n_n}")
    print(f"  Length after gap removal: {len(cleaned)}")

    if n_gap > 0:
        print(
            "\nWARNING: the Reference sequence itself contains alignment gaps.\n"
            "This usually means some other genome(s) in the alignment have an insertion "
            "the reference lacks. If your VCF's POS values were assigned as raw alignment "
            "COLUMN numbers (the typical snp-sites behavior), those positions include the "
            "gap columns - but this script strips gaps to reconstruct a clean reference "
            "FASTA, which shifts everything after each gap by however many bases were "
            "removed. In that case, positions in the extracted, gap-free reference will NOT "
            "line up 1-to-1 with your VCF POS values anymore, and you would need to keep "
            "track of the gap positions to re-derive the correct offset per site rather "
            "than relying on straightforward name-based coordinate matching downstream."
        )

    out_header = args.header if args.header is not None else args.sample
    with open(args.out, "w") as out:
        out.write(f">{out_header}\n")
        for i in range(0, len(cleaned), args.wrap):
            out.write(cleaned[i:i + args.wrap] + "\n")

    print(f"Saved extracted reference to {args.out} (header: '>{out_header}')")
    print("\nNext step: run Bakta on this file, e.g.:")
    print(f"  bakta --db <path-to-bakta-db> --output bakta_out {args.out}")
    print("Then use the resulting .gff3 with annotate_fst_genes.py - since the contig "
          f"name is now '{out_header}', it should match your Fst file's CHROM directly "
          "(no concatenated-coordinate mode needed).")


if __name__ == "__main__":
    main()