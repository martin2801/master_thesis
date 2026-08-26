#!/usr/bin/env python3
"""
plot_fst.py

Visualize Fst values across the genome from vcftools --weir-fst-pop output.
Auto-detects whether the input is windowed (.fst.windowed) or per-site (.weir.fst).

Usage:
    python plot_fst.py --input fst_results.windowed.weir.fst --out fst_plot.png
    python plot_fst.py --input fst_results.weir.fst --out fst_plot.png

Optional:
    --min-variants N     (windowed only) drop windows with fewer than N variants
    --title "My title"
    --fst-hit-threshold  draw a vertical marker line for every site/window at or
                         above this Fst value (default: 0.999)
    --floor-negative     clamp genuine negative Fst values to 0 (common convention;
                         negative Fst is sampling noise around true ~0 differentiation)
    --drop-undefined     drop sites where Fst was undefined (vcftools reports "-nan"/"nan")
                         instead of keeping them as missing/excluded from the plot (default
                         behavior already excludes them from plotting; this flag just makes
                         the exclusion and count explicit in the console output)

Note on vcftools' "-nan" values:
    vcftools writes "-nan" (or "nan") when the Weir & Cockerham Fst estimator is
    mathematically undefined at a site (e.g. no informative allele-frequency variance
    to divide by). These are NOT the same as Fst = 0, and should never be edited in the
    raw file with a text tool like sed (e.g. "sed 's/-/0/g'") — doing so corrupts the
    "-nan" into "0nan", which is not a valid number, and can also silently flip real
    negative Fst values (e.g. "-0.011" -> "00.011"). This script parses "-nan"/"nan"
    directly and safely; do not pre-process the file yourself.
"""

import argparse
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_fst(path):
    df = pd.read_csv(path, sep="\t")
    cols = set(df.columns)

    if {"BIN_START", "BIN_END", "WEIGHTED_FST", "MEAN_FST", "N_VARIANTS"}.issubset(cols):
        mode = "windowed"
        df["POS"] = (df["BIN_START"] + df["BIN_END"]) / 2
        fst_col = "WEIGHTED_FST"
    elif {"POS", "WEIR_AND_COCKERHAM_FST"}.issubset(cols):
        mode = "persite"
        fst_col = "WEIR_AND_COCKERHAM_FST"
    else:
        sys.exit(
            f"Could not recognize file format. Columns found: {list(df.columns)}\n"
            "Expected either windowed vcftools output (BIN_START, BIN_END, WEIGHTED_FST, ...) "
            "or per-site output (POS, WEIR_AND_COCKERHAM_FST)."
        )

    # Robustly coerce to numeric. vcftools writes "-nan"/"nan" for undefined sites;
    # pandas parses these natively as NaN as long as the file hasn't been hand-edited
    # (e.g. via sed) beforehand. errors="coerce" is a safety net for any other stray
    # non-numeric text, turning it into NaN rather than silently corrupting the column.
    df["FST"] = pd.to_numeric(df[fst_col], errors="coerce")

    n_total = len(df)
    n_undefined = df["FST"].isna().sum()
    if n_undefined:
        print(f"Note: {n_undefined} of {n_total} sites/windows had undefined Fst "
              f"('-nan'/'nan' in the vcftools output) and will be excluded from the plot.")

    return df, mode


def main():
    ap = argparse.ArgumentParser(description="Plot Fst distribution across the genome.")
    ap.add_argument("--input", required=True, help="vcftools Fst output file (windowed or per-site)")
    ap.add_argument("--out", default="fst_plot.png", help="Output image path (png)")
    ap.add_argument("--min-variants", type=int, default=1,
                     help="Minimum N_VARIANTS to keep a window (windowed mode only, default: 1)")
    ap.add_argument("--title", default=None, help="Custom plot title")
    ap.add_argument("--fst-hit-threshold", type=float, default=0.999,
                     help="Draw a vertical marker line for sites/windows at or above this "
                          "Fst value (default: 0.999)")
    ap.add_argument("--floor-negative", action="store_true",
                     help="Clamp genuine negative Fst values to 0 (does NOT affect "
                          "undefined '-nan' sites, which are always excluded, not zeroed)")
    args = ap.parse_args()

    df, mode = load_fst(args.input)

    if mode == "windowed" and "N_VARIANTS" in df.columns:
        before = len(df)
        df = df[df["N_VARIANTS"] >= args.min_variants]
        dropped = before - len(df)
        if dropped:
            print(f"Dropped {dropped} windows with fewer than {args.min_variants} variants.")

    df = df.dropna(subset=["FST"]).sort_values(["CHROM", "POS"])

    if args.floor_negative:
        n_neg = (df["FST"] < 0).sum()
        if n_neg:
            print(f"Flooring {n_neg} genuine negative Fst values to 0 (--floor-negative).")
            df.loc[df["FST"] < 0, "FST"] = 0.0

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [3, 1]})

    # --- Manhattan-style plot across genome position ---
    ax = axes[0]
    chroms = df["CHROM"].unique()
    multi_chrom = len(chroms) > 1
    offset = 0
    chrom_label_pos, chrom_labels = [], []
    all_x = []  # keep track of true genomic x-coordinate for every row, in df's row order

    for i, chrom in enumerate(chroms):
        sub = df[df["CHROM"] == chrom]
        # Only offset (concatenate) positions if there's more than one contig/chromosome.
        # With a single chromosome, x is just the real genomic position (POS) as-is.
        x = sub["POS"] + offset
        all_x.append(pd.Series(x.values, index=sub.index))
        ax.scatter(x, sub["FST"], s=8, alpha=0.5,
                   color="tab:blue" if i % 2 == 0 else "tab:orange",
                   zorder=2)
        if multi_chrom:
            chrom_label_pos.append(x.mean())
            chrom_labels.append(str(chrom))
            offset += sub["POS"].max() if len(sub) else 0

    x_all = pd.concat(all_x).sort_index()

    # --- mark sites/windows that hit (or exceed) the Fst threshold with vertical lines ---
    hits = df[df["FST"] >= args.fst_hit_threshold]
    hit_x = x_all.loc[hits.index]
    for xv in hit_x:
        ax.axvline(xv, color="red", linewidth=0.6, alpha=0.35, zorder=1)
    if len(hits):
        # dummy handle for the legend, since axvline doesn't auto-label well with many lines
        ax.axvline(hit_x.iloc[0], color="red", linewidth=0.6, alpha=0.35,
                   label=f"Fst \u2265 {args.fst_hit_threshold} (n={len(hits)})")
        ax.legend(loc="upper right", fontsize=8)

    ax.axhline(0, color="grey", linewidth=0.8, linestyle="--", zorder=1)

    if multi_chrom:
        # Numeric position ticks, plus chromosome name labels as a secondary annotation.
        ax.set_xlabel("Genomic position (concatenated across contigs)")
        # add faint vertical separators + labels at the top for each contig
        running = 0
        for i, chrom in enumerate(chroms):
            sub = df[df["CHROM"] == chrom]
            if i > 0:
                ax.axvline(running, color="black", linewidth=0.8, alpha=0.5, zorder=3)
            ax.text(running + (sub["POS"].max() if len(sub) else 0) / 2,
                    ax.get_ylim()[1] * 1.02, str(chrom),
                    ha="center", va="bottom", fontsize=8, rotation=0)
            running += sub["POS"].max() if len(sub) else 0
    else:
        ax.set_xlabel(f"Genomic position (bp) \u2014 {chroms[0]}")

    ax.set_ylabel("Fst" + (" (weighted, windowed)" if mode == "windowed" else " (per-site)"))
    ax.set_title(args.title or f"Fst across the genome ({mode})")
    ax.ticklabel_format(style="plain", axis="x")

    # --- Distribution histogram ---
    ax2 = axes[1]
    ax2.hist(df["FST"], bins=40, color="tab:blue", alpha=0.8)
    ax2.axvline(df["FST"].mean(), color="red", linestyle="--", linewidth=1,
                label=f"mean = {df['FST'].mean():.3f}")
    ax2.set_xlabel("Fst")
    ax2.set_ylabel("Count")
    ax2.set_title("Distribution")
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"Saved plot to {args.out}")
    print(f"Mode detected: {mode}")
    print(f"N points plotted: {len(df)}")
    print(f"Sites/windows with Fst >= {args.fst_hit_threshold}: {len(hits)}")
    print(f"Fst summary:\n{df['FST'].describe()}")


if __name__ == "__main__":
    main()