#!/usr/bin/env python3
"""
pangenome_heatmap.py

Tree + metadata strips + gene presence/absence heatmap for a Panaroo
pangenome, in the same visual style as amrfinder_summary.py (fig1).

Because a pangenome has thousands of clusters, genes must be SELECTED.
Default: the --top-n genes with the largest prevalence difference between
the levels of the main grouping (the first --color-by column), which is the
same ranking amrfinder_summary.py uses for its prevalence dot plot.
Columns are then grouped into blocks by which group they are enriched in,
the analogue of the AMR / Point / Stress / Virulence blocks.

Example:
    python3 pangenome_heatmap.py \
        --rtab      output/infantis/panaroo/gene_presence_absence.Rtab \
        --results   output/infantis/fisher_results_annotated.tsv \
        --tree      .../infantis_clean_tree_noambiguous.treefile \
        --midpoint \
        --metadata  .../09_AMRfinder_analysis/metadata.tsv \
        --color-by  pESI pESI_detail group \
        --top-n     60 \
        --outdir    output/infantis/pangenome_figs

Requirements: pandas, numpy, scipy, matplotlib, biopython (for --tree).
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list

# --------------------------------------------------------------------------
# Style (identical to amrfinder_summary.py)
# --------------------------------------------------------------------------
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER_GREY = "#a3a29c"
MISSING_GREY = "#e4e3df"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
ABSENT = "#f0efec"
PRESENT = "#1f2d3d"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": INK_2,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
})


def log(msg):
    print(msg, file=sys.stderr)


def clean_id(s):
    s = str(s).strip().strip("'\"")
    s = re.sub(r"\.(fna|fasta|fa|fas|gbk|gff|gff3)$", "", s, flags=re.I)
    return s


def group_colors(values):
    vals = pd.Series(values).dropna().astype(str)
    levels = sorted(vals.unique(), key=lambda x: str(x))
    if len(levels) > len(CATEGORICAL):
        top = vals.value_counts().index[:len(CATEGORICAL) - 1].tolist()
        levels = sorted(top)
    return dict((lv, CATEGORICAL[i]) for i, lv in enumerate(levels))


def map_group(values, cmap):
    out = []
    for v in values:
        if pd.isna(v) or str(v) in ("", "nan", "NA"):
            out.append("NA")
        elif str(v) in cmap:
            out.append(str(v))
        else:
            out.append("Other")
    return out


def color_of(level, cmap):
    if level == "NA":
        return MISSING_GREY
    if level == "Other":
        return OTHER_GREY
    return cmap[level]


# --------------------------------------------------------------------------
# Tree (identical logic to amrfinder_summary.py)
# --------------------------------------------------------------------------
def load_tree(path, midpoint):
    from Bio import Phylo
    tree = Phylo.read(path, "newick")
    if midpoint:
        tree.root_at_midpoint()
    tree.ladderize(reverse=True)
    for t in tree.get_terminals():
        t.name = clean_id(t.name)
    return tree


def tree_segments(tree):
    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}
    depths = tree.depths()
    if max(depths.values()) == 0:
        depths = tree.depths(unit_branch_lengths=True)

    def set_y(clade):
        if clade.is_terminal():
            return y[clade]
        ys = [set_y(c) for c in clade.clades]
        y[clade] = (min(ys) + max(ys)) / 2
        return y[clade]
    set_y(tree.root)

    segs = []
    for clade in tree.find_clades(order="preorder"):
        x0 = depths[clade]
        if not clade.is_terminal():
            ys = [y[c] for c in clade.clades]
            segs.append([(x0, min(ys)), (x0, max(ys))])
            for c in clade.clades:
                segs.append([(x0, y[c]), (depths[c], y[c])])
    return [t.name for t in tips], segs, max(depths.values())


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------
def read_rtab(path):
    """Panaroo Rtab: genes as rows, genomes as columns -> genomes x genes."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.columns = [clean_id(c) for c in df.columns]
    mat = (df.T > 0).astype(int)
    mat.index.name = "sample"
    return mat


def read_annotation(path, gene_col="gene"):
    """fisher_results_annotated.tsv -> {cluster_id: readable label}."""
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if gene_col not in df.columns:
        sys.exit(f"--results has no '{gene_col}' column; "
                 f"found {list(df.columns)[:8]}")
    ann = {}
    for _, r in df.iterrows():
        gid = r[gene_col]
        name = str(r.get("panaroo_non_unique_gene_name") or "").split(";")[0]
        prod = str(r.get("panaroo_annotation") or "").split(";")[0]
        ann[gid] = (name.strip(), prod.strip())
    return ann


def gene_label(gid, ann, maxlen=38):
    """Readable column label: real gene name if there is one, else a
    trimmed product description, else the cluster ID."""
    name, prod = ann.get(gid, ("", ""))
    if name and name.lower() not in ("nan", "none"):
        return name
    if gid and not gid.startswith("group_"):
        return gid
    if prod and prod.lower() not in ("", "nan", "hypothetical protein",
                                     "not_found_in_panaroo_csv"):
        p = prod if len(prod) <= maxlen else prod[:maxlen - 1] + "\u2026"
        return f"{p} [{gid}]"
    return gid


# --------------------------------------------------------------------------
# Gene selection
# --------------------------------------------------------------------------
def select_genes(mat, meta, main, top_n, min_genomes, max_prev,
                 genes_file=None):
    """Return (ordered gene list, block label per gene, prevalence table)."""
    n = len(mat)
    keep = mat.columns[(mat.sum() >= min_genomes) &
                       (mat.sum() <= max_prev * n)]
    if len(keep) == 0:
        sys.exit("no genes left after --min-genomes / --max-prevalence")
    sub = mat[keep]

    if genes_file:
        wanted = [l.strip() for l in open(genes_file) if l.strip()]
        missing = [g for g in wanted if g not in mat.columns]
        if missing:
            log(f"    {len(missing)} genes from --genes not in the Rtab "
                f"(e.g. {missing[:3]})")
        feats = [g for g in wanted if g in mat.columns]
        prev = pd.DataFrame({"n_genomes": mat[feats].sum()})
        return feats, pd.Series("selected", index=feats), prev

    if not main:
        # no grouping: just take the most variable genes
        p = sub.mean() * 100
        score = (p - 50).abs().sort_values()
        feats = list(score.index[:top_n])
        prev = pd.DataFrame({"pct_all": p[feats].round(2)})
        return feats, pd.Series("variable", index=feats), prev

    cmap = group_colors(meta[main])
    lv = pd.Series(map_group(meta[main].reindex(sub.index), cmap),
                   index=sub.index)
    levels = [l for l in cmap if (lv == l).any()]
    if len(levels) < 2:
        sys.exit(f"--color-by '{main}' has fewer than 2 populated levels")

    prev = pd.DataFrame({l: sub[lv == l].mean() * 100 for l in levels})
    prev["max_diff_pct"] = (prev.max(axis=1) - prev.min(axis=1))
    prev["enriched_in"] = prev[levels].idxmax(axis=1)
    prev = prev.sort_values("max_diff_pct", ascending=False)

    feats = list(prev.index[:top_n])
    blocks = prev.loc[feats, "enriched_in"]
    # order: block (by level order), then by difference within block
    order = sorted(feats, key=lambda g: (levels.index(blocks[g]),
                                         -prev.loc[g, "max_diff_pct"]))
    return order, blocks.reindex(order), prev.round(2)


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
def plot(mat, feats, blocks, labels, tree_path, midpoint, meta, color_by,
         main, outdir, title_extra=""):
    if tree_path:
        tree = load_tree(tree_path, midpoint)
        tips, segs, maxdepth = tree_segments(tree)
        no_data = [t for t in tips if t not in mat.index]
        no_tree = [s for s in mat.index if s not in set(tips)]
        if no_data:
            log(f"    {len(no_data)} tree tips absent from the Rtab "
                f"(e.g. {no_data[:3]}) - drawn as empty rows")
        if no_tree:
            log(f"    {len(no_tree)} genomes absent from the tree "
                f"(e.g. {no_tree[:3]}) - left out of this figure")
        rows = tips
    else:
        log("    no --tree: ordering genomes by Jaccard clustering")
        X = mat[feats].values.astype(bool)
        d = np.nan_to_num(pdist(X, "jaccard"))
        rows = list(mat.index[leaves_list(linkage(d, "average"))])
        segs, maxdepth = None, None

    M = mat.reindex(index=rows, columns=feats).fillna(0).values
    n_rows, n_cols = M.shape

    n_strips = len(color_by)
    w_tree = 3.2 if segs else 0.0
    w_strip = 0.22 * n_strips
    w_heat = max(4.0, 0.16 * n_cols)
    height = min(max(6.0, 0.012 * n_rows + 3), 16)
    widths = ([w_tree] if segs else []) + ([w_strip] if n_strips else []) \
        + [w_heat]
    fig = plt.figure(figsize=(sum(widths) + 3.5, height))
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=0.02)
    k = 0

    if segs:
        axt = fig.add_subplot(gs[0, k]); k += 1
        axt.add_collection(LineCollection(segs, colors=INK_2, linewidths=0.4))
        axt.set_xlim(0, maxdepth * 1.02)
        axt.set_ylim(n_rows - 0.5, -0.5)
        axt.axis("off")
        sb = float(f"{maxdepth / 5:.1g}")
        axt.plot([0, sb], [n_rows + n_rows * 0.01] * 2, color=INK_2, lw=1,
                 clip_on=False)
        axt.text(sb / 2, n_rows + n_rows * 0.025, f"{sb:g}", ha="center",
                 va="top", fontsize=7, color=INK_2)

    legend_handles = []
    if n_strips:
        axs = fig.add_subplot(gs[0, k]); k += 1
        rgb = np.ones((n_rows, n_strips, 3))
        for j, col in enumerate(color_by):
            cmap = group_colors(meta[col])
            lv = map_group(meta[col].reindex(rows), cmap)
            rgb[:, j, :] = [to_rgb(color_of(l, cmap)) for l in lv]
            legend_handles.append(Line2D(
                [], [], ls="",
                label=f"$\\bf{{{col.replace('_', chr(92) + '_')}}}$"))
            for l, c in cmap.items():
                legend_handles.append(Patch(fc=c, label=l))
            if "Other" in lv:
                legend_handles.append(Patch(fc=OTHER_GREY, label="Other"))
            if "NA" in lv:
                legend_handles.append(Patch(fc=MISSING_GREY, label="NA"))
        axs.imshow(rgb, aspect="auto", interpolation="nearest")
        axs.set_xticks(range(n_strips))
        axs.set_xticklabels(color_by, rotation=90, fontsize=8)
        axs.xaxis.tick_top()
        axs.set_yticks([])
        for s in axs.spines.values():
            s.set_visible(False)

    axh = fig.add_subplot(gs[0, k])
    rgb = np.empty((n_rows, n_cols, 3))
    rgb[:] = to_rgb(ABSENT)
    rgb[M.astype(bool)] = to_rgb(PRESENT)
    axh.imshow(rgb, aspect="auto", interpolation="nearest")
    axh.set_yticks([])
    axh.set_xticks(range(n_cols))
    axh.set_xticklabels([labels[g] for g in feats], rotation=90,
                        fontsize=6.5, fontstyle="italic")
    axh.xaxis.tick_top()
    axh.tick_params(axis="x", length=0, pad=2)
    for s in axh.spines.values():
        s.set_visible(False)
    for x in np.arange(-0.5, n_cols, 1):
        axh.axvline(x, color="white", lw=0.4)

    # block separators + labels (enriched-in group)
    bl = [str(blocks[g]) for g in feats]
    start = 0
    for i in range(1, n_cols + 1):
        if i == n_cols or bl[i] != bl[start]:
            if i < n_cols:
                axh.axvline(i - 0.5, color=INK, lw=1.2)
            axh.text((start + i - 1) / 2, n_rows - 0.5 + n_rows * 0.012,
                     f"enriched in {bl[start]}" if main else bl[start],
                     ha="center", va="top", fontsize=8, color=INK,
                     fontweight="bold")
            start = i

    legend_handles += [Line2D([], [], ls="", label=r"$\bf{Gene}$"),
                       Patch(fc=PRESENT, label="present"),
                       Patch(fc=ABSENT, ec=GRID, label="absent")]
    axh.legend(handles=legend_handles, loc="upper left",
               bbox_to_anchor=(1.01, 1), fontsize=8, handlelength=1.2)
    fig.suptitle(f"Pangenome gene presence across {n_rows} genomes"
                 f"{' (ordered by core-genome tree)' if segs else ''}",
                 x=0.01, ha="left", va="bottom", fontsize=12,
                 fontweight="bold", y=1.035)
    fig.text(0.01, 1.01, title_extra, ha="left", va="bottom",
             fontsize=8.5, color=INK_2)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"fig1_pangenome_heatmap.{ext}"))
    plt.close(fig)
    pd.Series(rows, name="sample").to_csv(
        os.path.join(outdir, "fig1_row_order.txt"), index=False)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--rtab", required=True,
                    help="Panaroo gene_presence_absence.Rtab")
    ap.add_argument("--results",
                    help="fisher_results_annotated.tsv, for readable labels")
    ap.add_argument("--gene-col", default="gene",
                    help="gene-ID column in --results (default: gene)")
    ap.add_argument("--tree", help="Newick core-genome tree")
    ap.add_argument("--midpoint", action="store_true")
    ap.add_argument("--metadata", help="TSV, header, first column = sample")
    ap.add_argument("--color-by", nargs="+", default=[])
    ap.add_argument("--genes", help="file with one gene/cluster ID per line; "
                                    "overrides the automatic selection")
    ap.add_argument("--top-n", type=int, default=60,
                    help="genes to draw (default 60)")
    ap.add_argument("--min-genomes", type=int, default=5,
                    help="skip genes in fewer genomes than this (default 5)")
    ap.add_argument("--max-prevalence", type=float, default=0.99,
                    help="skip genes above this fraction of genomes "
                         "(default 0.99, i.e. drop the core)")
    ap.add_argument("--outdir", default="pangenome_figs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    log("[1] reading Rtab")
    mat = read_rtab(args.rtab)
    log(f"    {mat.shape[0]} genomes x {mat.shape[1]} gene clusters")

    meta = pd.DataFrame(index=mat.index)
    if args.metadata:
        m = pd.read_csv(args.metadata, sep="\t", dtype=str)
        m.iloc[:, 0] = m.iloc[:, 0].map(clean_id)
        m = m.set_index(m.columns[0])
        miss = [c for c in args.color_by if c not in m.columns]
        if miss:
            sys.exit(f"--color-by columns not in metadata: {miss}; "
                     f"available: {list(m.columns)}")
        n_match = len(set(m.index) & set(mat.index))
        log(f"    metadata matched for {n_match}/{len(mat)} genomes")
        if n_match == 0:
            sys.exit("no metadata sample IDs match the Rtab column names")
        meta = m
    elif args.color_by:
        sys.exit("--color-by needs --metadata")
    main_grp = args.color_by[0] if args.color_by else None

    log("[2] selecting genes")
    feats, blocks, prev = select_genes(mat, meta, main_grp, args.top_n,
                                       args.min_genomes, args.max_prevalence,
                                       args.genes)
    log(f"    {len(feats)} genes selected")
    prev.to_csv(os.path.join(args.outdir, "gene_prevalence_by_group.tsv"),
                sep="\t")

    ann = read_annotation(args.results, args.gene_col)
    labels = {g: gene_label(g, ann) for g in feats}
    pd.DataFrame({"gene": feats,
                  "label": [labels[g] for g in feats],
                  "block": [blocks[g] for g in feats]}).to_csv(
        os.path.join(args.outdir, "fig1_gene_order.tsv"),
        sep="\t", index=False)

    note = (f"{len(feats)} genes with the largest prevalence difference "
            f"between {main_grp} groups" if main_grp and not args.genes
            else f"{len(feats)} selected genes")
    if args.max_prevalence < 1:
        note += f"; genes in >{100*args.max_prevalence:.0f}% of genomes removed"

    log("[3] drawing")
    plot(mat, feats, blocks, labels, args.tree, args.midpoint, meta,
         args.color_by, main_grp, args.outdir, note)
    log(f"done -> {args.outdir}/")


if __name__ == "__main__":
    main()
