#!/usr/bin/env python3
"""
amrfinder_summary.py

Summarise AMRFinderPlus results across many genomes:

  1) merge all per-genome output tables into one long table (+ `sample` column)
  2) build presence/absence matrices on `Element symbol`
       - Panaroo/Roary-style .Rtab (feature rows x genome columns, 0/1),
         directly usable with fisher_gene_test.py and pyseer --pres
       - one matrix for everything + one per element type (AMR, STRESS,
         VIRULENCE) + one for point mutations only
  3) tree + heatmap figure (core-genome tree, metadata strip(s), hit matrix)
  4) PCoA on Jaccard distances (+ PERMANOVA), prevalence-by-group dot plot,
     per-genome resistance burden plot, feature co-occurrence heatmap

Requirements: python>=3.8, pandas, numpy, scipy, matplotlib, biopython
(biopython only needed when --tree is given)

Example:
    python3 amrfinder_summary.py \
        --indir     amrfinder_out/ \
        --tree      core_gene_alignment.treefile \
        --metadata  metadata.tsv \
        --color-by  pESI clade \
        --outdir    amrfinder_summary/

    metadata.tsv: tab-separated, header row, FIRST column = sample ID
    (must match the AMRFinderPlus file names minus the suffix, and the tree
    tip labels), other columns = categorical groupings, e.g.
        sample      pESI    clade
        SRR123456   1       Sweep7
        ...

The first column given to --color-by is the "main grouping" used for the
PCoA colours, PERMANOVA, prevalence and burden plots. Every --color-by
column is drawn as a strip next to the tree.
"""

import argparse
import glob
import os
import re
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, to_rgb
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
# Categorical slots in fixed order (never cycled; >8 groups fold into "Other")
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER_GREY = "#a3a29c"
MISSING_GREY = "#e4e3df"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
ABSENT = "#f0efec"
PRESENT = "#1f2d3d"   # dark ink for "present" cells

# Element type -> display order and label
TYPE_ORDER = ["AMR", "POINT", "STRESS", "VIRULENCE"]
TYPE_LABEL = {"AMR": "AMR genes", "POINT": "Point mutations",
              "STRESS": "Stress (metal / biocide / heat / acid)",
              "VIRULENCE": "Virulence"}

SHORT_LABEL = {"AMR": "AMR", "POINT": "Point", "STRESS": "Stress",
               "VIRULENCE": "Virulence"}

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

# Column names differ between AMRFinderPlus v3 and v4 -> normalise to v4
COLUMN_ALIASES = {
    "Gene symbol": "Element symbol",
    "Sequence name": "Element name",
    "Element type": "Type",
    "Element subtype": "Subtype",
    "Contig id": "Contig id",
}

PARTIAL_METHODS = {"PARTIALP", "PARTIALX", "PARTIAL_CONTIG_ENDP",
                   "PARTIAL_CONTIG_ENDX", "INTERNAL_STOP"}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def log(msg):
    print(msg, file=sys.stderr)


def clean_id(s):
    """Normalise a sample / tip label so file names and tree tips match."""
    s = str(s).strip().strip("'\"")
    s = re.sub(r"\.(fna|fasta|fa|fas|gbk|gff)$", "", s, flags=re.I)
    return s


def group_colors(values):
    """Fixed-order categorical colours for the levels of a grouping.
    Levels are sorted (so colours are stable between runs); >8 levels: the
    8 most common keep a colour, the rest become 'Other' (grey)."""
    vals = pd.Series(values).dropna().astype(str)
    levels = sorted(vals.unique(), key=lambda x: (str(x)))
    if len(levels) > len(CATEGORICAL):
        top = vals.value_counts().index[:len(CATEGORICAL) - 1].tolist()
        levels = sorted(top)
    cmap = OrderedDict((lv, CATEGORICAL[i]) for i, lv in enumerate(levels))
    return cmap


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


def write_rtab(mat, path):
    """mat: genomes x features (0/1). Writes Panaroo-style Rtab
    (features as rows, genomes as columns, first header cell 'Gene')."""
    t = mat.T.copy()
    t.index.name = "Gene"
    t.to_csv(path, sep="\t")


# --------------------------------------------------------------------------
# 1) Merge
# --------------------------------------------------------------------------
def merge_tables(indir, pattern, suffix):
    files = sorted(glob.glob(os.path.join(indir, pattern)))
    if not files:
        sys.exit(f"No files matching {pattern} in {indir}")
    frames, samples = [], []
    for f in files:
        name = os.path.basename(f)
        sample = re.sub(re.escape(suffix) + r"(\.txt|\.tsv)?$", "", name)
        sample = clean_id(re.sub(r"\.(txt|tsv)$", "", sample))
        samples.append(sample)
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
        except pd.errors.EmptyDataError:
            log(f"  WARNING: {name} is completely empty (no header) - "
                f"kept as a genome with zero hits, but check that run")
            continue
        df = df.rename(columns=COLUMN_ALIASES)
        if "Name" in df.columns:          # from --name; replace with ours
            df = df.drop(columns="Name")
        df.insert(0, "sample", sample)
        frames.append(df)
    if len(set(samples)) != len(samples):
        dup = pd.Series(samples)[pd.Series(samples).duplicated()].tolist()
        sys.exit(f"Duplicate sample names after stripping suffix: {dup[:5]}")
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    need = {"Element symbol", "Type"}
    if not need.issubset(merged.columns):
        sys.exit(f"Expected columns {need}; found {list(merged.columns)}")
    return merged, samples


def point_mask(df):
    m = pd.Series(False, index=df.index)
    if "Subtype" in df.columns:
        m |= df["Subtype"].fillna("").str.upper().eq("POINT")
    if "Method" in df.columns:
        m |= df["Method"].fillna("").str.upper().str.startswith("POINT")
    return m


# --------------------------------------------------------------------------
# 2) Matrices
# --------------------------------------------------------------------------
def build_matrices(merged, samples, outdir):
    df = merged.copy()
    df["feature_class"] = np.where(point_mask(df), "POINT", df["Type"])

    mat = (pd.crosstab(df["sample"], df["Element symbol"])
             .reindex(samples, fill_value=0)
             .clip(upper=1).astype(int))
    mat.index.name = "sample"

    # feature annotation: the most common annotation per symbol
    def mode(s):
        s = s.dropna()
        return s.mode().iloc[0] if len(s) else ""
    cols = [c for c in ["feature_class", "Type", "Subtype", "Class",
                        "Subclass", "Element name"] if c in df.columns]
    info = df.groupby("Element symbol")[cols].agg(mode)
    info["n_genomes"] = mat.sum(axis=0).reindex(info.index)
    info["pct_genomes"] = (100 * info["n_genomes"] / len(samples)).round(2)
    info["n_hits_total"] = df.groupby("Element symbol").size()
    info = info.sort_values(["feature_class", "Class", "n_genomes"],
                            ascending=[True, True, False]) \
        if "Class" in info.columns else info
    info.index.name = "Element symbol"
    info.to_csv(os.path.join(outdir, "feature_info.tsv"), sep="\t")

    os.makedirs(os.path.join(outdir, "matrices"), exist_ok=True)
    write_rtab(mat, os.path.join(outdir, "matrices", "all_features.Rtab"))
    mat.to_csv(os.path.join(outdir, "matrices",
                            "all_features_genome_x_feature.tsv"), sep="\t")
    for fc in TYPE_ORDER:
        feats = info.index[info["feature_class"] == fc]
        if len(feats):
            write_rtab(mat[feats],
                       os.path.join(outdir, "matrices", f"{fc}.Rtab"))
    return mat, info, df


# --------------------------------------------------------------------------
# Tree utilities (Biopython)
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
    """Return tip order (top->bottom), and line segments for a rectangular
    phylogram with tips at y = 0..n-1."""
    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}
    depths = tree.depths()
    if max(depths.values()) == 0:            # cladogram without lengths
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
            segs.append([(x0, min(ys)), (x0, max(ys))])        # vertical
            for c in clade.clades:
                segs.append([(x0, y[c]), (depths[c], y[c])])  # horizontal
    return [t.name for t in tips], segs, max(depths.values())


# --------------------------------------------------------------------------
# 3) Tree + heatmap
# --------------------------------------------------------------------------
def heatmap_feature_order(info, mat, min_genomes, drop_core):
    n = len(mat)
    keep = info[(info["n_genomes"] >= min_genomes)]
    if drop_core:
        keep = keep[keep["n_genomes"] < 0.99 * n]
    order = []
    for fc in TYPE_ORDER:
        sub = keep[keep["feature_class"] == fc]
        sort_cols = ["Class", "n_genomes"] if "Class" in sub.columns \
            else ["n_genomes"]
        asc = [True, False] if len(sort_cols) == 2 else [False]
        order.extend(sub.sort_values(sort_cols, ascending=asc).index)
    return order


def plot_tree_heatmap(mat, info, tree_path, midpoint, meta, color_by,
                      min_genomes, drop_core, outdir):
    if tree_path:
        tree = load_tree(tree_path, midpoint)
        tips, segs, maxdepth = tree_segments(tree)
        missing_amr = [t for t in tips if t not in mat.index]
        missing_tree = [s for s in mat.index if s not in set(tips)]
        if missing_amr:
            log(f"  {len(missing_amr)} tree tips have no AMRFinderPlus file "
                f"(e.g. {missing_amr[:3]}) - drawn as empty rows")
        if missing_tree:
            log(f"  {len(missing_tree)} genomes are not in the tree "
                f"(e.g. {missing_tree[:3]}) - left out of this figure")
        rows = tips
    else:
        # no tree: order genomes by hierarchical clustering of their profile
        log("  no --tree given: ordering genomes by Jaccard clustering")
        X = mat.values.astype(bool)
        d = np.nan_to_num(pdist(X, "jaccard"))
        rows = list(mat.index[leaves_list(linkage(d, "average"))])
        segs, maxdepth = None, None

    feats = heatmap_feature_order(info, mat, min_genomes, drop_core)
    M = mat.reindex(index=rows, columns=feats).fillna(0).values
    n_rows, n_cols = M.shape

    # layout
    n_strips = len(color_by)
    w_tree = 3.2 if segs else 0.0
    w_strip = 0.22 * n_strips
    w_heat = max(4.0, 0.16 * n_cols)
    height = min(max(6.0, 0.012 * n_rows + 3), 16)
    widths = ([w_tree] if segs else []) + ([w_strip] if n_strips else []) \
        + [w_heat]
    fig = plt.figure(figsize=(sum(widths) + 2.5, height))
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=0.02)
    k = 0

    if segs:
        axt = fig.add_subplot(gs[0, k]); k += 1
        axt.add_collection(LineCollection(segs, colors=INK_2, linewidths=0.4))
        axt.set_xlim(0, maxdepth * 1.02)
        axt.set_ylim(n_rows - 0.5, -0.5)
        axt.axis("off")
        # scale bar
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
            legend_handles.append(Line2D([], [], ls="", label=f"$\\bf{{{col.replace('_', chr(92) + '_')}}}$"))
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
    axh.set_xticklabels(feats, rotation=90, fontsize=6.5,
                        fontstyle="italic")
    axh.xaxis.tick_top()
    axh.tick_params(axis="x", length=0, pad=2)
    for s in axh.spines.values():
        s.set_visible(False)
    # thin white separators between columns + type blocks
    for x in np.arange(-0.5, n_cols, 1):
        axh.axvline(x, color="white", lw=0.4)
    fcs = info.loc[feats, "feature_class"].tolist()
    start = 0
    for i in range(1, n_cols + 1):
        if i == n_cols or fcs[i] != fcs[start]:
            if i < n_cols:
                axh.axvline(i - 0.5, color=INK, lw=1.2)
            axh.text((start + i - 1) / 2, n_rows - 0.5 + n_rows * 0.012,
                     SHORT_LABEL.get(fcs[start], fcs[start]),
                     ha="center", va="top", fontsize=8, color=INK,
                     fontweight="bold")
            start = i

    legend_handles += [Line2D([], [], ls="", label=r"$\bf{Hit}$"),
                       Patch(fc=PRESENT, label="present"),
                       Patch(fc=ABSENT, ec=GRID, label="absent")]
    axh.legend(handles=legend_handles, loc="upper left",
               bbox_to_anchor=(1.01, 1), fontsize=8, handlelength=1.2)
    fig.suptitle(f"AMRFinderPlus hits across {n_rows} genomes"
                 f"{' (ordered by core-genome tree)' if segs else ''}",
                 x=0.01, ha="left", va="bottom", fontsize=12, fontweight="bold",
                 y=1.035)
    fig.text(0.01, 1.01, f"Features found in ≥{min_genomes} genomes"
             f"{'; features in ≥99% of genomes removed' if drop_core else ''}",
             ha="left", va="bottom", fontsize=8.5, color=INK_2)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"fig1_tree_heatmap.{ext}"))
    plt.close(fig)
    pd.Series(rows, name="sample").to_csv(
        os.path.join(outdir, "fig1_row_order.txt"), index=False)


# --------------------------------------------------------------------------
# 4a) PCoA + PERMANOVA
# --------------------------------------------------------------------------
def pcoa(D):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = vals > 1e-10
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    expl = vals[pos] / vals[pos].sum()
    return coords, expl


def permanova(D, groups, n_perm=999, seed=1):
    """One-way PERMANOVA (Anderson 2001). Returns pseudo-F, R2, p."""
    rng = np.random.default_rng(seed)
    groups = np.asarray(groups)
    n = len(groups)
    D2 = D ** 2
    ss_t = D2[np.triu_indices(n, 1)].sum() / n

    def ss_w(g):
        s = 0.0
        for lv in np.unique(g):
            idx = np.where(g == lv)[0]
            if len(idx) > 1:
                sub = D2[np.ix_(idx, idx)]
                s += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return s
    a = len(np.unique(groups))
    w = ss_w(groups)
    F = ((ss_t - w) / (a - 1)) / (w / (n - a))
    R2 = (ss_t - w) / ss_t
    perm_F = np.empty(n_perm)
    for i in range(n_perm):
        wp = ss_w(rng.permutation(groups))
        perm_F[i] = ((ss_t - wp) / (a - 1)) / (wp / (n - a))
    p = (np.sum(perm_F >= F) + 1) / (n_perm + 1)
    return F, R2, p


def plot_pcoa(mat, meta, main, outdir, n_perm):
    var = mat.loc[:, (mat.sum() > 0) & (mat.sum() < len(mat))]
    X = var.values.astype(bool)
    D = squareform(np.nan_to_num(pdist(X, "jaccard")))
    coords, expl = pcoa(D)
    if coords.shape[1] < 2:
        log("  PCoA: fewer than 2 positive axes - skipped")
        return
    pc = pd.DataFrame(coords[:, :min(5, coords.shape[1])], index=mat.index,
                      columns=[f"PCo{i+1}" for i in
                               range(min(5, coords.shape[1]))])

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    stats_txt = ""
    if main:
        cmap = group_colors(meta[main])
        lv = pd.Series(map_group(meta[main].reindex(mat.index), cmap),
                       index=mat.index)
        pc[main] = lv
        levels = list(cmap) + [l for l in ("Other", "NA") if l in set(lv)]
        ok = lv != "NA"
        if ok.sum() > 2 and lv[ok].nunique() > 1:
            idx = np.where(ok.values)[0]
            F, R2, p = permanova(D[np.ix_(idx, idx)], lv[ok].values, n_perm)
            stats_txt = (f"PERMANOVA ({main}): pseudo-F = {F:.1f}, "
                         f"R² = {R2:.2f}, p = {p:.3g} ({n_perm} perm.)")
    else:
        lv = pd.Series("all", index=mat.index)
        cmap = {"all": CATEGORICAL[0]}
        levels = ["all"]

    # many genomes share identical profiles -> one marker per unique
    # (position, group), area proportional to the number of genomes
    tmp = pc[["PCo1", "PCo2"]].round(6).assign(g=lv.values)
    agg = tmp.groupby(["PCo1", "PCo2", "g"]).size().reset_index(name="n")
    for l in levels:
        sub = agg[agg["g"] == l]
        if sub.empty:
            continue
        ax.scatter(sub["PCo1"], sub["PCo2"], s=18 + 14 * np.sqrt(sub["n"]),
                   color=color_of(l, cmap) if l in cmap or l in ("NA", "Other")
                   else OTHER_GREY,
                   alpha=0.75, edgecolors="white", linewidths=0.8,
                   label=f"{l} (n={int((lv == l).sum())})")
    ax.set_xlabel(f"PCo1 ({100*expl[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({100*expl[1]:.1f}%)")
    ax.axhline(0, color=GRID, lw=0.8, zorder=0)
    ax.axvline(0, color=GRID, lw=0.8, zorder=0)
    ax.set_title("PCoA of AMRFinderPlus profiles (Jaccard)")
    leg = ax.legend(title=main if main else None, loc="best", fontsize=8,
                    markerscale=0.8)
    if leg.get_title():
        leg.get_title().set_fontweight("bold")
    note = (f"{var.shape[1]} variable features; marker area ∝ number of "
            f"genomes with an identical profile")
    if stats_txt:
        note = stats_txt + "\n" + note
    ax.text(0.0, -0.13, note, fontsize=7.5, color=INK_2, ha="left",
             va="top", transform=ax.transAxes)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"fig2_pcoa_jaccard.{ext}"))
    plt.close(fig)
    pc.to_csv(os.path.join(outdir, "pcoa_coordinates.tsv"), sep="\t")
    with open(os.path.join(outdir, "pcoa_permanova.txt"), "w") as fh:
        fh.write("Variance explained (positive eigenvalues):\n")
        for i, e in enumerate(expl[:10]):
            fh.write(f"  PCo{i+1}\t{100*e:.2f}%\n")
        if stats_txt:
            fh.write("\n" + stats_txt + "\n")
            fh.write("NOTE: permutations assume genomes are exchangeable. "
                     "Clonal lineages violate this, so the p-value is "
                     "anti-conservative; treat R² as descriptive and use "
                     "lineage-aware tests (pyseer) for inference.\n")


# --------------------------------------------------------------------------
# 4b) Prevalence by group (Cleveland dot plot) + table
# --------------------------------------------------------------------------
def plot_prevalence(mat, info, meta, main, outdir, top_n):
    cmap = group_colors(meta[main])
    lv = pd.Series(map_group(meta[main].reindex(mat.index), cmap),
                   index=mat.index)
    levels = [l for l in cmap if (lv == l).any()]
    prev = pd.DataFrame({l: mat[lv == l].mean() * 100 for l in levels})
    counts = pd.DataFrame({f"n_{l}": mat[lv == l].sum() for l in levels})
    sizes = {l: int((lv == l).sum()) for l in levels}
    table = pd.concat([info[["feature_class"] +
                            (["Class"] if "Class" in info else [])],
                       counts, prev.add_prefix("pct_").round(2)], axis=1)
    table["max_diff_pct"] = (prev.max(axis=1) - prev.min(axis=1)).round(2)
    table = table.sort_values("max_diff_pct", ascending=False)
    table.index.name = "Element symbol"
    table.to_csv(os.path.join(outdir, f"prevalence_by_{main}.tsv"), sep="\t")

    top = table.head(top_n)
    top = top.sort_values(["feature_class", "max_diff_pct"],
                          ascending=[True, True])
    feats = list(top.index)
    if not feats:
        return
    fig, ax = plt.subplots(figsize=(6.8, 0.22 * len(feats) + 1.6))
    y = np.arange(len(feats))
    ax.hlines(y, prev.loc[feats].min(axis=1), prev.loc[feats].max(axis=1),
              color=GRID, lw=2, zorder=1)
    for l in levels:
        ax.scatter(prev.loc[feats, l], y, s=40, color=cmap[l],
                   edgecolors="white", linewidths=1, zorder=3,
                   label=f"{l} (n={sizes[l]})")
    ax.set_yticks(y)
    ax.set_yticklabels(feats, fontstyle="italic", fontsize=8)
    for t, f in zip(ax.get_yticklabels(), feats):
        t.set_color(INK)
    # type labels on the right
    fcs = top["feature_class"].tolist()
    start = 0
    for i in range(1, len(feats) + 1):
        if i == len(feats) or fcs[i] != fcs[start]:
            if i < len(feats):
                ax.axhline(i - 0.5, color=INK_2, lw=0.6, ls=":")
            ax.text(101, (start + i - 1) / 2, fcs[start], va="center",
                    fontsize=8, color=INK_2, fontweight="bold")
            start = i
    ax.set_xlim(-2, 100)
    ax.set_xlabel("% of genomes carrying the feature")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f"Feature prevalence by {main}", pad=20)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.6 / len(feats)
              - 0.06), ncol=min(4, len(levels)), fontsize=8)
    ax.text(0, 1.012, f"{len(feats)} features with the largest between-group "
            f"difference", transform=ax.transAxes, fontsize=8, color=INK_2)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"fig3_prevalence_by_{main}.{ext}"))
    plt.close(fig)


# --------------------------------------------------------------------------
# 4c) Resistance burden per genome
# --------------------------------------------------------------------------
def burden_table(long_df, samples, ignore=()):
    amr = long_df[(long_df["Type"] == "AMR") &
                  (~long_df["Element symbol"].isin(ignore))].copy()
    rows = []
    for s in samples:
        sub = amr[amr["sample"] == s]
        # intrinsic chromosomal efflux (mdsA/mdsB) is not an acquired gene
        acq = (sub["feature_class"] == "AMR")
        if "Class" in sub.columns:
            acq &= sub["Class"].fillna("").ne("EFFLUX")
        genes = sub.loc[acq, "Element symbol"]
        points = sub.loc[sub["feature_class"] == "POINT", "Element symbol"]
        classes = set()
        if "Class" in sub.columns:
            for c in sub["Class"].dropna():
                classes.update(x.strip() for x in c.split("/") if x.strip())
        classes.discard("EFFLUX")      # efflux pumps are not a drug class
        rows.append({"sample": s, "n_AMR_genes": genes.nunique(),
                     "n_point_mutations": points.nunique(),
                     "n_drug_classes": len(classes),
                     "MDR_3plus_classes": int(len(classes) >= 3),
                     "drug_classes": ";".join(sorted(classes))})
    return pd.DataFrame(rows).set_index("sample")


def plot_burden(burden, meta, main, outdir):
    cmap = group_colors(meta[main])
    lv = pd.Series(map_group(meta[main].reindex(burden.index), cmap),
                   index=burden.index)
    levels = [l for l in cmap if (lv == l).any()]
    metrics = [("n_AMR_genes", "Acquired AMR genes"),
               ("n_drug_classes", "Drug classes affected"),
               ("n_point_mutations", "Resistance point mutations")]
    fig, axes = plt.subplots(1, 3, figsize=(10, 4.2))
    rng = np.random.default_rng(0)
    for ax, (col, title) in zip(axes, metrics):
        for i, l in enumerate(levels):
            v = burden.loc[lv == l, col].values
            # integer counts: show distribution as a jittered strip + median
            ax.scatter(i + rng.uniform(-0.28, 0.28, len(v)),
                       v + rng.uniform(-0.15, 0.15, len(v)), s=8,
                       color=cmap[l], alpha=0.45, linewidths=0)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([i - 0.36, i + 0.36], [med, med], color=INK, lw=2)
            ax.plot([i, i], [q1, q3], color=INK, lw=1)
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels([f"{l}\n(n={(lv == l).sum()})" for l in levels],
                           fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.yaxis.get_major_locator().set_params(integer=True)
    mdr = burden.groupby(lv)["MDR_3plus_classes"].mean().mul(100)
    txt = "MDR (≥3 drug classes): " + ", ".join(
        f"{l} {mdr[l]:.0f}%" for l in levels if l in mdr)
    fig.suptitle(f"Resistance burden per genome by {main}", x=0.01,
                 ha="left", fontweight="bold", fontsize=12)
    fig.text(0.01, -0.02, txt + "   ·   black bar = median, line = IQR",
             fontsize=8, color=INK_2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"fig4_burden_by_{main}.{ext}"))
    plt.close(fig)


# --------------------------------------------------------------------------
# 4d) Co-occurrence (phi correlation) between variable features
# --------------------------------------------------------------------------
def plot_cooccurrence(mat, info, min_genomes, outdir, max_feats=80):
    n = len(mat)
    keep = info[(info["n_genomes"] >= min_genomes) &
                (info["n_genomes"] <= n - min_genomes)].index
    if len(keep) < 3:
        log("  co-occurrence: <3 variable features - skipped")
        return
    if len(keep) > max_feats:
        keep = info.loc[keep].sort_values("n_genomes",
                                          ascending=False).index[:max_feats]
    X = mat[keep].values.astype(float)
    r = np.corrcoef(X.T)               # Pearson on 0/1 = phi coefficient
    r = np.nan_to_num(r)
    order = leaves_list(linkage(squareform(1 - r, checks=False) .clip(0),
                                "average"))
    r = r[np.ix_(order, order)]
    labels = np.array(keep)[order]
    cm = LinearSegmentedColormap.from_list(
        "div", ["#2a78d6", "#f0efec", "#eb6834"])
    sz = max(5, 0.16 * len(labels) + 2)
    fig, ax = plt.subplots(figsize=(sz + 1.2, sz))
    im = ax.imshow(r, cmap=cm, vmin=-1, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    fs = 7 if len(labels) < 50 else 5.5
    ax.set_xticklabels(labels, rotation=90, fontsize=fs, fontstyle="italic")
    ax.set_yticklabels(labels, fontsize=fs, fontstyle="italic")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("phi correlation (co-occurrence across genomes)")
    cb.outline.set_visible(False)
    ax.set_title("Which features travel together?", pad=20)
    ax.text(0, 1.012, "Blocks of orange = genes almost always found together "
            "(e.g. one plasmid / resistance island)", transform=ax.transAxes,
            fontsize=8, color=INK_2)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"fig5_cooccurrence.{ext}"))
    plt.close(fig)
    pd.DataFrame(r, index=labels, columns=labels).round(3).to_csv(
        os.path.join(outdir, "cooccurrence_phi.tsv"), sep="\t")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--indir", required=True,
                    help="folder with per-genome AMRFinderPlus outputs")
    ap.add_argument("--pattern", default="*_amrfinder.txt",
                    help="glob for output files (default: *_amrfinder.txt)")
    ap.add_argument("--suffix", default="_amrfinder",
                    help="suffix stripped from file names to get the sample "
                         "ID (default: _amrfinder)")
    ap.add_argument("--outdir", default="amrfinder_summary")
    ap.add_argument("--tree", help="Newick core-genome tree (tip labels = "
                                   "sample IDs). Optional.")
    ap.add_argument("--midpoint", action="store_true",
                    help="midpoint-root the tree before plotting")
    ap.add_argument("--metadata", help="TSV, header, first column = sample")
    ap.add_argument("--color-by", nargs="+", default=[],
                    help="metadata column(s) for strips/colours; the first "
                         "is the main grouping for PCoA/prevalence/burden")
    ap.add_argument("--exclude-partial", action="store_true",
                    help="drop PARTIAL*/INTERNAL_STOP hits before building "
                         "matrices (they are always kept in the merged table)")
    ap.add_argument("--exclude-types", nargs="*", default=[],
                    help="feature classes to leave out of matrices/figures, "
                         "e.g. VIRULENCE")
    ap.add_argument("--min-genomes", type=int, default=2,
                    help="min. genomes carrying a feature for it to be shown "
                         "in the heatmap / co-occurrence plot (default 2)")
    ap.add_argument("--drop-core", action="store_true",
                    help="hide features present in >=99%% of genomes in fig 1")
    ap.add_argument("--top-n", type=int, default=40,
                    help="features shown in the prevalence plot (default 40)")
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--burden-ignore", nargs="*", default=["aac(6')-Iaa"],
                    help="AMR element symbols NOT counted in the per-genome "
                         "burden (default: aac(6')-Iaa, intrinsic/cryptic in "
                         "all S. enterica). Pass with no value to count all.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) merge ---------------------------------------------------------------
    log("[1] merging AMRFinderPlus tables")
    merged, samples = merge_tables(args.indir, args.pattern, args.suffix)
    merged.to_csv(os.path.join(args.outdir, "amrfinder_merged.tsv"),
                  sep="\t", index=False)
    log(f"    {len(samples)} genomes, {len(merged)} hits, "
        f"{merged['Element symbol'].nunique()} distinct element symbols")
    zero = sorted(set(samples) - set(merged["sample"]))
    if zero:
        log(f"    {len(zero)} genomes have no hits at all "
            f"(kept as all-zero rows)")

    # 2) matrices ------------------------------------------------------------
    log("[2] building presence/absence matrices")
    work = merged.copy()
    if args.exclude_partial and "Method" in work.columns:
        before = len(work)
        work = work[~work["Method"].fillna("").isin(PARTIAL_METHODS)]
        log(f"    dropped {before - len(work)} partial / internal-stop hits")
    work["feature_class"] = np.where(point_mask(work), "POINT", work["Type"])
    if args.exclude_types:
        work = work[~work["feature_class"].isin(args.exclude_types)]
    mat, info, long_df = build_matrices(work.drop(columns="feature_class"),
                                        samples, args.outdir)
    log(f"    matrix: {mat.shape[0]} genomes x {mat.shape[1]} features")
    for fc in TYPE_ORDER:
        log(f"      {fc:<10} {int((info['feature_class'] == fc).sum())}")

    burden = burden_table(long_df, samples, args.burden_ignore)
    burden.to_csv(os.path.join(args.outdir, "per_genome_burden.tsv"),
                  sep="\t")

    # metadata ---------------------------------------------------------------
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
            sys.exit("No metadata sample IDs match the AMRFinderPlus "
                     "sample names - check the first column / --suffix")
        meta = m
    elif args.color_by:
        sys.exit("--color-by needs --metadata")
    main_grp = args.color_by[0] if args.color_by else None

    # 3) tree + heatmap ------------------------------------------------------
    log("[3] tree + heatmap")
    plot_tree_heatmap(mat, info, args.tree, args.midpoint, meta,
                      args.color_by, args.min_genomes, args.drop_core,
                      args.outdir)

    # 4) other figures -------------------------------------------------------
    log("[4] PCoA / PERMANOVA")
    plot_pcoa(mat, meta, main_grp, args.outdir, args.permutations)
    if main_grp:
        log("[4] prevalence by group + burden")
        plot_prevalence(mat, info, meta, main_grp, args.outdir, args.top_n)
        plot_burden(burden, meta, main_grp, args.outdir)
    log("[4] co-occurrence")
    plot_cooccurrence(mat, info, args.min_genomes, args.outdir)
    log(f"done -> {args.outdir}/")


if __name__ == "__main__":
    main()