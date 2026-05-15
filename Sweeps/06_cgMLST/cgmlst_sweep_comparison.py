#!/usr/bin/env python3
"""
cgMLST vs Selective Sweep Comparison Pipeline
----------------------------------------------
1. Computes pairwise allelic distance matrix from cgMLST profiles
2. Builds a Neighbor-Joining tree from the distance matrix
3. Loads sweep group assignments
4. Produces:
   - Distance matrix TSV
   - cgMLST NJ tree (Newick)
   - Colored tree figure (cgMLST tree annotated with sweep groups)
   - Tanglegram comparing cgMLST tree vs sweep tree
   - Adjusted Rand Index between cgMLST clusters and sweep groups

Usage:
    python cgmlst_sweep_comparison.py \
        --profiles cgMLST_profiles.tsv \
        --sweeps sweeps.tsv \
        --sweep_tree sweep_tree.nwk \
        --outdir results/

Dependencies:
    pip install pandas numpy scipy scikit-learn biopython matplotlib
"""

import argparse
import os
import numpy as np
import pandas as pd


# ─────────────────────────────────────────────
# 1. ARGUMENT PARSING
# ─────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles",     required=True,
                   help="cgMLST_profiles.tsv from ExtractCgMLST")
    p.add_argument("--sweeps",       required=True,
                   help="TSV with columns: sweep_id, genome (and optionally others)")
    p.add_argument("--sweep_tree",   required=True,
                   help="Newick tree from selective sweep analysis")
    p.add_argument("--outdir",       default="cgmlst_sweep_results",
                   help="Output directory (created if needed)")
    p.add_argument("--threshold",    type=int, default=None,
                   help="Allelic distance threshold for cgMLST clustering "
                        "(default: auto-detect using 5th percentile of nonzero distances)")
    p.add_argument("--strip_prefix", default="out_",
                   help="Prefix to strip from genome IDs in sweep file (default: 'out_')")
    return p.parse_args()


# ─────────────────────────────────────────────
# 2. LOAD & CLEAN PROFILES
# ─────────────────────────────────────────────

def load_profiles(path):
    print(f"[1/7] Loading cgMLST profiles from {path} ...")
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    # Force sample names to strings
    df.index = [str(s) for s in df.index]
    print(f"      {df.shape[0]} samples x {df.shape[1]} loci")
    df = df.apply(pd.to_numeric, errors="coerce")
    missing_pct = df.isna().mean().mean() * 100
    print(f"      Missing data after masking invalid calls: {missing_pct:.1f}%")
    return df


# ─────────────────────────────────────────────
# 3. PAIRWISE ALLELIC DISTANCE
# ─────────────────────────────────────────────

def compute_distance_matrix(df):
    print("[2/7] Computing pairwise allelic distance matrix ...")
    arr = df.values.astype(float)
    n = arr.shape[0]
    dist = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            a, b = arr[i], arr[j]
            mask = ~(np.isnan(a) | np.isnan(b))
            if mask.sum() == 0:
                dist[i, j] = dist[j, i] = np.nan
            else:
                dist[i, j] = dist[j, i] = np.sum(a[mask] != b[mask])

    samples = [str(s) for s in df.index]
    print(f"      Done. Distance range: {np.nanmin(dist[dist>0]):.0f} - {np.nanmax(dist):.0f} alleles")
    return dist, samples


def save_distance_matrix(dist, samples, outdir):
    path = os.path.join(outdir, "cgmlst_distance_matrix.tsv")
    pd.DataFrame(dist, index=samples, columns=samples).to_csv(path, sep="\t")
    print(f"      Distance matrix saved -> {path}")


# ─────────────────────────────────────────────
# 4. NEIGHBOR-JOINING TREE
# ─────────────────────────────────────────────

def build_nj_tree(dist, samples, outdir):
    print("[3/7] Building Neighbor-Joining tree ...")
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio import Phylo

    # Ensure all sample names are strings (BioPython requires this)
    samples_str = [str(s) for s in samples]

    # Build lower triangle for BioPython
    n = len(samples_str)
    lower = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            val = dist[i][j]
            row.append(0.0 if np.isnan(val) else float(val))
        lower.append(row)

    dm = DistanceMatrix(samples_str, lower)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    nwk_path = os.path.join(outdir, "cgmlst_nj_tree.nwk")
    Phylo.write(tree, nwk_path, "newick")
    print(f"      NJ tree saved -> {nwk_path}")
    return nwk_path


# ─────────────────────────────────────────────
# 5. LOAD SWEEPS
# ─────────────────────────────────────────────

def load_sweeps(path, strip_prefix="out_"):
    print(f"[4/7] Loading sweep assignments from {path} ...")
    df = pd.read_csv(path, sep="\t")
    df["genome_clean"] = df["genome"].str.replace(f"^{strip_prefix}", "", regex=True)
    sweep_map = dict(zip(df["genome_clean"], df["sweep_id"].astype(str)))
    print(f"      {len(sweep_map)} genomes in {df['sweep_id'].nunique()} sweeps")
    return sweep_map


# ─────────────────────────────────────────────
# 6. DISTANCE-BASED CLUSTERING
# ─────────────────────────────────────────────

def cluster_by_distance(dist, samples, threshold, outdir):
    print(f"[5/7] Clustering cgMLST profiles at threshold <= {threshold} alleles ...")
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    from collections import Counter

    dist_filled = np.where(np.isnan(dist), np.nanmax(dist), dist)
    condensed = squareform(dist_filled)
    Z = linkage(condensed, method="single")
    labels = fcluster(Z, t=threshold, criterion="distance")

    counts = Counter(labels)
    cluster_map = {}
    cg_id = 1
    label_to_cg = {}
    for sample, lbl in zip(samples, labels):
        if counts[lbl] > 1:
            if lbl not in label_to_cg:
                label_to_cg[lbl] = f"cgC{cg_id}"
                cg_id += 1
            cluster_map[sample] = label_to_cg[lbl]
        else:
            cluster_map[sample] = "background"

    n_clusters = cg_id - 1
    n_bg = sum(1 for v in cluster_map.values() if v == "background")
    print(f"      {n_clusters} cgMLST clusters found, {n_bg} background singletons")

    out_path = os.path.join(outdir, "cgmlst_clusters.tsv")
    pd.DataFrame([{"sample": s, "cgmlst_cluster": cluster_map[s]} for s in samples]
                 ).to_csv(out_path, sep="\t", index=False)
    print(f"      Cluster assignments saved -> {out_path}")
    return cluster_map


# ─────────────────────────────────────────────
# 7. ADJUSTED RAND INDEX
# ─────────────────────────────────────────────

def compute_ari(samples, cgmlst_clusters, sweep_map):
    print("[6/7] Computing Adjusted Rand Index ...")
    from sklearn.metrics import adjusted_rand_score

    cg_labels = [cgmlst_clusters[s] for s in samples]
    sw_labels = [sweep_map.get(s, "no_sweep") for s in samples]
    ari = adjusted_rand_score(sw_labels, cg_labels)

    if ari > 0.7:
        interp = "strong agreement"
    elif ari > 0.4:
        interp = "moderate agreement"
    elif ari > 0.1:
        interp = "weak agreement"
    else:
        interp = "little to no agreement"

    print(f"      Adjusted Rand Index (cgMLST vs sweeps): {ari:.4f}  [{interp}]")
    return ari


# ─────────────────────────────────────────────
# 8. VISUALISATION
# ─────────────────────────────────────────────

SWEEP_COLORS = {
    "1":  "#E24B4A", "2":  "#3B8BD4", "3":  "#1D9E75",
    "4":  "#EF9F27", "5":  "#D4537E", "6":  "#7F77DD",
    "7":  "#D85A30", "8":  "#639922", "9":  "#BA7517",
    "10": "#0F6E56", "no_sweep": "#CCCCCC", "background": "#CCCCCC",
}

def get_color(label):
    return SWEEP_COLORS.get(str(label), "#999999")


def plot_annotated_tree(nwk_path, sweep_map, cgmlst_clusters, ari, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from Bio import Phylo

    tree = Phylo.read(nwk_path, "newick")
    tree.ladderize()
    terminals = tree.get_terminals()
    n = len(terminals)

    fig, ax = plt.subplots(figsize=(14, max(10, n * 0.12)))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda x: "", branch_labels=lambda x: "")

    for i, clade in enumerate(terminals):
        name = clade.name.strip("'") if clade.name else ""
        sweep = sweep_map.get(name, "no_sweep")
        color = get_color(sweep)
        x = tree.distance(clade)
        ax.plot(x, i + 1, "o", color=color, markersize=4, zorder=5)

    present_sweeps = sorted(set(sweep_map.values()), key=lambda x: int(x))
    handles = [mpatches.Patch(color=get_color(str(s)), label=f"Sweep {s}")
               for s in present_sweeps]
    handles.append(mpatches.Patch(color=get_color("no_sweep"), label="No sweep"))
    ax.legend(handles=handles, loc="lower right", framealpha=0.9, fontsize=9)
    ax.set_title(f"cgMLST NJ Tree — tips colored by sweep group\n"
                 f"ARI = {ari:.4f}", fontsize=12)
    ax.set_xlabel("Allelic distance")
    ax.set_yticks([])

    fig_path = os.path.join(outdir, "cgmlst_tree_sweep_colored.pdf")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    plt.close()
    print(f"      Annotated tree saved -> {fig_path}")


def plot_tanglegram(nwk_cgmlst, nwk_sweep, sweep_map, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo

    tree_cg = Phylo.read(nwk_cgmlst, "newick")
    tree_sw = Phylo.read(nwk_sweep,  "newick")
    tree_cg.ladderize()
    tree_sw.ladderize()

    cg_tips  = [c.name.strip("'") for c in tree_cg.get_terminals() if c.name]
    sw_tips_raw = [c.name.strip("'") for c in tree_sw.get_terminals() if c.name]
    sw_tips  = [t.replace("out_", "") for t in sw_tips_raw]
    common   = set(cg_tips) & set(sw_tips)
    print(f"      Tanglegram: {len(common)} samples present in both trees")

    fig, axes = plt.subplots(1, 2, figsize=(20, max(12, len(cg_tips) * 0.12)))
    Phylo.draw(tree_cg, axes=axes[0], do_show=False,
               label_func=lambda x: "", branch_labels=lambda x: "")
    Phylo.draw(tree_sw, axes=axes[1], do_show=False,
               label_func=lambda x: "", branch_labels=lambda x: "")
    axes[0].set_title("cgMLST NJ tree", fontsize=11)
    axes[1].set_title("Selective sweep tree", fontsize=11)
    axes[0].set_yticks([])
    axes[1].set_yticks([])

    cg_y = {c.name.strip("'"): i + 1 for i, c in enumerate(tree_cg.get_terminals())}
    sw_y = {t.replace("out_", ""): i + 1
            for i, t in enumerate(c.name.strip("'")
                                  for c in tree_sw.get_terminals() if c.name)}

    n_cg = len(cg_tips)
    n_sw = len(sw_tips)
    fig.canvas.draw()  # needed to get axis positions

    for sample in common:
        sweep  = sweep_map.get(sample, "no_sweep")
        color  = get_color(sweep)
        alpha  = 0.7 if sweep != "no_sweep" else 0.08
        lw     = 1.2 if sweep != "no_sweep" else 0.3

        y_cg = cg_y.get(sample)
        y_sw = sw_y.get(sample)
        if y_cg is None or y_sw is None:
            continue

        ax0_bbox = axes[0].get_position()
        ax1_bbox = axes[1].get_position()
        y0_fig = ax0_bbox.y0 + (y_cg / (n_cg + 1)) * ax0_bbox.height
        y1_fig = ax1_bbox.y0 + (y_sw / (n_sw + 1)) * ax1_bbox.height

        fig.add_artist(plt.Line2D(
            [ax0_bbox.x1, ax1_bbox.x0], [y0_fig, y1_fig],
            transform=fig.transFigure,
            color=color, alpha=alpha, linewidth=lw, zorder=0))

    fig.suptitle("Tanglegram: cgMLST NJ tree vs selective sweep tree\n"
                 "Colored lines = sweep members  |  grey = background", fontsize=12)
    tang_path = os.path.join(outdir, "tanglegram_cgmlst_vs_sweeps.pdf")
    fig.tight_layout()
    fig.savefig(tang_path, dpi=150)
    plt.close()
    print(f"      Tanglegram saved -> {tang_path}")


def plot_confusion(samples, cgmlst_clusters, sweep_map, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.DataFrame({
        "sweep":  [sweep_map.get(s, "no_sweep") for s in samples],
        "cgmlst": [cgmlst_clusters[s] for s in samples],
    })
    ct = pd.crosstab(df["sweep"], df["cgmlst"])

    fig, ax = plt.subplots(figsize=(max(8, len(ct.columns) * 0.6),
                                    max(4, len(ct.index) * 0.5)))
    im = ax.imshow(ct.values, aspect="auto", cmap="Blues")
    ax.set_xticks(range(len(ct.columns)))
    ax.set_xticklabels(ct.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(ct.index)))
    ax.set_yticklabels(ct.index, fontsize=9)
    ax.set_xlabel("cgMLST cluster")
    ax.set_ylabel("Sweep group")
    ax.set_title("Co-occurrence: sweep groups x cgMLST clusters")
    plt.colorbar(im, ax=ax, label="Number of samples")

    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            val = ct.values[i, j]
            if val > 0:
                ax.text(j, i, str(val), ha="center", va="center", fontsize=7,
                        color="white" if val > ct.values.max() * 0.6 else "black")

    conf_path = os.path.join(outdir, "sweep_vs_cgmlst_cooccurrence.pdf")
    fig.tight_layout()
    fig.savefig(conf_path, dpi=150)
    plt.close()
    print(f"      Co-occurrence heatmap saved -> {conf_path}")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    profiles = load_profiles(args.profiles)
    dist, samples = compute_distance_matrix(profiles)
    save_distance_matrix(dist, samples, args.outdir)

    nwk_cgmlst = build_nj_tree(dist, samples, args.outdir)

    sweep_map = load_sweeps(args.sweeps, strip_prefix=args.strip_prefix)

    threshold = args.threshold
    if threshold is None:
        valid = dist[(dist > 0) & ~np.isnan(dist)]
        threshold = max(1, int(np.percentile(valid, 5)))
        print(f"      Auto-threshold: {threshold} alleles (5th percentile of nonzero distances)")
        print(f"      Re-run with --threshold N to adjust")

    cgmlst_clusters = cluster_by_distance(dist, samples, threshold, args.outdir)
    ari = compute_ari(samples, cgmlst_clusters, sweep_map)

    summary_path = os.path.join(args.outdir, "summary.txt")
    n_cg = len(set(v for v in cgmlst_clusters.values() if v != "background"))
    with open(summary_path, "w") as f:
        f.write(f"Samples:               {len(samples)}\n")
        f.write(f"cgMLST loci:           {profiles.shape[1]}\n")
        f.write(f"Clustering threshold:  {threshold} alleles\n")
        f.write(f"cgMLST clusters:       {n_cg}\n")
        f.write(f"Sweep groups:          {len(set(sweep_map.values()))}\n")
        f.write(f"Adjusted Rand Index:   {ari:.4f}\n")
    print(f"      Summary saved -> {summary_path}")

    plot_annotated_tree(nwk_cgmlst, sweep_map, cgmlst_clusters, ari, args.outdir)
    plot_tanglegram(nwk_cgmlst, args.sweep_tree, sweep_map, args.outdir)
    plot_confusion(samples, cgmlst_clusters, sweep_map, args.outdir)

    print("\nDone. Output files:")
    for fname in sorted(os.listdir(args.outdir)):
        print(f"  {args.outdir}/{fname}")


if __name__ == "__main__":
    main()