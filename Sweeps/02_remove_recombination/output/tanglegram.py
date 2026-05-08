"""
tanglegram.py
=============
Draw a publication-quality tanglegram comparing two phylogenetic trees,
with bezier curves connecting matching leaves coloured by a top-to-bottom
gradient based on the left tree tip position.  This makes it easy to
visually follow where each tip moves between the two trees.
 
Usage
-----
    python3 tanglegram.py raw_tree.treefile clean_tree.treefile tanglegram.png
    Add the value for the Robinson-Foulds distance between the trees in the subtitle (see example).
 
Dependencies
------------
    pip install ete3 matplotlib numpy
"""
 
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.path import Path
import numpy as np
 
# ── Try to import a tree parser ───────────────────────────────────────────────
try:
    from ete3 import Tree
    PARSER = "ete3"
except ImportError:
    try:
        from Bio import Phylo
        PARSER = "biopython"
    except ImportError:
        raise ImportError("Install ete3 or biopython: pip install ete3")
 
 
# =============================================================================
# Constants
# =============================================================================
 
TREE_COLOR     = "#2E4057"    # dark blue-grey for branches
LABEL_FONTSIZE = 4.5          # small enough for 117 taxa
 
# Colormap for the gradient: top = one colour, bottom = another.
# "turbo" gives a vivid, perceptually distinct rainbow that's easy to follow.
# Good alternatives: "rainbow", "hsv", "nipy_spectral"
CMAP = cm.get_cmap("nipy_spectral")
 
 
# =============================================================================
# Tree loading
# =============================================================================
 
def load_tree(path):
    """Load a newick tree. Returns ete3 Tree or biopython Clade."""
    if PARSER == "ete3":
        return Tree(path, format=0)   # format=0: branch lengths + support
    else:
        return list(Phylo.parse(path, "newick"))[0]
 
 
def get_ordered_tips(tree):
    """Ladderize tree and return tip labels in top-to-bottom display order."""
    if PARSER == "ete3":
        tree.ladderize()
        return [leaf.name for leaf in tree.get_leaves()]
    else:
        tree.ladderize()
        return [c.name for c in tree.get_terminals()]
 
 
# =============================================================================
# Coordinate computation
# =============================================================================
 
def compute_coords_ete3(tree, x_origin, x_width, tips_order):
    """
    Assign (x, y) to every node.
    Leaves get y = their index in tips_order.
    Internal nodes get y = mean of children y.
    x is scaled from root-to-tip distance.
 
    Returns
    -------
    coords : dict  node -> (x, y)
    edges  : list of ((px,py),(px,cy),(cx,cy)) elbow segments
    tip_y  : dict  name -> y
    """
    tip_y    = {name: i for i, name in enumerate(tips_order)}
    max_dist = max(tree.get_distance(l) for l in tree.get_leaves()) or 1.0
 
    def sx(d):
        return x_origin + (d / max_dist) * x_width
 
    coords = {}
    for node in tree.traverse("postorder"):
        x = sx(tree.get_distance(node))
        y = tip_y[node.name] if node.is_leaf() else \
            sum(coords[c][1] for c in node.get_children()) / len(node.get_children())
        coords[node] = (x, y)
 
    edges = []
    for node in tree.traverse():
        if not node.is_root():
            px, py = coords[node.up]
            cx, cy = coords[node]
            edges.append(((px, py), (px, cy), (cx, cy)))
 
    return coords, edges, tip_y
 
 
def compute_coords_biopython(tree, x_origin, x_width, tips_order):
    """Equivalent for biopython trees."""
    tip_y    = {name: i for i, name in enumerate(tips_order)}
    max_dist = max(tree.distance(c) for c in tree.get_terminals()) or 1.0
 
    def sx(d):
        return x_origin + (d / max_dist) * x_width
 
    coords = {}
    edges  = []
 
    def assign(clade, parent_dist=0.0):
        dist = parent_dist + (clade.branch_length or 0.0)
        x    = sx(dist)
        if clade.is_terminal():
            y = tip_y[clade.name]
        else:
            for child in clade.clades:
                assign(child, dist)
            y = sum(coords[id(c)][1] for c in clade.clades) / len(clade.clades)
        coords[id(clade)] = (x, y)
 
    assign(tree.root)
 
    def build_edges(clade):
        px, py = coords[id(clade)]
        for child in clade.clades:
            cx, cy = coords[id(child)]
            edges.append(((px, py), (px, cy), (cx, cy)))
            build_edges(child)
 
    build_edges(tree.root)
    return coords, edges, tip_y
 
 
def compute_coords(tree, x_origin, x_width, tips_order):
    if PARSER == "ete3":
        return compute_coords_ete3(tree, x_origin, x_width, tips_order)
    else:
        return compute_coords_biopython(tree, x_origin, x_width, tips_order)
 
 
# =============================================================================
# Drawing helpers
# =============================================================================
 
def draw_tree_edges(ax, edges, color=TREE_COLOR, lw=0.6):
    """Draw rectangular elbow edges (vertical connector + horizontal branch)."""
    for (px, py), (px2, cy), (cx, cy2) in edges:
        ax.plot([px, px2], [py, cy],  color=color, lw=lw, solid_capstyle="round")
        ax.plot([px2, cx], [cy, cy2], color=color, lw=lw, solid_capstyle="round")
 
 
def draw_bezier(ax, y_left, y_right, x_left, x_right, color, alpha=0.75, lw=0.65):
    """
    Draw a cubic bezier curve from (x_left, y_left) to (x_right, y_right).
    Control points at 40/60% of the horizontal span produce a smooth S-curve.
    """
    cx1 = x_left  + (x_right - x_left) * 0.4
    cx2 = x_right - (x_right - x_left) * 0.4
 
    verts = [(x_left,  y_left),
             (cx1,     y_left),
             (cx2,     y_right),
             (x_right, y_right)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
 
    ax.add_patch(mpatches.PathPatch(
        Path(verts, codes),
        facecolor="none",
        edgecolor=color,
        lw=lw,
        alpha=alpha,
        zorder=2
    ))
 
 
# =============================================================================
# Main tanglegram
# =============================================================================
 
def make_tanglegram(tree1_path, tree2_path, out_path):
 
    # ── Load & order tips ─────────────────────────────────────────────────────
    t1 = load_tree(tree1_path)
    t2 = load_tree(tree2_path)
 
    order1 = get_ordered_tips(t1)   # left tree: top → bottom
    order2 = get_ordered_tips(t2)   # right tree: top → bottom
 
    n = len(order1)
 
    # ── Assign one gradient colour per tip (based on left-tree rank) ──────────
    # Tip at top (index n-1) gets the top of the colormap;
    # tip at bottom (index 0) gets the bottom.  We normalise rank to [0, 1].
    tip_color = {}
    for name in order1:
        rank      = order1.index(name)    # 0 = bottom, n-1 = top
        norm_rank = rank / (n - 1)        # normalise to [0.0, 1.0]
        tip_color[name] = CMAP(norm_rank) # RGBA tuple from colormap
 
    # ── Figure layout ─────────────────────────────────────────────────────────
    fig_h = max(20, n * 0.18)
    fig_w = 22
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off()
    ax.set_xlim(-1, fig_w + 1)
    ax.set_ylim(-1, n)
 
    # x layout constants (data units)
    LEFT_LABEL_X     = 0.0
    LEFT_TREE_START  = 0.8
    LEFT_TREE_END    = 7.5
    RIGHT_TREE_START = fig_w - 7.5
    RIGHT_TREE_END   = fig_w - 0.8
    RIGHT_LABEL_X    = fig_w
    GAP_LEFT         = LEFT_TREE_END
    GAP_RIGHT        = RIGHT_TREE_START
 
    # ── Compute node coordinates ───────────────────────────────────────────────
    _, edges1, tip_y1 = compute_coords(
        t1, LEFT_TREE_START, LEFT_TREE_END - LEFT_TREE_START, order1)
    _, edges2, tip_y2 = compute_coords(
        t2, RIGHT_TREE_START, RIGHT_TREE_END - RIGHT_TREE_START, order2)
 
    # Mirror the right tree so its root faces inward
    def mirror_x(x):
        return RIGHT_TREE_START + (RIGHT_TREE_END - x)
 
    edges2_mirrored = [
        ((mirror_x(px), py), (mirror_x(px2), cy), (mirror_x(cx), cy2))
        for (px, py), (px2, cy), (cx, cy2) in edges2
    ]
 
    # ── Draw tree branches ────────────────────────────────────────────────────
    draw_tree_edges(ax, edges1,          color=TREE_COLOR, lw=0.55)
    draw_tree_edges(ax, edges2_mirrored, color=TREE_COLOR, lw=0.55)
 
    # ── Draw gradient bezier curves ───────────────────────────────────────────
    # Iterate in bottom-to-top order so colours higher in the gradient
    # (warmer/brighter) render on top of lower ones when lines cross
    for name in order1:
        if name not in tip_y2:
            continue
        draw_bezier(ax,
                    y_left  = tip_y1[name],
                    y_right = tip_y2[name],
                    x_left  = GAP_LEFT,
                    x_right = GAP_RIGHT,
                    color   = tip_color[name],
                    alpha   = 0.70,
                    lw      = 0.65)
 
    # ── Tip labels — coloured to match their bezier line ─────────────────────
    for name in order1:
        ax.text(LEFT_LABEL_X, tip_y1[name], name,
                ha="right", va="center",
                fontsize=LABEL_FONTSIZE,
                color=tip_color[name],
                fontfamily="monospace")
 
    # Right labels get the same colour their tip has on the left tree,
    # so you can immediately match left label → line → right label by colour
    for name in order2:
        color = tip_color.get(name, "#333333")
        ax.text(RIGHT_LABEL_X, tip_y2[name], name,
                ha="left", va="center",
                fontsize=LABEL_FONTSIZE,
                color=color,
                fontfamily="monospace")
 
    # ── Titles & subtitle ────────────────────────────────────────────────────
    ax.text((LEFT_TREE_START + LEFT_TREE_END) / 2, n - 0.2,
            "Raw alignment tree",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
            color=TREE_COLOR)
    ax.text((RIGHT_TREE_START + RIGHT_TREE_END) / 2, n - 0.2,
            "Recombination-corrected tree",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
            color=TREE_COLOR)
    ax.text(fig_w / 2, n - 0.2,
            "Robinson-Foulds distance = 152  |  Lines & labels coloured by left-tree position (top → bottom)",
            ha="center", va="bottom", fontsize=8, color="#555555", style="italic")
 
    # ── Colorbar as gradient legend ───────────────────────────────────────────
    cbar_ax = fig.add_axes([0.92, 0.1, 0.008, 0.8])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(vmin=0, vmax=n - 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Left-tree position\n(top → bottom)", fontsize=7, labelpad=6)
    cbar.set_ticks([0, n - 1])
    cbar.set_ticklabels(["bottom", "top"], fontsize=6)
 
    # ── Save ─────────────────────────────────────────────────────────────────
    plt.savefig(out_path, dpi=200, bbox_inches="tight",
                format=out_path.split(".")[-1])
    print(f"Saved: {out_path}")
 
 
# =============================================================================
# Entry point
# =============================================================================
 
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 tanglegram.py <tree1> <tree2> <output.png>")
        sys.exit(1)
    make_tanglegram(sys.argv[1], sys.argv[2], sys.argv[3])
 