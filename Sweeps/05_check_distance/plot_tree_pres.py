#!/usr/bin/env python3
"""
plot_tree.py
Draws the phylogenetic tree and colors tips/branches by sweep membership.
Reads sweeps.txt output from detect_sweeps.py.

Usage:
    python3 plot_tree.py --tree full_tree.treefile --sweeps sweeps.txt --out tree.png [--labels]
"""

import argparse
import itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgba
import numpy as np
from Bio import Phylo
import io

parser = argparse.ArgumentParser()
parser.add_argument("--tree",   required=True)
parser.add_argument("--sweeps", required=True)
parser.add_argument("--out",    default="tree.png")
parser.add_argument("--dpi",    type=int, default=200)
parser.add_argument("--labels", action="store_true",
                    help="Show tip labels (only for sweep members and sisters by default; use --labels-all for every tip)")
parser.add_argument("--labels-all", action="store_true",
                    help="Show labels for all tips (very dense with 400+ tips)")
args = parser.parse_args()

# --- Load sweep assignments ---
genome_to_sweep = {}   # genome -> sweep_id (int)
sweep_meta = {}        # sweep_id -> {sister, within_mean, sister_mean, ratio}
with open(args.sweeps) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 7:
            continue
        sweep_id, genome, sister, within_mean, sister_mean, ratio, clade_size = parts
        sweep_id = int(sweep_id)
        genome_to_sweep[genome] = sweep_id
        if sweep_id not in sweep_meta:
            sweep_meta[sweep_id] = {
                "sister": sister,
                "within_mean": float(within_mean),
                "sister_mean": float(sister_mean),
                "ratio": float(ratio),
                "size": int(clade_size),
            }

n_sweeps = max(sweep_meta.keys()) if sweep_meta else 0

# collect all sister tips
sister_tips = set(meta["sister"] for meta in sweep_meta.values())

# --- Color palette (distinct, colorblind-friendly-ish) ---
SWEEP_COLORS = [
    "#E63946",  # red
    "#F4A261",  # orange
    "#2ECC71",  # green
    "#F1C40F",  # yellow
    "#0000FF",  # dark blue
    "#6A0572",  # purple
    "#457B9D",  # steel blue
    "#E67E22",  # dark orange
    "#1ABC9C",  # mint
    "#2A9D8F",  # teal
]
DEFAULT_COLOR = "#000000"   # was #CCCCCC (light gray, for dark bg)
SISTER_COLOR  = "#000000"   # was #FFFFFF (white, for dark bg)
SWEEP_COLOR_MAP = {i+1: SWEEP_COLORS[i % len(SWEEP_COLORS)] for i in range(n_sweeps)}

# --- Load and root tree ---
print("Loading tree...")
tree = Phylo.read(args.tree, "newick")
tree.root_at_midpoint()

# --- Compute layout (manual cladogram) ---
# We draw a rectangular cladogram: x = depth, y = tip order
terminals = tree.get_terminals()
n_tips = len(terminals)

# Assign y positions to tips in tree order
tip_y = {}
counter = [0]
def assign_y(clade):
    if clade.is_terminal():
        tip_y[clade.name] = counter[0]
        counter[0] += 1
    else:
        for child in clade.clades:
            assign_y(child)
assign_y(tree.root)

# Compute x positions (cumulative branch length from root)
node_x = {}
node_y = {}

def assign_x(clade, depth=0.0):
    bl = clade.branch_length if clade.branch_length else 0.0
    x = depth + bl
    node_x[id(clade)] = x
    if clade.is_terminal():
        node_y[id(clade)] = tip_y[clade.name]
    else:
        for child in clade.clades:
            assign_x(child, x)
        ys = [node_y[id(c)] for c in clade.clades]
        node_y[id(clade)] = sum(ys) / len(ys)

assign_x(tree.root)

max_x = max(node_x.values())

# Clip display to 99th percentile of tip x positions so outlier long branches
# do not compress the rest of the tree. Long branches still draw but run off-edge.
tip_xs = [node_x[id(c)] for c in tree.get_terminals()]
display_max_x = np.percentile(tip_xs, 99) * 1.15

# --- Draw tree ---
print("Drawing tree...")
fig_height = max(8, n_tips * 0.045)
fig_width  = 28 if (args.labels or args.labels_all) else 22
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
ax.set_facecolor("none")

def get_tip_color(name):
    if name in sister_tips:
        return SISTER_COLOR
    sid = genome_to_sweep.get(name)
    if sid:
        return SWEEP_COLOR_MAP[sid]
    return DEFAULT_COLOR

def get_clade_color(clade):
    """Color a branch by sweep if ALL its tips belong to the same sweep."""
    tips = [c.name for c in clade.get_terminals()]
    sweeps = set(genome_to_sweep.get(t) for t in tips)
    if len(sweeps) == 1 and None not in sweeps:
        sid = sweeps.pop()
        return SWEEP_COLOR_MAP.get(sid, DEFAULT_COLOR)
    return DEFAULT_COLOR

def draw_clade(clade, parent_x=None):
    x = node_x[id(clade)]
    y = node_y[id(clade)]
    color = get_clade_color(clade)
    lw = 0.8 if color == DEFAULT_COLOR else 1.6
    alpha = 0.85 if color == DEFAULT_COLOR else 0.95

    # horizontal line from parent to this node
    if parent_x is not None:
        ax.plot([parent_x, x], [y, y], color=color, lw=lw, alpha=alpha, solid_capstyle="round")

    if not clade.is_terminal():
        child_ys = [node_y[id(c)] for c in clade.clades]
        # vertical connector
        ax.plot([x, x], [min(child_ys), max(child_ys)],
                color=color, lw=lw, alpha=alpha, solid_capstyle="round")
        for child in clade.clades:
            draw_clade(child, x)
    else:
        # tip dot
        tip_color = get_tip_color(clade.name)
        if clade.name in sister_tips:
            ax.scatter(x, y, s=18, color=SISTER_COLOR, zorder=6,
                       marker="D", linewidths=0)
        elif tip_color != DEFAULT_COLOR:
            ax.scatter(x, y, s=6, color=tip_color, zorder=5, linewidths=0)

        # tip label
        show_label = (
            args.labels_all or
            (args.labels and (clade.name in genome_to_sweep or clade.name in sister_tips))
        )
        if show_label:
            label_color = get_tip_color(clade.name)
            fontsize = 4.5 if args.labels_all else 6.0
            ax.text(x + display_max_x * 0.005, y, clade.name,
                    va="center", ha="left",
                    fontsize=fontsize,
                    fontfamily="monospace",
                    color=label_color,
                    alpha=0.9,
                    zorder=7)

for child in tree.root.clades:
    draw_clade(child, node_x[id(tree.root)])

# --- Legend ---
legend_patches = []
for sid in sorted(sweep_meta.keys()):
    meta = sweep_meta[sid]
    label = (f"Sweep {sid}  "
             f"n={meta['size']}  "
             f"within={meta['within_mean']:.1f}  "
             f"sister={meta['sister_mean']:.1f}  "
             f"ratio={meta['ratio']:.1f}x")
    legend_patches.append(
        mpatches.Patch(color=SWEEP_COLOR_MAP[sid], label=label)
    )
legend_patches.append(mpatches.Patch(color=DEFAULT_COLOR, alpha=0.5, label="No sweep"))
legend_patches.append(
    mpatches.Patch(color=SISTER_COLOR, label=f"Sister tip ({len(sister_tips)} genomes)")
)

leg = ax.legend(
    handles=legend_patches,
    loc="lower right",
    fontsize=6.5,
    framealpha=0.7,
    facecolor="white",
    edgecolor="#888888",
    labelcolor="black",
    handlelength=1.2,
    borderpad=0.8,
    labelspacing=0.5,
)

# --- Axes styling ---
ax.set_xlim(-display_max_x * 0.01, display_max_x)
ax.set_ylim(-2, n_tips + 1)
ax.axis("off")
ax.set_title("Phylogenetic tree — selective sweeps highlighted",
             color="#222222", fontsize=11, fontfamily="monospace", pad=12)

# --- Scale bar (substitutions per site) ---
# Pick the largest round number (1, 2, or 5 * 10^n) that fits ~10-20% of display width
import math
raw = display_max_x * 0.12
magnitude = 10 ** math.floor(math.log10(raw))
for nice in [magnitude, 2 * magnitude, 5 * magnitude]:
    if nice <= display_max_x * 0.20:
        bar_len = nice
bar_x = display_max_x * 0.02
bar_y = -1
# line with end ticks
ax.plot([bar_x, bar_x + bar_len], [bar_y, bar_y], color="#444444", lw=1.5, solid_capstyle="butt")
ax.plot([bar_x, bar_x], [bar_y - 0.3, bar_y + 0.3], color="#444444", lw=1.5)
ax.plot([bar_x + bar_len, bar_x + bar_len], [bar_y - 0.3, bar_y + 0.3], color="#444444", lw=1.5)
# label with correct significant figures
ax.text(bar_x + bar_len / 2, bar_y - 1.0,
        f"{bar_len:.2e} substitutions/site",
        ha="center", va="top", color="#444444",
        fontsize=6.5, fontfamily="monospace")

plt.tight_layout(pad=0.5)
plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight",
            #transparent=True
            )
print(f"Saved: {args.out}")