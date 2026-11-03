#!/usr/bin/env python3
"""
plot_tree.py
Draws the phylogenetic tree and colors tips/branches by sweep membership.
Reads sweeps.txt output from detect_sweeps.py, plus any number of extra
highlight groups given as plain-text genome lists.

Usage:
    python3 plot_tree.py --tree full_tree.treefile --sweeps sweeps.txt --out tree.png [--labels]

    python3 plot_tree.py --tree full_tree.treefile --sweeps sweeps.txt \
        --highlight sweep7.txt "#1B7837" "Sweep 7" \
        --highlight Infantis_AGES.txt "#00BCD4" "Infantis AGES" \
        --out tree.png

--highlight FILE [COLOR] [LABEL]   (repeat as often as you like)
    FILE   txt file, one genome ID per line (a trailing ".fna" is stripped;
           blank lines and lines starting with "#" are ignored; only the first
           whitespace-separated column is used)
    COLOR  optional, any matplotlib color ("#E63946", "red", ...).
           Default: next unused color from the palette.
    LABEL  optional legend label. Default: file name without extension.

Priority: sweeps from --sweeps are assigned first, then --highlight groups in
the order given on the command line. By default a genome keeps the first group
it was assigned to, so list nested (smaller) clades BEFORE the larger clades
that contain them. Use --highlight-overwrite to make later groups win instead.
"""

import argparse
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import is_color_like
import numpy as np
from Bio import Phylo

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--tree",   required=True)
parser.add_argument("--sweeps", default=None,
                    help="sweeps.txt from detect_sweeps.py (optional)")
parser.add_argument("--out",    default="tree.png")
parser.add_argument("--dpi",    type=int, default=200)
parser.add_argument("--labels", action="store_true",
                    help="Show tip labels (only for highlighted tips and sisters; "
                         "use --labels-all for every tip)")
parser.add_argument("--labels-all", action="store_true",
                    help="Show labels for all tips (very dense with 400+ tips)")
parser.add_argument("--highlight", action="append", nargs="+", default=[],
                    metavar=("FILE", "COLOR LABEL"),
                    help="Genome list to highlight: FILE [COLOR] [LABEL]. Repeatable.")
parser.add_argument("--highlight-overwrite", action="store_true",
                    help="Later groups override earlier assignments (incl. --sweeps). "
                         "Default: a genome keeps the first group it was assigned to.")
args = parser.parse_args()

# --- Color palette (distinct, colorblind-friendly-ish) ---
PALETTE = [
    "#E63946",  # red
    "#F4A261",  # orange
    "#2A9D8F",  # teal
    "#457B9D",  # steel blue
    "#6A0572",  # purple
    "#1B7837",  # forest green
    "#00BCD4",  # cyan
    "#E67E22",  # dark orange
    "#8E44AD",  # violet
    "#D81B60",  # magenta
    "#827717",  # olive
    "#5D4037",  # brown
]
DEFAULT_COLOR = "#000000"
SISTER_COLOR  = "#000000"

genome_to_group = {}   # genome -> group_id (int)
groups = {}            # group_id -> metadata dict (insertion order = legend order)
_palette_i = [0]

def next_palette_color():
    used = {g["color"] for g in groups.values()}
    for _ in range(len(PALETTE)):
        c = PALETTE[_palette_i[0] % len(PALETTE)]
        _palette_i[0] += 1
        if c not in used:
            return c
    c = PALETTE[_palette_i[0] % len(PALETTE)]
    _palette_i[0] += 1
    return c

# --- Load sweep assignments (optional) ---
if args.sweeps:
    with open(args.sweeps) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            sweep_id, genome, sister, within_mean, sister_mean, ratio, clade_size = parts[:7]
            sweep_id = int(sweep_id)
            genome_to_group[genome] = sweep_id
            if sweep_id not in groups:
                groups[sweep_id] = {
                    "label": f"Sweep {sweep_id}",
                    "sister": sister,
                    "within_mean": float(within_mean),
                    "sister_mean": float(sister_mean),
                    "ratio": float(ratio),
                    "size": int(clade_size),
                    "color": None,
                }
    # sweeps get palette colors in sweep-id order
    for sid in sorted(groups):
        groups[sid]["color"] = next_palette_color()

# --- Load highlight groups from txt files ---
def read_genome_list(path):
    ids = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            g = line.split()[0]
            if g.endswith(".fna"):
                g = g[:-len(".fna")]
            ids.append(g)
    return ids

highlight_ids = []
for spec in args.highlight:
    if len(spec) > 3:
        parser.error(f"--highlight takes FILE [COLOR] [LABEL], got {len(spec)} values: {spec} "
                     f"(quote labels that contain spaces)")
    path = spec[0]
    color = spec[1] if len(spec) >= 2 else None
    label = spec[2] if len(spec) == 3 else os.path.splitext(os.path.basename(path))[0]
    if color is not None and not is_color_like(color):
        parser.error(f"--highlight {path}: '{color}' is not a valid color "
                     f"(remember to quote hex colors in the shell: \"#1B7837\")")
    if color is None:
        color = next_palette_color()

    gid = (max(groups) if groups else 0) + 1
    genomes = read_genome_list(path)
    assigned = 0
    for g in genomes:
        if args.highlight_overwrite or g not in genome_to_group:
            genome_to_group[g] = gid
            assigned += 1
    groups[gid] = {
        "label": label,
        "sister": None,
        "within_mean": None,
        "size": assigned,
        "n_in_list": len(genomes),
        "genomes": set(genomes),
        "color": color,
    }
    highlight_ids.append(gid)
    print(f"Highlight '{label}' ({path}): {len(genomes)} genomes read, "
          f"{assigned} assigned, {len(genomes) - assigned} already in an earlier group.")

# With --highlight-overwrite, earlier groups may have lost members: recount.
if args.highlight_overwrite:
    for gid, meta in groups.items():
        meta["size"] = sum(1 for v in genome_to_group.values() if v == gid)

GROUP_COLOR = {gid: meta["color"] for gid, meta in groups.items()}

# collect all sister tips
sister_tips = set(meta["sister"] for meta in groups.values() if meta.get("sister"))

# --- Load and root tree ---
print("Loading tree...")
tree = Phylo.read(args.tree, "newick")
tree.root_at_midpoint()

terminals = tree.get_terminals()
n_tips = len(terminals)
tip_names = {c.name for c in terminals}

# Sanity check: how many listed genomes are actually tips in the tree
for gid in highlight_ids:
    meta = groups[gid]
    found = len(meta["genomes"] & tip_names)
    msg = f"  '{meta['label']}': {found}/{meta['n_in_list']} listed genomes found in tree"
    if found < meta["n_in_list"]:
        missing = sorted(meta["genomes"] - tip_names)[:5]
        msg += f"  (e.g. missing: {', '.join(missing)})"
    print(msg)

# --- Compute layout (manual cladogram) ---
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

# Clip display to 99th percentile of tip x positions so outlier long branches
# do not compress the rest of the tree.
tip_xs = [node_x[id(c)] for c in terminals]
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
    gid = genome_to_group.get(name)
    if gid is not None:
        return GROUP_COLOR[gid]
    return DEFAULT_COLOR

def get_clade_color(clade):
    """Color a branch by group if ALL its tips belong to the same group."""
    ids = set(genome_to_group.get(c.name) for c in clade.get_terminals())
    if len(ids) == 1 and None not in ids:
        return GROUP_COLOR.get(ids.pop(), DEFAULT_COLOR)
    return DEFAULT_COLOR

def draw_clade(clade, parent_x=None):
    x = node_x[id(clade)]
    y = node_y[id(clade)]
    color = get_clade_color(clade)
    lw = 0.8 if color == DEFAULT_COLOR else 1.6
    alpha = 0.85 if color == DEFAULT_COLOR else 0.95

    if parent_x is not None:
        ax.plot([parent_x, x], [y, y], color=color, lw=lw, alpha=alpha, solid_capstyle="round")

    if not clade.is_terminal():
        child_ys = [node_y[id(c)] for c in clade.clades]
        ax.plot([x, x], [min(child_ys), max(child_ys)],
                color=color, lw=lw, alpha=alpha, solid_capstyle="round")
        for child in clade.clades:
            draw_clade(child, x)
    else:
        tip_color = get_tip_color(clade.name)
        if clade.name in sister_tips:
            ax.scatter(x, y, s=18, color=SISTER_COLOR, zorder=6, marker="D", linewidths=0)
        elif tip_color != DEFAULT_COLOR:
            ax.scatter(x, y, s=6, color=tip_color, zorder=5, linewidths=0)

        show_label = (
            args.labels_all or
            (args.labels and (clade.name in genome_to_group or clade.name in sister_tips))
        )
        if show_label:
            fontsize = 4.5 if args.labels_all else 6.0
            ax.text(x + display_max_x * 0.005, y, clade.name,
                    va="center", ha="left", fontsize=fontsize,
                    fontfamily="monospace", color=tip_color, alpha=0.9, zorder=7)

for child in tree.root.clades:
    draw_clade(child, node_x[id(tree.root)])

# --- Legend (sweeps first, then highlight groups in command-line order) ---
legend_patches = []
for gid, meta in groups.items():
    if meta["within_mean"] is None:
        label = f"{meta['label']}  n={meta['size']}"
    else:
        label = (f"{meta['label']}  n={meta['size']}  "
                 f"within={meta['within_mean']:.1f}  "
                 f"sister={meta['sister_mean']:.1f}  "
                 f"ratio={meta['ratio']:.1f}x")
    legend_patches.append(mpatches.Patch(color=meta["color"], label=label))
legend_patches.append(mpatches.Patch(color=DEFAULT_COLOR, alpha=0.5, label="Not highlighted"))
if sister_tips:
    legend_patches.append(
        mpatches.Patch(color=SISTER_COLOR, label=f"Sister tip ({len(sister_tips)} genomes)"))

ax.legend(
    handles=legend_patches, loc="lower right", fontsize=6.5, framealpha=0.7,
    facecolor="white", edgecolor="#888888", labelcolor="black",
    handlelength=1.2, borderpad=0.8, labelspacing=0.5,
)

# --- Axes styling ---
ax.set_xlim(-display_max_x * 0.01, display_max_x)
ax.set_ylim(-2, n_tips + 1)
ax.axis("off")
ax.set_title("Phylogenetic tree — highlighted clades",
             color="#222222", fontsize=11, fontfamily="monospace", pad=12)

# --- Scale bar (substitutions per site) ---
raw = display_max_x * 0.12
magnitude = 10 ** math.floor(math.log10(raw))
bar_len = magnitude
for nice in [magnitude, 2 * magnitude, 5 * magnitude]:
    if nice <= display_max_x * 0.20:
        bar_len = nice
bar_x = display_max_x * 0.02
bar_y = -1
ax.plot([bar_x, bar_x + bar_len], [bar_y, bar_y], color="#444444", lw=1.5, solid_capstyle="butt")
ax.plot([bar_x, bar_x], [bar_y - 0.3, bar_y + 0.3], color="#444444", lw=1.5)
ax.plot([bar_x + bar_len, bar_x + bar_len], [bar_y - 0.3, bar_y + 0.3], color="#444444", lw=1.5)
ax.text(bar_x + bar_len / 2, bar_y - 1.0, f"{bar_len:.2e} substitutions/site",
        ha="center", va="top", color="#444444", fontsize=6.5, fontfamily="monospace")

plt.tight_layout(pad=0.5)
plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
print(f"Saved: {args.out}")