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

# --- Hard-coded sweep 7 ---
# Not present in the sweeps file; defined manually here. Names keep the "out_"
# prefix so they match the tree tip names directly (same as the file-based sweeps).
SWEEP7_ID = 7
SWEEP7_GENOMES = [
    "out_11452375", "out_11452785", "out_11453035", "out_8046645", "out_8091485", "out_11452255",
    "out_11452675", "out_45176625", "out_9318515", "out_31904705", "out_29413815", "out_27606475",
    "out_11168375", "out_11595885", "out_11588445", "out_32316375", "out_25671015", "out_28916195",
    "out_28916085", "out_31765335", "out_24242935", "out_24000255", "out_32498505", "out_32498225",
    "out_31895175", "out_14051005", "out_24910175", "out_8713855", "out_8813075", "out_22363375",
    "out_11740395", "out_32773835", "out_17181415", "out_32858445", "out_16626965", "out_41292945",
    "out_7823845", "out_8778775", "out_14599635", "out_10805445", "out_7794975", "out_17337285",
    "out_10134355", "out_9164925", "out_19205865", "out_27171365", "out_9536415", "out_18277005",
    "out_18276945", "out_44464205", "out_23514945", "out_32498565", "out_14768535", "out_14075525",
    "out_33637285", "out_31725495", "out_23416215", "out_44468645", "out_31765255", "out_10895095",
    "out_33877195", "out_7792295", "out_21308675", "out_9207225", "out_27606455", "out_16459765",
    "out_8559715", "out_7824425", "out_17689515", "out_15165115", "out_26044135", "out_14600175",
    "out_20607955", "out_19562115", "out_8765535", "out_21761305", "out_10141235", "out_9167145",
    "out_20516485", "out_11519185", "out_29769425", "out_31725395", "out_20081485", "out_16626985",
    "out_7883985", "out_33632345", "out_23837735", "out_22905815", "out_32643285", "out_21180765",
    "out_7792575", "out_33015095", "out_8222005", "out_15376085", "out_8461385", "out_10142775",
    "out_10142075", "out_9046505", "out_31904685", "out_9474185", "out_26252495", "out_11518205",
    "out_9342185", "out_16058255", "out_9523365", "out_9026405", "out_31059355", "out_14550555",
    "out_20964875", "out_11156655", "out_8666235", "out_9579385", "out_16976055", "out_9124805",
    "out_33688975", "out_33015535", "out_31895395", "out_32858425", "out_20653895", "out_18276585",
    "out_11528025", "out_8715535", "out_22617875", "out_9401885", "out_11169795", "out_8461765",
    "out_20808085", "out_46616815", "out_32450575", "out_19357875", "out_32323275", "out_10902475",
    "out_25806255", "out_19267755", "out_14599955", "out_10833965", "out_9206385", "out_8649995",
    "out_8769955", "out_10142695", "out_31895675", "out_15165015", "out_12310975", "out_40741165",
    "out_8989965", "out_20487425", "out_7918635", "out_29323955", "out_19233445", "out_16058695",
    "out_9167845", "out_26370835", "out_11582185", "out_21007735", "out_14581275", "out_8444725",
    "out_9476285", "out_14053045", "out_9123725", "out_9365255", "out_24851145", "out_8770075",
    "out_24001045", "out_8714695", "out_9131505", "out_33015275", "out_31059995", "out_32498545",
    "out_31770105", "out_31184075", "out_31027665", "out_26044315", "out_31722655", "out_24688825",
    "out_17659675", "out_31059475", "out_11600345", "out_22330065", "out_11592475", "out_11551395",
    "out_33757895", "out_32759455", "out_11464235", "out_16349885", "out_22331385", "out_15176205",
    "out_8547575", "out_7795135", "out_32858385", "out_33687925", "out_25597515", "out_15165415",
    "out_9217935", "out_10892735", "out_14045005", "out_10133815", "out_14581255", "out_33014935",
    "out_33631805", "out_31895515", "out_25671005", "out_28547265", "out_20658025", "out_19562235",
    "out_15169455", "out_7736995", "out_9317375", "out_8461165", "out_20864745", "out_16121955",
    "out_14738275", "out_11163695", "out_10133415",
]

# Assign sweep 7 ONLY to genomes not already in another sweep, so the tighter
# sweep 3 (which is a subset of sweep 7) keeps its own color.
# To make sweep 7 take priority instead, change the condition to overwrite
# unconditionally: genome_to_sweep[g] = SWEEP7_ID
_assigned_to_7 = 0
for g in SWEEP7_GENOMES:
    if g not in genome_to_sweep:
        genome_to_sweep[g] = SWEEP7_ID
        _assigned_to_7 += 1
sweep_meta[SWEEP7_ID] = {
    "sister": None,            # no sister/stats for the hard-coded sweep
    "within_mean": None,
    "sister_mean": None,
    "ratio": None,
    "size": _assigned_to_7,    # genomes actually colored as sweep 7
}

n_sweeps = max(sweep_meta.keys()) if sweep_meta else 0

# collect all sister tips (ignore None from the hard-coded sweep)
sister_tips = set(meta["sister"] for meta in sweep_meta.values() if meta["sister"])

# --- Color palette (distinct, colorblind-friendly-ish) ---
SWEEP_COLORS = [
    "#E63946",  # red
    "#F4A261",  # orange
    "#2A9D8F",  # teal
    "#457B9D",  # steel blue
    "#A8DADC",  # light blue
    "#6A0572",  # purple
    "#F1C40F",  # yellow
    "#2ECC71",  # green
    "#E67E22",  # dark orange
    "#1ABC9C",  # mint
]
DEFAULT_COLOR = "#000000"   # was #CCCCCC (light gray, for dark bg)
SISTER_COLOR  = "#000000"   # was #FFFFFF (white, for dark bg)
SWEEP_COLOR_MAP = {i+1: SWEEP_COLORS[i % len(SWEEP_COLORS)] for i in range(n_sweeps)}
# Sweep 7 would land on the pale yellow slot, which is weak on a light background.
# Override it with a strong, distinct color (not used by sweeps 1-6).
SWEEP_COLOR_MAP[7] = "#1B7837"  # forest green

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
    if meta["within_mean"] is None:
        # hard-coded sweep with no computed statistics
        label = f"Sweep {sid}  n={meta['size']}"
    else:
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