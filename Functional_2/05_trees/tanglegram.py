#!/usr/bin/env python3
"""
tanglegram.py
Draws 2 or 3 phylogenetic trees side by side and connects the same sample in
each tree with a line, so you can see where a genome sits in each topology.
Trees may contain different tip sets (one tree being a subset of the others is
fine) - only shared tips get connectors.

Colouring works exactly like plot_tree.py: sweeps.txt from detect_sweeps.py
plus any number of --highlight genome lists.

Usage:
    python3 tanglegram.py \
        --left  core_tree.treefile   "Core SNP" \
        --mid   accessory.treefile   "Accessory" \
        --right pesi_tree.treefile   "pESI plasmid" \
        --sweeps sweeps.txt \
        --highlight sweep7.txt "#1B7837" "Sweep 7" \
        --highlight Infantis_AGES.txt "#00BCD4" "Infantis AGES" \
        --out tanglegram.png

    # two trees only, with labels
    python3 tanglegram.py --left a.treefile --right b.treefile --labels --out t.png

--left/--mid/--right  FILE [LABEL]
    FILE   newick tree. LABEL is the panel title (default: file name).
    --mid is optional; with two trees you get a classic tanglegram.

--highlight FILE [COLOR] [LABEL]   (repeat as often as you like)
    Same semantics as plot_tree.py: one genome ID per line, trailing ".fna"
    stripped, "#" comments and blank lines ignored, first column only.

Layout
    The left tree is drawn root-left / tips-right, the right tree mirrored
    (root-right / tips-left), so their tip columns face each other across the
    gutter. The middle tree can only face one way: by default its tips face
    right (towards the right tree) and the left gutter's connectors land on an
    anchor column of dots at the middle panel's left edge, one dot per row.
    Use --mid-face left to flip which gutter is the "true" tip-to-tip one, and
    --mid-leaders to draw faint guide lines across the middle panel.

Untangling
    Tip order is fixed for one tree (--order-ref, default left) and the others
    are rotated at their internal nodes to minimise crossings (barycentre
    heuristic, repeated --untangle-passes times). Crossing counts before and
    after are printed. --no-untangle keeps each tree's own ladderized order.
"""

import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import is_color_like
import numpy as np
from Bio import Phylo

# --------------------------------------------------------------------------
# Arguments
# --------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--left",  required=True, nargs="+", metavar=("FILE", "LABEL"))
parser.add_argument("--mid",   default=None,  nargs="+", metavar=("FILE", "LABEL"))
parser.add_argument("--right", required=True, nargs="+", metavar=("FILE", "LABEL"))
parser.add_argument("--sweeps", default=None,
                    help="sweeps.txt from detect_sweeps.py (optional)")
parser.add_argument("--highlight", action="append", nargs="+", default=[],
                    metavar=("FILE", "COLOR LABEL"),
                    help="Genome list to highlight: FILE [COLOR] [LABEL]. Repeatable.")
parser.add_argument("--highlight-overwrite", action="store_true",
                    help="Later groups override earlier assignments (incl. --sweeps).")
parser.add_argument("--out", default="tanglegram.png")
parser.add_argument("--dpi", type=int, default=200)
parser.add_argument("--labels", action="store_true",
                    help="Show tip labels for highlighted tips only")
parser.add_argument("--labels-all", action="store_true",
                    help="Show labels for every tip (dense above ~150 tips)")
parser.add_argument("--label-size", type=float, default=0.0,
                    help="Label font size (default: auto from tip count)")
parser.add_argument("--order-ref", choices=["left", "mid", "right"], default="left",
                    help="Tree whose tip order is held fixed (default: left)")
parser.add_argument("--no-untangle", action="store_true")
parser.add_argument("--untangle-passes", type=int, default=12)
parser.add_argument("--mid-face", choices=["left", "right"], default="right",
                    help="Which gutter the middle tree's tips point into (default: right)")
parser.add_argument("--mid-leaders", action="store_true",
                    help="Draw faint guide lines across the middle panel to its anchor column")
parser.add_argument("--gutter", type=float, default=0.13,
                    help="Gutter width as a fraction of total figure width (default 0.13)")
parser.add_argument("--clip-pct", type=float, default=99.0,
                    help="Percentile of tip depth used as panel width (default 99)")
parser.add_argument("--connect", choices=["all", "highlighted"], default="all",
                    help="Draw connectors for all shared tips or only highlighted ones")
parser.add_argument("--connector-alpha", type=float, default=0.0,
                    help="Alpha of unhighlighted connectors (default: auto from tip count)")
parser.add_argument("--title", default="Tanglegram")
parser.add_argument("--width", type=float, default=0.0, help="Figure width in inches (auto)")
parser.add_argument("--height", type=float, default=0.0, help="Figure height in inches (auto)")
args = parser.parse_args()

# --------------------------------------------------------------------------
# Groups / colours (same machinery as plot_tree.py)
# --------------------------------------------------------------------------
PALETTE = [
    "#E63946", "#F4A261", "#2A9D8F", "#457B9D", "#6A0572", "#1B7837",
    "#00BCD4", "#E67E22", "#8E44AD", "#D81B60", "#827717", "#5D4037",
]
DEFAULT_COLOR = "#333333"

genome_to_group = {}
groups = {}
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


if args.sweeps:
    with open(args.sweeps) as f:
        f.readline()
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
                    "within_mean": float(within_mean),
                    "sister_mean": float(sister_mean),
                    "ratio": float(ratio),
                    "size": int(clade_size),
                    "color": None,
                }
    for sid in sorted(groups):
        groups[sid]["color"] = next_palette_color()


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
        parser.error(f"--highlight takes FILE [COLOR] [LABEL], got {len(spec)}: {spec} "
                     f"(quote labels that contain spaces)")
    path = spec[0]
    color = spec[1] if len(spec) >= 2 else None
    label = spec[2] if len(spec) == 3 else os.path.splitext(os.path.basename(path))[0]
    if color is not None and not is_color_like(color):
        parser.error(f"--highlight {path}: '{color}' is not a valid color "
                     f"(quote hex colours in the shell: \"#1B7837\")")
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
        "label": label, "within_mean": None, "size": assigned,
        "n_in_list": len(genomes), "genomes": set(genomes), "color": color,
    }
    highlight_ids.append(gid)
    print(f"Highlight '{label}' ({path}): {len(genomes)} read, {assigned} assigned, "
          f"{len(genomes) - assigned} already in an earlier group.")

if args.highlight_overwrite:
    for gid, meta in groups.items():
        meta["size"] = sum(1 for v in genome_to_group.values() if v == gid)

GROUP_COLOR = {gid: meta["color"] for gid, meta in groups.items()}


def tip_color(name):
    gid = genome_to_group.get(name)
    return GROUP_COLOR[gid] if gid is not None else DEFAULT_COLOR


def is_highlighted(name):
    return name in genome_to_group


# --------------------------------------------------------------------------
# Panels
# --------------------------------------------------------------------------
class Panel:
    def __init__(self, spec, slot):
        if len(spec) > 2:
            parser.error(f"--{slot} takes FILE [LABEL], got {len(spec)} values "
                         f"(quote labels that contain spaces)")
        self.path = spec[0]
        self.title = spec[1] if len(spec) == 2 else os.path.splitext(os.path.basename(spec[0]))[0]
        self.slot = slot
        print(f"Loading {slot} tree: {self.path}")
        self.tree = Phylo.read(self.path, "newick")
        try:
            self.tree.root_at_midpoint()
        except Exception as e:
            print(f"  midpoint rooting failed ({e}); using the file's rooting")
        self.tree.ladderize()
        self.tips = [t.name for t in self.tree.get_terminals()]
        if any(n is None for n in self.tips):
            parser.error(f"{self.path}: tree has unnamed tips")
        self.n = len(self.tips)
        print(f"  {self.n} tips")
        self.face = 1          # +1 root-left/tips-right, -1 mirrored
        self.x0 = 0.0
        self.width = 1.0

    def order(self):
        return [t.name for t in self.tree.get_terminals()]

    def ypos(self, H):
        o = self.order()
        if len(o) == 1:
            return {o[0]: H / 2.0}
        return {name: i * H / (len(o) - 1) for i, name in enumerate(o)}


panels = [Panel(args.left, "left")]
if args.mid:
    panels.append(Panel(args.mid, "mid"))
panels.append(Panel(args.right, "right"))

slot_index = {p.slot: i for i, p in enumerate(panels)}
if args.order_ref not in slot_index:
    parser.error(f"--order-ref {args.order_ref} but no {args.order_ref} tree was given")
ref_i = slot_index[args.order_ref]

# orientation: left faces right, right faces left, middle as requested
panels[0].face = 1
panels[-1].face = -1
if len(panels) == 3:
    panels[1].face = 1 if args.mid_face == "right" else -1

H = max(p.n for p in panels) - 1.0
if H <= 0:
    H = 1.0

all_tips = set().union(*[set(p.tips) for p in panels])
for i in range(len(panels) - 1):
    a, b = panels[i], panels[i + 1]
    shared = set(a.tips) & set(b.tips)
    print(f"{a.slot} vs {b.slot}: {len(shared)} shared tips, "
          f"{len(set(a.tips) - set(b.tips))} only in {a.slot}, "
          f"{len(set(b.tips) - set(a.tips))} only in {b.slot}")


# --------------------------------------------------------------------------
# Untangling
# --------------------------------------------------------------------------
def count_inversions(seq):
    """Number of out-of-order pairs == number of line crossings."""
    seq = list(seq)
    if len(seq) < 2:
        return 0
    tmp = [0] * len(seq)

    def sort_count(lo, hi):
        if hi - lo < 2:
            return 0
        mid = (lo + hi) // 2
        c = sort_count(lo, mid) + sort_count(mid, hi)
        i, j, k = lo, mid, lo
        while i < mid and j < hi:
            if seq[i] <= seq[j]:
                tmp[k] = seq[i]; i += 1
            else:
                tmp[k] = seq[j]; j += 1
                c += mid - i
            k += 1
        while i < mid:
            tmp[k] = seq[i]; i += 1; k += 1
        while j < hi:
            tmp[k] = seq[j]; j += 1; k += 1
        seq[lo:hi] = tmp[lo:hi]
        return c

    return sort_count(0, len(seq))


def crossings(ya, yb):
    shared = sorted(set(ya) & set(yb), key=lambda n: ya[n])
    return count_inversions([yb[n] for n in shared])


def untangle(panel, ref_y, passes):
    """Rotate internal nodes so tips line up with ref_y as well as possible."""
    last = None
    for _ in range(max(1, passes)):
        cur = panel.ypos(H)
        val = {}

        def rec(clade):
            if clade.is_terminal():
                v = ref_y.get(clade.name, cur[clade.name])
                val[id(clade)] = v
                return v
            vs = [rec(c) for c in clade.clades]
            clade.clades.sort(key=lambda c: val[id(c)])
            v = sum(vs) / len(vs)
            val[id(clade)] = v
            return v

        rec(panel.tree.root)
        now = tuple(panel.order())
        if now == last:
            break
        last = now


ypos = {}
ypos[ref_i] = panels[ref_i].ypos(H)

if not args.no_untangle:
    # propagate outwards from the reference panel
    for i in range(ref_i - 1, -1, -1):
        before = crossings(panels[i].ypos(H), ypos[i + 1])
        untangle(panels[i], ypos[i + 1], args.untangle_passes)
        ypos[i] = panels[i].ypos(H)
        print(f"untangle {panels[i].slot} vs {panels[i+1].slot}: "
              f"{before} -> {crossings(ypos[i], ypos[i+1])} crossings")
    for i in range(ref_i + 1, len(panels)):
        before = crossings(ypos[i - 1], panels[i].ypos(H))
        untangle(panels[i], ypos[i - 1], args.untangle_passes)
        ypos[i] = panels[i].ypos(H)
        print(f"untangle {panels[i-1].slot} vs {panels[i].slot}: "
              f"{before} -> {crossings(ypos[i-1], ypos[i])} crossings")
else:
    for i, p in enumerate(panels):
        ypos[i] = p.ypos(H)
    for i in range(len(panels) - 1):
        print(f"{panels[i].slot} vs {panels[i+1].slot}: "
              f"{crossings(ypos[i], ypos[i+1])} crossings (untangling off)")


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
n_panels = len(panels)
g = min(max(args.gutter, 0.02), 0.45 / max(1, n_panels - 1))
panel_w = (1.0 - g * (n_panels - 1)) / n_panels
for i, p in enumerate(panels):
    p.x0 = i * (panel_w + g)
    p.width = panel_w


def depths(panel):
    d = {}

    def rec(clade, acc):
        x = acc + (clade.branch_length or 0.0)
        d[id(clade)] = x
        for c in clade.clades:
            rec(c, x)

    rec(panel.tree.root, 0.0)
    return d


for p in panels:
    p.depth = depths(p)
    tipd = [p.depth[id(t)] for t in p.tree.get_terminals()]
    if args.clip_pct >= 100:
        maxd = max(tipd)
    else:
        maxd = float(np.percentile(tipd, args.clip_pct))
    p.maxd = maxd if maxd > 0 else 1.0
    p.scale = p.width / (p.maxd * 1.04)

    def X(d, p=p):
        if p.face > 0:
            return min(p.x0 + d * p.scale, p.x0 + p.width)
        return max(p.x0 + p.width - d * p.scale, p.x0)

    p.X = X
    p.tip_edge = p.x0 + p.width if p.face > 0 else p.x0       # side the tips point to
    p.root_edge = p.x0 if p.face > 0 else p.x0 + p.width      # the other side


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
n_rows = int(H) + 1
show_any_label = args.labels_all or (args.labels and bool(genome_to_group))
row_in = 0.10 if show_any_label else 0.06
fig_h = args.height if args.height > 0 else max(7.0, n_rows * row_in)
fig_w = args.width if args.width > 0 else 18.0
# Strip reserved inside each gutter for tip labels, so labels never sit on top
# of the connector lines.
label_strip = g * 0.38 if show_any_label else 0.0
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_facecolor("none")

if args.label_size > 0:
    label_size = args.label_size
else:
    label_size = 6.0 if n_rows <= 120 else max(2.6, 6.0 * 120.0 / n_rows)

if args.connector_alpha > 0:
    base_alpha = args.connector_alpha
else:
    base_alpha = 0.35 if n_rows <= 150 else max(0.08, 0.35 * 150.0 / n_rows)
base_lw = 0.6 if n_rows <= 150 else 0.4


def clade_color(clade):
    ids = {genome_to_group.get(c.name) for c in clade.get_terminals()}
    if len(ids) == 1 and None not in ids:
        return GROUP_COLOR.get(ids.pop(), DEFAULT_COLOR)
    return DEFAULT_COLOR


def draw_panel(p, y):
    nodey = {}

    def ycalc(clade):
        if clade.is_terminal():
            nodey[id(clade)] = y[clade.name]
        else:
            for c in clade.clades:
                ycalc(c)
            ys = [nodey[id(c)] for c in clade.clades]
            nodey[id(clade)] = sum(ys) / len(ys)

    ycalc(p.tree.root)

    def draw(clade, parent_x):
        x = p.X(p.depth[id(clade)])
        yy = nodey[id(clade)]
        col = clade_color(clade)
        lw = 0.7 if col == DEFAULT_COLOR else 1.3
        alpha = 0.8 if col == DEFAULT_COLOR else 0.95
        if parent_x is not None:
            ax.plot([parent_x, x], [yy, yy], color=col, lw=lw, alpha=alpha,
                    solid_capstyle="round", zorder=3)
        if not clade.is_terminal():
            cys = [nodey[id(c)] for c in clade.clades]
            ax.plot([x, x], [min(cys), max(cys)], color=col, lw=lw, alpha=alpha,
                    solid_capstyle="round", zorder=3)
            for c in clade.clades:
                draw(c, x)
        else:
            tc = tip_color(clade.name)
            if tc != DEFAULT_COLOR:
                ax.scatter(x, yy, s=5, color=tc, zorder=5, linewidths=0)
            # stub out to the panel edge so connectors start from a straight column
            if abs(x - p.tip_edge) > 1e-9:
                ax.plot([x, p.tip_edge], [yy, yy], color=tc,
                        lw=0.4, alpha=0.18 if tc == DEFAULT_COLOR else 0.50,
                        zorder=2)

    draw(p.tree.root, None)
    return nodey


print("Drawing...")
for i, p in enumerate(panels):
    draw_panel(p, ypos[i])

# anchor column on the root side of the middle tree, where connectors land
mid_i = 1 if n_panels == 3 else None
if mid_i is not None:
    p = panels[mid_i]
    neighbour = mid_i - 1 if p.face > 0 else mid_i + 1   # panel on the root side
    shared = set(p.tips) & set(panels[neighbour].tips)
    for name in p.tips:
        if name not in shared:
            continue
        yy = ypos[mid_i][name]
        tc = tip_color(name)
        ax.scatter(p.root_edge, yy, s=4, color=tc, zorder=5, linewidths=0,
                   alpha=0.9 if tc != DEFAULT_COLOR else 0.6)
        if args.mid_leaders:
            ax.plot([p.root_edge, p.tip_edge],
                    [yy, yy], color=tc, lw=0.3,
                    alpha=0.20 if tc == DEFAULT_COLOR else 0.45,
                    zorder=1, linestyle=(0, (1, 3)))


def connector_anchor(p, side):
    """x where connectors attach on the 'left'/'right' side of a panel, backing
    off by the label strip when that side is the one carrying tip labels."""
    tip_side = "right" if p.face > 0 else "left"
    pad = label_strip if side == tip_side else 0.0
    return (p.x0 - pad) if side == "left" else (p.x0 + p.width + pad)


n_drawn = 0
for i in range(n_panels - 1):
    a, b = panels[i], panels[i + 1]
    xa = connector_anchor(a, "right")
    xb = connector_anchor(b, "left")
    shared = set(a.tips) & set(b.tips)
    for name in sorted(shared):
        if args.connect == "highlighted" and not is_highlighted(name):
            continue
        ya, yb = ypos[i][name], ypos[i + 1][name]
        hl = is_highlighted(name)
        col = tip_color(name)
        ax.plot([xa, xb], [ya, yb], color=col,
                lw=base_lw * (1.8 if hl else 1.0),
                alpha=min(0.95, base_alpha * 2.2) if hl else base_alpha,
                zorder=4 if hl else 1, solid_capstyle="round")
        n_drawn += 1
print(f"{n_drawn} connectors drawn")

# --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------
if show_any_label:
    for i, p in enumerate(panels):
        for name in p.tips:
            if not (args.labels_all or is_highlighted(name)):
                continue
            yy = ypos[i][name]
            # label just outside the tip column, pointing away from the gutter
            if p.face > 0:
                x = p.tip_edge + label_strip * 0.08
                ha = "left"
            else:
                x = p.tip_edge - label_strip * 0.08
                ha = "right"
            ax.text(x, yy, name, va="center", ha=ha, fontsize=label_size,
                    fontfamily="monospace", color=tip_color(name), alpha=0.9, zorder=7)

# --------------------------------------------------------------------------
# Panel titles, scale bars, legend
# --------------------------------------------------------------------------
for i, p in enumerate(panels):
    ax.text(p.x0 + p.width / 2, H + max(2.0, n_rows * 0.02), p.title,
            ha="center", va="bottom", fontsize=10, fontfamily="monospace",
            color="#222222")
    ax.text(p.x0 + p.width / 2, H + max(0.8, n_rows * 0.008),
            f"n={p.n}" + ("  (mirrored)" if p.face < 0 else ""),
            ha="center", va="bottom", fontsize=7, fontfamily="monospace",
            color="#777777")

    # scale bar in substitutions/site
    raw = p.maxd * 0.2
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    bar = mag
    for nice in (mag, 2 * mag, 5 * mag):
        if nice <= p.maxd * 0.35:
            bar = nice
    blen = bar * p.scale
    bx = p.x0 + p.width * 0.03
    by = -max(1.5, n_rows * 0.015)
    ax.plot([bx, bx + blen], [by, by], color="#444444", lw=1.4, solid_capstyle="butt")
    for xx in (bx, bx + blen):
        ax.plot([xx, xx], [by - n_rows * 0.006, by + n_rows * 0.006],
                color="#444444", lw=1.4)
    ax.text(bx + blen / 2, by - n_rows * 0.012, f"{bar:.1e} subs/site",
            ha="center", va="top", color="#444444", fontsize=6.5,
            fontfamily="monospace")

legend_patches = []
for gid, meta in groups.items():
    if meta["within_mean"] is None:
        lab = f"{meta['label']}  n={meta['size']}"
    else:
        lab = (f"{meta['label']}  n={meta['size']}  within={meta['within_mean']:.1f}  "
               f"sister={meta['sister_mean']:.1f}  ratio={meta['ratio']:.1f}x")
    legend_patches.append(mpatches.Patch(color=meta["color"], label=lab))
legend_patches.append(mpatches.Patch(color=DEFAULT_COLOR, alpha=0.5, label="Not highlighted"))
if legend_patches:
    ax.legend(handles=legend_patches, loc="lower right", fontsize=6.5, framealpha=0.7,
              facecolor="white", edgecolor="#888888", labelcolor="black",
              handlelength=1.2, borderpad=0.8, labelspacing=0.5)

ax.set_xlim(-0.03, 1.03)
ax.set_ylim(-max(4.0, n_rows * 0.05), H + max(5.0, n_rows * 0.05))
ax.axis("off")
ax.set_title(args.title, color="#222222", fontsize=12, fontfamily="monospace", pad=18)

plt.tight_layout(pad=0.5)
plt.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
print(f"Saved: {args.out}")
