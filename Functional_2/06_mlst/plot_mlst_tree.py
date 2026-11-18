#!/usr/bin/env python3
r"""
plot_mlst_tree.py — map 7-gene MLST sequence types, lineages and AGES sampling
metadata (year, facility) onto a core-genome tree.

What it does
------------
1. Reads a Newick tree, roots it at the midpoint and ladderizes it.
2. Reads `mlst` (tseemann/mlst) output and labels each tip by sequence type (ST).
   ST32 is treated as the background; every other ST is highlighted. Genomes
   with ST "-" are labelled "novel ST" (all alleles exact, unregistered
   combination) or "novel/partial (...)" (listing the alleles flagged with
   ~ or ?).
3. Optionally reads lineage lists (AGES, Sweep7); all remaining genomes are
   "other". With neither list given, every tip is "other": branches are drawn
   grey and the lineage strip and legend are omitted.
4. Optionally reads sampling metadata (isolation date, facility). In the AGES
   project these exist for AGES genomes only; tips without metadata are left
   blank in those strips.
5. For every non-ST32 type, checks whether its genomes form one clade (MRCA test).
6. Draws two figures:
   a) the full tree, and
   b) the tree pruned to AGES genomes only (larger, easier to read the metadata).
   Each figure has:
     - branches coloured by lineage (a branch takes a lineage colour only when
       all tips below it belong to that lineage; mixed branches are grey)
     - strips, left to right: lineage | year | facility | MLST
     - dots + labels at tips for non-ST32 sequence types

Inputs
------
--tree     Newick tree (e.g. IQ-TREE .treefile). Tip names must match genome IDs.
--mlst     Tab-separated `mlst` output (default format, no header). Column 1 is the
           FASTA path; the genome ID is the file name with `.fna` stripped.
--ages     (optional) Text file, one AGES genome ID per line.
--sweep7   (optional) Text file, one Sweep7 genome ID per line.
--meta     (optional) sampling metadata, semicolon-separated, with columns
           `facility`, `date` (DD.MM.YY) and `Seq_number` (= genome ID).
           Facility names garbled by a Mac-Roman/UTF-8 mix-up (e.g. "Gefl√ºgel")
           are repaired automatically.

Outputs (written to --outdir, default: current directory)
-------
<prefix>.pdf / .png               full tree
<prefix>_AGES.pdf / .png          AGES-only subtree (only if both --meta and --ages
                                  are given and some, not all, tips are AGES)
<prefix>_midpoint.nwk             the midpoint-rooted, ladderized tree
stdout                            tip counts, ST counts, monophyly check per ST,
                                  ST x lineage table, and (with --meta) ST x year
                                  and ST x facility tables for AGES

Requirements
------------
Python 3, biopython, matplotlib   (pip install biopython matplotlib)

Example
-------
python3 plot_mlst_tree.py \
    --tree   infantis_clean_tree_noambiguous.treefile \
    --mlst   mlst_results.tsv \
    --ages   AGES_genomes.txt \
    --sweep7 sweep7_genomes.txt \
    --meta   metadata_infantis.csv \
    --outdir mlst_tree_plot \
    --prefix mlst_on_tree
"""
import argparse, copy, csv, os, re
from collections import Counter, defaultdict
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tree", required=True, help="Newick tree file")
ap.add_argument("--mlst", required=True, help="mlst output (tsv)")
ap.add_argument("--ages", help="(optional) AGES genome IDs, one per line")
ap.add_argument("--sweep7", help="(optional) Sweep7 genome IDs, one per line")
ap.add_argument("--meta", help="(optional) sampling metadata csv (';'-separated: facility, date, Seq_number)")
ap.add_argument("--outdir", default=".", help="output directory (default: .)")
ap.add_argument("--prefix", default="mlst_on_tree", help="output file prefix (default: mlst_on_tree)")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)
out = os.path.join(args.outdir, args.prefix)

# ---------------------------------------------------------------- inputs
st, alleles = {}, {}
for line in open(args.mlst):
    f = line.rstrip("\n").split("\t")
    gid = re.sub(r"\.fna$", "", f[0].split("/")[-1])
    st[gid] = f[2]
    alleles[gid] = f[3:]

lineage = {}
for fn, name in ((args.ages, "AGES"), (args.sweep7, "Sweep7")):
    if not fn: continue
    for l in open(fn):
        if l.strip(): lineage[l.strip()] = name
LIN_COL = {"AGES": "#e08a00", "Sweep7": "#0f8a8a", "other": "0.6"}
def lin(g): return lineage.get(g, "other")
# metadata is AGES-only when an AGES list is given; without one it may cover any tip
META_TAG = "AGES" if args.ages else "metadata"
META_SUFFIX = " (AGES only)" if args.ages else ""

def fix_text(s):
    """Repair UTF-8 text that was decoded as Mac-Roman (e.g. 'Gefl√ºgel' -> 'Geflügel')."""
    try:
        return s.encode("mac_roman").decode("utf-8")
    except (UnicodeEncodeError, UnicodeDecodeError):
        return s

year, facility = {}, {}
if args.meta:
    with open(args.meta, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter=";"):
            gid = r["Seq_number"].strip()
            facility[gid] = fix_text(r["facility"].strip())
            d = r["date"].strip()
            if re.fullmatch(r"\d{2}\.\d{2}\.\d{2,4}", d):
                y = d.split(".")[2]
                year[gid] = int(y) + 2000 if len(y) == 2 else int(y)

tree = Phylo.read(args.tree, "newick")
tree.root_at_midpoint()
tree.ladderize()
names = [t.name for t in tree.get_terminals()]
print("tips:", len(names), "| in MLST:", sum(n in st for n in names),
      "| missing from MLST:", [n for n in names if n not in st])
if args.meta:
    print("tips with metadata:", sum(n in facility or n in year for n in names), "/", len(names),
          "| metadata IDs not in tree:", sorted(set(facility) - set(names)))

def label(g):
    s = st.get(g)
    if s is None: return "no MLST"
    if s == "32": return "ST32"
    if s == "-":
        if all(re.fullmatch(r"\w+\(\d+\)", x) for x in alleles[g]): return "novel ST"
        a = [x for x in alleles[g] if not re.fullmatch(r"\w+\(\d+\)", x)]
        return "novel/partial (" + ",".join(a) + ")"
    return "ST" + s

lab = {g: label(g) for g in names}
print(Counter(lab.values()))

groups = defaultdict(list)
for g, l in lab.items():
    if l != "ST32": groups[l].append(g)
for l, gs in sorted(groups.items()):
    ca = tree.common_ancestor(gs) if len(gs) > 1 else None
    n_clade = len(ca.get_terminals()) if ca else 1
    print(f"{l}: n={len(gs)} -> MRCA clade size {n_clade} "
          f"({'monophyletic' if n_clade == len(gs) else 'NOT monophyletic'}) : {sorted(gs)}")

# ---------------------------------------------------------------- colours
palette = {"ST2283": "#d62728", "ST11062": "#1f77b4", "ST11750": "#2ca02c",
           "ST8662": "#9467bd", "ST3815": "#e377c2"}
NOVEL_COL = "#8c564b"
def st_colour(l):
    if l in palette: return palette[l]
    if l.startswith("novel"): return NOVEL_COL
    return None

years = sorted(set(year.values()))
ycmap = plt.get_cmap("viridis")
YEAR_COL = {y: ycmap(i / max(1, len(years) - 1)) for i, y in enumerate(years)}
fac_order = [f for f, _ in Counter(facility.values()).most_common()]
# distinct colours, most-sampled facilities get the most distinct ones
_fac_pal = list(plt.get_cmap("tab10").colors) + [plt.get_cmap("tab20").colors[i] for i in (1, 3, 5, 9, 13, 17, 19, 15, 11, 7)]
FAC_COL = {f: _fac_pal[i % len(_fac_pal)] for i, f in enumerate(fac_order)}

# ---------------------------------------------------------------- drawing
def draw(tr, outbase, title, height, show_meta, label_singletons=True):
    tips = tr.get_terminals()
    tnames = [t.name for t in tips]
    depths = tr.depths()
    if not max(depths.values()):
        depths = tr.depths(unit_branch_lengths=True)
    ypos = {t: i for i, t in enumerate(tips)}
    def set_y(c):
        if c.is_terminal(): return ypos[c]
        ys = [set_y(ch) for ch in c.clades]
        ypos[c] = (min(ys) + max(ys)) / 2
        return ypos[c]
    set_y(tr.root)
    yof = {t.name: ypos[t] for t in tips}

    clade_lin = {}
    for c in tr.find_clades(order="postorder"):
        if c.is_terminal():
            clade_lin[c] = lin(c.name)
        else:
            ls = {clade_lin[ch] for ch in c.clades}
            clade_lin[c] = ls.pop() if len(ls) == 1 else "mixed"
    bcol = lambda c: LIN_COL.get(clade_lin[c], "0.6")

    fig = plt.figure(figsize=(17 if show_meta else 13, height))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.2, 1.7 if show_meta else 1], wspace=0.01)
    ax = fig.add_subplot(gs[0]); lax = fig.add_subplot(gs[1]); lax.axis("off")
    for c in tr.find_clades(order="preorder"):
        x, y = depths[c], ypos[c]
        for ch in c.clades:
            ax.plot([x, x], [y, ypos[ch]], color=bcol(c), lw=0.5)
            ax.plot([x, depths[ch]], [ypos[ch], ypos[ch]], color=bcol(ch), lw=0.5)

    xmax = max(depths[t] for t in tips)
    w = xmax * 0.03
    multi_lin = len({lin(g) for g in tnames}) > 1
    cols = [("lineage", lambda g: LIN_COL[lin(g)])] if multi_lin else []
    if show_meta:
        cols += [("year", lambda g: YEAR_COL.get(year.get(g))),
                 ("facility", lambda g: FAC_COL.get(facility.get(g)))]
    cols += [("MLST", lambda g: st_colour(lab[g]))]
    x0 = xmax * 1.03
    for j, (hdr, fn) in enumerate(cols):
        xs = x0 + j * (w * 1.4)
        for g in tnames:
            c = fn(g)
            if c is not None:
                ax.add_patch(plt.Rectangle((xs, yof[g] - 0.5), w, 1, color=c, lw=0))
        ax.text(xs + w / 2, -2, hdr, ha="center", va="bottom", fontsize=7, rotation=90)
    xt = x0 + len(cols) * (w * 1.4) + w * 0.3

    # tip dots + labels for non-ST32
    in_tree = defaultdict(list)
    for g in tnames:
        in_tree[lab[g]].append(g)
    done = set()
    for g in tnames:
        l = lab[g]; col = st_colour(l)
        if not col: continue
        ax.scatter(depths[tr.find_any(name=g)], yof[g], s=22, color=col, zorder=3,
                   edgecolor="k", lw=0.3)
        members = in_tree[l]
        if len(members) > 1 and l in palette:
            if l not in done:
                done.add(l)
                ys = [yof[m] for m in members]
                ax.text(xt, sum(ys) / len(ys), f"{l} (n={len(members)})",
                        fontsize=8, va="center", color=col, weight="bold")
        elif label_singletons:
            ax.text(xt, yof[g], f"{g}  {l}", fontsize=6.5, va="center")

    # legends
    leg_kw = dict(loc="upper left", fontsize=8, title_fontsize=9, frameon=False, alignment="left")
    present_lins = [k for k in ("AGES", "Sweep7", "other") if any(lin(g) == k for g in tnames)]
    legends = []
    if len(present_lins) > 1:
        legends.append(("Lineage", [Line2D([], [], color=LIN_COL[k], lw=3,
                         label=f"{k} (n={sum(lin(g)==k for g in tnames)})") for k in present_lins]))
    st_h = [Line2D([], [], marker="o", ls="", color=c, mec="k", mew=0.3,
                   label=f"{k} (n={len(in_tree[k])})") for k, c in palette.items() if in_tree[k]]
    nn = sum(len(v) for k, v in in_tree.items() if k.startswith("novel"))
    if nn:
        st_h.append(Line2D([], [], marker="o", ls="", color=NOVEL_COL, mec="k", mew=0.3,
                           label=f"ST unassigned / novel (n={nn})"))
    st_h.append(Line2D([], [], marker="s", ls="", color="white", mec="0.7",
                       label=f"ST32 (n={len(in_tree['ST32'])})"))
    legends.append(("7-gene MLST", st_h))
    if show_meta:
        yc = Counter(year.get(g) for g in tnames if g in year)
        legends.append(("Year" + META_SUFFIX, [Patch(color=YEAR_COL[y], label=f"{y} (n={yc[y]})")
                                             for y in years if yc[y]]))
        fc = Counter(facility.get(g) for g in tnames if g in facility)
        legends.append(("Facility" + META_SUFFIX, [Patch(color=FAC_COL[f], label=f"{f} (n={fc[f]})")
                                                 for f in fac_order if fc[f]]))
    # stack legends top-down in the right-hand legend panel
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ytop = 1.0
    for t, h in legends:
        lg = lax.legend(handles=h, title=t, bbox_to_anchor=(0, ytop), **leg_kw)
        lax.add_artist(lg)
        bb = lg.get_window_extent(renderer).transformed(lax.transAxes.inverted())
        ytop = bb.y0 - 0.01

    ax.set_ylim(-6, len(tips) + 1); ax.invert_yaxis()
    ax.set_xlim(0, xt + xmax * 0.45)
    ax.set_yticks([]); ax.set_xlabel("substitutions / site (midpoint-rooted)")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=11, loc="left")
    plt.savefig(outbase + ".pdf", bbox_inches="tight"); plt.savefig(outbase + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", outbase + ".pdf/.png")

_shown = ([ "lineage"] if lineage else []) + ([f"{META_TAG}"] if args.meta else [])
_shown = (", ".join(_shown) + " and " if _shown else "") + "non-ST32 7-gene MLST types"
draw(tree, out, f"S. Infantis core tree (n={len(names)}) — {_shown}",
     height=30, show_meta=bool(args.meta))

ages_names = [n for n in names if lin(n) == "AGES"]
if args.meta and 1 < len(ages_names) < len(names):
    sub = copy.deepcopy(tree)
    for t in list(sub.get_terminals()):
        if lin(t.name) != "AGES":
            sub.prune(t)
    sub.ladderize()
    n_sub = sub.count_terminals()
    draw(sub, out + "_AGES",
         f"AGES genomes only (n={n_sub}, pruned from the midpoint-rooted tree) — year, facility, MLST",
         height=max(12, n_sub * 0.075), show_meta=True)

Phylo.write(tree, out + "_midpoint.nwk", "newick")
print("wrote", out + "_midpoint.nwk")

# ---------------------------------------------------------------- tables
if lineage:
    print("\nST x lineage:")
    for l in sorted(set(lab.values())):
        c = Counter(lin(g) for g in names if lab[g] == l)
        print(f"  {l:32s} AGES={c['AGES']:3d}  Sweep7={c['Sweep7']:3d}  other={c['other']:3d}")

if args.meta:
    # restrict to AGES only when a lineage list was given, else use every tip
    meta_tips = [g for g in names if lin(g) == "AGES"] if ages_names else names
    print(f"\n{META_TAG}: ST x year")
    print("  " + " " * 32 + "".join(f"{y:>6d}" for y in years))
    for l in sorted({lab[g] for g in meta_tips}):
        c = Counter(year.get(g) for g in meta_tips if lab[g] == l)
        print(f"  {l:32s}" + "".join(f"{c[y]:6d}" for y in years))
    print(f"\n{META_TAG}: non-ST32 genomes with year and facility")
    for g in sorted((g for g in meta_tips if lab[g] != "ST32"), key=lambda g: (lab[g], year.get(g, 0))):
        print(f"  {g}  {lab[g]:28s} {year.get(g, 'NA')}  {facility.get(g, 'NA')}")