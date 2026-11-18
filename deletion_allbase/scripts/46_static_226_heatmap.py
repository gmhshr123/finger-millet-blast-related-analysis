#!/usr/bin/env python3
"""Static copy of the interactive 226-strain deletion heat map (script 36), run from manuscripts_all/.

View: all chromosomes, 20 kb windows, order "Strain name",
"Show sample IDs, compact". The figure is drawn from the data embedded in the interactive page, with the
page's own rules, so the two match:
  value = deleted / uniquely mappable bases per window (capped at 1); grey if < 20% uniquely mappable;
  populations: ET, admixed (Data S1 "neither" with ET membership < 0.85; >= 0.85 counts as ET), KU;
  15 empty windows (300 kb) between chromosomes; strain names in grey italic if median depth < 20x;
  bar on the right = deleted Mb in the shown windows.
Output: static_figures_from_interactive_pages/heatmap_226_strains_all_chr_20kb_by_name.{png,pdf,svg}
Usage: 46_static_226_heatmap.py [window_bp]   (default 20000)
"""
import base64, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

SRC = "deletion_allbase/interactive_226_strain_heatmap/interactive_226_strains_heatmap.html"
OUTD = "deletion_allbase/static_figures_from_interactive_pages/"
W = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
os.makedirs(OUTD, exist_ok=True)

html = open(SRC, encoding="utf-8").read()
start = html.index("const D = ") + len("const D = ")
D, _ = json.JSONDecoder().raw_decode(html[start:])
SC = np.frombuffer(base64.b64decode(D["sc"]), "<u2").astype(np.int64)
DEL = [np.frombuffer(base64.b64decode(x), np.uint8).astype(np.int64) for x in D["dels"]]

# windows, as in build(): k bins of 10 kb per window, a gap of 300 kb between chromosomes
k, gap = W // D["bin"], max(1, round(300000 / W))
cols, x, starts = [], 0, {}
for c in D["chroms"]:
    if cols:
        x += gap
    starts[c] = x
    n = -(-D["nb"][c] // k)
    for i in range(n):
        lo, hi = D["off"][c] + i * k, D["off"][c] + min((i + 1) * k, D["nb"][c])
        s, e = i * W, min((i + 1) * W, D["sizes"][c])
        sc = int(SC[lo:hi].sum())
        cols.append((x, lo, hi, sc, sc < 0.2 * (e - s)))
        x += 1
N = x


def pop_of(s):
    if s["p"] == "neither":
        try:
            return "ET" if float(s["m"]) >= 0.85 else "mix"
        except (TypeError, ValueError):
            return "mix"
    return s["p"]


# rows: ET, admixed, KU; strain-name order within each
groups = []
for p, lab in [("ET", "ET"), ("mix", "admixed"), ("KU", "KU")]:
    g = [i for i, s in enumerate(D["strains"]) if pop_of(s) == p]
    g.sort(key=lambda i: [(0, int(t)) if t.isdigit() else (1, t) for t in __import__("re").split(r"(\d+)", D["strains"][i]["n"])])
    if g:
        groups.append((lab, g))

NAN = -1.0
rows, names, low, totals, sep = [], [], [], [], []
for gi, (lab, g) in enumerate(groups):
    if gi:
        sep.append(len(rows))
        rows.append(np.full(N, np.nan)); names.append(""); low.append(False); totals.append(np.nan)
    for i in g:
        a, v, tot = DEL[i], np.full(N, np.nan), 0
        for (cx, lo, hi, sc, lw) in cols:
            d = int(a[lo:hi].sum()) * D["unit"]
            tot += d
            v[cx] = NAN if lw or not sc else min(d / sc, 1)
        rows.append(v); names.append(D["strains"][i]["n"]); low.append(D["strains"][i]["d"] < D["mindepth"])
        totals.append(tot / 1e6)
M = np.array(rows)

# figure: ~0.09 inch per strain row, so the IDs are readable (like "compact" in the page)
nr = len(rows)
fig_h = 1.6 + 0.09 * nr
fig = plt.figure(figsize=(16, fig_h))
top, bot = 1 - 0.55 / fig_h, 1.0 / fig_h
ax = fig.add_axes([0.11, bot, 0.775, top - bot])
axb = fig.add_axes([0.895, bot, 0.09, top - bot], sharey=ax)

cmap = ListedColormap(D["lut"])
rgb = np.ones((nr, N, 3))                                  # white = gap between chromosomes / groups
val = np.nan_to_num(M, nan=-2)
idx = np.clip(np.round(np.where(val >= 0, val, 0) * 100), 0, 100).astype(int)
lut = np.array([[int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)] for h in D["lut"]])
rgb[val >= 0] = lut[idx[val >= 0]]
rgb[val == NAN] = (0.784, 0.784, 0.784)                    # grey: < 20% uniquely mappable
ax.imshow(rgb, aspect="auto", interpolation="none", extent=(0, N, nr, 0))   # "none": every window kept in the SVG/PDF

ax.set_yticks(np.arange(nr) + 0.5)
ax.set_yticklabels(names, fontsize=5.2)
for t, lw in zip(ax.get_yticklabels(), low):
    if lw:
        t.set_color("#80868f"); t.set_fontstyle("italic")
ax.tick_params(axis="y", length=0, pad=1.5)
chrom_mid = [starts[c] + (-(-D["nb"][c] // k)) / 2 for c in D["chroms"]]
for a_ in (ax,):
    a_.set_xticks(chrom_mid); a_.set_xticklabels(D["chroms"], fontsize=9)
    a_.tick_params(axis="x", length=0, top=True, labeltop=True)
for sp in ax.spines.values():
    sp.set_visible(False)
# population labels on the left
g0 = 0
for gi, (lab, g) in enumerate(groups):
    if gi:
        g0 += 1
    y = ax.transData.transform((0, g0 + len(g) / 2))[1]
    fig.text(0.012, fig.transFigure.inverted().transform((0, y))[1], f"{lab}\nn = {len(g)}", ha="left", va="center",
             fontsize=10, fontweight="bold")
    g0 += len(g)

t = np.array(totals)
axb.barh(np.arange(nr) + 0.5, np.nan_to_num(t), height=0.8, color="#55595f")
axb.set_xlim(0, np.nanmax(t) * 1.05)
axb.tick_params(axis="y", left=False, labelleft=False)
axb.tick_params(axis="x", labelsize=8, top=True, labeltop=True)
axb.set_xlabel("Deleted, Mb", fontsize=9)
for sp in ("right", "left"):
    axb.spines[sp].set_visible(False)
axb.set_ylim(nr, 0)

# colour bar and grey key below the map, under the bottom chromosome names
cax = fig.add_axes([0.11, bot - 0.62 / fig_h, 0.16, 0.10 / fig_h])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
cb.ax.tick_params(labelsize=7)
fig.text(0.28, bot - 0.57 / fig_h, f"fraction of the window's uniquely mappable bases deleted ({W // 1000} kb windows);  "
         "grey = < 20% uniquely mappable;  grey italic name = median depth < 20x", fontsize=8, va="center")

base = OUTD + f"heatmap_226_strains_all_chr_{W // 1000}kb_by_name"
for ext, kw in [("png", dict(dpi=300)), ("pdf", {}), ("svg", {})]:
    fig.savefig(f"{base}.{ext}", **kw)
print(nr - len(sep), "strains;", len(cols), "windows; written", base + ".{png,pdf,svg}")
