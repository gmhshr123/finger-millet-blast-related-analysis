#!/usr/bin/env python3
"""SNP-based deletion heat map of the 226 strains (one row per strain), run from manuscripts_all/.

SNP-based version of the 226-strain heat map, because the main text
reports SNP-based deletions. Same layout as script 46 (all chromosomes, strain-name order, ET / admixed /
KU as in the interactive page), but the deletions follow the manuscript's definition: a run of at least
five consecutive SNPs with a missing call ('-') in a strain, in the genome-wide SNP file (Data S3,
198,454 SNPs). Value per window = fraction of the window's SNPs that lie in such a run in that strain.
Windows without SNPs are grey. Bar on the right = number of SNP-based deletions (runs) per strain.
Output: static_figures_from_interactive_pages/heatmap_226_strains_SNP_based_all_chr_<W>kb_by_name.svg
Usage: 48_snp_based_226_heatmap.py [window_bp]   (default 20000)
"""
import json, os, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

PAGE = "deletion_allbase/interactive_226_strain_heatmap/interactive_226_strains_heatmap.html"
GENO = "deletion_allbase/inputs/DataS3_genotypes.tsv.gz"
OUTD = "deletion_allbase/static_figures_from_interactive_pages/"
W = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
MIN_RUN = 5
os.makedirs(OUTD, exist_ok=True)

# strains, populations, chromosome sizes and colours from the interactive page, so both heat maps match
html = open(PAGE, encoding="utf-8").read()
D, _ = json.JSONDecoder().raw_decode(html[html.index("const D = ") + len("const D = "):])
g = pd.read_csv(GENO, sep="\t", skiprows=1, dtype=str, keep_default_na=False)
g["POS"] = g.POS.astype(int)

# windows: W bp, 300 kb of empty space between chromosomes
gap = max(1, round(300000 / W))
starts, x = {}, 0
for c in D["chroms"]:
    if starts:
        x += gap
    starts[c] = x
    x += -(-D["sizes"][c] // W)
N = x
col_of = {c: starts[c] + g.POS[g.CHROM == c].values // W for c in D["chroms"]}
nsnp = np.zeros(N)
for c in D["chroms"]:
    np.add.at(nsnp, col_of[c], 1)


def pop_of(s):
    if s["p"] == "neither":
        try:
            return "ET" if float(s["m"]) >= 0.85 else "mix"
        except (TypeError, ValueError):
            return "mix"
    return s["p"]


def natural(name):
    return [(0, int(t)) if t.isdigit() else (1, t) for t in re.split(r"(\d+)", name)]


groups = []
for p, lab in [("ET", "ET"), ("mix", "admixed"), ("KU", "KU")]:
    idx = sorted((i for i, s in enumerate(D["strains"]) if pop_of(s) == p), key=lambda i: natural(D["strains"][i]["n"]))
    if idx:
        groups.append((lab, idx))


def strain_row(name):
    """Fraction of each window's SNPs in runs of >= MIN_RUN missing SNPs, and the number of runs."""
    inrun, runs = np.zeros(N), 0
    for c in D["chroms"]:
        miss = (g.loc[g.CHROM == c, name].values == "-").astype(np.int8)
        d = np.diff(np.concatenate(([0], miss, [0])))
        mark = np.zeros(len(miss), bool)
        for i, j in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
            if j - i >= MIN_RUN:
                mark[i:j] = True
                runs += 1
        np.add.at(inrun, col_of[c][mark], 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(nsnp > 0, inrun / nsnp, -1.0), runs


rows, names, low, counts = [], [], [], []
for gi, (lab, idx) in enumerate(groups):
    if gi:
        rows.append(np.full(N, np.nan)); names.append(""); low.append(False); counts.append(np.nan)
    for i in idx:
        s = D["strains"][i]
        v, n = strain_row(s["n"])
        rows.append(v); names.append(s["n"]); low.append(s["d"] < D["mindepth"]); counts.append(n)
M = np.array(rows)
is_gap = np.ones(N, bool)   # the empty space between chromosomes stays white, not grey
for c in D["chroms"]:
    is_gap[starts[c]:starts[c] + -(-D["sizes"][c] // W)] = False

# figure (layout as script 46)
nr = len(rows)
fig_h = 1.6 + 0.09 * nr
fig = plt.figure(figsize=(16, fig_h))
top, bot = 1 - 0.55 / fig_h, 1.0 / fig_h
ax = fig.add_axes([0.11, bot, 0.775, top - bot])
axb = fig.add_axes([0.895, bot, 0.09, top - bot], sharey=ax)

lut = np.array([[int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)] for h in D["lut"]])
rgb = np.ones((nr, N, 3))
val = np.nan_to_num(M, nan=-2)
ok = val >= 0
rgb[ok] = lut[np.clip(np.round(val[ok] * 100), 0, 100).astype(int)]
rgb[(val == -1) & ~is_gap[None, :]] = (0.784, 0.784, 0.784)   # grey: window without SNPs
ax.imshow(rgb, aspect="auto", interpolation="none", extent=(0, N, nr, 0))   # "none": every window kept in the SVG

ax.set_yticks(np.arange(nr) + 0.5)
ax.set_yticklabels(names, fontsize=5.2)
for t, lw in zip(ax.get_yticklabels(), low):
    if lw:
        t.set_color("#80868f"); t.set_fontstyle("italic")
ax.tick_params(axis="y", length=0, pad=1.5)
ax.set_xticks([starts[c] + (-(-D["sizes"][c] // W)) / 2 for c in D["chroms"]])
ax.set_xticklabels(D["chroms"], fontsize=9)
ax.tick_params(axis="x", length=0, top=True, labeltop=True)
for sp in ax.spines.values():
    sp.set_visible(False)
g0 = 0
for gi, (lab, idx) in enumerate(groups):
    if gi:
        g0 += 1
    y = ax.transData.transform((0, g0 + len(idx) / 2))[1]
    fig.text(0.012, fig.transFigure.inverted().transform((0, y))[1], f"{lab}\nn = {len(idx)}", ha="left", va="center",
             fontsize=10, fontweight="bold")
    g0 += len(idx)

cnt = np.array(counts, float)
axb.barh(np.arange(nr) + 0.5, np.nan_to_num(cnt), height=0.8, color="#55595f")
axb.set_xlim(0, np.nanmax(cnt) * 1.05)
axb.tick_params(axis="y", left=False, labelleft=False)
axb.tick_params(axis="x", labelsize=8, top=True, labeltop=True)
axb.set_xlabel("SNP-based\ndeletions (n)", fontsize=9)
for sp in ("right", "left"):
    axb.spines[sp].set_visible(False)
axb.set_ylim(nr, 0)

cax = fig.add_axes([0.11, bot - 0.62 / fig_h, 0.16, 0.10 / fig_h])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=ListedColormap(D["lut"])), cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
cb.ax.tick_params(labelsize=7)
fig.text(0.28, bot - 0.57 / fig_h, f"fraction of the window's SNPs in a SNP-based deletion (>= {MIN_RUN} consecutive missing SNPs; "
         f"{W // 1000} kb windows);  grey = no SNPs in the window;  grey italic name = median depth < 20x", fontsize=8, va="center")

out = OUTD + f"heatmap_226_strains_SNP_based_all_chr_{W // 1000}kb_by_name.svg"
fig.savefig(out)
tot = int(np.nansum(cnt))
print(f"{nr - len(groups) + 1} strains; {N} columns; {int((nsnp > 0).sum())} windows with SNPs; {tot} SNP-based deletions "
      f"(runs of >= {MIN_RUN}); written {out}")
