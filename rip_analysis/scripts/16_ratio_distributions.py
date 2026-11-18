#!/usr/bin/env python3
"""Distribution of the three product-to-substrate ratios and the matching substrate indices over all 1,000-bp windows, and where The RIPper's RIP-affected windows fall in it.

Where the RIP-affected windows fall is described without choosing a cut-off between groups: every window is ranked by each ratio
(percentile among all windows of its genome), and the table gives the median rank of the RIP-affected windows and how many of them
are among the 10% of windows with the highest and the lowest values.

Output: results/dinucleotide/ratio_distribution_figure.(png|pdf), ratio_rip_overlap.tsv
"""
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "dinucleotide"
GENOMES = ["E2", "MZ5", "B71", "PoP131"]
STRAIN = {"E2": "E2", "MZ5": "MZ5-1-6", "B71": "B71", "PoP131": "P131"}
RATIOS = {"r1": "Ratio 1", "r2": "Ratio 2", "r3": "Ratio 3"}
# rows of the figure: each ratio followed by the substrate index of the same CpN dinucleotide
PANELS = [("r1", "Ratio 1"), ("S_CA", "CpA index"), ("r2", "Ratio 2"), ("S_CT", "CpT index"), ("r3", "Ratio 3"), ("S_CC", "CpC index")]
DI = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]


def windows(g):
    w = pd.read_csv(D / f"{g}.dinuc_w1000.tsv.gz", sep="\t")
    w = w[w.acgt >= 900].reset_index(drop=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        p, s = w.TA / w.AT, (w.CA + w.TG) / (w.AC + w.GT)
        w["rip"] = (p >= 1.1) & (s <= 0.75) & (p - s > 0)          # The RIPper's three criteria, as in 14_product_ratios.py
        w["r1"] = w.TA / (w.CA + w.TG)
        w["r2"] = (w.TT + w.AA) / (w.CT + w.AG)
        w["r3"] = (w.TC + w.GA) / (w.CC + w.GG)
        # substrate indices: the CpN dinucleotide against the same two bases in reverse order (The RIPper's design).
        # CC reversed is CC, so CpC has no such partner; its index is CpC+GpG observed / expected from base composition.
        n = w[DI].sum(axis=1)
        fc, fg = w[["C" + x for x in "ACGT"]].sum(axis=1) / n, w[["G" + x for x in "ACGT"]].sum(axis=1) / n
        w["S_CA"] = (w.CA + w.TG) / (w.AC + w.GT)
        w["S_CT"] = (w.CT + w.AG) / (w.TC + w.GA)
        w["S_CC"] = ((w.CC + w.GG) / n) / (fc ** 2 + fg ** 2)
    return w.replace([np.inf, -np.inf], np.nan)


rows, W = [], {}
for g in GENOMES:
    w = W[g] = windows(g)
    for k in RATIOS:
        pct = w[k].rank(pct=True) * 100                          # percentile of each window among all windows of the genome
        top, bottom = pct > 90, pct <= 10
        rows.append(dict(genome=g, ratio=k, windows=len(w), rip_windows=int(w.rip.sum()),
                         rip_median_percentile=float(pct[w.rip].median()),
                         rip_in_top10=int((top & w.rip).sum()), rip_in_bottom10=int((bottom & w.rip).sum()),
                         top10_value=float(w[k].quantile(0.9)), top10_windows=int(top.sum()),
                         top10_share_rip=float((top & w.rip).sum() / top.sum())))
T = pd.DataFrame(rows)
T.to_csv(D / "ratio_rip_overlap.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(T.round(3).to_string())

# ---- figure: rows = ratios and substrate indices, columns = genomes; all windows in grey, RIP-affected windows in orange ----
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "pdf.fonttype": 42})
fig, ax = plt.subplots(len(PANELS), 4, figsize=(7.2, 10.4), sharey=True, gridspec_kw={"hspace": 0.62, "wspace": 0.08})
XR = {"r1": (0.03, 10), "r2": (0.1, 100), "r3": (0.1, 10), "S_CA": (0.1, 10), "S_CT": (0.03, 10), "S_CC": (0.1, 10)}
TICKS = {"r1": [0.1, 0.3, 1, 3], "r2": [0.3, 1, 3, 10, 30], "r3": [0.3, 1, 3],     # range ends left unlabelled so neighbouring panels do not collide
         "S_CA": [0.3, 1, 3], "S_CT": [0.1, 0.3, 1, 3], "S_CC": [0.3, 1, 3]}
for i, (k, name) in enumerate(PANELS):
    lo, hi = XR[k]
    bins = np.logspace(math.log10(lo), math.log10(hi), 81)
    for j, g in enumerate(GENOMES):
        a, w = ax[i, j], W[g]
        v = w[k].clip(lo, hi)
        a.hist(v, bins=bins, color="#BBBBBB", edgecolor="none", label="All windows")
        a.hist(v[w.rip], bins=bins, color="#D55E00", edgecolor="none", label="RIP-affected windows (The RIPper)")
        a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lo, hi); a.set_ylim(0.8, 2e4)
        a.set_xticks(TICKS[k]); a.set_xticklabels([f"{t:g}" for t in TICKS[k]]); a.minorticks_off()
        a.spines[["top", "right"]].set_visible(False)
        if i == 0:
            a.set_title(STRAIN[g], fontsize=8.5, fontweight="bold")
        if j == 0:
            a.set_ylabel("Windows")
        a.set_xlabel(name, fontsize=7.5)
ax[0, 0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.25), ncol=2, fontsize=7.5, handlelength=1.2)
fig.savefig(D / "ratio_distribution_figure.pdf", bbox_inches="tight", dpi=300)
fig.savefig(D / "ratio_distribution_figure.png", dpi=300, bbox_inches="tight", facecolor="white")
