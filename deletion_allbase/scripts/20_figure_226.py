#!/usr/bin/env python3
"""Deletion distribution along the E2 genome for all 226 strains.

Input: the 100-kb window table written by 19_call_226.sbatch (download export_226/ into
deletion_allbase/results/cluster_226/ first). One row of the heat map per strain, one column per
100-kb window; colour = fraction of the uniquely mappable, non-repeat bases of the window that lie in an
all-base deletion call in that strain. Strains are grouped by the population label of Data S1
(ET, KU, or neither = below 90% membership) and sorted by their total deleted sequence.

Usage: python 20_figure_226.py [window table] [output prefix]
The defaults make the 226-strain figure; pass the 38-strain table to test the script.
"""
import sys, warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd

R = "deletion_allbase/results/"
src = sys.argv[1] if len(sys.argv) > 1 else R + "cluster_226/windows100kb.bwa.uniq.all226.tsv.gz"
out = sys.argv[2] if len(sys.argv) > 2 else R + "figure_226_strains"

warnings.simplefilter("ignore")
ws = openpyxl.load_workbook("DataS1_Finger Millet Blast Samples_Final_Ed_July26.xlsx", data_only=True, read_only=True).active
pop = {}
for row in ws.iter_rows(values_only=True):
    if row[0] is not None and str(row[0]).strip().isdigit():
        p = str(row[4]).strip()
        pop[str(row[2]).strip()] = p if p in ("ET", "KU") else "neither (<90%)"

w = pd.read_csv(src, sep="\t")
w = w[w.strain.isin(pop)]
w["frac"] = w.deleted_bases_in_calls / w.scorable_bases.where(w.scorable_bases >= 20000)      # windows with little scorable sequence are blank
chroms = sorted(w.chrom.unique())
cols = w[["chrom", "start"]].drop_duplicates().sort_values(["chrom", "start"]).reset_index(drop=True)
M = w.pivot_table(index="strain", columns=["chrom", "start"], values="frac", dropna=False).reindex(columns=pd.MultiIndex.from_frame(cols))
tot = w.groupby("strain").deleted_bases_in_calls.sum() / 1e6
order = []
for g in ("ET", "neither (<90%)", "KU"):
    order += sorted([s for s in M.index if pop[s] == g], key=lambda s: tot[s])
M = M.loc[order]

fig, (ax, bx) = plt.subplots(1, 2, figsize=(14, max(4, 0.045 * len(order) + 2)), gridspec_kw={"width_ratios": [12, 1.3], "wspace": 0.02})
cmap = plt.get_cmap("magma_r").copy(); cmap.set_bad("#dddddd")
im = ax.imshow(np.ma.masked_invalid(M.values), aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
edges = np.flatnonzero(cols.chrom.values[1:] != cols.chrom.values[:-1]) + 0.5
for e in edges:
    ax.axvline(e, color="white", lw=1.2)
mids = [np.flatnonzero(cols.chrom.values == c).mean() for c in chroms]
ax.set_xticks(mids); ax.set_xticklabels(chroms, fontsize=9)
y0 = 0
ticks, labels = [], []
for g in ("ET", "neither (<90%)", "KU"):
    n = sum(pop[s] == g for s in order)
    if n:
        ax.axhline(y0 + n - 0.5, color="white", lw=1.2); bx.axhline(y0 + n - 0.5, color="#888888", lw=0.5)
        ticks.append(y0 + n / 2 - 0.5); labels.append(f"{g.replace('neither (<90%)', 'neither')}\n(n = {n})"); y0 += n
ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=9)
bx.barh(np.arange(len(order)), tot[order].values, height=1.0, color="#555555")
bx.set_ylim(len(order) - 0.5, -0.5); bx.set_yticks([]); bx.set_xlabel("deleted, Mb", fontsize=8); bx.tick_params(labelsize=7)
bx.spines[["top", "right"]].set_visible(False)
cb = fig.colorbar(im, ax=[ax, bx], orientation="horizontal", fraction=0.03, pad=0.07, aspect=50)
cb.set_label("fraction of the uniquely mappable, non-repeat bases in a 100-kb window that are deleted in the strain (grey: fewer than 20 kb of such bases)", fontsize=8)
fig.suptitle(f"All-base deletions relative to E2 in {len(order)} strains, per 100-kb window", fontsize=11)
fig.savefig(out + ".png", dpi=200, bbox_inches="tight"); fig.savefig(out + ".pdf", bbox_inches="tight")
summ = pd.DataFrame({"population": pd.Series(pop), "deleted_Mb": tot}).dropna()
print(summ.groupby("population").deleted_Mb.describe().round(3).to_string())
summ.to_csv(out + ".per_strain.tsv", sep="\t", index_label="strain")
print("saved", out + ".png")
