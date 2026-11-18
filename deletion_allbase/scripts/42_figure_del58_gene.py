#!/usr/bin/env python3
"""Figure: why EcMO2g00043960 is 'located in' Data S10 deletion Del58 but scored present in Data S2 (strain E13).
Run from manuscripts_all/.

Tracks (Chr02, E2 coordinates): SNP genotypes of E13 (Data S3; '-' = missing), the Data S10 interval of Del58,
E2 genes, uniquely mappable bases, read depth of E13 (uniquely mapped reads / strain median; script 34 profiles),
and the all-base deletions of E13 (deletion_extent, as in the browser).
Output: deletion_allbase/present_genes_and_Del58/Del58_EcMO2g00043960.png/.pdf
"""
import json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from depth_tracks import deletion_extent

R, I = "deletion_allbase/results_edta/", "deletion_allbase/inputs/"
OUT = "deletion_allbase/present_genes_and_Del58/Del58_EcMO2g00043960"
ST, CH, GENE, DEL = "E13", "Chr02", "EcMO2g00043960", "Del58"
A, B = 7_970_300, 7_985_400
INK, INK2, GRID = "#1d1f23", "#5b6068", "#e3e1db"
MISS, CALL, S10C, GENEC, GENEH, DELC, MAPC, LOWC = "#d9822b", "#8a9099", "#d9822b", "#b9bec5", "#2a78d6", "#c00000", "#2f8f5b", "#e3e1db"
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "svg.fonttype": "none", "pdf.fonttype": 42})

# ---- data
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t").set_index("del_id").loc[DEL]
snp = pd.read_csv(I + "DataS3_genotypes.tsv.gz", sep="\t", header=1, usecols=lambda c: c in ["CHROM", "POS", ST], dtype=str)
snp["POS"] = snp.POS.astype(int)
snp = snp[(snp.CHROM == CH) & (snp.POS > A) & (snp.POS < B)]
genes = []
for l in open(R + "cluster_uniq/E2.genes.bed", "rb").read().decode().split("\n"):
    f = l.replace("\r", "").split("\t")
    if len(f) >= 4 and f[0] == CH and int(f[2]) > A and int(f[1]) < B:
        genes.append((int(f[1]), int(f[2]), f[3]))
ga, gb = [(a, b) for a, b, g in genes if g == GENE][0]
scor = pd.read_csv(R + "cluster_uniq/scorable.bwa.uniq.bed", sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
scor = scor[(scor.c == CH) & (scor.e > A) & (scor.s < B)]
med = pd.read_csv(R + "cluster_uniq/depth_calls.bwa.uniq.primary.tsv", sep="\t").groupby("strain").strain_median_depth.first()[ST]
rw = pd.read_csv(R + "region_windows.tsv", sep="\t")
rw = rw[(rw.strain == ST) & (rw.chrom == CH) & (rw.win_end > A) & (rw.win_start < B)].sort_values("win_start")
prof = {}
for line in open(R + "cluster_uniq/region_depth_profiles.jsonl"):
    x = json.loads(line)
    if x["idx"] in set(rw.idx):
        prof[x["idx"]] = x
# stitch the overview profiles of the overlapping windows: each window is used up to the middle of its overlap with the next
ws = list(rw.itertuples())
cuts = [A] + [(ws[i].win_end + ws[i + 1].win_start) / 2 for i in range(len(ws) - 1)] + [B]
xs, ys = [], []
for (w, lo, hi) in zip(ws, cuts[:-1], cuts[1:]):
    u = np.array(prof[w.idx]["unique"]["overview"], float)
    pos = w.win_start + (np.arange(len(u)) + 0.5) * (w.win_end - w.win_start) / len(u)
    k = (pos >= lo) & (pos < hi)
    xs += list(pos[k]); ys += list(u[k] / med)
dels = [(s, e, ev) for s, e, ev, _, _ in deletion_extent(ST, CH, A, B)]
L = pd.read_excel("deletion_allbase/present_genes_and_Del58/present_genes_position.xlsx", sheet_name="Listings")   # script 41
call = L[(L.gene == GENE) & (L.strain == ST) & (L.DataS10_deletion == DEL)].DataS2_call.iloc[0]
gene_rel = np.mean([y for x, y in zip(xs, ys) if ga <= x <= gb])
present_in_gap = snp[(snp[ST] != "-") & (snp.POS > s10.start) & (snp.POS < s10.end)].POS.tolist()

# ---- figure
fig = plt.figure(figsize=(11, 5.4))
ax = fig.add_axes([0.13, 0.11, 0.84, 0.82])
ax.set_xlim(A, B); ax.set_ylim(0, 10.6)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.set_yticks([])
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.3f}"))
ax.set_xlabel(f"Position on {CH} of the E2 reference (Mb)")
lab = lambda y, t: ax.text(A - 250, y, t, ha="right", va="center", fontsize=10, color=INK2)

# gene column highlighted across all tracks
ax.axvspan(ga, gb, color=GENEH, alpha=0.10, lw=0)

# 1. SNPs
y = 9.6
lab(y, f"SNPs in {ST}\n(Data S3)")
for p, v in zip(snp.POS, snp[ST]):
    m = v == "-"
    ax.plot([p, p], [y - 0.35, y + 0.35], color=MISS if m else CALL, lw=1.4 if m else 0.9)
for p in present_in_gap:
    ax.annotate("", xy=(p, y + 0.42), xytext=(p, y + 0.95), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1))
ax.text(np.mean(present_in_gap), y + 1.05, f"only {len(present_in_gap)} SNPs present between the deletions", ha="center", va="bottom", fontsize=9.5)

# 2. Data S10 interval
y = 8.1
lab(y, f"Data S10 {DEL}\n(SNP-based)")
ax.add_patch(Rectangle((s10.start, y - 0.28), s10.end - s10.start, 0.56, color=S10C, alpha=0.85, lw=0))
ax.text((s10.start + s10.end) / 2, y, f"{DEL}: {s10.start:,}–{s10.end:,} ({(s10.end - s10.start + 1) / 1000:.1f} kb), first to last missing SNP",
        ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")

# 3. genes
y = 6.6
lab(y, "E2 genes")
for a, b, g in genes:
    hl = g == GENE
    ax.add_patch(Rectangle((a, y - 0.25), b - a, 0.5, color=GENEH if hl else GENEC, lw=0))
ax.text(gb, y - 0.4, f"{GENE} (no SNPs inside)", ha="right", va="top", fontsize=9.5, color=GENEH, fontweight="bold")
others = [(a, b, g) for a, b, g in genes if g != GENE]
if others:
    ax.text(min(a for a, _, _ in others), y - 0.4,
            ", ".join(g.replace("EcMO2g000", "-") if i else g for i, (_, _, g) in enumerate(others)),
            ha="left", va="top", fontsize=9, color=INK2)

# 4. uniquely mappable
y = 5.4
lab(y, "Uniquely\nmappable")
ax.add_patch(Rectangle((A, y - 0.18), B - A, 0.36, color=LOWC, lw=0))
for s, e in zip(scor.s, scor.e):
    ax.add_patch(Rectangle((max(s, A), y - 0.18), min(e, B) - max(s, A), 0.36, color=MAPC, lw=0))

# 5. read depth
y0, h = 2.1, 2.6
lab(y0 + h / 2, f"Read depth, {ST}\n(÷ strain median)")
for v in (0, 0.5, 1.0):
    ax.plot([A, B], [y0 + v / 1.4 * h] * 2, color=GRID, lw=0.8, zorder=0)
    ax.text(B + 120, y0 + v / 1.4 * h, f"{v:g}", va="center", fontsize=8.5, color=INK2)
ax.plot([A, B], [y0 + 0.05 / 1.4 * h] * 2, color=DELC, lw=0.9, ls=(0, (3, 3)))
ax.text(A + 100, y0 + 0.05 / 1.4 * h + 0.1, "5% cut-off", fontsize=8.5, color=DELC, va="bottom")
ax.plot(xs, [y0 + min(v, 1.4) / 1.4 * h for v in ys], color="#0a6fc2", lw=1.3)
ax.text((ga + gb) / 2, y0 + h + 0.05, f"gene: {gene_rel:.2f} × median", ha="center", va="bottom", fontsize=9.5, color=GENEH, fontweight="bold")

# 6. all-base deletions
y = 1.0
lab(y, f"All-base\ndeletions, {ST}")
for s, e, ev in dels:
    ax.add_patch(Rectangle((s, y - 0.28), e - s, 0.56, color=DELC, lw=0))
    ax.text((s + e) / 2, y - 0.45, ev, ha="center", va="top", fontsize=9, color=DELC)

for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200)
print("written", OUT + ".png/.pdf", "| present SNPs in interval:", present_in_gap, "| gene depth", round(gene_rel, 2), "| call", call)
