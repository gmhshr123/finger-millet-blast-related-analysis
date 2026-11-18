#!/usr/bin/env python3
"""Two figures for the report, in the same journal style as Figure 1.

1. Alignment_figure: the RIPCaL input alignment of E2 LTR25 (results/ripcal/alignment_plots/E2.LTR25.all.input_alignment.fa,
   copied from ripcal/groups/ on the cluster), redrawn. Reference = majority base of each column (gaps not counted; ties left
   undetermined), as in RIPCaL's consensus model. A C-to-T change is classed by the base after the C in the reference;
   a G-to-A change is the same event on the other strand and is classed by the base before the G (TpG = CpA, GpG = CpC, CpG = CpG, ApG = CpT).
2. Genome_overview_figure: the seven chromosomes of each genome, with windows that have a CpT index of 0.5 or below
   (12_alternative_indices.py) and windows that The RIPper website marked as RIP-affected (09_ripper_website.py).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
D = RES / "dinucleotide"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "pdf.fonttype": 42})
SITE_FILL = {"ct": "#0072B2", "ca": "#E69F00", "cc": "#999999", "cg": "#CC79A7"}   # same colours as Figure 1b
SITE_NAME = {"ct": "C followed by T", "ca": "C followed by A", "cc": "C followed by C", "cg": "C followed by G"}

# ---- 1. alignment ---------------------------------------------------------------------------------------------------
aln = {}
for l in open(RES / "ripcal" / "alignment_plots" / "E2.LTR25.all.input_alignment.fa"):
    if l.startswith(">"):
        n = l[1:].strip(); aln[n] = []
    else:
        aln[n].append(l.strip().upper())
names = sorted(aln, key=lambda n: (n.split("_")[1], int(n.split("_")[2].split("-")[0])))   # by chromosome, then position
A = np.array([list("".join(aln[n])) for n in names])
ncol = A.shape[1]
ref = np.full(ncol, "N")
for j in range(ncol):
    b, c = np.unique(A[:, j][np.isin(A[:, j], list("ACGT"))], return_counts=True)
    if len(c) and (c == c.max()).sum() == 1:
        ref[j] = b[c.argmax()]


def neighbour(j, step):
    """The nearest determined, non-gap reference base to the right (step 1) or left (step -1) of column j."""
    k = j + step
    while 0 <= k < ncol and ref[k] not in "ACGT":
        k += step
    return ref[k] if 0 <= k < ncol else "N"


SITE_AFTER_C = {"T": "ct", "A": "ca", "C": "cc", "G": "cg"}
SITE_BEFORE_G = {"A": "ct", "T": "ca", "G": "cc", "C": "cg"}
events = []           # (row, column, kind)
for i in range(len(names)):
    for j in range(ncol):
        q, r = A[i, j], ref[j]
        if q == "-" or r not in "ACGT" or q == r or q not in "ACGT":
            continue
        if r == "C" and q == "T":
            kind = SITE_AFTER_C.get(neighbour(j, 1), "other")
        elif r == "G" and q == "A":
            kind = SITE_BEFORE_G.get(neighbour(j, -1), "other")
        else:
            kind = "other"
        events.append((i, j, kind))
E = pd.DataFrame(events, columns=["row", "col", "kind"])
print("alignment events:", E.kind.value_counts().to_dict())

fig, (a, b) = plt.subplots(2, 1, figsize=(7.2, 4.2), sharex=True, gridspec_kw={"height_ratios": [1.15, 1], "hspace": 0.12})
for i in range(len(names)):
    present = A[i] != "-"
    e = np.diff(np.concatenate(([0], present.astype(int), [0])))
    for s, t in zip(np.flatnonzero(e == 1), np.flatnonzero(e == -1)):
        a.plot([s, t], [i, i], color="#d9d9d9", linewidth=3.2, solid_capstyle="butt", zorder=1)
oth = E[E.kind == "other"]
a.scatter(oth.col, oth.row, marker="|", s=34, linewidths=0.6, color="black", zorder=2)
for k in SITE_FILL:
    x = E[E.kind == k]
    a.scatter(x.col, x.row, marker="|", s=48, linewidths=1.0, color=SITE_FILL[k], zorder=3)
a.set_yticks(range(len(names)))
a.set_yticklabels([n.replace("E2_", "").replace("_", ":").replace("-", "–") for n in names], fontsize=6.5)
a.set_ylim(len(names) - 0.4, -0.6)
a.tick_params(axis="y", length=0)
a.spines[["top", "right", "left"]].set_visible(False)
a.tick_params(axis="x", bottom=False)
# number of C-to-T changes per site class in a sliding window, summed over copies (window as in RIPCaL: 1/50 of the alignment)
win = max(10, round(ncol / 50))
pos = np.arange(ncol)
bottom = np.zeros(ncol)
for k in SITE_FILL:
    v = np.zeros(ncol)
    np.add.at(v, E[E.kind == k].col.values, 1)
    sm = np.convolve(v, np.ones(win), mode="same")
    b.fill_between(pos, bottom, bottom + sm, color=SITE_FILL[k], linewidth=0, step="mid")
    bottom += sm
b.set_xlim(0, ncol)
b.set_ylabel(f"C-to-T changes\nper {win} bp")
b.set_xlabel("Alignment position (bp)")
b.spines[["top", "right"]].set_visible(False)
handles = [Line2D([], [], marker="|", linestyle="none", markersize=8, markeredgewidth=1.2, color=SITE_FILL[k], label=SITE_NAME[k]) for k in SITE_FILL]
handles += [Line2D([], [], marker="|", linestyle="none", markersize=7, markeredgewidth=0.8, color="black", label="Other difference"),
            Patch(color="#d9d9d9", label="Copy sequence (no line: gap)")]
a.legend(handles=handles, ncol=3, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=7, handletextpad=0.3, columnspacing=1.2)
for ax, t in ((a, "a"), (b, "b")):
    ax.text(-0.135, 1.0, t, transform=ax.transAxes, fontweight="bold", fontsize=11, ha="left", va="bottom")
fig.savefig(RES / "ripcal" / "Alignment_figure_E2_LTR25.pdf", bbox_inches="tight")
fig.savefig(RES / "ripcal" / "Alignment_figure_E2_LTR25.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---- 2. genome overview -----------------------------------------------------------------------------------------------
GENOMES = ["E2", "MZ5", "B71", "PoP131"]
STRAIN = {"E2": "E2", "MZ5": "MZ5-1-6", "B71": "B71", "PoP131": "P131"}
CHRS = [f"Chr0{i}" for i in range(1, 8)]
LOWC, RIPC = "#0072B2", "#D55E00"
omitted = {}
fig, axes = plt.subplots(4, 1, figsize=(7.2, 8.4), sharex=True, gridspec_kw={"hspace": 0.32})
for ax, g in zip(axes, GENOMES):
    w = pd.read_csv(D / f"{g}.dinuc_w1000.tsv.gz", sep="\t")
    w = w[w.acgt >= 900]
    low = w[(w.CT + w.AG) / (w.TC + w.GA) <= 0.5]
    web = pd.read_csv(RES / "ripper_website" / f"{g}.website_windows.tsv.gz", sep="\t")
    web["chrom"] = web.SequenceName.str.split().str[0]
    rip = web[web.rip]
    lens = w.groupby("chrom").end.max()
    omitted[g] = (int((~low.chrom.isin(CHRS)).sum()), int((~rip.chrom.isin(CHRS)).sum()))
    for i, c in enumerate(CHRS):
        L = lens[c] / 1e6
        ax.add_patch(plt.Rectangle((0, i - 0.09), L, 0.18, facecolor="#e6e6e6", edgecolor="black", linewidth=0.4))
        x = low[low.chrom == c]
        ax.vlines((x.start + 500) / 1e6, i - 0.42, i - 0.12, color=LOWC, linewidth=0.5)
        x = rip[rip.chrom == c]
        ax.vlines((x.Start + 500) / 1e6, i + 0.12, i + 0.42, color=RIPC, linewidth=0.5)
    ax.set_yticks(range(7)); ax.set_yticklabels([c.replace("Chr0", "Chr ") for c in CHRS])
    ax.set_ylim(6.6, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title(f"{STRAIN[g]}   ({len(low[low.chrom.isin(CHRS)]):,} low-CpT windows, {len(rip[rip.chrom.isin(CHRS)]):,} RIP-affected windows)",
                 loc="left", fontsize=8.5, fontweight="bold")
axes[-1].set_xlabel("Position (Mb)")
axes[-1].set_xlim(0, max(9, axes[-1].get_xlim()[1]))
axes[0].legend(handles=[Line2D([], [], color=LOWC, linewidth=1.2, label="CpT index 0.5 or below (above the chromosome)"),
                        Line2D([], [], color=RIPC, linewidth=1.2, label="RIP-affected, The RIPper website (below the chromosome)")],
               frameon=False, fontsize=7, loc="lower right", bbox_to_anchor=(1.0, 1.08), ncol=1)
fig.savefig(D / "Genome_overview_figure.pdf", bbox_inches="tight")
fig.savefig(D / "Genome_overview_figure.png", dpi=300, bbox_inches="tight", facecolor="white")
pd.DataFrame([(g, *v) for g, v in omitted.items()], columns=["genome", "low_ct_windows_not_on_chr1_7", "ripper_windows_not_on_chr1_7"]) \
    .to_csv(D / "overview_omitted.tsv", sep="\t", index=False)
print("omitted (not on Chr01-07):", omitted)
