#!/usr/bin/env python3
"""100-kb tracks for the two deletion groups of script 35 (run from manuscripts_all/, after
scripts 10, 16, 29 and 35). Same layout as the 100-kb figure of script 13.

Tracks, one bar per 100-kb window of the seven E2 chromosomes:
  0a all-base deletions: fraction of uniquely mappable bases deleted, mean of the 33 strains (script 10)
  0b SNP-based deletions: fraction of SNP genotype calls inside runs of 5 or more missing SNPs (script 10)
  1  events of group A: read evidence "simple deletion" and "flanks joined" at one or both ends
  2  events of group B: any read evidence except "DNA still present", with a TE joined to the flank at one or both ends
     (groups as in 35_deletion_groups.py; an event counts in every window it overlaps; an event in both groups is
     counted in both tracks)
  3  SNP density: SNPs (Data S3) per kb of uniquely mappable sequence
  4  repeat fraction: share of the window annotated as repeat (RepeatMasker and TE library)
  5  effector genes (Data S2) per window
Grey background: windows with fewer than 20 kb of uniquely mappable sequence (deletions cannot be called well there).

Output: deletion_allbase/DataS10_check_and_deletion_tracks/Fig4_window_tracks_two_groups.png / .svg / .pdf
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "deletion_allbase/results_edta/"
OUT = "deletion_allbase/DataS10_check_and_deletion_tracks/"
CA, CB, INK, INK2, GREY, SEP = "#2a78d6", "#eb6834", "#1d1f23", "#5b6068", "#e6e4de", "#c9c6bd"
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": SEP, "svg.fonttype": "none", "pdf.fonttype": 42})

# ---- the two groups, exactly as in script 35
g = pd.read_csv(R + "regions_for_browser.tsv", sep="\t", keep_default_na=False)
DIR, TE = "opposite flank (direct deletion junction)", "transposable element (TE)"
g["A"] = (g["Read evidence"] == "simple deletion") & ((g["Left origin"] == DIR) | (g["Right origin"] == DIR))
g["B"] = (g["Read evidence"] != "DNA still present") & ((g["Left origin"] == TE) | (g["Right origin"] == TE))
# --seven: group B only counts a TE end when the TE is one of the seven families selected in the manuscript
# (as in script 40); outputs get the suffix _sevenTE
SEVEN_ONLY = "--seven" in sys.argv
SEVEN = {"Pot2", "MAGGY/Fosbury", "Grasshopper", "Pyret", "MZ5.TE00000001", "MZ5.TE00000027", "Inago2"}
SUF = "_sevenTE" if SEVEN_ONLY else ""
if SEVEN_ONLY:
    bl = pd.read_csv(R + "breakpoints_long.tsv", sep="\t", keep_default_na=False)
    s7 = set(map(tuple, bl[(bl.origin == TE) & bl.family.isin(SEVEN)][["strain", "event", "side"]].values))
    left7 = np.array([(s, e, "left") in s7 for s, e in zip(g.Strain, g["Deletion event ID"])])
    right7 = np.array([(s, e, "right") in s7 for s, e in zip(g.Strain, g["Deletion event ID"])])
    g["B"] = (g["Read evidence"] != "DNA still present") & (((g["Left origin"] == TE) & left7) | ((g["Right origin"] == TE) & right7))
ev = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t").set_index("event_id")
E = g.groupby("Deletion event ID").agg(A=("A", "any"), B=("B", "any"))
E = E[E.A | E.B].join(ev[["chrom", "start", "end"]])

# ---- windows
W = pd.read_csv(R + "windows100kb_tracks.tsv", sep="\t")
for grp in "AB":
    sub = E[E[grp]]
    W[f"n_{grp}"] = [int(((sub.chrom == c) & (sub.start < e) & (sub.end > s)).sum()) for c, s, e in zip(W.chrom, W.start, W.end)]
chroms = sorted(W.chrom.unique())
off, x0 = {}, 0
for c in chroms:
    off[c] = x0; x0 += W[W.chrom == c].end.max() + 300_000
W["x"] = W.chrom.map(off) + W.start
W["w"] = W.end - W.start
low = W.scorable_bases < 20000

nA, nB = int(E.A.sum()), int(E.B.sum())
# short track labels; units and definitions go in the figure legend (see the docstring)
tracks = [("allbase_deleted_frac", "All-base\ndeletions", "#3b6fb6"),
          ("snpbased_deleted_frac", "SNP-based\ndeletions", "#c2571a"),
          ("n_A", "Simple deletion,\nflanks joined", CA),
          ("n_B", "Seven TEs at one\nor both ends" if SEVEN_ONLY else "TE at one or\nboth ends", CB),
          ("snps_per_scorable_kb", "SNP density", "#4a4f57"),
          ("repeat_frac", "Repeat fraction", "#7a4fa3"),
          ("n_effectors", "Effector genes", "#2f8f6b")]

fig, axes = plt.subplots(len(tracks), 1, figsize=(13, 11.2), sharex=True)
for ax, (col, lab, color) in zip(axes, tracks):
    v = W[col].fillna(0)
    top = v.max() * 1.08 if v.max() > 0 else 1
    ax.bar(W.x[low], [top] * int(low.sum()), width=W.w[low], align="edge", color=GREY, linewidth=0, zorder=0)
    ax.bar(W.x, v, width=W.w, align="edge", color=color, linewidth=0, zorder=2)
    ax.set_ylim(0, top)
    ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=8.5, color=INK2, labelpad=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=7.5, length=2)
    if col in ("n_A", "n_B", "n_effectors"):
        ax.yaxis.get_major_locator().set_params(integer=True)
    for c in chroms[1:]:
        ax.axvline(off[c] - 150_000, color=SEP, lw=0.6, zorder=1)
axes[-1].set_xticks([off[c] + W[W.chrom == c].end.max() / 2 for c in chroms], chroms, fontsize=9, color=INK)
axes[-1].tick_params(axis="x", length=0)
axes[-1].set_xlim(-200_000, x0 - 100_000)
fig.tight_layout()
for ext in ("png", "svg", "pdf"):
    fig.savefig(OUT + f"Fig4_window_tracks_two_groups{SUF}.{ext}", dpi=200)
plt.close(fig)
print("events A", nA, "B", nB, "| windows", len(W), "| low-mappability windows", int(low.sum()))
print("windows with A events:", int((W.n_A > 0).sum()), " with B events:", int((W.n_B > 0).sum()))
