#!/usr/bin/env python3
"""Benchmark figure: 100-kb tracks along the seven E2 chromosomes (33-strain cohort)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = __import__("os").environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")    # same result set as scripts 10, 16 and 18
W = pd.read_csv(R + "windows100kb_tracks.tsv", sep="\t")
chroms = sorted(W.chrom.unique())
off, x0 = {}, 0
for c in chroms:
    off[c] = x0; x0 += W[W.chrom == c].end.max() + 300_000
W["x"] = W.chrom.map(off) + W.start
W["w"] = W.end - W.start
tracks = [("allbase_deleted_frac", "All-base deletions\n(fraction of uniquely mappable bases\ndeleted, mean of 33 strains)", "#3b6fb6"),
          ("snpbased_deleted_frac", "SNP-based deletions\n(fraction of SNP calls inside\nruns of 5 or more missing SNPs)", "#c2571a"),
          ("snps_per_scorable_kb", "SNP density\n(SNPs per uniquely mappable kb)", "#555555"),
          ("repeat_frac", "Repeat fraction\n(RepeatMasker and TE library)", "#7a4fa3"),
          ("n_effectors", "Effector genes\nper window", "#2f8f6b")]
fig, axes = plt.subplots(len(tracks), 1, figsize=(13, 9), sharex=True)
for ax, (col, lab, color) in zip(axes, tracks):
    ax.bar(W.x, W[col].fillna(0), width=W.w, align="edge", color=color, linewidth=0)
    low = W[W.scorable_bases < 20000]
    ax.bar(low.x, [ax.get_ylim()[1]] * len(low), width=low.w, align="edge", color="#dddddd", linewidth=0, zorder=0)
    ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=7)
    for c in chroms[1:]:
        ax.axvline(off[c] - 150_000, color="#bbbbbb", lw=0.5)
axes[-1].set_xticks([off[c] + W[W.chrom == c].end.max() / 2 for c in chroms]); axes[-1].set_xticklabels(chroms, fontsize=9)
ok = W[(W.scorable_bases >= 20000) & W.snpbased_deleted_frac.notna()]
rho = ok[["allbase_deleted_frac", "snpbased_deleted_frac"]].corr("spearman").iloc[0, 1]
fig.suptitle(f"Deletions in 33 Ethiopian strains per 100-kb window: all-base versus SNP-based calls "
             f"(Spearman rho = {rho:.2f}, {len(ok)} windows; grey = fewer than 20 kb uniquely mappable)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(R + "figure_windows100kb.png", dpi=200)
fig.savefig(R + "figure_windows100kb.pdf")
print("saved; rho", round(rho, 3))
