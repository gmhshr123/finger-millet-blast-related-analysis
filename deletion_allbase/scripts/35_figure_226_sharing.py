#!/usr/bin/env python3
"""How many of the 226 strains share each deleted base.
Run from manuscripts_all/.

Panel A: for every uniquely mappable base of the E2 chromosomes, the share of ET strains (up) and KU strains (down)
         in which the base lies in an all-base deletion call; shown as the maximum per 10-kb bin.
Panel B: deleted sequence (Mb of uniquely mappable bases) by the number of strains in which a base is deleted,
         split by the populations of those strains.

Inputs: results/cluster_226/depth_calls.bwa.uniq.all226.primary.tsv (all deletion calls, 5% cut-off, >= 100 deleted
        bases), results/cluster_226/depth_summary_all.tsv (median depth of uniquely mapped reads per strain),
        results/cluster_uniq/scorable.bwa.uniq.bed (uniquely mappable bases), Data S1 (population labels).
Strains with a median depth of uniquely mapped reads < MIN_DEPTH are excluded (low depth inflates the calls).
A base inside a call counts as deleted (calls can contain short runs of < 25 non-deleted bases).
Outputs: results/figure_226_sharing.{png,pdf}, .bins10kb.tsv, .spectrum.tsv, .summary.json
"""
import json, warnings
import numpy as np
import pandas as pd
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "deletion_allbase/results/"
OUT = R + "figure_226_sharing"
MIN_DEPTH, BIN = 20, 10_000
CHROMS = [f"Chr0{i}" for i in range(1, 8)]
COL = {"ET": "#1D7A3A", "KU": "#E9A04A", "both": "#7A4FB5", "admixed": "#D04F7E"}   # validated (dataviz validator)
INK, MUTED, GRID = "#1d1f23", "#6b7078", "#e3e1db"

# ---------------------------------------------------------------- strains: population and depth
warnings.simplefilter("ignore")
ws = openpyxl.load_workbook("DataS1_Finger Millet Blast Samples_Final_Ed_July26.xlsx", data_only=True, read_only=True).active
pop = {}
for row in ws.iter_rows(values_only=True):
    if row[0] is not None and str(row[0]).strip().isdigit():
        p = str(row[4]).strip()
        pop[str(row[2]).strip()] = p if p in ("ET", "KU") else "admixed"
dep = pd.read_csv(R + "cluster_226/depth_summary_all.tsv", sep="\t", header=None, names=["strain", "kind", "glen", "mean", "median", "zero"])
dep = dep[dep.kind == "bwa.q20"].set_index("strain")["median"]
strains = sorted(s for s in pop if s in dep.index)
keep = [s for s in strains if dep[s] >= MIN_DEPTH]
dropped = sorted((s for s in strains if dep[s] < MIN_DEPTH), key=lambda s: dep[s])
n = {p: sum(pop[s] == p for s in keep) for p in ("ET", "KU", "admixed")}
print(f"strains: {len(strains)}; kept {len(keep)} (ET {n['ET']}, KU {n['KU']}, admixed {n['admixed']}); "
      f"excluded < {MIN_DEPTH}x: {len(dropped)} {[f'{s} {dep[s]:.0f}x' for s in dropped]}")

# ---------------------------------------------------------------- per-base counts of deleted strains, per population
sizes = {l.split()[0]: int(l.split()[1]) for l in open(R + "cluster_uniq/E2.chr.genome")}
scor = {c: np.zeros(sizes[c], bool) for c in CHROMS}
for l in open(R + "cluster_uniq/scorable.bwa.uniq.bed"):
    c, s, e = l.split("\t")[:3]
    if c in scor:
        scor[c][int(s):int(e)] = True
calls = pd.read_csv(R + "cluster_226/depth_calls.bwa.uniq.all226.primary.tsv", sep="\t")
calls = calls[calls.strain.isin(keep) & calls.chrom.isin(CHROMS)]
cnt = {p: {c: np.zeros(sizes[c] + 1, np.int32) for c in CHROMS} for p in ("ET", "KU", "admixed")}
for (s, c), g in calls.groupby(["strain", "chrom"]):
    d = cnt[pop[s]][c]
    np.add.at(d, g.start.values, 1); np.add.at(d, g.end.values, -1)
for p in cnt:
    for c in CHROMS:
        cnt[p][c] = np.cumsum(cnt[p][c])[:-1] * scor[c]

# ---------------------------------------------------------------- panel A data: max share per 10-kb bin
bins = []
for c in CHROMS:
    nb = -(-sizes[c] // BIN)
    for p in ("ET", "KU", "admixed"):
        a = np.zeros(nb * BIN, np.int32); a[:sizes[c]] = cnt[p][c]
        cnt_bin = a.reshape(nb, BIN).max(1)
        sc = np.zeros(nb * BIN, bool); sc[:sizes[c]] = scor[c]
        bins.append(pd.DataFrame({"chrom": c, "start": np.arange(nb) * BIN, "pop": p, "max_strains": cnt_bin,
                                  "share": cnt_bin / n[p], "scorable_bp": sc.reshape(nb, BIN).sum(1)}))
B = pd.concat(bins)
B.to_csv(OUT + ".bins10kb.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- panel B data: sharing spectrum of deleted bases
EDGES = [1, 2, 3, 6, 11, 26, 51, 101, 10_000]
LAB = ["1", "2", "3–5", "6–10", "11–25", "26–50", "51–100", ">100"]
rows = []
for c in CHROMS:
    et, ku, ad = cnt["ET"][c], cnt["KU"][c], cnt["admixed"][c]
    tot = et + ku + ad
    m = tot > 0
    grp = np.where((et > 0) & (ku > 0), "both", np.where(et > 0, "ET", np.where(ku > 0, "KU", "admixed")))
    rows.append(pd.DataFrame({"n": tot[m], "grp": grp[m]}))
D = pd.concat(rows)
D["cls"] = pd.cut(D.n, EDGES, right=False, labels=LAB)
S = D.groupby(["cls", "grp"], observed=False).size().unstack(fill_value=0).reindex(columns=["ET", "KU", "both", "admixed"], fill_value=0) / 1e6
S.to_csv(OUT + ".spectrum.tsv", sep="\t")
tot_mb = S.values.sum()
summary = dict(min_depth=MIN_DEPTH, strains_total=len(strains), strains_kept=len(keep), kept_by_pop=n,
               excluded=[(s, float(dep[s])) for s in dropped], deleted_mb_any_strain=round(float(tot_mb), 3),
               share_in_one_strain=round(float(S.loc["1"].sum() / tot_mb), 3),
               share_in_2_to_10=round(float(S.loc[["2", "3–5", "6–10"]].values.sum() / tot_mb), 3),
               share_in_more_than_10=round(float(S.loc[["11–25", "26–50", "51–100", ">100"]].values.sum() / tot_mb), 3),
               share_by_group={g: round(float(S[g].sum() / tot_mb), 3) for g in S.columns},
               scorable_mb=round(sum(int(v.sum()) for v in scor.values()) / 1e6, 2))
json.dump(summary, open(OUT + ".summary.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(summary, indent=1, ensure_ascii=False))

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
fig = plt.figure(figsize=(7.2, 10.0))
gs = fig.add_gridspec(len(CHROMS) + 2, 1, height_ratios=[1] * len(CHROMS) + [0.55, 2.4], hspace=0.6)   # row before B: spacer
maxlen = max(sizes[c] for c in CHROMS) / 1e6
for i, c in enumerate(CHROMS):
    ax = fig.add_subplot(gs[i])
    for p, sign in (("ET", 1), ("KU", -1)):
        b = B[(B.chrom == c) & (B["pop"] == p)]
        x = (b.start.values + BIN / 2) / 1e6
        y = np.where(b.scorable_bp.values > 0, b.share.values, np.nan) * sign
        ax.fill_between(x, 0, y, step="mid", color=COL[p], lw=0, alpha=0.9)
    ax.axhline(0, color=MUTED, lw=0.5)
    ax.set_xlim(0, maxlen); ax.set_ylim(-1, 1)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_yticklabels(["1", "", "0", "", "1"])
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.text(0.0, 1.02, c, transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=INK, fontweight="bold")
    if i < len(CHROMS) - 1:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("Position on the E2 chromosome (Mb)")
    if i == 0:
        ax.set_title("A   Share of strains in which each base is deleted (maximum per 10-kb bin)", loc="left", fontsize=9, color=INK, pad=16)
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=COL["ET"], label=f"ET strains (n = {n['ET']}), above the axis"),
                           Patch(color=COL["KU"], label=f"KU strains (n = {n['KU']}), below the axis")],
                  loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, frameon=False, fontsize=7.5, handlelength=1.2, borderaxespad=0)
fig.text(0.04, 0.64, "Share of strains deleted", rotation=90, va="center", ha="center", color=INK)
fig.add_subplot(gs[len(CHROMS)]).axis("off")

ax = fig.add_subplot(gs[-1])
x = np.arange(len(LAB)); bottom = np.zeros(len(LAB))
names = {"ET": "ET strains only", "KU": "KU strains only", "both": "ET and KU strains", "admixed": "admixed strains only"}
for g in ["ET", "KU", "both", "admixed"]:
    v = S[g].values
    ax.bar(x, v, bottom=bottom, color=COL[g], width=0.72, edgecolor="white", linewidth=1.0, label=names[g])
    bottom += v
for xi, t in zip(x, bottom):
    ax.text(xi, t + bottom.max() * 0.015, f"{t:.2f}", ha="center", va="bottom", fontsize=7, color=INK)
ax.set_xticks(x); ax.set_xticklabels(LAB)
ax.set_xlabel(f"Number of strains in which a base is deleted (of {len(keep)})")
ax.set_ylabel("Deleted sequence (Mb)")
ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.set_ylim(0, bottom.max() * 1.12)
ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.55, 1.0), ncol=2)
ax.set_title("B   Deleted sequence by the number of strains sharing it", loc="left", fontsize=9, color=INK, pad=8)
fig.savefig(OUT + ".png", dpi=300, bbox_inches="tight"); fig.savefig(OUT + ".pdf", bbox_inches="tight")
print("written", OUT + ".png/.pdf")
