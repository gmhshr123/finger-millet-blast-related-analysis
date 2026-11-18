#!/usr/bin/env python3
"""Static manuscript versions of the two interactive pages, in 100-kb windows, without titles. Run from manuscripts_all/.

The values are read from the data stored in the two pages and summed to 100 kb by the same rules as the pages, so the
static and the interactive figures are the same:
  Heat map: one row per strain, deleted / uniquely mappable bases per window; populations by the 85% ET membership rule;
            within a population, rows sorted by total deleted sequence; right panel: total deleted sequence.
  Tracks:   as interactive_all_tracks_version3.html, without the simple-deletion and low-mappability tracks; the all-base track leaves out
            the strain-level regions classed as "DNA still present". Windows with < 20% uniquely mappable sequence are left out of
            the all-base and SNP density tracks.
  Grey in the heat map: windows with < 20% uniquely mappable sequence.
Outputs: deletion_allbase/supplement_figures_and_NoteS5/Figure_226_strains_heatmap_100kb.pdf/.png, Figure_genome_tracks_100kb.pdf/.png,
         Figure_226_strains_row_order.tsv
"""
import base64, json, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps

OUT = "deletion_allbase/supplement_figures_and_NoteS5/"
WIN = 100_000
GAP = 3                                                       # windows between chromosomes = 300 kb, as in the pages
INK, INK2 = "#1d1f23", "#5b6068"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": "#9a9fa6", "axes.linewidth": 0.6, "xtick.color": INK2,
                     "ytick.color": INK2, "pdf.fonttype": 42})


def page_data(path):
    return json.loads(re.search(r"const D = (\{.*?\});?\n", open(path, encoding="utf-8").read()).group(1))


# ---------------------------------------------------------------- heat map of the 226 strains
D = page_data("deletion_allbase/interactive_226_strain_heatmap/interactive_226_strains_heatmap.html")
k = WIN // D["bin"]
sc = np.frombuffer(base64.b64decode(D["sc"]), "<u2").astype(float)
dels = np.array([np.frombuffer(base64.b64decode(x), np.uint8) for x in D["dels"]], float) * D["unit"]
Wn = []                                                       # windows: chromosome, first and last bin, x position
x = 0
for c in D["chroms"]:
    if Wn:
        x += GAP
    for i in range(int(np.ceil(D["nb"][c] / k))):
        lo = D["off"][c] + i * k
        hi = D["off"][c] + min((i + 1) * k, D["nb"][c])
        Wn.append((c, lo, hi, x, min((i + 1) * WIN, D["sizes"][c]) - i * WIN)); x += 1
ncol = x
S = pd.DataFrame(D["strains"])
S["m"] = pd.to_numeric(S.m, errors="coerce")
S["pop"] = [("ET" if m >= 0.85 else "Admixed") if p == "neither" else p for p, m in zip(S.p, S.m)]
V = np.full((len(S), ncol), np.nan)
for c, lo, hi, xp, length in Wn:
    b = sc[lo:hi].sum()
    if b >= 0.2 * length and b:
        V[:, xp] = np.minimum(dels[:, lo:hi].sum(1) / b, 1)
S["total"] = dels.sum(1) / 1e6
GROUPS = [("ET", "ET"), ("Admixed", "Admixed"), ("KU", "KU")]
ROWGAP = 4
rows, ticks, y = [], [], 0
for g, lab in GROUPS:
    idx = S[S["pop"] == g].sort_values("total", kind="stable").index.tolist()
    ticks.append((y, y + len(idx), f"{lab}\n(n = {len(idx)})"))
    rows += [(i, y + j) for j, i in enumerate(idx)]
    y += len(idx) + ROWGAP
nrow = y - ROWGAP
cmap = colormaps["magma_r"]
GREY = np.array([0.82, 0.82, 0.82])
img = np.ones((nrow, ncol, 3))                                # white between chromosomes and populations
inmap = sorted(xp for _, _, _, xp, _ in Wn)
for i, yp in rows:
    v = V[i, inmap]
    col = cmap(np.nan_to_num(v))[:, :3]
    col[np.isnan(v)] = GREY
    img[yp, inmap] = col

fig = plt.figure(figsize=(7.2, 7.4))
ax = fig.add_axes([0.095, 0.065, 0.765, 0.925])
bx = fig.add_axes([0.872, 0.065, 0.122, 0.925])
ax.imshow(img, aspect="auto", interpolation="nearest")
starts = {}
for c, _, _, xp, _ in Wn:
    starts.setdefault(c, []).append(xp)
ax.set_xticks([np.mean(v) for v in starts.values()], [c for c in starts], color=INK, fontsize=7)
ax.set_yticks([(a + b - 1) / 2 for a, b, _ in ticks], [t for _, _, t in ticks], color=INK, fontsize=7)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
for i, yp in rows:
    bx.barh(yp, S.total[i], height=1, color="#555555", linewidth=0)
bx.set_ylim(nrow - 0.5, -0.5); bx.set_yticks([])
bx.set_xlabel("Deleted (Mb)", fontsize=7, color=INK, labelpad=1)
bx.tick_params(axis="x", labelsize=6, length=2, width=0.5, pad=1)
bx.spines[["top", "right", "left"]].set_visible(False)
cax = fig.add_axes([0.095, 0.022, 0.2, 0.011])
cax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=cmap, extent=[0, 1, 0, 1])
cax.set_yticks([]); cax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"], fontsize=6)
cax.tick_params(length=2, width=0.5, pad=1)
for sp in cax.spines.values():
    sp.set_linewidth(0.4)
for ext in ("pdf", "png"):
    fig.savefig(OUT + f"Figure_226_strains_heatmap_100kb.{ext}", dpi=600)
plt.close(fig)
o = S.loc[[i for i, _ in rows], ["n", "pop", "m", "d", "total"]].round({"total": 3})
o.columns = ["Strain", "Population", "ET membership", "Median depth", "Deleted (Mb)"]
o.insert(0, "Row, top to bottom", range(1, len(o) + 1))
o.to_csv(OUT + "Figure_226_strains_row_order.tsv", sep="\t", index=False)
print("heat map:", S["pop"].value_counts().to_dict(), "windows", len(Wn), "grey", int(np.isnan(V[0, inmap]).sum()), "below 20x", int((S.d < D["mindepth"]).sum()))
print(o.groupby("Population")["Deleted (Mb)"].describe().round(2).to_string())

# ---------------------------------------------------------------- tracks
T = page_data("deletion_allbase/DataS10_check_and_deletion_tracks/interactive_all_tracks_version3.html")
k = WIN // T["bin"]
nst = T["n_strains"]
# all-base track without the "DNA still present" regions: same weighting as script 37, for all fragments (to check against the
# page) and for the fragments of those regions, which are then subtracted
R = "deletion_allbase/results_edta/"
coh = [l.strip() for l in open("deletion_allbase/inputs/cohort33.txt")]
frag = pd.read_csv(R + "allbase_fragments_strain_level.tsv", sep="\t", usecols=["strain", "chrom", "start", "end", "deleted_bases", "scorable_bases"])
frag = frag[frag.strain.isin(coh)].reset_index(drop=True)
reg = pd.read_csv(R + "regions_for_browser.tsv", sep="\t", keep_default_na=False)
dna = set(zip(*[reg.loc[reg["Read evidence"] == "DNA still present", q] for q in ("Strain", "Deletion event ID")]))
calls = pd.read_csv(R + "allbase_calls_strain_level.tsv", sep="\t", usecols=["strain", "chrom", "start", "end", "event_id"])
drop = calls[[(a, b) in dna for a, b in zip(calls.strain, calls.event_id)]]
isdna = np.zeros(len(frag), bool)
for r in drop.itertuples():
    isdna |= (frag.strain.values == r.strain) & (frag.chrom.values == r.chrom) & (frag.start.values >= r.start) & (frag.end.values <= r.end)
print("regions classed DNA still present:", len(dna), "| matched to calls:", len(drop), "| fragments left out:", int(isdna.sum()), "of", len(frag),
      "| deleted bases left out:", int(frag.deleted_bases[isdna].sum()), "of", int(frag.deleted_bases.sum()))
scor = {c: np.zeros(n, bool) for c, n in T["sizes"].items()}
for l in open(R + "cluster_uniq/scorable.bwa.uniq.bed"):
    c, a, b = l.split("\t")[:3]
    if c in scor:
        scor[c][int(a):int(b)] = True
def binned(fr, c):
    v = np.zeros(T["sizes"][c])
    for r in fr[fr.chrom == c].itertuples():
        v[r.start:r.end + 1] += r.deleted_bases / max(r.scorable_bases, 1)
    n = int(np.ceil(T["sizes"][c] / T["bin"]))
    w = np.zeros(n * T["bin"]); w[:T["sizes"][c]] = v * scor[c]
    return np.round(w.reshape(n, T["bin"]).sum(1))
for c in T["chroms"]:
    assert np.array_equal(binned(frag, c), np.array(T["bins"][c]["dl"], float)), c            # same as the page before the removal
    T["bins"][c]["dl"] = binned(frag[~isdna], c).tolist()
L, x = [], 0
for c in T["chroms"]:
    if L:
        x += GAP
    b = T["bins"][c]
    ev = [e for e in T["events"] if e[1] == c]
    for i in range(int(np.ceil(T["sizes"][c] / WIN))):
        s, e = i * WIN, min((i + 1) * WIN, T["sizes"][c])
        lo, hi = i * k, min((i + 1) * k, len(b["sc"]))
        n = {q: sum(b[q][lo:hi]) for q in ("sc", "rp", "sn", "ef", "dl", "sd")}
        low = n["sc"] < 0.2 * (e - s)
        inw = [z for z in ev if z[2] < e and z[3] > s]
        L.append(dict(chrom=c, x=x, low=low, all=0 if low or not n["sc"] else n["dl"] / (n["sc"] * nst),
                      snpd=n["sd"] / (n["sn"] * nst) if n["sn"] else 0, A=sum(z[4] for z in inw), B=sum(z[5] for z in inw),
                      snp=0 if low or not n["sc"] else n["sn"] / (n["sc"] / 1000), rep=n["rp"] / (e - s), eff=n["ef"]))
        x += 1
L = pd.DataFrame(L)
TR = [("all", "All-base\ndeletions", "#3b6fb6"), ("snpd", "SNP-based\ndeletions", "#c2571a"), ("B", "TE at one or\nboth ends", "#eb6834"),
      ("snp", "SNP density", "#4a4f57"), ("rep", "Repeat fraction", "#7a4fa3"), ("eff", "Effector genes", "#2f8f6b")]
fig, axes = plt.subplots(len(TR), 1, figsize=(7.2, 6.0), sharex=True, gridspec_kw={"hspace": 0.18})
seps = [L[L.chrom == c].x.min() - GAP / 2 - 0.5 for c in T["chroms"][1:]]
for ax, (q, lab, col) in zip(axes, TR):
    v = L[q].astype(float)
    ax.bar(L.x, v, width=1, align="center", color=col, linewidth=0)
    ax.set_ylim(0, v.max() * 1.05)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, integer=q in ("B", "eff")))
    ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=7, color=INK, labelpad=4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=6, length=2, width=0.5, pad=1)
    ax.tick_params(axis="x", length=0)
    for z in seps:
        ax.axvline(z, color=INK2, lw=0.6)
g = L.groupby("chrom").x.mean()
axes[-1].set_xticks(g.values, g.index, fontsize=7, color=INK)
axes[-1].set_xlim(-1, L.x.max() + 1)
fig.subplots_adjust(left=0.15, right=0.995, top=0.995, bottom=0.04)
fig.align_ylabels(axes)
for ext in ("pdf", "png"):
    fig.savefig(OUT + f"Figure_genome_tracks_100kb.{ext}", dpi=300)
plt.close(fig)
print("tracks: windows", len(L), "| grey", int(L.low.sum()), "| B max", L.B.max(), "| all max", round(L["all"].max(), 3),
      "| effectors", int(L.eff.sum()), "| SNP density max", round(L.snp.max(), 1))
