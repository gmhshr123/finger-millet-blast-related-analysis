#!/usr/bin/env python3
"""Product-to-substrate ratios for three CpN dinucleotides, and the association of the CpT and CpA indices.

Ratios per 1,000-bp window (counts on one strand; the substrate counts both strands through its reverse complement):
  CpA   TpA / (CpA + TpG)
  CpT   (TpT + ApA) / (CpT + ApG)
  CpC   (TpC + GpA) / (CpC + GpG)
Each ratio is also given as observed / expected, where the expected value comes from the base composition of the window
alone (bases independent of their neighbours). A ratio that is high only because the window is AT-rich has observed/expected
near 1; a dinucleotide-specific loss of substrate and gain of product gives observed/expected above 1.

Association of the CpT index (CpT+ApG)/(TpC+GpA) with the CpA index (CpA+TpG)/(ApC+GpT): Spearman rank correlation over
non-overlapping windows (every second window), overall, after removing the effect of GC content (partial rank correlation),
and within GC classes.

Output results/dinucleotide/: product_ratios_summary.tsv, index_association.tsv, Product_ratio_figure.(png|pdf)
"""
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
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)


def prep(g):
    w = pd.read_csv(D / f"{g}.dinuc_w1000.tsv.gz", sep="\t")
    w = w[w.acgt >= 900].reset_index(drop=True)
    n = w[["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]].sum(axis=1)
    # base frequencies from the first base of each dinucleotide
    f = {b: w[[b + x for x in "ACGT"]].sum(axis=1) / n for b in "ACGT"}
    with np.errstate(divide="ignore", invalid="ignore"):
        w["R_CA"] = w.TA / (w.CA + w.TG)
        w["R_CT"] = (w.TT + w.AA) / (w.CT + w.AG)
        w["R_CC"] = (w.TC + w.GA) / (w.CC + w.GG)
        w["E_CA"] = (f["T"] * f["A"]) / (2 * f["C"] * f["A"])          # TpA / (CpA + TpG); with fT=fA... uses fC*fA + fT*fG
        w["E_CA"] = (f["T"] * f["A"]) / (f["C"] * f["A"] + f["T"] * f["G"])
        w["E_CT"] = (f["T"] ** 2 + f["A"] ** 2) / (f["C"] * f["T"] + f["A"] * f["G"])
        w["E_CC"] = (f["T"] * f["C"] + f["G"] * f["A"]) / (f["C"] ** 2 + f["G"] ** 2)
        for k in ("CA", "CT", "CC"):
            w["OE_" + k] = w["R_" + k] / w["E_" + k]
        w["S_CA"] = (w.CA + w.TG) / (w.AC + w.GT)
        w["S_CT"] = (w.CT + w.AG) / (w.TC + w.GA)
        w["P_CA"] = w.TA / w.AT
    w["gcf"] = w.gc / w.acgt
    w["ripper"] = (w.P_CA >= 1.1) & (w.S_CA <= 0.75) & (w.P_CA - w.S_CA > 0)
    w["low_ct"] = w.S_CT <= 0.5
    return w.replace([np.inf, -np.inf], np.nan)


def spearman(a, b):
    m = a.notna() & b.notna()
    return float(np.corrcoef(a[m].rank(), b[m].rank())[0, 1])


def partial_spearman(a, b, z):
    """Rank correlation of a and b after removing the linear effect of the rank of z from both ranks."""
    m = a.notna() & b.notna() & z.notna()
    ra, rb, rz = a[m].rank().values, b[m].rank().values, z[m].rank().values
    X = np.column_stack([np.ones_like(rz), rz])
    ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    eb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
    return float(np.corrcoef(ea, eb)[0, 1])


W = {g: prep(g) for g in GENOMES}
rows = []
for g, w in W.items():
    for lab, x in (("all windows", w), ("CpT index <= 0.5", w[w.low_ct]), ("RIP-affected (The RIPper criteria)", w[w.ripper])):
        r = dict(genome=g, windows_class=lab, windows=len(x), gc=x.gcf.median() * 100)
        for k in ("CA", "CT", "CC"):
            r["ratio_" + k] = x["R_" + k].median(); r["obs_exp_" + k] = x["OE_" + k].median()
        rows.append(r)
P = pd.DataFrame(rows)
P.to_csv(D / "product_ratios_summary.tsv", sep="\t", index=False, float_format="%.3f")
print(P.round(2).to_string())

assoc = []
for g, w in W.items():
    x = w.iloc[::2]                                     # non-overlapping windows
    r = dict(genome=g, windows=len(x), rho_all=spearman(x.S_CT, x.S_CA), rho_gc_removed=partial_spearman(x.S_CT, x.S_CA, x.gcf),
             rho_CT_gc=spearman(x.S_CT, x.gcf), rho_CA_gc=spearman(x.S_CA, x.gcf),
             rho_oe_CT_CA=spearman(x.OE_CT, x.OE_CA), rho_oe_CT_CC=spearman(x.OE_CT, x.OE_CC))
    for lo, hi in ((0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 1.0)):
        y = x[(x.gcf >= lo) & (x.gcf < hi)]
        r[f"rho_gc_{int(lo * 100)}_{int(hi * 100)}"] = spearman(y.S_CT, y.S_CA); r[f"n_gc_{int(lo * 100)}_{int(hi * 100)}"] = len(y)
    # how often a low CpT window also has a low CpA index, compared with all windows
    r["share_CA_le075_in_low_ct"] = float((x[x.low_ct].S_CA <= 0.75).mean()); r["share_CA_le075_all"] = float((x.S_CA <= 0.75).mean())
    r["share_low_ct_in_CA_le075"] = float((x[x.S_CA <= 0.75].low_ct).mean()); r["share_low_ct_all"] = float(x.low_ct.mean())
    assoc.append(r)
A = pd.DataFrame(assoc)
A.to_csv(D / "index_association.tsv", sep="\t", index=False, float_format="%.3f")
print(A.set_index("genome").T.round(3).to_string())

# figure: CpA index against CpT index, E2, non-overlapping windows, coloured by GC content
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6, "pdf.fonttype": 42})
fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.2), sharex=True, sharey=True, gridspec_kw={"wspace": 0.1})
for a, g in zip(ax, GENOMES):
    x = W[g].iloc[::2]
    sc = a.scatter(x.S_CA, x.S_CT, c=x.gcf * 100, s=1.2, cmap="viridis", vmin=20, vmax=65, linewidths=0, rasterized=True)
    a.set_xlim(0.2, 2.0); a.set_ylim(0, 2.0)
    a.axvline(0.75, color="black", linewidth=0.5, linestyle=(0, (4, 3))); a.axhline(0.5, color="black", linewidth=0.5, linestyle=(0, (4, 3)))
    a.set_title(STRAIN[g], fontsize=8, fontweight="bold", loc="left")
    a.set_xlabel("CpA index")
    a.spines[["top", "right"]].set_visible(False)
ax[0].set_ylabel("CpT index")
cb = fig.colorbar(sc, ax=ax, fraction=0.02, pad=0.01); cb.set_label("GC content (%)"); cb.outline.set_linewidth(0.4)
fig.savefig(D / "Index_association_figure.pdf", bbox_inches="tight", dpi=300)
fig.savefig(D / "Index_association_figure.png", dpi=300, bbox_inches="tight", facecolor="white")
print("figure written")
