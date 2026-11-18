#!/usr/bin/env python3
"""Dinucleotide indices other than The RIPper's: genome-wide windows and the analysed TE groups.

Input  results/dinucleotide/<genome>.dinuc_w1000.tsv.gz   all 16 dinucleotides per 1,000-bp window (11_dinucleotide_windows.py)
       results/four_genomes/copies, consensus, subfamily_membership.tsv
       results/ripcal/ripcal_per_copy.tsv, results/dinucleotide/ripcal_dinuc_groups.tsv
       inputs/E2.repeat.out                                RepeatMasker annotation of E2

Indices (counts on one strand; each numerator holds a dinucleotide and its reverse complement)
  The RIPper       product   TpA/ApT                 substrate (CpA+TpG)/(ApC+GpT)      for CpA -> TpA
  CpT version      substrate (CpT+ApG)/(TpC+GpA)     built the same way: the target dinucleotide over the same two bases in reverse order
  CpA+CpT version  substrate (CpA+TpG+CpT+ApG)/(ApC+GpT+TpC+GpA)
  RIPCaL           (CpA+TpG)/TpA and (CpT+ApG)/(TpT+ApA)   substrate over product (Hane & Oliver 2008, equations 4 and 7)

Output results/dinucleotide/: genome_index_summary.tsv, window_class_composition.tsv, E2_repeat_share.tsv, te_group_indices.tsv,
       ripcal_other_output.tsv, low_CpT_regions.tsv, Dinucleotide_figure.(png|pdf)
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
D = RES / "dinucleotide"
GENOMES = ["E2", "MZ5", "B71", "PoP131"]
STRAIN = {"E2": "E2", "MZ5": "MZ5-1-6", "B71": "B71", "PoP131": "P131"}
COLOR = {"E2": "#0072B2", "MZ5": "#E69F00", "B71": "#009E73", "PoP131": "#CC79A7"}   # Okabe-Ito
DI = [x + y for x in "ACGT" for y in "ACGT"]
LOW = 0.5          # cut-off used to describe windows with a low CpT index
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)


def indices(t):
    """Add every index to a table that has the 16 dinucleotide counts as columns."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t["P_CA"] = t.TA / t.AT
        t["S_CA"] = (t.CA + t.TG) / (t.AC + t.GT)
        t["comp"] = t.P_CA - t.S_CA
        t["S_CT"] = (t.CT + t.AG) / (t.TC + t.GA)
        t["S_both"] = (t.CA + t.TG + t.CT + t.AG) / (t.AC + t.GT + t.TC + t.GA)
        t["P_CT"] = (t.TT + t.AA) / (2 * t.AT)
        t["R_CA"] = (t.CA + t.TG) / t.TA
        t["R_CT"] = (t.CT + t.AG) / (t.TT + t.AA)
        t["R_both"] = (t.CA + t.TG + t.CT + t.AG) / (t.TA + t.TT + t.AA)
    return t.replace([np.inf, -np.inf], np.nan)


def rank_corr(a, b):
    m = a.notna() & b.notna()
    return float(np.corrcoef(a[m].rank(), b[m].rank())[0, 1])


def runs(w, flag, n=7):
    """Runs of n or more consecutive flagged windows (the definition of a large RIP-affected region in The RIPper)."""
    out = []
    for c, d in w.groupby("chrom", sort=False):
        v = d[flag].values.astype(int)
        e = np.diff(np.concatenate(([0], v, [0])))
        for i, j in zip(np.flatnonzero(e == 1), np.flatnonzero(e == -1)):
            if j - i >= n:
                out.append((c, int(d.start.values[i]), int(d.end.values[j - 1]), j - i))
    return pd.DataFrame(out, columns=["chrom", "start", "end", "windows"])


# ---- 1. genome-wide windows ---------------------------------------------------------------------------------------
copies = pd.read_csv(ROOT / "inputs" / "full_length_copies.tsv", sep="\t")
mem = pd.read_csv(RES / "four_genomes" / "subfamily_membership.tsv", sep="\t")
WIN, summ, comp, regions = {}, [], [], []
for g in GENOMES:
    w = pd.read_csv(D / f"{g}.dinuc_w1000.tsv.gz", sep="\t")
    w = indices(w[w.acgt >= 900].reset_index(drop=True))            # windows that are at least 90% sequence, not N
    w["gcf"] = w.gc / w.acgt
    w["ripper"] = (w.P_CA >= 1.1) & (w.S_CA <= 0.75) & (w.comp > 0)
    w["low_ct"] = w.S_CT <= LOW
    WIN[g] = w
    r = runs(w, "low_ct")
    r.insert(0, "genome", g)
    regions.append(r)
    lo = w[w.low_ct]
    row = dict(genome=g, windows=len(w), ripper_rule=int(w.ripper.sum()), ripper_rule_pct=w.ripper.mean() * 100)
    for k in ("S_CA", "S_CT", "S_both"):
        row[f"{k}_median"] = w[k].median()
        row[f"{k}_le075"] = int((w[k] <= 0.75).sum()); row[f"{k}_le075_pct"] = (w[k] <= 0.75).mean() * 100
        row[f"{k}_le05"] = int((w[k] <= 0.5).sum()); row[f"{k}_le05_pct"] = (w[k] <= 0.5).mean() * 100
    row.update(low_ct_also_ripper=int((w.low_ct & w.ripper).sum()), low_ct_regions=len(r), low_ct_regions_bp=int((r.end - r.start).sum()),
               low_ct_longest_bp=int((r.end - r.start).max()) if len(r) else 0,
               gc_genome=w.gc.sum() / w.acgt.sum() * 100, gc_low_ct=lo.gc.sum() / lo.acgt.sum() * 100,
               gc_ripper=w[w.ripper].gc.sum() / w[w.ripper].acgt.sum() * 100,
               corr_gc_composite=rank_corr(w.gcf, w.comp), corr_gc_product=rank_corr(w.gcf, w.P_CA), corr_gc_S_CA=rank_corr(w.gcf, w.S_CA),
               corr_gc_S_CT=rank_corr(w.gcf, w.S_CT))
    for k in ("P_CA", "S_CA", "P_CT", "R_CA", "R_CT", "R_both"):
        row[f"{k}_median_all"] = w[k].median(); row[f"{k}_median_low_ct"] = lo[k].median()
    # analysed full-length copies that overlap a low-CpT window
    cg = copies[copies.copy_id.isin(mem.copy_id) & (copies.genome == g)]
    hit = 0
    for c, d in lo.groupby("chrom"):
        x = cg[cg.chrom == c]
        hit += sum(((d.start < e) & (d.end > s)).any() for s, e in zip(x.start, x.end))
    row.update(analysed_copies=len(cg), analysed_copies_in_low_ct=int(hit))
    summ.append(row)
    for lab, x in (("low CpT index", lo), ("other windows", w[~w.low_ct])):
        f = x[DI].sum() / x[DI].sum().sum() * 100
        comp.append(dict(genome=g, windows_class=lab, windows=len(x), gc=x.gc.sum() / x.acgt.sum() * 100,
                         CpT=(f.CT + f.AG) / 2, CpA=(f.CA + f.TG) / 2, CpC=(f.CC + f.GG) / 2, CpG=f.CG,
                         TpC=(f.TC + f.GA) / 2, ApC=(f.AC + f.GT) / 2, TpT=(f.TT + f.AA) / 2, TpA=f.TA, ApT=f.AT))
S = pd.DataFrame(summ)
S.to_csv(D / "genome_index_summary.tsv", sep="\t", index=False, float_format="%.4f")
C = pd.DataFrame(comp)
C.to_csv(D / "window_class_composition.tsv", sep="\t", index=False, float_format="%.3f")
pd.concat(regions).to_csv(D / "low_CpT_regions.tsv", sep="\t", index=False)
print(S.set_index("genome").T.round(3).to_string())
print(C.round(2).to_string())

# ---- 2. E2: repeat annotation of the low-CpT windows ----------------------------------------------------------------
rm = pd.read_csv(ROOT / "inputs" / "E2.repeat.out", sep=r"\s+", skiprows=3, header=None, usecols=[4, 5, 6, 9, 10],
                 names=["chrom", "s", "e", "name", "cls"])
rm["k"] = rm.cls.str.split("/").str[0].replace({"Unknown": "Unclassified", "RC": "DNA", "Simple_repeat": "Simple or low complexity",
                                                 "Low_complexity": "Simple or low complexity"})
w = WIN["E2"]
cov = {}
for c, d in w.groupby("chrom", sort=False):
    L = int(d.end.max())
    anyr = np.zeros(L + 2, np.int8)
    for k, r in rm[rm.chrom == c].groupby("k"):
        a = np.zeros(L + 2, np.int8)
        for s, e in zip(r.s.values, r.e.values):
            a[s - 1:e] = 1; anyr[s - 1:e] = 1
        cs = np.concatenate(([0], np.cumsum(a, dtype=np.int64)))
        cov.setdefault(k, pd.Series(0.0, index=w.index)).loc[d.index] = (cs[d.end.values] - cs[d.start.values]) / 1000
    cs = np.concatenate(([0], np.cumsum(anyr, dtype=np.int64)))
    cov.setdefault("Any repeat", pd.Series(0.0, index=w.index)).loc[d.index] = (cs[d.end.values] - cs[d.start.values]) / 1000
RC = pd.DataFrame(cov)
E2rep = pd.DataFrame({"low_CpT_windows": RC[w.low_ct].mean() * 100, "other_windows": RC[~w.low_ct].mean() * 100,
                      "ripper_windows": RC[w.ripper].mean() * 100}).rename_axis("repeat_class").reset_index()
E2rep.to_csv(D / "E2_repeat_share.tsv", sep="\t", index=False, float_format="%.2f")
print(E2rep.round(1).to_string())
print("E2 low-CpT windows at least half repeat:", round(float((RC[w.low_ct]["Any repeat"] >= 0.5).mean()) * 100, 1), "%")

# ---- 3. the analysed TE groups: indices of the copies, of the group reference, and of the genome --------------------
def read_fasta(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l.strip().upper().replace("-", ""))
    return {k: "".join(v) for k, v in d.items()}


def count(seqs):
    c = dict.fromkeys(DI, 0)
    for s in seqs:
        for i in range(len(s) - 1):
            k = s[i:i + 2]
            if k in c:
                c[k] += 1
    return c


meta = pd.read_csv(RES / "ripcal" / "ripcal_vs_ours.tsv", sep="\t")[["group", "Strain", "Family", "Group", "Copies"]]
rows = []
for r in meta.itertuples():
    g, fam, sub = r.group.split(".")
    seqs = read_fasta(RES / "four_genomes" / "copies" / f"{g}.{fam}.fa")
    m = mem[(mem.genome == g) & (mem.family == fam)]
    ids = set(m.copy_id) if sub == "all" else set(m[m.subfamily == sub].copy_id)
    cp = count([s for k, s in seqs.items() if k in ids])
    ref = count(read_fasta(RES / "four_genomes" / "consensus" / f"{r.group}.consensus.fa").values())
    rows.append(dict(group=r.group, what="copies", **cp)); rows.append(dict(group=r.group, what="reference", **ref))
TE = indices(pd.DataFrame(rows))
TE["gc"] = (TE[[k for k in DI if k[0] in "CG"]].sum(axis=1)) / TE[DI].sum(axis=1) * 100
TEw = TE.pivot(index="group", columns="what", values=["P_CA", "S_CA", "comp", "S_CT", "S_both", "R_CA", "R_CT", "gc"])
TEw.columns = [f"{a}_{b}" for a, b in TEw.columns]
TEw = meta.merge(TEw.reset_index(), on="group")
# RIPCaL's own dinucleotide-frequency output for the same copies (it counts AA, CC, GG, TT without overlap)
rc = pd.read_csv(D / "ripcal_dinuc_groups.tsv", sep="\t")
rc.columns = [c.upper() if c != "group" else c for c in rc.columns]
rc = indices(rc)
TEw["ripcal_P_CA"] = TEw.group.map(rc.set_index("group").P_CA); TEw["ripcal_S_CA"] = TEw.group.map(rc.set_index("group").S_CA)
print("largest difference between RIPCaL's dinucleotide mode and the count of script 11, TpA/ApT and (CpA+TpG)/(ApC+GpT):",
      float((TEw.ripcal_P_CA - TEw.P_CA_copies).abs().max()), float((TEw.ripcal_S_CA - TEw.S_CA_copies).abs().max()))
TEw.to_csv(D / "te_group_indices.tsv", sep="\t", index=False, float_format="%.4f")
print(TEw[["group", "P_CA_copies", "P_CA_reference", "S_CA_copies", "S_CA_reference", "S_CT_copies", "S_CT_reference"]].round(2).to_string())

# other RIPCaL output per group: transitions, transversions, dominance
PC = pd.read_csv(RES / "ripcal" / "ripcal_per_copy.tsv", sep="\t")
SM = pd.read_csv(RES / "ripcal" / "ripcal_summary.tsv", sep="\t")
o = PC[PC.reference == "consensus"].groupby("group")[["Transitions", "Transversions", "Mismatches"]].sum().reset_index()
o["ti_tv"] = o.Transitions / o.Transversions
o = o.merge(SM.loc[SM.reference == "consensus", ["group", "dominance ct->tt", "dominance ca->ta", "dominance cc->tc", "dominance cg->tg"]], on="group")
o = meta.merge(o, on="group")
o.to_csv(D / "ripcal_other_output.tsv", sep="\t", index=False, float_format="%.4f")
print(o.round(2).to_string())

# ---- 4. figures (journal style) -------------------------------------------------------------------------------------
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "pdf.fonttype": 42})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True, gridspec_kw={"wspace": 0.08})
bins = np.arange(0, 2.0001, 0.05)
for a, k, lab in ((ax[0], "S_CA", "(CpA + TpG) / (ApC + GpT)"), (ax[1], "S_CT", "(CpT + ApG) / (TpC + GpA)")):
    for g in GENOMES:
        h, _ = np.histogram(WIN[g][k].dropna(), bins=bins)             # the few windows above 2 are not drawn
        a.step(bins[:-1], h / len(WIN[g]) * 100, where="post", color=COLOR[g], linewidth=0.9, label=STRAIN[g])
    a.set_yscale("log"); a.set_ylim(0.005, 30); a.set_xlim(0, 2)
    a.set_xlabel(lab); a.set_xticks([0, 0.5, 1, 1.5, 2]); a.set_xticklabels(["0", "0.5", "1.0", "1.5", "2.0"])
    a.spines[["top", "right"]].set_visible(False)
ax[0].set_ylabel("Windows (%)")
ax[0].set_yticks([0.01, 0.1, 1, 10]); ax[0].set_yticklabels(["0.01", "0.1", "1", "10"])
ax[0].legend(frameon=False, fontsize=7, loc="upper left", handlelength=1.4)
for a, t in zip(ax, "ab"):
    a.text(0, 1.04, t, transform=a.transAxes, fontweight="bold", fontsize=11, ha="left", va="bottom")
fig.savefig(D / "Dinucleotide_figure.pdf", bbox_inches="tight")
fig.savefig(D / "Dinucleotide_figure.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

# the genome overview (Figure 4 of the report) is drawn by 13_figures.py
print("figures written")
