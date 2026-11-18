#!/usr/bin/env python3
"""How fragmented are the depth calls, and what separates neighbouring calls?"""
import numpy as np, pandas as pd
R = "deletion_allbase/results/cluster/"
pd.set_option("display.width", 250)
cohort = [l.strip() for l in open("deletion_allbase/inputs/cohort33.txt")]
calls = pd.read_csv("deletion_allbase/results/depth_calls_with_junction.tsv", sep="\t")
calls = calls[calls.strain.isin(cohort)].sort_values(["strain", "chrom", "start"]).reset_index(drop=True)

def bed_cov(bed, sizes):
    cov = {c: np.zeros(n + 1, np.int32) for c, n in sizes.items()}
    for c, s, e in bed.itertuples(index=False):
        if c in cov: cov[c][s:e] = 1
    return {c: np.concatenate(([0], np.cumsum(v[:-1]))) for c, v in cov.items()}
sizes = {l.split()[0]: int(l.split()[1]) for l in open(R + "E2.chr.genome")}
rep = pd.read_csv(R + "E2.repeats.bed", sep="\t", header=None).iloc[:, :3]
te = pd.read_csv(R + "E2.named_TEs.bed", sep="\t", header=None).iloc[:, :3]
cal = pd.read_csv(R + "callable.bwa.bed", sep="\t", header=None)
c_rep, c_te, c_cal = bed_cov(rep, sizes), bed_cov(te, sizes), bed_cov(cal, sizes)
frac = lambda cs, c, s, e: (cs[c][e] - cs[c][s]) / max(e - s, 1)

# gap between a call and the next call in the same strain / chromosome
nx = calls.shift(-1)
same = (nx.strain == calls.strain) & (nx.chrom == calls.chrom)
g = pd.DataFrame({"chrom": calls.chrom[same], "g0": calls.end[same], "g1": nx.start[same].astype(int)})
g["gap"] = g.g1 - g.g0
g["rep"] = [frac(c_rep, c, s, e) for c, s, e in g[["chrom", "g0", "g1"]].itertuples(index=False)]
g["te"] = [frac(c_te, c, s, e) for c, s, e in g[["chrom", "g0", "g1"]].itertuples(index=False)]
g["callable"] = [frac(c_cal, c, s, e) for c, s, e in g[["chrom", "g0", "g1"]].itertuples(index=False)]
g["gapbin"] = pd.cut(g.gap, [0, 20, 50, 100, 200, 500, 1000, 5000, 1e9])
print("gaps between neighbouring calls (present bases in between):")
print(g.groupby("gapbin", observed=True).agg(n=("gap", "size"), repeat_frac=("rep", "mean"), namedTE_frac=("te", "mean"),
                                             callable_frac=("callable", "mean")).round(3))
for mg in [0, 20, 50, 100, 200, 500, 1000]:
    newev = (~same.shift(1, fill_value=False)) | (calls.start - calls.end.shift(1) > mg)
    grp = newev.cumsum()
    m = calls.groupby(grp).agg(strain=("strain", "first"), bp=("deleted_bases", "sum"), s=("start", "min"), e=("end", "max"))
    print(f"merge gap<={mg:5d}: {len(m):6d} calls, {len(m)/33:7.1f} per strain; >=1kb: {(m.bp>=1000).sum()/33:.1f} per strain")
print("\ncalls: repeat content of the call itself")
calls["rep"] = [frac(c_rep, c, s, e) for c, s, e in calls[["chrom", "start", "end"]].itertuples(index=False)]
calls["te"] = [frac(c_te, c, s, e) for c, s, e in calls[["chrom", "start", "end"]].itertuples(index=False)]
print(calls.groupby("sizebin", observed=True).agg(n=("rep", "size"), repeat_frac=("rep", "mean"), namedTE_frac=("te", "mean"),
                                                  junction=("junction", "mean")).round(3))
print("genome-wide repeat fraction:", round(sum(v[-1] for v in c_rep.values()) / sum(sizes.values()), 3),
      "named TE:", round(sum(v[-1] for v in c_te.values()) / sum(sizes.values()), 3))
