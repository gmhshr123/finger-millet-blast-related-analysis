#!/usr/bin/env python3
"""Are the genes listed in Data S10 ('Genes located in deletion') deleted by the all-base method? Run from manuscripts_all/.

For every Data S10 deletion, every strain listed for it and every gene listed for it: the gene's uniquely mappable
(scored) bases on E2, and how many of them lie in an all-base deletion call of that strain (primary: EDTA, 5% cut-off).
  whole gene deleted   >= 90% of the gene's uniquely mappable bases are in a deletion call of the strain
  partly deleted       some, but < 90%
  not deleted          none
  not testable         < 30 uniquely mappable bases in the gene (same rule as the Data S10 interval check)
A gene counts once per Data S10 deletion and strain; 'unique genes' takes the best result over all its listings.
Input:  inputs/DataS10_events.tsv, inputs/DataS10_coloured_genes.tsv, results_edta/allbase_fragments_strain_level.tsv,
        results_edta/cluster_uniq/scorable.bwa.uniq.bed, results_edta/cluster_uniq/E2.genes.bed
Output: deleted_genes_and_effectors/DataS10_genes_allbase.xlsx (sheets Summary, Genes, Gene_by_strain)
"""
import os
import numpy as np
import pandas as pd

I, R = "deletion_allbase/inputs/", "deletion_allbase/results_edta/"
OUTD = "deletion_allbase/deleted_genes_and_effectors/"
os.makedirs(OUTD, exist_ok=True)
WHOLE, MINB = 0.9, 30

s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")
eff = set(pd.read_csv(I + "DataS10_coloured_genes.tsv", sep="\t").gene)
frag = pd.read_csv(R + "allbase_fragments_strain_level.tsv", sep="\t")          # 0-based, end exclusive
genes = {}
for l in open(R + "cluster_uniq/E2.genes.bed", "rb").read().decode().split("\n"):   # lines carry a stray \r before the gene ID
    f = l.replace("\r", "").split("\t")
    if len(f) < 4:
        continue
    genes[f[3]] = (f[0], int(f[1]), int(f[2]))
sizes = {l.split()[0]: int(l.split()[1]) for l in open(R + "cluster_uniq/E2.chr.genome")}
scor = {c: np.zeros(n, bool) for c, n in sizes.items()}
for l in open(R + "cluster_uniq/scorable.bwa.uniq.bed"):
    c, a, b = l.split("\t")[:3]
    scor[c][int(a):int(b)] = True
by = {k: g for k, g in frag.groupby(["strain", "chrom"])}

rows = []
print("Data S10 deletions with no listed gene:", int(s10.genes.isna().sum()), "of", len(s10))
for x in s10[s10.genes.notna()].itertuples():
    for st in x.strains.split(","):
        for gid in x.genes.split(","):
            c, a, b = genes[gid]
            sc = scor[c][a:b]
            d = np.zeros(b - a, bool)
            g = by.get((st, c))
            if g is not None:
                for p, q in zip(g.start.values, g.end.values):
                    if q > a and p < b:
                        d[max(p, a) - a:min(q, b) - a] = True
            n, k = int(sc.sum()), int((d & sc).sum())
            res = ("not testable" if n < MINB else "whole gene deleted" if k >= WHOLE * n else "partly deleted" if k else "not deleted")
            rows.append(dict(DataS10_deletion=x.del_id, strain=st, gene=gid, effector=gid in eff, chrom=c, gene_start=a + 1,
                             gene_end=b, gene_length=b - a, mappable_bases=n, deleted_bases=k,
                             deleted_share=round(k / n, 3) if n else None, result=res))
L = pd.DataFrame(rows)
rank = {"whole gene deleted": 0, "partly deleted": 1, "not deleted": 2, "not testable": 3}
L["rk"] = L.result.map(rank)
G = (L.groupby("gene").agg(effector=("effector", "first"), chrom=("chrom", "first"), gene_start=("gene_start", "first"),
                           gene_end=("gene_end", "first"), DataS10_deletions=("DataS10_deletion", lambda s: ",".join(dict.fromkeys(s))),
                           listings=("strain", "size"), whole=("result", lambda s: (s == "whole gene deleted").sum()),
                           partly=("result", lambda s: (s == "partly deleted").sum()), none=("result", lambda s: (s == "not deleted").sum()),
                           untestable=("result", lambda s: (s == "not testable").sum()), best=("rk", "min"))
     .reset_index())
inv = {v: k for k, v in rank.items()}
G["best_result"] = G.best.map(inv)
testable = G.listings - G.untestable
G["whole_in_all_testable_strains"] = (G.whole == testable) & (testable > 0)
G = G.drop(columns="best").rename(columns={"whole": "strains_whole_deleted", "partly": "strains_partly_deleted",
                                           "none": "strains_not_deleted", "untestable": "strains_not_testable"})

order = list(rank)
summ = pd.DataFrame({
    "Result": order,
    "Gene listings (gene x strain)": [int((L.result == r).sum()) for r in order],
    "Unique genes (best result over their strains)": [int((G.best_result == r).sum()) for r in order],
    "Effector genes (best result)": [int((G[G.effector].best_result == r).sum()) for r in order]})
summ.loc[len(summ)] = ["Total", len(L), len(G), int(G.effector.sum())]
with pd.ExcelWriter(OUTD + "DataS10_genes_allbase.xlsx") as w:
    summ.to_excel(w, sheet_name="Summary", index=False)
    G.to_excel(w, sheet_name="Genes", index=False)
    L.drop(columns="rk").to_excel(w, sheet_name="Gene_by_strain", index=False)
print(summ.to_string(index=False))
print("unique genes whole-deleted in every testable listed strain:", int(G.whole_in_all_testable_strains.sum()))
print("genes listed under more than one Data S10 deletion:", int(G.DataS10_deletions.str.contains(",").sum()))
