#!/usr/bin/env python3
"""Deleted-gene statistics for the TE-associated all-base deletions, counting only WHOLE deleted genes, as in the
manuscript (a partly deleted gene is not counted). Run from manuscripts_all/.

Event sets (TE seen in soft-clipped reads at a breakpoint; regions_for_browser.tsv):
  245  every event with a TE at one or both ends
  191  one TE family only: at one end (121) or at both ends of the same strain's region (70)  [rule: memory te-both-ends-rule]
Deleted bases of a strain: inside one of its deletion regions of the event set (between the read-defined breakpoints
where they exist, otherwise from the first to the last deleted uniquely mappable base; deletion_extent() in
depth_tracks.py, as in the browsers) and called deleted (bases that are not uniquely mappable count as deleted when
they lie inside the region).
Effectors: Data S2 effector sections. Genome: E2 genes on the seven chromosomes.
Output: printed, deletion_allbase/deleted_genes_and_effectors/whole_gene_stats.xlsx and .tsv (input of 39_gene_stats_page.py)
Thresholds: a gene counts as deleted in a strain when at least 50%, 80% or 100% of its uniquely mappable bases are deleted
by that strain's deletions of the event set (100% = whole gene, the manuscript's rule).
"""
import os, sys
from math import comb
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from depth_tracks import deletion_extent

R, I = "deletion_allbase/results_edta/", "deletion_allbase/inputs/"
OUT = "deletion_allbase/deleted_genes_and_effectors/whole_gene_stats.xlsx"
cohort = [l.strip() for l in open(I + "cohort33.txt")]
TE = "transposable element (TE)"

# ---- event sets
g = pd.read_csv(R + "regions_for_browser.tsv", sep="\t", keep_default_na=False)
B = g[(g["Left origin"] == TE) | (g["Right origin"] == TE)].copy()


def fams(s, side):
    out = set()
    for part in s.split(";"):
        part = part.strip()
        if part and part.split(":", 1)[0].strip() == side:
            out |= {x.strip() for x in part.split(":", 1)[1].split(",")}
    return out


B["L"] = B["TE at a breakpoint, from reads"].map(lambda s: fams(s, "left"))
B["R"] = B["TE at a breakpoint, from reads"].map(lambda s: fams(s, "right"))
B["dif"] = [bool(l and r and not (l & r)) for l, r in zip(B.L, B.R)]
B["both"] = [bool(l & r) for l, r in zip(B.L, B.R)]
ev = B.groupby("Deletion event ID").agg(L=("L", lambda x: set().union(*x)), R=("R", lambda x: set().union(*x)),
                                         dif=("dif", "any"), both=("both", "any"))
ev["fam"] = [l | r for l, r in zip(ev.L, ev.R)]
ev["cls"] = np.select([ev.dif, ev.fam.map(len) > 1, ev.both], ["two families at the two ends", "different families in different strains",
                                                                 "one family, both ends"], "one family, one end")
print(ev.cls.value_counts().to_string())
SETS = {"245 events": set(ev.index), "191 events (one TE family)": set(ev[ev.cls.str.startswith("one family")].index)}

# ---- genes, effectors, uniquely mappable bases, deletion calls
genes = []
for l in open(R + "cluster_uniq/E2.genes.bed", "rb").read().decode().split("\n"):
    f = l.replace("\r", "").split("\t")
    if len(f) >= 4 and f[0].startswith("Chr"):
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))
G = pd.DataFrame(genes, columns=["chrom", "a", "b", "gene"]).drop_duplicates("gene")
s2 = pd.read_csv(I + "DataS2_gene_sections.tsv", sep="\t")
EFF = set(s2[s2.section.str.startswith("Effector", na=False)].gene) | {"EcMO3g00052690"}   # marked as effector in Data S10, missing from Data S2
G["eff"] = G.gene.isin(EFF)
sizes = {l.split()[0]: int(l.split()[1]) for l in open(R + "cluster_uniq/E2.chr.genome")}
scor = {c: np.zeros(n, bool) for c, n in sizes.items()}
for l in open(R + "cluster_uniq/scorable.bwa.uniq.bed"):
    c, a, b = l.split("\t")[:3]; scor[c][int(a):int(b)] = True
frag = pd.read_csv(R + "allbase_fragments_strain_level.tsv", sep="\t")
dele = {}
for (st, c), x in frag.groupby(["strain", "chrom"]):
    m = np.zeros(sizes[c], bool)
    for a, b in zip(x.start.values, x.end.values):
        m[a:b] = True
    dele[(st, c)] = m
Gc = {c: x for c, x in G.groupby("chrom")}


THRESH = [0.5, 0.8, 1.0]          # share of a gene's uniquely mappable bases that must be deleted
eff = set(G[G.eff].gene)
K, N = int(G.eff.sum()), len(G)
calls = pd.read_csv(R + "allbase_calls_strain_level.tsv", sep="	")


def gene_shares(ids):
    """per strain: {gene: share of its uniquely mappable bases deleted by this strain's regions of the event set};
    genes without uniquely mappable bases count as fully deleted when they lie entirely inside a region."""
    sub = calls[calls.event_id.isin(ids) & calls.strain.isin(cohort)]
    out, owner = {st: {} for st in cohort}, {}
    for (st, c), x in sub.groupby(["strain", "chrom"]):
        ext = [(s, e, eid) for s, e, eid, _, _ in deletion_extent(st, c, 0, sizes[c]) if eid in set(x.event_id)]
        m = np.zeros(sizes[c], bool)
        for s, e, _ in ext:
            m[s:e] = True
        m &= dele[(st, c)] | ~scor[c]            # deleted = in a deletion call, or not scorable but inside the region
        gx = Gc.get(c)
        for s, e, eid in ext:
            for a, b, gene in zip(gx.a, gx.b, gx.gene):
                if b <= s or a >= e:
                    continue
                sc = scor[c][a:b]; n = int(sc.sum())
                share = (m[a:b] & sc).sum() / n if n else (1.0 if (a >= s and b <= e) else 0.0)
                if share > out[st].get(gene, 0):
                    out[st][gene] = share; owner[(st, gene)] = eid
    return out, owner


rows, per_strain_tabs = {}, {}
for name, ids in SETS.items():
    shares, owner = gene_shares(ids)
    for th in THRESH:
        per = {st: {g for g, v in d.items() if v >= th - 1e-9} for st, d in shares.items()}
        evg = {owner[(st, g)] for st, gs in per.items() for g in gs}
        P = pd.DataFrame([(st, len(v), len(v & eff)) for st, v in per.items()], columns=["strain", "genes", "effectors"])
        P["effector_share_%"] = (P.effectors / P.genes.replace(0, np.nan) * 100).round(1)
        allg = set().union(*per.values()); k, n = len(allg & eff), len(allg)
        p = sum(comb(K, i) * comb(N - K, n - i) for i in range(k, min(n, K) + 1)) / comb(N, n) if n else float("nan")
        sh = P["effector_share_%"].dropna()
        key = f"{name}, {int(th * 100)}%"
        rows[key] = dict(set=name, threshold=f"{int(th * 100)}%", events=len(ids), events_with_gene=len(evg),
                         pct_events_with_gene=round(len(evg) / len(ids) * 100, 1),
                         genes_per_strain_mean=round(P.genes.mean(), 1), genes_per_strain_min=int(P.genes.min()), genes_per_strain_max=int(P.genes.max()),
                         effector_share_per_strain_mean=round(sh.mean(), 1), effector_share_per_strain_min=round(sh.min(), 1),
                         effector_share_per_strain_max=round(sh.max(), 1), unique_genes=n, unique_effectors=k,
                         unique_effector_share=round(k / n * 100, 1) if n else None, genome_effector_share=round(K / N * 100, 1),
                         hypergeometric_p=float(f"{p:.2g}"))
        per_strain_tabs[key] = P
S = pd.DataFrame(rows).T
print(S.drop(columns=["set", "threshold"]).to_string())
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with pd.ExcelWriter(OUT) as w:
    S.to_excel(w, sheet_name="Summary")
    for k, v in per_strain_tabs.items():
        v.to_excel(w, sheet_name=k.replace("events (one TE family)", "ev").replace(" events", " ev").replace(",", "")[:31], index=False)
S.to_csv(OUT.replace(".xlsx", ".tsv"), sep="	")
print("written", OUT)
