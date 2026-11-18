#!/usr/bin/env python3
"""Genes listed in a Data S10 deletion but scored present (Y) in Data S2 for a strain that carries the deletion:
where do they lie in the deletion? Run from manuscripts_all/.

Unit: one gene in one strain of one Data S10 deletion (as in script 37).
Data S2 call (first sheet, strain columns): 'Y-<covered bp>-<depth>-<% covered>' present, 'N-...' absent, 'A-...' ambiguous.
A gene with several mRNAs takes the call of its first mRNA.
Position: the gene crosses an end of the Data S10 interval, or lies inside it. For genes inside: distance of the gene
midpoint to the nearer end, as a share of the deletion length (0 = at an end, 0.5 = in the middle; random position
gives a uniform spread of 0-0.5, mean 0.25).
Also: the all-base result of the same listing (script 37: whole / partly / not deleted).
Output: printed, deletion_allbase/present_genes_and_Del58/present_genes_position.xlsx
"""
import os, re, warnings
import numpy as np
import pandas as pd
import openpyxl
from scipy.stats import mannwhitneyu, fisher_exact, kstest

warnings.filterwarnings("ignore")
I, R = "deletion_allbase/inputs/", "deletion_allbase/results_edta/"
OUTD = "deletion_allbase/present_genes_and_Del58/"
os.makedirs(OUTD, exist_ok=True)

ws = openpyxl.load_workbook("DataS2_AnnotatedGenes_PresenceAbsence_fixed.xlsx", read_only=True, data_only=True).worksheets[0]
rows = list(ws.iter_rows(values_only=True))
hdr = rows[3]
col = {h: i for i, h in enumerate(hdr) if isinstance(h, str) and re.fullmatch(r"[A-Z]+\d+", h)}
call = {}
for r in rows[4:]:
    if r[1] and str(r[1]).startswith("EcMO") and str(r[1]) not in call:
        call[str(r[1])] = {s: r[i] for s, i in col.items()}

genes = {}
for l in open(R + "cluster_uniq/E2.genes.bed", "rb").read().decode().split("\n"):
    f = l.replace("\r", "").split("\t")
    if len(f) >= 4:
        genes[f[3]] = (f[0], int(f[1]) + 1, int(f[2]))          # 1-based, inclusive
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")
ab = pd.read_excel("deletion_allbase/deleted_genes_and_effectors/DataS10_genes_allbase.xlsx", sheet_name="Gene_by_strain")
ab = ab.set_index(["DataS10_deletion", "strain", "gene"]).result

out = []
for x in s10[s10.genes.notna()].itertuples():
    L = x.end - x.start + 1
    for st in x.strains.split(","):
        for gid in x.genes.split(","):
            c, a, b = genes[gid]
            v = (call.get(gid) or {}).get(st)
            pa = str(v)[0] if v and str(v)[0] in "YNA" else "no call"
            inside = a >= x.start and b <= x.end
            mid = (a + b) / 2
            rel = min(mid - x.start, x.end - mid) / L if inside else np.nan
            out.append(dict(DataS10_deletion=x.del_id, strain=st, gene=gid, DataS2_call=v, call=pa, deletion_length=L,
                            gene_start=a, gene_end=b, gene_length=b - a + 1,
                            position="inside the deletion" if inside else "crosses an end of the deletion",
                            midpoint_to_nearer_end_bp=round(min(mid - x.start, x.end - mid)) if inside else np.nan,
                            midpoint_to_nearer_end_share=round(rel, 3) if inside else np.nan,
                            allbase=ab.get((x.del_id, st, gid))))
D = pd.DataFrame(out)
print("listings:", len(D)); print(D.call.value_counts().to_string()); print()
D["present"] = D.call == "Y"

# 1. crosses an end vs inside
T1 = pd.crosstab(D.position, D.present.map({True: "present (Y)", False: "not present (N/A)"}), margins=True)
T1["% present"] = (T1["present (Y)"] / T1["All"] * 100).round(1)
f1 = fisher_exact(pd.crosstab(D.position, D.present).values)
# 2. inside: position from the nearer end
In = D[D.position == "inside the deletion"]
y, n = In[In.present].midpoint_to_nearer_end_share, In[~In.present].midpoint_to_nearer_end_share
mw = mannwhitneyu(y, n)
ks_y = kstest(y, "uniform", args=(0, 0.5))
bins = pd.cut(In.midpoint_to_nearer_end_share, [0, .1, .2, .3, .4, .5001], right=False,
              labels=["0-10% (near an end)", "10-20%", "20-30%", "30-40%", "40-50% (middle)"])
T2 = pd.crosstab(bins, In.present.map({True: "present (Y)", False: "not present (N/A)"}))
T2["% present"] = (T2["present (Y)"] / T2.sum(1) * 100).round(1)
# 3. size of the deletion
sz = pd.cut(D.deletion_length, [0, 1000, 5000, 20000, 10**9], labels=["< 1 kb", "1-5 kb", "5-20 kb", "> 20 kb"])
T3 = pd.crosstab(sz, D.present.map({True: "present (Y)", False: "not present (N/A)"}))
T3["% present"] = (T3["present (Y)"] / T3.sum(1) * 100).round(1)
# 4. all-base result
T4 = pd.crosstab(D.allbase, D.call, margins=True)
# 5. per Data S10 deletion
T5 = (D[D.present].groupby("DataS10_deletion").agg(present_listings=("gene", "size"), genes=("gene", "nunique"))
      .join(D.groupby("DataS10_deletion").size().rename("listings")).sort_values("present_listings", ascending=False))

TESTS = pd.DataFrame([
    {"Test": "Fisher: share present, genes crossing an end vs genes inside", "P": float(f"{f1.pvalue:.2g}"), "Odds ratio": round(f1.statistic, 2)},
    {"Test": "Mann-Whitney: position inside (share from nearer end), present vs not present",
     "P": float(f"{mw.pvalue:.2g}"), "Median present": round(y.median(), 3), "Median not present": round(n.median(), 3), "N present": len(y), "N not present": len(n)},
    {"Test": "Kolmogorov-Smirnov: position of present genes inside vs uniform 0-0.5 (random)", "P": float(f"{ks_y.pvalue:.2g}"), "Mean present": round(y.mean(), 3)}])
for name, t in (("Ends vs inside", T1), ("Position inside", T2), ("Deletion size", T3), ("All-base result", T4), ("Tests", TESTS), ("Per deletion", T5.head(15))):
    print("==", name); print(t.to_string()); print()
with pd.ExcelWriter(OUTD + "present_genes_position.xlsx") as w:
    for name, t in (("Ends_vs_inside", T1), ("Position_inside", T2), ("Deletion_size", T3), ("Allbase_result", T4),
                    ("Tests", TESTS), ("Per_deletion", T5)):
        t.to_excel(w, sheet_name=name, index=name not in ("Tests",))
    D.drop(columns="present").to_excel(w, sheet_name="Listings", index=False)
print("written", OUTD + "present_genes_position.xlsx")
