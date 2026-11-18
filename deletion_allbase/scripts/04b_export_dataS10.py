#!/usr/bin/env python3
"""Export Data S10 to plain tables and test whether the SNP-gap runs reproduce it."""
import collections, warnings
import openpyxl
import pandas as pd
warnings.filterwarnings("ignore")

wb = openpyxl.load_workbook("DataS10_Deletions_WholeGenomeLevel.xlsx")   # full mode: font colour marks effectors
ws = wb.worksheets[0]
vals = openpyxl.load_workbook("DataS10_Deletions_WholeGenomeLevel.xlsx", data_only=True).worksheets[0]   # cached formula results
rows, coloured, allg = [], set(), set()
for r in ws.iter_rows(min_row=4):
    if not (r[3].value and str(r[3].value).startswith("Del")):
        continue
    genes = []
    for c in r[12:]:
        if c.value:
            g = str(c.value).strip(); genes.append(g); allg.add(g)
            col = c.font.color
            tag = None
            if col is not None:
                tag = col.rgb if col.type == "rgb" else f"theme{col.theme}"
            if tag not in (None, "FF000000", "theme1"):
                coloured.add((g, tag))
    rows.append(dict(del_id=r[3].value, strains=str(r[0].value).replace(" ", ""), te=r[1].value,
                     reads=str(r[2].value).replace(" ", ""), dist_above=int(r[4].value), start=int(r[5].value),
                     end=int(r[6].value), dist_below=int(r[7].value), min_size=int(vals[r[8].coordinate].value),
                     max_size=int(vals[r[9].coordinate].value), snp_count=int(r[10].value), chrom=r[11].value, genes=",".join(genes)))
s10 = pd.DataFrame(rows)
s10.to_csv("deletion_allbase/inputs/DataS10_events.tsv", sep="\t", index=False)
print(len(s10), "events; unique genes", len(allg))
print("font colours on genes:", collections.Counter(c for g, c in coloured))
pd.DataFrame(sorted(coloured), columns=["gene", "font"]).to_csv(
    "deletion_allbase/inputs/DataS10_coloured_genes.tsv", sep="\t", index=False)

iv = s10.assign(id=s10.del_id, s0=s10.start - 1)[["id", "chrom", "s0", "end"]].rename(columns={"s0": "start"})
iv.to_csv("deletion_allbase/inputs/DataS10.intervals.tsv", sep="\t", index=False)

g = pd.read_csv("deletion_allbase/results/snp_gap_runs.tsv", sep="\t")
u = g.drop_duplicates(["chrom", "start", "end"]).copy()
u["id"] = u.chrom + ":" + u.start.astype(str) + "-" + u.end.astype(str)
u.assign(s0=u.start - 1)[["id", "chrom", "s0", "end"]].rename(columns={"s0": "start"}).to_csv(
    "deletion_allbase/inputs/snp_gap_runs.intervals.tsv", sep="\t", index=False)
print("unique SNP-gap intervals:", len(u))

key = set(zip(g.strain, g.chrom, g.start, g.end))
hit = tot = 0; bad = []
for _, r in s10.iterrows():
    for st in r.strains.split(","):
        tot += 1
        if (st, r.chrom, r.start, r.end) in key:
            hit += 1
        else:
            bad.append((r.del_id, st))
print(f"S10 strain-events reproduced exactly by SNP-gap runs: {hit}/{tot}")
print("not reproduced (first 15):", bad[:15])
m = s10.merge(u[["chrom", "start", "end", "snp_count", "dist_above", "dist_below"]],
              on=["chrom", "start", "end"], how="left", suffixes=("", "_mine"))
print("events with identical interval:", int(m.snp_count_mine.notna().sum()),
      "; same SNP count:", int((m.snp_count == m.snp_count_mine).sum()),
      "; same flank distances:", int(((m.dist_above == m.dist_above_mine) & (m.dist_below == m.dist_below_mine)).sum()))
print("S10 SNP-count distribution (smallest):", s10.snp_count.value_counts().sort_index().head(6).to_dict())
