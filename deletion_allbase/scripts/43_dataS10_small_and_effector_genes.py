#!/usr/bin/env python3
"""Data S10 genes: which are smaller than 200 bp, and which are predicted effectors, by TE association of the deletion.
Run from manuscripts_all/. Reads the current Data S10 workbook directly (blue font = effector, footnote 7).

Length: mRNA length from Data S2 (fixed), and the gene span on E2 from E2.genes.bed.
Output: printed; table written by the caller only after the numbers are checked.
"""
import collections, warnings
import openpyxl
import pandas as pd
warnings.filterwarnings("ignore")

X = "DataS10_Deletions_WholeGenomeLevel.xlsx"
ws = openpyxl.load_workbook(X).worksheets[0]
rows = []
fonts = collections.Counter()
for r in ws.iter_rows(min_row=4):
    if not (r[3].value and str(r[3].value).startswith("Del")):
        continue
    for c in r[12:]:
        if c.value:
            col = c.font.color
            tag = None if col is None else (col.rgb if col.type == "rgb" else f"theme{col.theme}")
            blue = tag not in (None, "FF000000", "theme1")
            fonts[tag] += 1
            rows.append(dict(gene=str(c.value).strip(), del_id=r[3].value, te=r[1].value, strains=str(r[0].value).replace(" ", ""), blue=blue))
L = pd.DataFrame(rows)
print("font colours:", fonts)
print(len(L), "gene listings;", L.gene.nunique(), "unique genes; deletions with genes:", L.del_id.nunique())

s2 = pd.read_csv("deletion_allbase/inputs/DataS2_gene_sections.tsv", sep="\t")
s2 = s2[s2.gene.astype(str).str.startswith("EcMO")]
ml = s2.drop_duplicates("gene").set_index("gene").mrna_len
eff2 = set(s2[s2.section.astype(str).str.startswith("Effector")].gene)
bed = {}
for l in open("deletion_allbase/results_edta/cluster_uniq/E2.genes.bed", "rb").read().decode().split("\n"):
    f = l.replace("\r", "").split("\t")
    if len(f) >= 4:
        bed[f[3]] = (f[0], int(f[1]) + 1, int(f[2]), f[5] if len(f) > 5 else "")
G = L.groupby("gene").agg(blue=("blue", "any"), deletions=("del_id", lambda s: ",".join(dict.fromkeys(s))),
                          te=("te", lambda s: ",".join(dict.fromkeys(s))), noTE=("te", lambda s: (s == "No_TE").any()),
                          withTE=("te", lambda s: (s != "No_TE").any())).reset_index()
G["mrna_len"] = G.gene.map(ml)
G["span"] = G.gene.map(lambda g: bed[g][2] - bed[g][1] + 1 if g in bed else None)
G["in_S2"] = G.gene.isin(ml.index)
G["effector_S2"] = G.gene.isin(eff2)
print("in Data S2:", int(G.in_S2.sum()), "not in Data S2:", int((~G.in_S2).sum()), "; blue:", int(G.blue.sum()), "; effector in S2:", int(G.effector_S2.sum()))
print(G[~G.in_S2][["gene", "blue", "span", "deletions", "te"]].to_string())
print("mRNA < 200:", int((G.mrna_len < 200).sum()), " span < 200:", int((G.span < 200).sum()), " span missing:", int(G.span.isna().sum()))
print(G[G.blue != G.effector_S2][["gene", "blue", "effector_S2", "mrna_len", "span", "deletions"]].to_string())
print(G.groupby(["noTE", "withTE"]).agg(n=("gene", "size"), blue=("blue", "sum")))
G.to_csv(__import__("tempfile").gettempdir() + "/s10genes.tsv", sep="\t", index=False)

# ---- mRNA length from the exons of the E2 annotation (inputs/E2.gff3.gz); the workbook
import gzip, re, os
from openpyxl.styles import Alignment, Font, PatternFill
ex, info = collections.Counter(), {}
for l in gzip.open("deletion_allbase/inputs/E2.gff3.gz", "rt"):
    f = l.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    if f[2] == "gene":
        gid = re.search(r"ID=([^;]+)", f[8]).group(1)
        d = re.search(r"description=([^;]*)", f[8])
        info[gid] = (f[0], int(f[3]), int(f[4]), f[6], d.group(1) if d else "")
    elif f[2] == "exon":
        p = re.search(r"Parent=([^;,]+)", f[8]).group(1)
        if p.endswith("-mRNA-1"):
            ex[p[:-7]] += int(f[4]) - int(f[3]) + 1
G["mrna_gff"] = G.gene.map(ex)
print("mRNA length, annotation versus Data S2, differences:", int(((G.mrna_gff != G.mrna_len) & G.in_S2).sum()))
print("mRNA < 200 bp by annotation:", int((G.mrna_gff < 200).sum()), "; all of them absent from Data S2:", bool(set(G[G.mrna_gff < 200].gene) == set(G[~G.in_S2].gene)))
print(G[G.mrna_gff < 200][["gene", "mrna_gff", "span", "blue", "deletions", "te"]].to_string())

out = pd.DataFrame({"Gene": G.gene, "Chromosome": G.gene.map(lambda g: info[g][0]), "Start": G.gene.map(lambda g: info[g][1]), "End": G.gene.map(lambda g: info[g][2]),
                    "Strand": G.gene.map(lambda g: info[g][3]), "mRNA length (bp)": G.mrna_gff.astype(int), "Gene length on E2 (bp)": G.span.astype(int),
                    "Predicted effector": G.blue.map({True: "yes", False: "no"}), "Data S10 deletions": G.deletions, "TE label in Data S10": G.te,
                    "In Data S2": G.in_S2.map({True: "yes", False: "no"}), "Description": G.gene.map(lambda g: info[g][4])})
small, eff_ = G.mrna_gff < 200, G.blue
sheets = {"Under_200bp_no_TE": out[small & G.noTE], "Under_200bp_with_TE": out[small & G.withTE],
          "Effectors_no_TE": out[eff_ & G.noTE], "Effectors_with_TE": out[eff_ & G.withTE]}
summ = pd.DataFrame([
    ("Unique genes in Data S10", len(G)), ("mRNA shorter than 200 bp", int(small.sum())), ("Same, in a deletion without a TE", int((small & G.noTE).sum())),
    ("Same, in a deletion with a TE", int((small & G.withTE).sum())), ("Same, listed under both kinds of deletion", int((small & G.noTE & G.withTE).sum())),
    ("Genes of 200 bp or longer (as in the manuscript)", int((~small).sum())), ("Predicted effectors (blue in Data S10)", int(eff_.sum())),
    ("Same, in a deletion without a TE", int((eff_ & G.noTE).sum())), ("Same, in a deletion with a TE", int((eff_ & G.withTE).sum())),
    ("Same, listed under both kinds of deletion", int((eff_ & G.noTE & G.withTE).sum())), ("Effectors of 200 bp or longer (as in the manuscript)", int((eff_ & ~small).sum())),
    ("Effectors shorter than 200 bp", int((eff_ & small).sum()))], columns=["Item", "Genes"])
P = "deletion_allbase/supplement_figures_and_NoteS5/DataS10_small_genes_and_effectors.xlsx"
os.makedirs(os.path.dirname(P), exist_ok=True)
with pd.ExcelWriter(P, engine="openpyxl") as w:
    summ.to_excel(w, sheet_name="Summary", index=False)
    for k, v in sheets.items():
        v.sort_values(["Chromosome", "Start"]).to_excel(w, sheet_name=k, index=False)
wb = openpyxl.load_workbook(P)
for s in wb.worksheets:
    for row in s.iter_rows():
        for c in row:
            c.font = Font(name="Arial", size=10, bold=c.row == 1)
            if c.row == 1:
                c.fill = PatternFill("solid", start_color="DDE5EE"); c.alignment = Alignment(wrap_text=True, vertical="top")
    for col in s.columns:
        s.column_dimensions[col[0].column_letter].width = min(max(len(str(c.value or "")) for c in col[1:] or col) + 3, 48)
    s.freeze_panes = "A2"
wb.save(P)
print(summ.to_string(index=False))
print("written", P)
