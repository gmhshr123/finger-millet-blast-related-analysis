#!/usr/bin/env python3
"""Export the full-length TE copies listed in the two workbooks to one plain table (run from manuscripts_all/).

Sources
  Full_length_TE_copies_B71_E2_MZ_P131.xlsx   LTR retrotransposons; sheet "(2)" = all full-length copies,
                                              sheet "(3)" = families with 10 or more copies in a genome
  Full_length_Pot2_copies_B71_E2_MZ_P131.xlsx Pot2; sheet "(2)" = all RepeatMasker hits, sheet "(3)" = full-length hits
Output: rip_analysis/inputs/full_length_copies.tsv
  genome, family, chrom, start, end (1-based, inclusive), strand ('.' = not given), tsd, selected (in sheet 3)
"""
import warnings
import openpyxl
import pandas as pd
warnings.filterwarnings("ignore")


def rows_of(path, sheet_tag):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = [w for w in wb.worksheets if sheet_tag in w.title][0]
    return [r for r in ws.iter_rows(values_only=True) if any(c is not None for c in r)]


out = []
sel = set()
for tag in ("(3)", "(2)"):
    for r in rows_of("Full_length_TE_copies_B71_E2_MZ_P131.xlsx", tag)[1:]:
        if r[11] is None:
            continue                                            # note row
        g, c, span = str(r[6]).split(":")
        s, e = (int(x) for x in span.split(".."))
        key = (g, c, s, e)
        if tag == "(3)":
            sel.add(key)
        else:
            cls = str(r[14]).split("#")[1] if "#" in str(r[14]) else ""
            out.append(dict(genome=g, family=f"LTR{int(r[12]):02d}", te_class=cls, chrom=c, start=s, end=e, strand=".",
                            tsd=r[7] or "", ltr_identity=r[9], selected=key in sel))
for tag in ("(3)",):
    for r in rows_of("Full_length_Pot2_copies_B71_E2_MZ_P131.xlsx", tag)[2:]:
        g, c = str(r[4]).split(":")
        out.append(dict(genome=g, family="Pot2", te_class="DNA/Pot2", chrom=c, start=int(r[5]), end=int(r[6]),
                        strand="-" if r[8] == "C" else "+", tsd="", ltr_identity="", selected=True))
d = pd.DataFrame(out)
d["length"] = d.end - d.start + 1
d["copy_id"] = d.genome + ":" + d.chrom + ":" + d.start.astype(str) + ".." + d.end.astype(str)
assert not d.copy_id.duplicated().any()
d.to_csv("rip_analysis/inputs/full_length_copies.tsv", sep="\t", index=False)
t = d[d.selected].pivot_table(index="family", columns="genome", values="length", aggfunc="size", fill_value=0)
print(t.to_string()); print("all copies:", len(d), "| selected:", int(d.selected.sum()))
print("selected (family x genome) sets with >= 9 copies:", int((t.values >= 9).sum()))
