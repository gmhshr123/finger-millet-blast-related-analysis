#!/usr/bin/env python3
"""List the section headings of Data S2 and export gene IDs with their section (effector / non-effector)."""
import collections, re, warnings
import openpyxl
import pandas as pd
warnings.filterwarnings("ignore")
ws = openpyxl.load_workbook("DataS2_AnnotatedGenes_PresenceAbsence_fixed.xlsx", read_only=True, data_only=True).worksheets[0]
section, rows = None, []
for i, r in enumerate(ws.iter_rows(max_col=6, values_only=True)):
    a = r[0]
    if a is None:
        continue
    a = str(a)
    if not re.match(r"EcMO", a):
        if re.match(r"(Non-)?Effector [Gg]enes", a):
            section = a.strip()
            continue
        if section is None or i > 14084:      # header rows / footnotes
            continue
    rows.append((a, r[1], r[2], r[3], r[4], section))
df = pd.DataFrame(rows, columns=["mrna", "gene", "col_c", "col_d", "mrna_len", "section"])
print(df.groupby("section").size())
print(df.gene.astype(str).str.extract(r"^(EcMO\d?g?)")[0].value_counts())
print(df.head(3).to_string()); print(df[df.gene.astype(str).str.startswith("EcMO1g")].head(3).to_string())
df.to_csv("deletion_allbase/inputs/DataS2_gene_sections.tsv", sep="\t", index=False)
