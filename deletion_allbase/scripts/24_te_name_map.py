#!/usr/bin/env python3
"""Give every entry of the combined EDTA TE library a family label.

EDTA names its families with numbers (for example E2__TE_00000012_LTR). Where an entry matches
a TE with a published name (Pot2, Grasshopper, Pyret, ...), it gets that name, so results can be
read and compared with the manuscript. All other entries keep their EDTA number, written as
<genome>.TE<number>; the _LTR and _INT parts of one LTR retrotransposon share one label.

An entry gets a published name when its best BLASTN hit against the named sequences has at least
80% identity over at least 200 bp, or over at least half of the entry if the entry is shorter.

Usage: 24_te_name_map.py library.fa blast_outfmt6.tsv out.tsv
(blast columns: qseqid sseqid pident length qlen slen evalue bitscore; query = EDTA library)
"""
import re, sys
import pandas as pd

lib, blast, out = sys.argv[1:4]
rows = []
name, n = None, 0
for l in open(lib):
    if l.startswith(">"):
        if name:
            rows.append((name, n))
        name, n = l[1:].strip(), 0
    else:
        n += len(l.strip())
rows.append((name, n))
L = pd.DataFrame(rows, columns=["header", "lib_len"])
L["lib_id"] = L.header.str.split("#").str[0]
L["te_class"] = L.header.str.split("#").str[1].fillna("")


def edta_label(i):
    m = re.match(r"(.+?)__TE_(\d+)(?:_(LTR|INT))?$", i)
    return f"{m.group(1)}.TE{m.group(2)}" if m else i


L["edta_label"] = L.lib_id.map(edta_label)
b = pd.read_csv(blast, sep="\t", header=None, names="header sseqid pident length qlen slen evalue bitscore".split())
b["known"] = b.sseqid.str.replace(r"(_\d+)?#.*", "", regex=True)
# MAGGY (GenBank L35053) and the Fosbury elements of this project are 98.4% identical over their
# full length, i.e. one family; the library holds it as a single entry.
b["known"] = b.known.replace({"MAGGY": "MAGGY_Fosbury", "Fosbury": "MAGGY_Fosbury"})
b = b.sort_values("bitscore", ascending=False).drop_duplicates("header")
b = b[(b.pident >= 80) & (b.length >= (b.qlen * 0.5).clip(upper=200))]
L = L.merge(b[["header", "known", "pident", "length"]], on="header", how="left")
# the two parts of one LTR element get the same name if either part matched
fam_known = L.dropna(subset=["known"]).groupby("edta_label").known.first()
L["label"] = L.edta_label.map(fam_known).fillna(L.edta_label)
L[["lib_id", "label", "te_class", "lib_len", "known", "pident", "length"]].to_csv(out, sep="\t", index=False)
print(len(L), "library entries;", L.known.notna().sum(), "match a published TE;", L.label.nunique(), "family labels")
print(L[L.known.notna()].groupby("label").agg(entries=("lib_id", "size"), ids=("lib_id", lambda x: ",".join(x)[:120])).to_string())
found = set(L.known.dropna())
print("published TEs with no entry in the library:", sorted({"MAGGY_Fosbury" if x in ("MAGGY", "Fosbury") else x for x in open(sys.argv[4]).read().split()} - found) if len(sys.argv) > 4 else "")
