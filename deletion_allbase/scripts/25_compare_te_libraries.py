#!/usr/bin/env python3
"""Side-by-side numbers: stand-in TE library (results/) versus combined EDTA library (results_edta/)."""
import json
import pandas as pd

A, B = "deletion_allbase/results/", "deletion_allbase/results_edta/"
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 30); pd.set_option("display.max_rows", 200)
sa, sb = json.load(open(A + "summary.json")), json.load(open(B + "summary.json"))
ea, eb = pd.read_csv(A + "allbase_events.tsv", sep="\t"), pd.read_csv(B + "allbase_events.tsv", sep="\t")
assert len(ea) == len(eb) and (ea.event_id == eb.event_id).all() and (ea.start == eb.start).all()
print("events:", len(ea), "(identical call set)")
SZ = ["100-500 bp", "0.5-1 kb", "1-5 kb", "5-20 kb", ">20 kb"]
for lab, e in (("stand-in", ea), ("EDTA", eb)):
    e["te"] = e.te_class != "no TE evidence"; e["big"] = e.deleted_bases >= 1000
    print(f"\n== {lab}: TE at a breakpoint: all {e.te.mean():.3f}; >=1kb {e[e.big].te.mean():.3f} ({e[e.big].te.sum()}/{e.big.sum()}); "
          f"same family both ends {(e.te_class == 'TE at both breakpoints, same named family').sum()}")
    print(e.te_class.value_counts().to_dict()); print(e.te_evidence.value_counts().to_dict())
    print(e.groupby("size_class").te.agg(["size", "mean"]).reindex(SZ).round(3).T.to_string())
    print("found/missed by SNP method:", e.groupby("found_by_snp_method").te.mean().round(3).to_dict())
    print(e.groupby(["size_class", "found_by_snp_method"]).te.agg(["size", "mean"]).round(3).unstack().reindex(SZ).to_string())
print("\nchanged class:"); print(pd.crosstab(ea.te_class, eb.te_class))
for k in ("te_enrichment_reference_annotation", "genes"):
    print("\n", k); print(" stand-in:", sa[k]); print(" EDTA:    ", sb[k])
fa = pd.Series(sa["te_family_counts"]); fb = pd.Series(sb["te_family_counts"])
print("\nfamilies at breakpoints (events), EDTA top 25:"); print(fb.head(25).to_string())
print("stand-in:", fa.to_dict())
d = sb["DataS10"]; d0 = sa["DataS10"]
for k in ("te_rows", "te_rows_with_allbase_te", "te_rows_same_family", "no_te_rows", "no_te_rows_with_allbase_te"):
    print(k, d0[k], "->", d[k])
xa, xb = pd.read_csv(A + "DataS10_crosswalk.tsv", sep="\t"), pd.read_csv(B + "DataS10_crosswalk.tsv", sep="\t")
print("\nData S10 TE label vs EDTA family found at the overlapping all-base event(s):")
print(xb[xb.s10_te != "No_TE"].groupby("s10_te").allbase_te_family.agg(lambda s: s.fillna("-").value_counts().head(4).to_dict()).to_string())
ra, rb = json.load(open(A + "read_evidence_summary.json")), json.load(open(B + "read_evidence_summary.json"))


def flat(d, p=""):
    for k, v in d.items():
        if isinstance(v, dict):
            yield from flat(v, p + str(k) + ".")
        elif not isinstance(v, list):
            yield p + str(k), v


fa_, fb_ = dict(flat(ra)), dict(flat(rb))
print("\nread_evidence_summary values that changed:")
for k in fa_:
    if k in fb_ and fa_[k] != fb_[k]:
        print(f"  {k}: {fa_[k]} -> {fb_[k]}")
print("keys only in one:", sorted(set(fa_) ^ set(fb_))[:20])
s9a, s9b = pd.read_csv(A + "DataS9_crosswalk.tsv", sep="\t"), pd.read_csv(B + "DataS9_crosswalk.tsv", sep="\t")
print("\nData S9 crosswalk (EDTA):"); print(s9b.to_string())
