#!/usr/bin/env python3
"""Crosswalk of the 14 manually inspected breakpoints in Data S9 against the all-base calls.

Data S9 gives, for each example, a range for the proximal breakpoint and a range for the distal
breakpoint (each range lies between two neighbouring SNPs), the strain(s), and the TE family that
was seen at the breakpoint. For every example and strain, look up the all-base region(s) of
that strain between the two ranges and report:
  - where the all-base deletion starts and ends, and whether each edge falls inside the S9 range
    (the all-base edge may sit anywhere in its own breakpoint zone [outer_start, start] / [end, outer_end]);
  - the exact breakpoint position from clipped reads, when there is one;
  - the TE family seen in the clipped reads or annotated in E2 at the edges;
  - the read-evidence label from script 16.
E17 is not in the 33-strain cohort and was not processed.
Writes results/DataS9_crosswalk.tsv.
"""
import re
import openpyxl
import pandas as pd

R = __import__("os").environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")    # primary result set: combined EDTA TE library; set DEL_RESULTS=deletion_allbase/results/ for the earlier stand-in library
FAM = {"Pot2": "Pot2", "Grh": "Grasshopper", "Fos": "Fosbury", "Pyret": "Pyret"}
calls = pd.read_csv(R + "allbase_calls_strain_level_read_evidence.tsv", sep="\t")
ev = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t").set_index("event_id")
processed = set(calls.strain) | {"E2", "E5", "E6", "E19", "E36"}


def rng(txt):
    a, b = [int(x.replace(",", "")) for x in re.findall("[0-9][0-9,]*", str(txt))]
    return a, b


ws = openpyxl.load_workbook("DataS9_DeletionBP_FlankingSeq.xlsx", data_only=True).active
rows, chrom, k = [], None, 0
for r in ws.iter_rows(min_row=3, values_only=True):
    if r[0] and str(r[0]).startswith("Chr"):
        chrom = "Chr0" + str(r[0])[3:]
    if r[3] is None or not re.match("E[0-9]", str(r[3])) or str(r[3]).startswith("E2 ") or r[1] in (None, "-"):
        continue
    k += 1
    (p0, p1), (d0, d1) = rng(r[1]), rng(r[2])
    for s in str(r[3]).split("/"):
        base = dict(example=k, chrom=chrom, s9_proximal=f"{p0}-{p1}", s9_distal=f"{d0}-{d1}", strain=s, s9_te=r[4])
        if s not in processed:
            rows.append({**base, "result": "strain not processed (outside the 33-strain cohort)"})
            continue
        m = calls[(calls.strain == s) & (calls.chrom == chrom) & (calls.outer_end > p0) & (calls.outer_start < d1)]
        if m.empty:
            rows.append({**base, "result": "no all-base call here"})
            continue
        first, last = m.iloc[0], m.iloc[-1]
        # an edge agrees when its breakpoint zone overlaps the S9 range (20 bp tolerance)
        ok_p = first.outer_start <= p1 + 20 and first.start >= p0 - 20
        ok_d = last.end <= d1 + 20 and last.outer_end >= d0 - 20
        e = ev.loc[sorted(set(m.event_id))]
        fams = sorted({f for col in ("L_family", "R_family") for v in e[col].dropna() for f in str(v).split(",") if f})
        want = FAM.get(str(r[4]))
        rows.append({**base, "result": "found", "n_allbase_regions": len(m),
                     "allbase_start_zone": f"{first.outer_start}-{first.start}", "allbase_end_zone": f"{last.end}-{last.outer_end}",
                     "proximal_edge_agrees": ok_p, "distal_edge_agrees": ok_d,
                     "clipped_read_breakpoint_left": int(first.L_bp) if first.L_bp_reads >= 3 else None,
                     "clipped_read_breakpoint_right": int(last.R_bp) if last.R_bp_reads >= 3 else None,
                     "left_breakpoint_in_s9_range": bool(first.L_bp_reads >= 3 and p0 - 20 <= first.L_bp <= p1 + 20),
                     "right_breakpoint_in_s9_range": bool(last.R_bp_reads >= 3 and d0 - 20 <= last.R_bp <= d1 + 20),
                     "allbase_te_families": ",".join(fams), "te_class": ";".join(sorted(set(e.te_class))),
                     "same_family_as_s9": (want in [x.replace("MAGGY_Fosbury", "Fosbury") for x in fams]) if want else None,
                     "read_class": ";".join(sorted(set(m.read_class))), "events": ",".join(e.index)})
X = pd.DataFrame(rows)
X.to_csv(R + "DataS9_crosswalk.tsv", sep="\t", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 30)
print(X.drop(columns=["events"]).to_string())
f = X[X.result == "found"]
print("\nstrain-level examples:", len(X), "| processed:", int((X.result != "strain not processed (outside the 33-strain cohort)").sum()),
      "| found:", len(f), "| both edges agree:", int((f.proximal_edge_agrees & f.distal_edge_agrees).sum()),
      "| TE-labelled and same family:", int((f.same_family_as_s9 == True).sum()), "of", int(f.same_family_as_s9.notna().sum()))
