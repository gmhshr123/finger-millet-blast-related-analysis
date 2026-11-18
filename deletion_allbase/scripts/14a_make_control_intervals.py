#!/usr/bin/env python3
"""Matched control intervals for the read-evidence checks (run locally).

For every depth-based deletion call, one "fake" interval in the same strain, with the
same length and the same breakpoint-zone widths, at a random position where that strain has
no deletion call. The read-evidence scripts are then run on the real and on the fake
intervals in exactly the same way, so the fake ones show what the evidence looks like when
nothing is deleted.

A control is first tried within 100 kb of the real call (same local repeat / SNP context);
if no free place is found there it is put anywhere on the same chromosome.
"""
import numpy as np
import pandas as pd

C = "deletion_allbase/results/cluster_uniq/"
rng = np.random.default_rng(20261008)
CHROMS = [f"Chr0{i}" for i in range(1, 8)]
sizes = {l.split()[0]: int(l.split()[1]) for l in open(C + "E2.chr.genome")}
scor = {c: np.zeros(sizes[c], bool) for c in CHROMS}
for l in open(C + "scorable.bwa.uniq.bed"):
    c, s, e = l.split("\t")[:3]
    scor[c][int(s):int(e)] = True

calls = pd.read_csv(C + "depth_calls.bwa.uniq.primary.tsv", sep="\t")
FLANK = 1500          # a control (with its zones) must stay this far from any real call of the strain
out = []
for (s, c), sub in calls.groupby(["strain", "chrom"]):
    busy = np.zeros(sizes[c], bool)
    for r in sub.itertuples():
        busy[max(0, r.outer_start - FLANK):r.outer_end + FLANK] = True
    bcs = np.concatenate(([0], np.cumsum(busy)))
    for r in sub.itertuples():
        zl, zr, span = r.start - r.outer_start, r.outer_end - r.end, r.end - r.start
        lo, hi = zl + 2000, sizes[c] - span - zr - 2000
        placed = None
        for where, n_try in (("within 100 kb", 300), ("same chromosome", 3000)):
            for _ in range(n_try):
                st = int(rng.integers(max(lo, r.start - 100_000), min(hi, r.start + 100_000))) if where == "within 100 kb" \
                    else int(rng.integers(lo, hi))
                en = st + span
                if not (scor[c][st] and scor[c][en - 1]):
                    continue
                if bcs[en + zr] - bcs[st - zl] > 0:
                    continue
                placed = (st, en, where)
                break
            if placed:
                break
        if placed is None:
            continue
        st, en, where = placed
        busy[st - zl:en + zr] = True                      # controls of one strain do not overlap each other
        bcs = np.concatenate(([0], np.cumsum(busy)))
        out.append((s, c, st, en, st - zl, en + zr, span, f"{c}:{r.start}-{r.end}", where))
ctrl = pd.DataFrame(out, columns=["strain", "chrom", "start", "end", "outer_start", "outer_end", "span", "real_id", "placement"])
ctrl.to_csv("deletion_allbase/inputs/control_intervals.tsv", sep="\t", index=False)
print(len(calls), "calls;", len(ctrl), "controls placed")
print(ctrl.placement.value_counts().to_string())
