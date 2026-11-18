#!/usr/bin/env python3
"""Sequencing-coverage context of every strain-level deletion region.

Some stretches of a genome are simply hard to sequence (for example very GC-rich or very
AT-rich DNA), so every library has few reads there. Such a stretch can pass the "callable"
filter in E2 and still fall below the deletion cut-off in a strain with a weaker library.
To recognise these, the following is recorded for every region:
  gc                      GC content of the interval in E2
  E2_rel_depth            E2's own unique-read depth in the interval, relative to E2's median
  nearref_rel_depth_min   the lowest relative depth among E5, E6, E19 and E36 (almost identical to E2)
  flank_rel_depth_L / _R  the carrier strain's relative depth in the 400 bp outside each edge
Also writes E2's relative depth by GC content in 200-bp windows (the coverage-vs-GC curve).
"""
import sys
from multiprocessing import Pool
import numpy as np
import pandas as pd

CHROMS = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07"]
REFLIKE = ["E2", "E5", "E6", "E19", "E36"]
FLANK = 400
sizes = {l.split()[0]: int(l.split()[1]) for l in open("ref/E2.chr.genome")}
reg = pd.read_csv("inputs/regions_strain_level.tsv", sep="\t")
reg["id"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
scor = {c: np.zeros(sizes[c], bool) for c in CHROMS}
for l in open("calls_uniq/scorable.bwa.uniq.bed"):
    c, s, e = l.split("\t")[:3]
    scor[c][int(s):int(e)] = True


def load_depth(s):
    df = pd.read_csv(f"depth/{s}.bwa.q20.bga.gz", sep="\t", header=None, names=["c", "s", "e", "d"],
                     dtype={"c": "category", "s": np.int32, "e": np.int32, "d": np.int32})
    out = {}
    for c in CHROMS:
        sub = df[df.c == c]
        out[c] = np.repeat(np.minimum(sub.d.values, 65535).astype(np.uint16), (sub.e - sub.s).values)
    med = float(np.median(np.concatenate([out[c][scor[c]] for c in CHROMS])))
    return out, med


def scor_mean(d, c, s, e):
    s, e = max(0, s), min(sizes[c], e)
    m = scor[c][s:e]
    return float(d[c][s:e][m].mean()) if m.any() else np.nan


def reflike(s):
    d, med = load_depth(s)
    iv = reg[["chrom", "start", "end"]].drop_duplicates()
    iv[f"{s}_rel_depth"] = [scor_mean(d, c, a, b) / med for c, a, b in iv.itertuples(index=False)]
    curve = None
    if s == "E2":
        rows = []
        for c in CHROMS:
            n = sizes[c] // 200
            rows.append(pd.DataFrame({"gc": GC[c][:n * 200].reshape(n, 200).mean(1), "rel": d[c][:n * 200].reshape(n, 200).mean(1) / med,
                                      "ok": scor[c][:n * 200].reshape(n, 200).mean(1) >= 0.9}))
        w = pd.concat(rows); w = w[w.ok]
        curve = w.groupby((w.gc * 20).round() / 20).rel.agg(["size", "mean", "median", lambda x: (x < 0.25).mean()]).reset_index()
        curve.columns = ["gc_bin", "windows", "mean_rel_depth", "median_rel_depth", "share_below_25pct"]
    return iv, curve


def carrier(s):
    d, med = load_depth(s)
    sub = reg[reg.strain == s]
    return pd.DataFrame({"strain": s, "id": sub.id.values,
                         "flank_rel_depth_L": [scor_mean(d, c, a - FLANK, a) / med for c, a in zip(sub.chrom, sub.outer_start)],
                         "flank_rel_depth_R": [scor_mean(d, c, b, b + FLANK) / med for c, b in zip(sub.chrom, sub.outer_end)]})


# GC as a per-base 0/1 array
GC = {}
name = None
for line in open("ref/E2.fasta"):
    if line.startswith(">"):
        name = line[1:].split()[0]; GC[name] = []
    else:
        GC[name].append(line.strip())
GC = {c: np.isin(np.frombuffer("".join(GC[c]).upper().encode(), dtype="S1"), [b"G", b"C"]) for c in CHROMS}

if __name__ == "__main__":
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool:
        ref = pool.map(reflike, REFLIKE, chunksize=1)
        car = pool.map(carrier, sorted(reg.strain.unique()), chunksize=1)
    out = reg[["strain", "chrom", "start", "end", "id"]].copy()
    out["gc"] = [float(GC[c][a:b].mean()) for c, a, b in zip(out.chrom, out.start, out.end)]
    for iv, _ in ref:
        out = out.merge(iv, on=["chrom", "start", "end"])
    out["nearref_rel_depth_min"] = out[[f"{s}_rel_depth" for s in REFLIKE[1:]]].min(axis=1)
    out = out.merge(pd.concat(car), on=["strain", "id"])
    out.round(4).to_csv("export_uniq/region_context.tsv", sep="\t", index=False)
    ref[0][1].round(4).to_csv("export_uniq/E2_depth_by_gc.tsv", sep="\t", index=False)
    print(ref[0][1].to_string())
