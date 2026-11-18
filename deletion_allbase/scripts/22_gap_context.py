#!/usr/bin/env python3
"""Is the "present" DNA between two neighbouring deletion calls real, or a repeat artefact?

If a strain has lost a stretch of DNA and E2 has a repeat inside that stretch, reads from other
copies of the repeat still map there. The middle then looks present and one deletion is split
in two. Annotated repeats are already bridged by the caller (03). This script looks at what is
left: for every pair of neighbouring calls of one strain (less than 20 kb apart) it describes
the stretch between them (the "gap"):

  scorable_bases           unique, callable, non-repeat bases in the gap (the only ones that count as evidence)
  repeat_frac              share of the gap that is annotated repeat
  present_bases            scorable bases with unique depth above the deletion cut-off
  rel_depth_unique         the strain's unique-read depth on the scorable bases, relative to its median
  unique_share             unique reads / all reads on the scorable bases (low = reads there are multi-mapping)
  rel_depth_all_gap        the strain's all-read depth over the whole gap, repeats included
  E2_rel_depth_unique      the same for E2's own reads (low = hard to map uniquely even in E2)
  E2_rel_depth_all         E2 all-read depth relative to its median (high = more copies than the assembly shows)
  nearref_rel_depth_min    lowest relative unique depth among E5, E6, E19, E36

Second table: for every strain-level region, the repeats that the call was bridged over, and how
many reads map inside them (a situation the caller ignores by design).
"""
import sys
from multiprocessing import Pool
import numpy as np
import pandas as pd

CHROMS = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07"]
REFLIKE = ["E2", "E5", "E6", "E19", "E36"]
MAXGAP = 20_000
sizes = {l.split()[0]: int(l.split()[1]) for l in open("ref/E2.chr.genome")}


def mask(bed):
    m = {c: np.zeros(sizes[c], bool) for c in CHROMS}
    for l in open(bed):
        f = l.split("\t")
        if f[0] in m:
            m[f[0]][int(f[1]):int(f[2])] = True
    return m


scor = mask("calls_uniq/scorable.bwa.uniq.bed")
rep = mask("ref/E2.all_repeats.bed")
frag = pd.read_csv("calls_uniq/depth_calls.bwa.uniq.primary.tsv", sep="\t").sort_values(["strain", "chrom", "start"])
cohort = [l.strip() for l in open("meta/cohort33.txt")]
frag = frag[frag.strain.isin(cohort)]
nx = frag.groupby(["strain", "chrom"]).shift(-1)
gaps = pd.DataFrame({"strain": frag.strain, "chrom": frag.chrom, "gap_start": frag.end, "gap_end": nx.start,
                     "left_id": frag.chrom + ":" + frag.start.astype(str) + "-" + frag.end.astype(str)}).dropna()
gaps["gap_end"] = gaps.gap_end.astype(int)
gaps = gaps[(gaps.gap_end - gaps.gap_start <= MAXGAP) & (gaps.gap_end > gaps.gap_start)].reset_index(drop=True)
reg = pd.read_csv("inputs/regions_strain_level.tsv", sep="\t")


def load(s, q):
    df = pd.read_csv(f"depth/{s}.bwa.{q}.bga.gz", sep="\t", header=None, names=["c", "s", "e", "d"],
                     dtype={"c": "category", "s": np.int32, "e": np.int32, "d": np.int32})
    out = {}
    for c in CHROMS:
        sub = df[df.c == c]
        out[c] = np.repeat(np.minimum(sub.d.values, 65535).astype(np.uint16), (sub.e - sub.s).values)
    return out


def med(d):
    return float(np.median(np.concatenate([d[c][scor[c]] for c in CHROMS])))


def reflike(s):
    q20, q0 = load(s, "q20"), load(s, "q0")
    m20, m0 = med(q20), med(q0)
    iv = gaps[["chrom", "gap_start", "gap_end"]].drop_duplicates()
    u, a = [], []
    for c, x, y in iv.itertuples(index=False):
        k = scor[c][x:y]
        u.append(q20[c][x:y][k].mean() / m20 if k.any() else np.nan)
        a.append(q0[c][x:y].mean() / m0)
    iv[f"{s}_rel_depth_unique"] = u; iv[f"{s}_rel_depth_all"] = a
    return iv


def carrier(s):
    q20, q0 = load(s, "q20"), load(s, "q0")
    m20 = med(q20)
    rows = []
    for r in gaps[gaps.strain == s].itertuples():
        c, x, y = r.chrom, r.gap_start, r.gap_end
        k = scor[c][x:y]
        d20, d0 = q20[c][x:y], q0[c][x:y]
        rows.append((s, c, x, y, y - x, int(k.sum()), float(rep[c][x:y].mean()), int((d20[k] > 0.05 * m20).sum()),
                     d20[k].mean() / m20 if k.any() else np.nan, d20[k].sum() / max(d0[k].sum(), 1) if k.any() else np.nan,
                     d0.mean() / m20))
    g = pd.DataFrame(rows, columns=["strain", "chrom", "gap_start", "gap_end", "gap_len", "scorable_bases", "repeat_frac", "present_bases",
                                    "rel_depth_unique", "unique_share", "rel_depth_all_gap"])
    rows = []
    for r in reg[reg.strain == s].itertuples():
        c, x, y = r.chrom, r.start, r.end
        k = ~scor[c][x:y]                                   # bases the call was bridged over
        kr = k & rep[c][x:y]
        d0 = q0[c][x:y]
        rows.append((s, c, x, y, int(k.sum()), int(kr.sum()), d0[kr].mean() / m20 if kr.any() else np.nan,
                     int((d0[kr] >= 0.5 * m20).sum()), int((q20[c][x:y][kr] >= 0.5 * m20).sum())))
    b = pd.DataFrame(rows, columns=["strain", "chrom", "start", "end", "bridged_bases", "bridged_repeat_bases", "repeat_rel_depth_all",
                                    "repeat_bases_with_reads", "repeat_bases_with_unique_reads"])
    return g, b


if __name__ == "__main__":
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool:
        ref = pool.map(reflike, REFLIKE, chunksize=1)
        car = pool.map(carrier, cohort, chunksize=1)
    g = pd.concat([x[0] for x in car])
    for iv in ref:
        g = g.merge(iv, on=["chrom", "gap_start", "gap_end"])
    g["nearref_rel_depth_min"] = g[[f"{s}_rel_depth_unique" for s in REFLIKE[1:]]].min(axis=1)
    g.round(4).to_csv("export_uniq/gap_context.tsv", sep="\t", index=False)
    pd.concat([x[1] for x in car]).round(4).to_csv("export_uniq/region_bridged_repeats.tsv", sep="\t", index=False)
    print(len(g), "gaps written")
