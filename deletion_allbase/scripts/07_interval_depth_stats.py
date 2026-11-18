#!/usr/bin/env python3
"""Score any set of genome intervals in every strain from the all-base depth tracks.

For each interval x strain: number of scorable bases, number of those that are deleted
(unique-read depth <= the primary cut-off of the caller, 5% of the strain median), and mean depth relative to the strain
median for unique reads and for all reads. This is used to check SNP-based candidates and
Data S10 rows against the all-base data, and to genotype every depth call in every strain.

Intervals file: tab-separated with header, columns id, chrom, start, end (0-based, half-open).
"""
import argparse, sys
from multiprocessing import Pool
import numpy as np
import pandas as pd
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from importlib import import_module
caller = import_module("03_call_depth_deletions")
CHROMS, load_depth, median_depth = caller.CHROMS, caller.load_depth, caller.median_depth
THR = caller.PRIMARY_THR


def score(args):
    s, a = args
    q20 = load_depth(f"{a['depth']}/{s}.{a['aligner']}.q20.bga.gz", a["sizes"])
    q0 = load_depth(f"{a['depth']}/{s}.{a['aligner']}.q0.bga.gz", a["sizes"])
    med = median_depth(q20, a["callable"])
    out = []
    for c, sub in a["iv"].groupby("chrom"):
        if c not in CHROMS:
            continue
        cal = a["callable"][c]
        cs = lambda x: np.concatenate(([0], np.cumsum(x, dtype=np.int64)))
        c_cal, c_del = cs(cal), cs(cal & (q20[c] <= THR * med))
        c_q20, c_q0 = cs(np.where(cal, q20[c], 0)), cs(q0[c])
        st = sub.start.values.clip(0, a["sizes"][c]); en = sub.end.values.clip(0, a["sizes"][c])
        ncal = c_cal[en] - c_cal[st]; ln = np.maximum(en - st, 1)
        out.append(pd.DataFrame({
            "id": sub.id.values, "strain": s, "scorable_bases": ncal, "deleted_bases": c_del[en] - c_del[st],
            "rel_depth_unique_scorable": np.where(ncal > 0, (c_q20[en] - c_q20[st]) / np.maximum(ncal, 1) / med, np.nan),
            "rel_depth_allreads": (c_q0[en] - c_q0[st]) / ln / med}))
    return pd.concat(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intervals", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--depth", default="depth"); ap.add_argument("--ref", default="ref")
    ap.add_argument("--calls", default="calls"); ap.add_argument("--aligner", default="bwa")
    ap.add_argument("--scorable-bed", default=None, help="bases used as evidence (default: callable bed)")
    ap.add_argument("--samples", default="meta/copy38.txt"); ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(f"{a.ref}/E2.chr.genome")}
    callable_ = {c: np.zeros(sizes[c], bool) for c in CHROMS}
    for l in open(a.scorable_bed or f"{a.calls}/callable.{a.aligner}.bed"):
        c, s, e = l.split(); callable_[c][int(s):int(e)] = True
    iv = pd.read_csv(a.intervals, sep="\t", dtype={"id": str})
    job = dict(sizes=sizes, callable=callable_, depth=a.depth, aligner=a.aligner, iv=iv)
    samples = [l.strip() for l in open(a.samples) if l.strip()]
    with Pool(a.threads) as pool:
        res = pool.map(score, [(s, job) for s in samples], chunksize=1)
    pd.concat(res).to_csv(a.out, sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
