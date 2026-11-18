#!/usr/bin/env python3
"""Reproduce the manuscript's SNP-based deletion candidates from Data S3.

A candidate is a run of consecutive SNP positions with a missing call ('-') in one strain.
Columns follow Data S10: start = first missing SNP, end = last missing SNP,
dist_above = start minus the last SNP that is present, dist_below = first SNP that is
present again minus end. min_size = end - start; max_size = distance between the two
flanking present SNPs.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--geno", default="deletion_allbase/inputs/DataS3_genotypes.tsv.gz")
ap.add_argument("--strains", default=None, help="file with strain IDs; default all")
ap.add_argument("--min-snps", type=int, default=4, help="smallest run kept (filter to 5 later)")
ap.add_argument("--out", default="deletion_allbase/results/snp_gap_runs.tsv")
a = ap.parse_args()

g = pd.read_csv(a.geno, sep="\t", skiprows=1, dtype=str, keep_default_na=False)
g["POS"] = g.POS.astype(int)
strains = [l.strip() for l in open(a.strains)] if a.strains else list(g.columns[4:])
rows = []
for chrom, sub in g.groupby("CHROM", sort=False):
    pos = sub.POS.values
    assert (np.diff(pos) >= 0).all(), chrom
    for s in strains:
        miss = (sub[s].values == "-")
        d = np.diff(np.concatenate(([0], miss.astype(np.int8), [0])))
        for i, j in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):   # SNP index run [i, j)
            n = j - i
            if n < a.min_snps:
                continue
            above = pos[i - 1] if i > 0 else None          # last present SNP before the run
            below = pos[j] if j < len(pos) else None       # first present SNP after the run
            rows.append((s, chrom, pos[i], pos[j - 1], n,
                         pos[i] - above if above is not None else -1,
                         below - pos[j - 1] if below is not None else -1,
                         above if above is not None else -1, below if below is not None else -1))
out = pd.DataFrame(rows, columns=["strain", "chrom", "start", "end", "snp_count", "dist_above",
                                  "dist_below", "snp_above", "snp_below"])
out["min_size"] = out.end - out.start
out["max_size"] = np.where((out.snp_above >= 0) & (out.snp_below >= 0), out.snp_below - out.snp_above, -1)
out.to_csv(a.out, sep="\t", index=False)
print(len(out), "runs;", out.strain.nunique(), "strains")
print(out.groupby("snp_count").size().head(8))
