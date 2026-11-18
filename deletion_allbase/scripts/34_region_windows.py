#!/usr/bin/env python3
"""Plotting windows for the read-depth plots of the breakpoint browser (run from manuscripts_all/, after script 29).

Input:  deletion_allbase/results_edta/regions_for_browser.tsv (script 29), E2.chr.genome
Output: deletion_allbase/results_edta/region_windows.tsv
        Copy to the cluster and run 32_snp_unconfirmed_depth.sbatch with this file (see that script), then script 31.

One row per deleted region. The window is the region plus a margin on both sides (25% of its length, at least
1.5 kb). The two zoom plots are centred on the read-defined breakpoints where they exist, otherwise on the
first and last deleted base of the region (columns start / end, 1-based, as the depth script expects).
"""
import os
import pandas as pd

R = os.environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")
sizes = {l.split()[0]: int(l.split()[1]) for l in open(R + "cluster_uniq/E2.chr.genome")}
g = pd.read_csv(R + "regions_for_browser.tsv", sep="\t", keep_default_na=False)
d = pd.DataFrame({"strain": g.Strain, "event": g["Deletion event ID"], "chrom": g.Chromosome,
                  "reg_start": g["Region start"].astype(int), "reg_end": g["Region end"].astype(int)})
lb = pd.to_numeric(g["Left breakpoint"], errors="coerce"); rb = pd.to_numeric(g["Right breakpoint"], errors="coerce")
d["start"] = lb.fillna(d.reg_start).astype(int)
d["end"] = rb.fillna(d.reg_end).astype(int)
lo, hi = d[["reg_start", "start"]].min(axis=1), d[["reg_end", "end"]].max(axis=1)
pad = (0.25 * (hi - lo)).clip(lower=1500).astype(int)
d["win_start"] = (lo - pad).clip(lower=0)
d["win_end"] = pd.concat([hi + pad, d.chrom.map(sizes)], axis=1).min(axis=1)
d.insert(0, "idx", range(len(d)))
d.to_csv(R + "region_windows.tsv", sep="\t", index=False)
print(len(d), "regions; window length: median", int((d.win_end - d.win_start).median()), "max", int((d.win_end - d.win_start).max()))
