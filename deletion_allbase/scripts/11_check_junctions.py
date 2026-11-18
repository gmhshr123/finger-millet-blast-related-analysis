#!/usr/bin/env python3
"""Why do few depth calls have a Delly junction? Check overlap in both directions."""
import numpy as np, pandas as pd
R = "deletion_allbase/results/"; C = R + "cluster_uniq/"
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
cohort = [l.strip() for l in open("deletion_allbase/inputs/cohort33.txt")]
calls = pd.read_csv(R + "allbase_calls_strain_level.tsv", sep="\t")
dl = pd.read_csv(C + "delly_deletions.tsv", sep="\t", header=None, na_values=".",
                 names="strain chrom pos end filt pe sr mapq precise homlen gt dr dv rr rv".split())
dl["alt"] = dl.dv.fillna(0) + dl.rv.fillna(0); dl["ref"] = dl.dr.fillna(0) + dl.rr.fillna(0)
dl = dl[(dl.alt >= 3) & (dl.end - dl.pos >= 50) & dl.strain.isin(cohort) & dl.chrom.str.startswith("Chr")].copy()
dl["id"] = dl.chrom + ":" + dl.pos.astype(str) + "-" + dl.end.astype(str)
ps = pd.read_csv(C + "delly.per_strain_depth.tsv.gz", sep="\t")
dl = dl.merge(ps, on=["id", "strain"], how="left")
dl["len"] = dl.end - dl.pos
dl["frac_del"] = dl.deleted_bases / dl.scorable_bases.replace(0, np.nan)
dl["depth_status"] = np.select([dl.scorable_bases < 30, dl.frac_del >= 0.8, dl.frac_del >= 0.2], ["no scorable sequence (repeat)", "deleted", "partly"], "not deleted")
dl["vaf"] = dl.alt / (dl.alt + dl.ref)
dl["lenbin"] = pd.cut(dl.len, [0, 100, 500, 1000, 5000, 20000, 1e9])
print("Delly deletions (>=3 reads) in cohort strains:", len(dl))
print(pd.crosstab(dl.lenbin, dl.depth_status))
print("variant read fraction by depth status:\n", dl.groupby("depth_status").vaf.describe()[["count", "25%", "50%", "75%"]].round(2))

# depth calls: any Delly call overlapping with >=50 % reciprocal overlap (looser than breakpoint matching)
any50 = np.zeros(len(calls), bool); anyov = np.zeros(len(calls), bool)
for (s, ch), sub in calls.groupby(["strain", "chrom"]):
    dd = dl[(dl.strain == s) & (dl.chrom == ch)]
    for i, r in zip(sub.index, sub.itertuples()):
        ov = np.minimum(dd.end, r.outer_end) - np.maximum(dd.pos, r.outer_start)
        anyov[i] = (ov > 0).any()
        ov2 = np.minimum(dd.end, r.end) - np.maximum(dd.pos, r.start)
        any50[i] = (ov2 >= 0.5 * np.maximum(dd.len, r.end - r.start)).any()
calls["delly_any_overlap"], calls["delly_recip50"] = anyov, any50
calls["boundary_fuzzy"] = (calls.start - calls.outer_start > 20) | (calls.outer_end - calls.end > 20)
print(calls.groupby("size_class", observed=True).agg(n=("junction", "size"), junction=("junction", "mean"), recip50=("delly_recip50", "mean"),
                                                     any_overlap=("delly_any_overlap", "mean"), fuzzy_boundary=("boundary_fuzzy", "mean"),
                                                     scor_frac=("scorable_bases", lambda x: 0)).round(3))
print("mean unique depth inside calls relative to median:", (calls.mean_depth_unique_scorable / calls.strain_median_depth).describe().round(3).to_dict())
print("mean all-read depth inside calls relative to median:", (calls.mean_depth_allreads / calls.strain_median_depth).describe().round(3).to_dict())
print("deleted / scorable bases inside span:", (calls.deleted_bases / calls.scorable_bases).describe().round(3).to_dict())
ex = calls[(calls.strain == "E13") & ~calls.junction & ~calls.boundary_fuzzy]
print(ex.groupby("size_class", observed=True).head(3)[["strain", "chrom", "start", "end", "outer_start", "outer_end", "deleted_bases", "mean_depth_allreads", "snps_inside", "found_by_snp_method"]].to_string())
