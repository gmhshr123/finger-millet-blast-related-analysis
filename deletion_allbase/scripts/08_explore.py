#!/usr/bin/env python3
"""First look at the all-base calls: junction support, sizes, Data S10 crosswalk."""
import numpy as np, pandas as pd
R = "deletion_allbase/results/cluster/"
cohort = [l.strip() for l in open("deletion_allbase/inputs/cohort33.txt")]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)

calls = pd.read_csv(R + "depth_calls.bwa.primary.tsv", sep="\t")
dl = pd.read_csv(R + "delly_deletions.tsv", sep="\t", header=None,
                 names="strain chrom pos end filt pe sr mapq precise homlen gt dr dv rr rv".split(), na_values=".")
dl["start"] = dl.pos                       # 0-based first deleted base
dl["len"] = dl.end - dl.start
dl["alt"] = dl.dv.fillna(0) + dl.rv.fillna(0); dl["refr"] = dl.dr.fillna(0) + dl.rr.fillna(0)
dl["vaf"] = dl.alt / (dl.alt + dl.refr).replace(0, np.nan)
print("Delly deletions:", len(dl), "PASS:", (dl.filt == "PASS").sum(), "precise:", dl.precise.notna().sum())
good = dl[(dl.alt >= 3) & (dl.vaf >= 0.8) & (dl.len >= 50)]
print("Delly with >=3 supporting reads and >=80% variant reads, >=50bp:", len(good))
print(good.groupby("strain").size().describe()[["min", "50%", "max"]].to_dict())
print("size quantiles depth calls:", calls.deleted_bases.quantile([.1, .5, .9, .99]).to_dict())

# junction match: Delly breakpoints within tolerance of the depth-call edges
def match(c, d, tol=150):
    out = np.zeros(len(c), bool); alt = np.zeros(len(c))
    for (s, ch), sub in c.groupby(["strain", "chrom"]):
        dd = d[(d.strain == s) & (d.chrom == ch)]
        if dd.empty: continue
        for i, r in zip(sub.index, sub.itertuples()):
            m = dd[(dd.start >= r.outer_start - tol) & (dd.start <= r.start + tol) &
                   (dd.end >= r.end - tol) & (dd.end <= r.outer_end + tol)]
            if len(m): out[c.index.get_loc(i)] = True; alt[c.index.get_loc(i)] = m.alt.max()
    return out, alt
calls = calls.reset_index(drop=True)
calls["junction"], calls["junction_reads"] = match(calls, dl[dl.alt >= 3])
calls["sizebin"] = pd.cut(calls.deleted_bases, [0, 200, 500, 1000, 5000, 20000, 1e9],
                          labels=["100-200", "200-500", "500-1k", "1-5k", "5-20k", ">20k"])
calls["group"] = np.where(calls.strain.isin(cohort), "cohort33", "near_ref")
print(calls.groupby(["group", "sizebin"], observed=True).agg(n=("junction", "size"), with_junction=("junction", "mean"),
                                                             bp=("deleted_bases", "sum")).round(3))
calls["frac_callable"] = calls.callable_bases / calls.span
print("fraction of span callable, quantiles:", calls.frac_callable.quantile([.05, .25, .5]).round(3).to_dict())
print("uncertain boundary (unknown bases at an end) share:", ((calls.start > calls.outer_start) | (calls.end < calls.outer_end)).mean().round(3))

# Data S10 rows scored against the all-base data
s10 = pd.read_csv("deletion_allbase/inputs/DataS10_events.tsv", sep="\t")
ps = pd.read_csv(R + "DataS10.per_strain_depth.tsv.gz", sep="\t")
rows = []
for r in s10.itertuples():
    for st in r.strains.split(","):
        x = ps[(ps.id == r.del_id) & (ps.strain == st)].iloc[0]
        rows.append((r.del_id, st, r.te, r.min_size, r.snp_count, x.callable_bases, x.deleted_bases, x.rel_depth_allreads))
x = pd.DataFrame(rows, columns="del_id strain te min_size snp_count callable deleted rel_all".split())
x["frac_del"] = x.deleted / x.callable.replace(0, np.nan)
x["status"] = np.select([x.callable < 50, x.frac_del >= 0.8, x.frac_del >= 0.2], ["uncallable", "deleted", "partly"], "not_deleted")
print("\nData S10 strain-by-event (402) vs all-base depth:\n", x.status.value_counts())
ev = x.groupby("del_id").agg(te=("te", "first"), n=("strain", "size"), n_del=("status", lambda s: (s == "deleted").sum()),
                             n_part=("status", lambda s: (s == "partly").sum()), n_unc=("status", lambda s: (s == "uncallable").sum()),
                             min_size=("min_size", "first"))
ev["event_status"] = np.select([ev.n_del == ev.n, ev.n_del + ev.n_part == ev.n - ev.n_unc, ev.n_del + ev.n_part > 0], ["all_deleted", "all_deleted_or_partly", "mixed"], "none_deleted")
ev.loc[ev.n_unc == ev.n, "event_status"] = "uncallable"
print(ev.event_status.value_counts()); print(pd.crosstab(ev.te != "No_TE", ev.event_status))
print(x[x.status.isin(["not_deleted"])].head(20))
x.to_csv("deletion_allbase/results/DataS10_strain_event_depth_check.tsv", sep="\t", index=False)
calls.to_csv("deletion_allbase/results/depth_calls_with_junction.tsv", sep="\t", index=False)
