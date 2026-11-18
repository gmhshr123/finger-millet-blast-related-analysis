import pandas as pd, numpy as np
d = pd.read_csv("deletion_allbase/results/cluster/delly_deletions.tsv", sep="\t", header=None, na_values=".",
                names="strain chrom pos end filt pe sr mapq precise homlen gt dr dv rr rv".split())
d["alt"] = d.dv.fillna(0) + d.rv.fillna(0)
d = d[(d.alt >= 3) & (d.end - d.pos >= 50)].drop_duplicates(["chrom", "pos", "end"])
d = d[d.chrom.str.startswith("Chr")]
d.assign(id=d.chrom + ":" + d.pos.astype(str) + "-" + d.end.astype(str))[["id", "chrom", "pos", "end"]].rename(
    columns={"pos": "start"}).to_csv("deletion_allbase/inputs/delly.intervals.tsv", sep="\t", index=False)
print(len(d), "unique Delly deletion intervals")