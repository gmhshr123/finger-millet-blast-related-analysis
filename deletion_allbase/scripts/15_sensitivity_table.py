#!/usr/bin/env python3
"""Do the conclusions depend on the depth cut-off, the aligner, or strain E4?

Collects the headline numbers of the primary run and of the sensitivity runs made with
  python 10_analyze.py 1000 <depth cut-off> <aligner>
into one table, and compares each sensitivity call set with the primary one base by base.
Writes results/sensitivity_summary.tsv.
"""
import glob, json
import numpy as np
import pandas as pd

R = "deletion_allbase/results/"
C = R + "cluster_uniq/"
CHROMS = [f"Chr0{i}" for i in range(1, 8)]
sizes = {l.split()[0]: int(l.split()[1]) for l in open(C + "E2.chr.genome")}


def masks(calls):
    """strain -> chrom -> boolean array of bases inside a deletion region"""
    out = {}
    for (s, c), sub in calls.groupby(["strain", "chrom"]):
        m = out.setdefault(s, {}).setdefault(c, np.zeros(sizes[c], bool))
        for a, b in zip(sub.start, sub.end):
            m[a:b] = True
    return out


def overlap(m1, m2):
    inter = union = 0
    for s in set(m1) | set(m2):
        for c in CHROMS:
            a = m1.get(s, {}).get(c); b = m2.get(s, {}).get(c)
            if a is None and b is None:
                continue
            a = a if a is not None else np.zeros(sizes[c], bool); b = b if b is not None else np.zeros(sizes[c], bool)
            inter += int((a & b).sum()); union += int((a | b).sum())
    return inter, union


prim = pd.read_csv(R + "allbase_calls_strain_level.tsv", sep="\t")
pm = masks(prim)
runs = [("primary: BWA, depth <= 5% of median", R + "summary.json", R + "allbase_calls_strain_level.tsv")]
# the runs to compare with the primary one are named here, so that left-over files from earlier
# primary settings in results/sensitivity/ are not picked up
for tag in ("bwa.thr0", "bwa.thr0.1", "bwa.thr0.2", "bowtie2.thr0.05"):
    f = R + f"sensitivity/summary.{tag}.gap1000.json"
    aln, thr = tag.split(".thr")
    label = f"{'BWA' if aln == 'bwa' else 'Bowtie2 (original alignments)'}, depth <= {float(thr) * 100:g}% of median"
    runs.append((label, f, R + f"sensitivity/allbase_calls_strain_level.{tag}.gap1000.tsv"))

rows = []
for label, js, cf in runs:
    S = json.load(open(js)); calls = pd.read_csv(cf, sep="\t")
    inter, union = overlap(pm, masks(calls))
    ev = pd.read_csv(cf.replace("allbase_calls_strain_level", "allbase_events"), sep="\t")
    big = ev[ev.deleted_bases >= 1000]
    rows.append({
        "run": label, "regions": S["regions_strain_level"], "events": S["events"],
        "events_100_500bp": S["events_by_size"]["100-500 bp"], "events_ge1kb": len(big),
        "median_deleted_bp_per_strain": S["deleted_scorable_bp_per_strain"]["median"],
        "bases_shared_with_primary_pct": round(100 * inter / max(union, 1), 1),
        "primary_bases_recovered_pct": round(100 * inter / max(sum(int(v.sum()) for d in pm.values() for v in d.values()), 1), 1),
        "snp_method_recall_regions_pct": round(100 * S["snp_method_recall"]["regions"], 1),
        "snp_method_recall_bp_pct": round(100 * S["snp_method_recall"]["deleted_bp"], 1),
        "refTE_at_breakpoint_all_pct": round(100 * S["te_enrichment_reference_annotation"]["observed_te_at_either_breakpoint"], 1),
        "refTE_expected_same_windows_pct": round(100 * S["te_enrichment_reference_annotation"]["expected_random_same_100kb_windows"], 2),
        "refTE_at_breakpoint_100_500bp_pct": round(100 * S["ref_te_at_breakpoint_by_size"]["100-500 bp"]["mean"], 1),
        "refTE_at_breakpoint_ge1kb_pct": round(100 * S["ref_te_at_breakpoint_events_ge1kb"], 1),
        "refTE_found_by_snp_pct": round(100 * S["ref_te_by_snp_detection"]["true"]["mean"], 1),
        "refTE_missed_by_snp_pct": round(100 * S["ref_te_by_snp_detection"]["false"]["mean"], 1),
        "genes_deleted": S["genes"]["genes_at_least_half_deleted"], "effector_share_pct": round(100 * S["genes"]["effector_share"], 1),
        "effector_p": S["genes"]["hypergeometric_p"],
    })

# primary run without E4 (lowest depth, 26x): events that do not depend on E4
ev = pd.read_csv(R + "allbase_events.tsv", sep="\t")
e4only = ev.strains == "E4"
k = ev[~e4only]
ref = (k.L_te_ref.notna()) | k.L_rm_ref | (k.R_te_ref.notna()) | k.R_rm_ref
rows.append({"run": "primary, events carried only by E4 removed", "regions": int((prim.strain != "E4").sum()), "events": len(k),
             "events_100_500bp": int((k.size_class == "100-500 bp").sum()), "events_ge1kb": int((k.deleted_bases >= 1000).sum()),
             "refTE_at_breakpoint_100_500bp_pct": round(100 * ref[k.size_class == "100-500 bp"].mean(), 1),
             "refTE_at_breakpoint_ge1kb_pct": round(100 * ref[k.deleted_bases >= 1000].mean(), 1),
             "refTE_found_by_snp_pct": round(100 * ref[k.found_by_snp_method].mean(), 1),
             "refTE_missed_by_snp_pct": round(100 * ref[~k.found_by_snp_method].mean(), 1)})
T = pd.DataFrame(rows)
T.to_csv(R + "sensitivity_summary.tsv", sep="\t", index=False)
print(T.set_index("run").T.to_string())
print("\nE4: regions", int((prim.strain == "E4").sum()), "; events carried only by E4:", int(e4only.sum()),
      "of which 100-500 bp:", int((ev[e4only].size_class == "100-500 bp").sum()))
print("private events per strain (median over strains):", ev[ev.n_strains == 1].strains.value_counts().median())
print(ev[ev.n_strains == 1].strains.value_counts().head(6).to_string())
