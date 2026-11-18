#!/usr/bin/env python3
"""What kind of event is each all-base "deletion"? (run locally, after 10_analyze.py and job 14)

Read depth only says that E2 sequence gets no uniquely mapped reads in a strain. Here every
strain-level region gets one of four labels from the reads around it:

  simple deletion        the two sides of the interval are joined in the strain: a Delly junction, split
                         reads that continue at the other edge, or read pairs that straddle the interval
                         and look too far apart by the deletion size (script 14, LLR >= 10).
  other sequence at edge no join across the interval, but reads at an edge continue into something else:
                         a TE (clipped part matches the TE library), another place in the genome (mate maps
                         far away), or sequence that is not in E2 at all (mate unmapped). The E2 sequence is
                         absent, but this is not a clean two-sided deletion: it can be a TE-mediated deletion
                         (one TE copy is left behind, too long for a read pair to cross), a replacement,
                         or a larger rearrangement.
  DNA still present      read pairs straddle the interval at the normal distance, so DNA of about the same
                         length is still there. Reads just do not map to it (diverged sequence). Not a deletion.
                         This test only has power for intervals shorter than about 300 bp.
  depth only             none of the above. Nothing but missing depth supports the call.

Each region also gets a "poorly mappable" flag from script 17: E2's own unique-read depth in the interval
is below 60% of E2's median, or one of the four near-reference strains (E5, E6, E19, E36) is below 50% of
its median. These strains carry the E2 sequence, so low depth there means the place is hard to map
uniquely, not that DNA is missing.

The thresholds were set with the matched control intervals (same strain, same size, no deletion called),
and the script prints how often each kind of evidence turns up in those controls.
An event (locus) takes the label carried by at least half of its strains.
"""
import json
import numpy as np
import pandas as pd

R = __import__("os").environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")    # primary result set: combined EDTA TE library; set DEL_RESULTS=deletion_allbase/results/ for the earlier stand-in library
C = R + "cluster_uniq/"
cohort = [l.strip() for l in open("deletion_allbase/inputs/cohort33.txt")]
T_DEL, T_PRES, MIN_READS, MIN_SHARE = 10, 5, 3, 0.10
SIZES = ["100-500 bp", "0.5-1 kb", "1-5 kb", "5-20 kb", ">20 kb"]
S = {"thresholds": dict(llr_deletion=T_DEL, llr_present=-T_PRES, min_reads=MIN_READS, min_share_of_edge_reads=MIN_SHARE)}


def pair_flags(d):
    d = d.copy()
    d["pairs_deletion"] = (d.llr_deletion_vs_present >= T_DEL) & (d.n_straddling_pairs >= MIN_READS)
    d["pairs_present"] = ((d.llr_deletion_vs_present <= -T_PRES) & (d.n_straddling_pairs >= MIN_READS)
                          & (d.n_straddling_pairs >= 0.3 * d.expected_pairs_if_present))
    for side in "LR":
        n = d[f"{side}_anchor_reads"].clip(lower=1)
        d[f"{side}_unmapped"] = (d[f"{side}_mate_unmapped"] >= MIN_READS) & (d[f"{side}_mate_unmapped"] / n >= MIN_SHARE)
        d[f"{side}_far"] = (d[f"{side}_mate_far"] >= MIN_READS) & (d[f"{side}_mate_far"] / n >= MIN_SHARE)
    d["edge_mate_unmapped"] = d.L_unmapped | d.R_unmapped
    d["edge_mate_far"] = d.L_far | d.R_far
    return d


# ------------------------------------------------------------------ how often does each signal appear by chance?
ctrl = pair_flags(pd.read_csv(C + "pairs_controls.tsv", sep="\t")); ctrl = ctrl[ctrl.strain.isin(cohort)]
frg = pair_flags(pd.read_csv(C + "pairs_fragments.tsv", sep="\t")); frg = frg[frg.strain.isin(cohort)]
sc_c = pd.read_csv(C + "softclip_edges_controls.tsv", sep="\t"); sc_c = sc_c[sc_c.strain.isin(cohort)]
sc_r = pd.read_csv(C + "softclip_edges_all.tsv", sep="\t"); sc_r = sc_r[sc_r.strain.isin(cohort)]


def top_count(x):
    return int(x.split(",")[0].split(":")[1]) if isinstance(x, str) and x else 0


def clip_flags(d):
    te = np.array([max(top_count(a), top_count(b)) >= MIN_READS for a, b in zip(d.L_te, d.R_te)])
    return dict(clipped_reads_continue_into_TE=float(te.mean()),
                clipped_reads_continue_at_other_edge=float(((d.L_sa_other_edge >= MIN_READS) | (d.R_sa_other_edge >= MIN_READS)).mean()))


bins = pd.IntervalIndex.from_breaks([0, 150, 250, 500, 1000, 10**9])
chance = {}
for name, d, sc in (("control_intervals", ctrl, sc_c), ("real_fragments", frg, sc_r)):
    x = dict(n=len(d), pairs_say_deletion=float(d.pairs_deletion.mean()), pairs_say_present=float(d.pairs_present.mean()),
             edge_mate_unmapped=float(d.edge_mate_unmapped.mean()), edge_mate_far=float(d.edge_mate_far.mean()), **clip_flags(sc))
    g = d.groupby(pd.cut(d.span, bins), observed=True)
    x["by_span"] = {str(k): dict(n=int(v.shape[0]), pairs_say_deletion=round(float(v.pairs_deletion.mean()), 4),
                                 pairs_say_present=round(float(v.pairs_present.mean()), 4)) for k, v in g}
    chance[name] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in x.items()}
S["evidence_frequency_real_vs_control"] = chance

# ------------------------------------------------------------------ label every strain-level region
calls = pd.read_csv(R + "allbase_calls_strain_level.tsv", sep="\t")
calls["id"] = calls.chrom + ":" + calls.start.astype(str) + "-" + calls.end.astype(str)
pr = pair_flags(pd.read_csv(C + "pairs_regions.tsv", sep="\t"))
keep = ["strain", "id", "n_straddling_pairs", "median_pair_distance", "lib_median", "expected_pairs_if_present", "best_deletion_size",
        "llr_deletion_vs_present", "pairs_deletion", "pairs_present", "edge_mate_unmapped", "edge_mate_far", "mismatches_per_100bp_inside"]
calls = calls.merge(pr[keep], on=["strain", "id"], how="left", validate="1:1")
ctx = pd.read_csv(C + "region_context.tsv", sep="\t").drop(columns=["id"])
calls = calls.merge(ctx, on=["strain", "chrom", "start", "end"], how="left", validate="1:1")
calls["poorly_mappable"] = (calls.E2_rel_depth < 0.6) | (calls.nearref_rel_depth_min < 0.5)
calls["edge_te_reads"] = (calls.L_te_read_count >= MIN_READS) | (calls.R_te_read_count >= MIN_READS)
joined = calls.junction | calls.pairs_deletion
other = calls.edge_te_reads | calls.edge_mate_far | calls.edge_mate_unmapped
calls["read_class"] = np.select([joined, calls.pairs_present, other], ["simple deletion", "DNA still present", "other sequence at edge"], "depth only")
calls["edge_detail"] = np.select([calls.edge_te_reads, calls.edge_mate_far, calls.edge_mate_unmapped],
                                 ["TE (clipped reads)", "elsewhere in genome (mate far)", "sequence not in E2 (mate unmapped)"], "")

# ------------------------------------------------------------------ events
ev = pd.read_csv(R + "allbase_events.tsv", sep="\t")
sh = calls.groupby("event_id").read_class.value_counts(normalize=True).unstack(fill_value=0)
for k in ["simple deletion", "DNA still present", "other sequence at edge", "depth only"]:
    if k not in sh:
        sh[k] = 0.0
cls = np.select([sh["simple deletion"] >= 0.5, sh["DNA still present"] >= 0.5,
                 sh["simple deletion"] + sh["other sequence at edge"] >= 0.5], ["simple deletion", "DNA still present", "other sequence at edge"], "depth only")
ev = ev.merge(pd.DataFrame({"event_id": sh.index, "read_class": cls, "share_simple_deletion": sh["simple deletion"].round(3).values,
                            "share_other_sequence": sh["other sequence at edge"].round(3).values,
                            "share_present": sh["DNA still present"].round(3).values}), on="event_id")
ev = ev.merge(calls.groupby("event_id").agg(poorly_mappable=("poorly_mappable", lambda v: v.mean() >= 0.5), gc=("gc", "median"),
                                            E2_rel_depth=("E2_rel_depth", "median"), nearref_rel_depth_min=("nearref_rel_depth_min", "median"),
                                            carrier_flank_rel_depth=("flank_rel_depth_L", "median")).reset_index(), on="event_id")
ev["size_class"] = pd.Categorical(ev.size_class, SIZES, ordered=True)
ev["te"] = ev.te_class != "no TE evidence"
ev["ref_te_edge"] = ev.L_te_ref.notna() | ev.L_rm_ref | ev.R_te_ref.notna() | ev.R_rm_ref

S["regions_by_read_class"] = calls.read_class.value_counts().to_dict()
S["events_by_read_class"] = ev.read_class.value_counts().to_dict()
S["events_read_class_by_size"] = pd.crosstab(ev.size_class, ev.read_class).to_dict("index")
S["events_read_class_by_te_class"] = pd.crosstab(ev.te_class, ev.read_class).to_dict("index")
S["regions_poorly_mappable_by_read_class"] = calls.groupby("read_class").poorly_mappable.agg(["size", "mean"]).round(3).to_dict("index")
S["events_poorly_mappable_by_read_class"] = ev.groupby("read_class").poorly_mappable.agg(["size", "mean"]).round(3).to_dict("index")
S["events_poorly_mappable_by_read_class_and_size"] = (ev.groupby(["size_class", "read_class"], observed=True).poorly_mappable.agg(["size", "sum"])
                                                      .reset_index().to_dict("records"))
S["depth_only_events"] = dict(
    poorly_mappable=int(((ev.read_class == "depth only") & ev.poorly_mappable).sum()),
    n=int((ev.read_class == "depth only").sum()),
    with_reference_repeat_at_an_edge=int(((ev.read_class == "depth only") & ev.ref_te_edge).sum()),
    telomeric=int(((ev.read_class == "depth only") & ev.telomeric).sum()),
    private_to_one_strain=int(((ev.read_class == "depth only") & (ev.n_strains == 1)).sum()))


def headline(e):
    big = e[e.deleted_bases >= 1000]
    g = e.groupby(["size_class", "found_by_snp_method"], observed=True).te.agg(["size", "mean"]).round(3).reset_index()
    return dict(events=len(e), by_size=e.size_class.value_counts().sort_index().to_dict(),
                te_share=round(float(e.te.mean()), 3), te_share_by_size=e.groupby("size_class", observed=True).te.mean().round(3).to_dict(),
                te_share_ge1kb=round(float(big.te.mean()), 3), n_ge1kb=len(big),
                same_family_both_edges=int((e.te_class == "TE at both breakpoints, same named family").sum()),
                found_by_snp_method=round(float(e.found_by_snp_method.mean()), 3),
                te_share_found_by_snp=round(float(e[e.found_by_snp_method].te.mean()), 3),
                te_share_missed_by_snp=round(float(e[~e.found_by_snp_method].te.mean()), 3),
                te_share_by_size_and_snp_detection=g.to_dict("records"))


S["all_events"] = headline(ev)
S["without_DNA_still_present"] = headline(ev[ev.read_class != "DNA still present"])
S["read_supported_only"] = headline(ev[ev.read_class.isin(["simple deletion", "other sequence at edge"])])
S["simple_deletions_only"] = headline(ev[ev.read_class == "simple deletion"])
S["read_supported_and_mappable"] = headline(ev[ev.read_class.isin(["simple deletion", "other sequence at edge"]) & ~ev.poorly_mappable])
S["mappable_not_contradicted"] = headline(ev[(ev.read_class != "DNA still present") & ~ev.poorly_mappable])

# effectors among deleted genes, for the read-supported events only
import io, math
# E2.genes.bed has Windows line ends, which left a stray carriage return before the gene name
genes = pd.read_csv(io.StringIO(open(C + "E2.genes.bed", newline="").read().replace("\r\t", "\t").replace("\r", "")),
                    sep="\t", header=None, names=["chrom", "start", "end", "gene"])
genes = genes[genes.chrom.str.startswith("Chr")]
s2 = pd.read_csv("deletion_allbase/inputs/DataS2_gene_sections.tsv", sep="\t")
eff = set(s2[s2.section.str.startswith("Effector")].gene)


def effector_test(e):
    g = sorted(set(",".join(e.genes.dropna()).split(",")) - {""})
    k = sum(x in eff for x in g); n = len(g); K = int(genes.gene.isin(eff).sum()); N = len(genes)
    p = sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(k, min(K, n) + 1)) / math.comb(N, n)
    return dict(genes=n, effectors=k, effector_share=round(k / max(n, 1), 4), genome_share=round(K / N, 4), hypergeometric_p=float(f"{p:.3g}"))


S["effectors_all_events"] = effector_test(ev)
S["effectors_read_supported"] = effector_test(ev[ev.read_class.isin(["simple deletion", "other sequence at edge"])])
S["effectors_read_supported_and_mappable"] = effector_test(ev[ev.read_class.isin(["simple deletion", "other sequence at edge"]) & ~ev.poorly_mappable])

# Data S10 rows: read class of the all-base events they correspond to
x = pd.read_csv(R + "DataS10_crosswalk.tsv", sep="\t")
rc = dict(zip(ev.event_id, ev.read_class))
order = ["simple deletion", "other sequence at edge", "depth only", "DNA still present"]
x["read_class"] = [next((o for o in order if o in {rc[e] for e in str(a).split(",") if e in rc}), "no all-base event")
                   for a in x.allbase_events.fillna("")]
S["DataS10_rows_by_read_class"] = pd.crosstab(x.s10_te != "No_TE", x.read_class).rename(index={True: "TE", False: "No_TE"}).to_dict("index")

calls.drop(columns=["id"]).to_csv(R + "allbase_calls_strain_level_read_evidence.tsv", sep="\t", index=False)
ev.drop(columns=["te"]).to_csv(R + "allbase_events_read_evidence.tsv", sep="\t", index=False)
x.to_csv(R + "DataS10_crosswalk_read_evidence.tsv", sep="\t", index=False)
json.dump(S, open(R + "read_evidence_summary.json", "w"), indent=1, default=str)
print(json.dumps(S, indent=1, default=str))
