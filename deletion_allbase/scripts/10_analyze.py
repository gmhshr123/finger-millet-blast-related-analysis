#!/usr/bin/env python3
"""Main analysis of the all-base deletion calls for the 33-strain cohort.

Produces, in deletion_allbase/results/:
  allbase_fragments_strain_level.tsv  raw depth calls (one row per strain x uninterrupted deleted stretch)
  allbase_calls_strain_level.tsv      one row per strain x deletion region (neighbouring fragments merged)
  allbase_events.tsv                  one row per locus (regions merged across strains) with TE / gene annotation
  snp_candidates_checked.tsv          every SNP-gap candidate scored against the all-base data
  DataS10_crosswalk.tsv               every Data S10 row with its all-base disposition
  windows100kb_tracks.tsv             100-kb tracks: all-base deletions, SNP-based deletions, SNP density, repeats
  summary.json                        the headline numbers
"""
import io, json, math, sys
import numpy as np
import pandas as pd

R = __import__("os").environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")    # primary result set: combined EDTA TE library; set DEL_RESULTS=deletion_allbase/results/ for the earlier stand-in library
C = R + "cluster_uniq/"
I = "deletion_allbase/inputs/"
MERGE_GAP = int(sys.argv[1]) if len(sys.argv) > 1 else 1000     # fragments of one strain closer than this are one region
SUFFIX = "" if MERGE_GAP == 1000 else f".gap{MERGE_GAP}"
# optional sensitivity run: depth cut-off (fraction of the strain's median) and aligner. Tables go to results/sensitivity/.
# The split-read TE evidence was only collected for the primary calls, so in a sensitivity run only the
# reference-annotation TE numbers are comparable with the primary analysis.
THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
ALN = sys.argv[3] if len(sys.argv) > 3 else "bwa"
# optional 4th argument "gapaware": also join two neighbouring calls of a strain (up to 20 kb apart) when the DNA
# between them does not look truly present. A real single-copy stretch has about normal unique-read depth. If the
# depth there is below half of the strain's median, or more than 30% of the reads there are multi-mapping, the
# reads most likely come from other copies of a repeat (or from a few stray reads), and the two calls are one
# deletion. The stretches between calls are described by 22_gap_context.py.
GAPAWARE = len(sys.argv) > 4 and sys.argv[4] == "gapaware"
PRIMARY = (THR == 0.05 and ALN == "bwa" and not GAPAWARE)
OUT = R if PRIMARY else R + "sensitivity/"
if not PRIMARY:
    import os
    os.makedirs(OUT, exist_ok=True)
    SUFFIX = f".{ALN}.thr{THR:g}.gap{MERGE_GAP}" + (".gapaware" if GAPAWARE else "")
MIN_READS = 3
rng = np.random.default_rng(20261007)
cohort = [l.strip() for l in open(I + "cohort33.txt")]
CHROMS = [f"Chr0{i}" for i in range(1, 8)]
sizes = {l.split()[0]: int(l.split()[1]) for l in open(C + "E2.chr.genome")}
S = {"merge_gap_bp": MERGE_GAP, "depth_threshold_frac": THR, "aligner": ALN}


def cover(bed3):
    m = {c: np.zeros(sizes[c], bool) for c in CHROMS}
    for c, s, e in bed3.itertuples(index=False):
        if c in m:
            m[c][s:e] = True
    return m, {c: np.concatenate(([0], np.cumsum(v))) for c, v in m.items()}


def n_in(cs, c, s, e):
    s = max(0, min(s, sizes[c])); e = max(s, min(e, sizes[c]))
    return int(cs[c][e] - cs[c][s])


def top_te(txt):
    """'Pot2:7,Pyret:1' -> ('Pot2', 7)"""
    if not isinstance(txt, str) or not txt:
        return "", 0
    k, v = txt.split(",")[0].split(":")
    return k, int(v)


# ------------------------------------------------------------------ inputs
frag = pd.read_csv(C + f"depth_calls.{ALN}.uniq.all_thresholds.tsv.gz", sep="\t")
frag = frag[np.isclose(frag.thr_frac, THR) & (frag.deleted_bases >= 100)].reset_index(drop=True)
frag["id"] = frag.chrom + ":" + frag.start.astype(str) + "-" + frag.end.astype(str)
frag = frag.merge(pd.read_csv(C + "softclip_edges_all.tsv", sep="\t"), on=["strain", "id"], how="left")
near = frag[~frag.strain.isin(cohort)]
frag = frag[frag.strain.isin(cohort)].sort_values(["strain", "chrom", "start"]).reset_index(drop=True)

scor_m, scor_cs = cover(pd.read_csv(C + f"scorable.{ALN}.uniq.bed", sep="\t", header=None))
rm = pd.read_csv(C + "E2.repeats.bed", sep="\t", header=None, names=["chrom", "start", "end", "fam", "cls", "strand"])
te = pd.read_csv(C + "E2.named_TEs.bed", sep="\t", header=None).iloc[:, [0, 1, 2, 3, 5]]
te.columns = ["chrom", "start", "end", "fam", "strand"]
te["fam"] = te.fam.str.replace(r"_\d+$", "", regex=True)          # Fosbury_1..7 -> Fosbury
te = te[te.end - te.start >= 50]
interspersed = rm[~rm.cls.isin(["Simple_repeat", "Low_complexity", "Satellite"]) & (rm.end - rm.start >= 50)]
_, rep_cs = cover(pd.concat([rm.iloc[:, :3], te.iloc[:, :3]]))
_, te_cs = cover(te.iloc[:, :3])
# E2-gene-pos.txt has Windows line ends, which left a stray carriage return before the gene name
genes = pd.read_csv(io.StringIO(open(C + "E2.genes.bed", newline="").read().replace("\r\t", "\t").replace("\r", "")),
                    sep="\t", header=None, names=["chrom", "start", "end", "gene"])
genes = genes[genes.chrom.isin(CHROMS)]
s2 = pd.read_csv(I + "DataS2_gene_sections.tsv", sep="\t")
effectors = set(s2[s2.section.str.startswith("Effector")].gene)
genes["effector"] = genes.gene.isin(effectors)
S["cohort_strains"] = len(cohort)
S["scorable_bp"] = int(sum(v[-1] for v in scor_cs.values())); S["genome_bp"] = int(sum(sizes.values()))
S["E2_genes_on_chromosomes"] = len(genes); S["E2_effector_genes"] = int(genes.effector.sum())
S["fragments_strain_level"] = len(frag)

# ------------------------------------------------------------------ merge neighbouring fragments within a strain
prev = frag.shift(1)
new = (frag.strain != prev.strain) | (frag.chrom != prev.chrom) | (frag.start - prev.end > MERGE_GAP)
if GAPAWARE:
    gc_ = pd.read_csv(C + "gap_context.tsv", sep="\t")
    art = gc_[(gc_.rel_depth_unique < 0.5) | (gc_.unique_share < 0.7) | (gc_.scorable_bases == 0)]
    art = set(zip(art.strain, art.chrom, art.gap_start, art.gap_end))
    is_art = np.array([(s, c, a, b) in art for s, c, a, b in zip(frag.strain, frag.chrom, prev.end.fillna(-1).astype(int), frag.start)])
    S["gapaware_extra_joins"] = int((new & is_art & (frag.strain == prev.strain) & (frag.chrom == prev.chrom)).sum())
    new = new & ~(is_art & (frag.strain == prev.strain) & (frag.chrom == prev.chrom))
frag["region"] = new.cumsum()
first = frag.groupby("region").first(); last = frag.groupby("region").last()
calls = pd.DataFrame({
    "strain": first.strain, "chrom": first.chrom, "start": first.start, "end": last.end,
    "outer_start": first.outer_start, "outer_end": last.outer_end,
    "deleted_bases": frag.groupby("region").deleted_bases.sum(), "scorable_bases_in_fragments": frag.groupby("region").scorable_bases.sum(),
    "n_fragments": frag.groupby("region").size(), "strain_median_depth": first.strain_median_depth,
    "L_bp": first.L_bp, "L_bp_reads": first.L_bp_reads, "L_te": first.L_te, "L_sa_other_edge": first.L_sa_other_edge,
    "R_bp": last.R_bp, "R_bp_reads": last.R_bp_reads, "R_te": last.R_te, "R_sa_other_edge": last.R_sa_other_edge,
}).sort_values(["chrom", "start", "end"]).reset_index(drop=True)
calls["span"] = calls.end - calls.start
calls["scorable_bases_in_span"] = [n_in(scor_cs, c, s, e) for c, s, e in calls[["chrom", "start", "end"]].itertuples(index=False)]
calls["deleted_frac_of_scorable_span"] = calls.deleted_bases / calls.scorable_bases_in_span
for side in "LR":
    t = [top_te(x) for x in calls[f"{side}_te"]]
    calls[f"{side}_te_family_reads"] = [k if v >= MIN_READS else "" for k, v in t]
    calls[f"{side}_te_read_count"] = [v for k, v in t]

# ------------------------------------------------------------------ junction evidence
# (a) Delly deletion whose two breakpoints fall in the two breakpoint zones; (b) >= 3 soft-clipped reads at an
# edge whose other part aligns at the other edge of the same region.
dl = pd.read_csv(C + "delly_deletions.tsv", sep="\t", header=None, na_values=".",
                 names="strain chrom pos end filt pe sr mapq precise homlen gt dr dv rr rv".split())
dl["alt"] = dl.dv.fillna(0) + dl.rv.fillna(0)
dl = dl[(dl.alt >= MIN_READS) & (dl.end - dl.pos >= 50)]
TOL = 50
jr = np.zeros(len(calls)); js = np.full(len(calls), -1); je = np.full(len(calls), -1)
for (s, ch), sub in calls.groupby(["strain", "chrom"]):
    dd = dl[(dl.strain == s) & (dl.chrom == ch)]
    for i, r in zip(sub.index, sub.itertuples()):
        m = dd[(dd.pos >= r.outer_start - TOL) & (dd.pos <= r.start + TOL) & (dd.end >= r.end - TOL) & (dd.end <= r.outer_end + TOL)]
        if len(m):
            b = m.sort_values("alt").iloc[-1]
            jr[i], js[i], je[i] = b.alt, b.pos, b.end
calls["delly_reads"], calls["delly_start"], calls["delly_end"] = jr, js, je
calls["junction_delly"] = calls.delly_reads >= MIN_READS
calls["junction_splitread"] = (calls.L_sa_other_edge.fillna(0) >= MIN_READS) | (calls.R_sa_other_edge.fillna(0) >= MIN_READS)
calls["junction"] = calls.junction_delly | calls.junction_splitread

# ------------------------------------------------------------------ merge regions across strains into loci
ev = np.full(len(calls), -1); nev = 0
for ch, sub in calls.groupby("chrom"):
    idx = sub.index.values; st = sub.start.values; en = sub.end.values
    parent = list(range(len(idx)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a in range(len(idx)):
        b = a + 1
        while b < len(idx) and st[b] < en[a]:
            ov = min(en[a], en[b]) - st[b]
            if ov >= 0.8 * max(en[a] - st[a], en[b] - st[b]):      # reciprocal overlap >= 80 %
                parent[find(b)] = find(a)
            b += 1
    roots = {}
    for a in range(len(idx)):
        r = find(a)
        if r not in roots:
            roots[r] = nev; nev += 1
        ev[idx[a]] = roots[r]
calls["event"] = ev


def share(x):       # share of carrier strains with read evidence
    return float((x != "").mean())


def modal(x):
    x = x[x != ""]
    return x.value_counts().index[0] if len(x) else ""


E = calls.groupby("event").agg(
    chrom=("chrom", "first"), start=("start", "median"), end=("end", "median"),
    outer_start=("outer_start", "median"), outer_end=("outer_end", "median"),
    n_strains=("strain", "nunique"), strains=("strain", lambda x: ",".join(sorted(set(x), key=cohort.index))),
    deleted_bases=("deleted_bases", "median"), n_fragments=("n_fragments", "median"),
    n_junction=("junction", "sum"),
    L_te_reads_family=("L_te_family_reads", modal), L_te_reads_share=("L_te_family_reads", share),
    R_te_reads_family=("R_te_family_reads", modal), R_te_reads_share=("R_te_family_reads", share)).reset_index()
for k in ["start", "end", "outer_start", "outer_end", "deleted_bases"]:
    E[k] = E[k].astype(int)
E["span"] = E.end - E.start
E["event_id"] = ["AB%04d" % (i + 1) for i in range(len(E))]
E["junction"] = E.n_junction > 0
calls = calls.merge(E[["event", "event_id"]], on="event")

# ------------------------------------------------------------------ TE context of the two breakpoints
# Evidence 1 (reference): an annotated TE in E2 touches the breakpoint zone (zone padded by 20 bp).
# Evidence 2 (reads): in at least half of the carrier strains, >= 3 soft-clipped reads at the edge continue into a TE.
PAD = 20
def zone_hits(df, c, s, e):
    return df[(df.chrom == c) & (df.end > s - PAD) & (df.start < e + PAD)]
ann = []
for r in E.itertuples():
    L = zone_hits(te, r.chrom, r.outer_start, r.start); Rr = zone_hits(te, r.chrom, r.end, r.outer_end)
    L2 = zone_hits(interspersed, r.chrom, r.outer_start, r.start); R2 = zone_hits(interspersed, r.chrom, r.end, r.outer_end)
    both = set(L.fam) & set(Rr.fam)
    same_or = any(set(L[L.fam == f].strand) & set(Rr[Rr.fam == f].strand) for f in both)
    ann.append((",".join(sorted(set(L.fam))), ",".join(sorted(set(Rr.fam))), len(L2) > 0, len(R2) > 0, ",".join(sorted(both)), same_or))
E["L_te_ref"], E["R_te_ref"], E["L_rm_ref"], E["R_rm_ref"], E["ref_same_family_both"], E["ref_same_orientation"] = zip(*ann)
E["L_te_reads"] = np.where(E.L_te_reads_share >= 0.5, E.L_te_reads_family, "")
E["R_te_reads"] = np.where(E.R_te_reads_share >= 0.5, E.R_te_reads_family, "")
E["L_te_any"] = (E.L_te_ref != "") | E.L_rm_ref | (E.L_te_reads != "")
E["R_te_any"] = (E.R_te_ref != "") | E.R_rm_ref | (E.R_te_reads != "")
E["L_family"] = np.where(E.L_te_reads != "", E.L_te_reads, E.L_te_ref)
E["R_family"] = np.where(E.R_te_reads != "", E.R_te_reads, E.R_te_ref)
samefam = np.array([bool((set(a.split(",")) & set(b.split(","))) - {""}) for a, b in zip(E.L_family, E.R_family)])
E["te_class"] = np.select([E.L_te_any & E.R_te_any & samefam, E.L_te_any & E.R_te_any, E.L_te_any | E.R_te_any],
                          ["TE at both breakpoints, same named family", "TE at both breakpoints, other", "TE at one breakpoint"],
                          "no TE evidence")
E["te_evidence"] = np.select(
    [((E.L_te_reads != "") | (E.R_te_reads != "")) & ((E.L_te_ref != "") | (E.R_te_ref != "") | E.L_rm_ref | E.R_rm_ref),
     (E.L_te_reads != "") | (E.R_te_reads != ""), E.L_te_any | E.R_te_any],
    ["reference annotation and reads", "reads only (TE not in E2 at this place)", "reference annotation only"], "none")
E["te_family"] = [",".join(sorted((set(a.split(",")) | set(b.split(","))) - {""})) for a, b in zip(E.L_family, E.R_family)]
E["repeat_frac_inside"] = [n_in(rep_cs, c, s, e) / max(e - s, 1) for c, s, e in E[["chrom", "start", "end"]].itertuples(index=False)]
E["telomeric"] = (E.outer_start < 20000) | (E.outer_end > E.chrom.map(sizes) - 20000)
E["size_class"] = pd.cut(E.deleted_bases, [0, 500, 1000, 5000, 20000, 1e9], labels=["100-500 bp", "0.5-1 kb", "1-5 kb", "5-20 kb", ">20 kb"])

# genes: counted when at least half of the gene lies inside the deleted interval
gl = []
for r in E.itertuples():
    g = genes[(genes.chrom == r.chrom) & (genes.end > r.start) & (genes.start < r.end)]
    ov = np.minimum(g.end, r.end) - np.maximum(g.start, r.start)
    g = g[ov >= 0.5 * (g.end - g.start)]
    gl.append((",".join(g.gene), ",".join(g[g.effector].gene), len(g)))
E["genes"], E["effector_genes"], E["n_genes"] = zip(*gl)

# ------------------------------------------------------------------ SNP-based candidates scored with all-base depth
geno = pd.read_csv(I + "DataS3_genotypes.tsv.gz", sep="\t", skiprows=1, usecols=["CHROM", "POS"])
geno = geno[geno.CHROM.isin(CHROMS)]
snp_pos = {c: np.sort(g.POS.values) for c, g in geno.groupby("CHROM")}
runs = pd.read_csv(R + "snp_gap_runs.tsv", sep="\t")
runs = runs[runs.chrom.isin(CHROMS)].copy()
runs["id"] = runs.chrom + ":" + runs.start.astype(str) + "-" + runs.end.astype(str)
runs = runs.merge(pd.read_csv(C + "snp_gap_runs.per_strain_depth.tsv.gz", sep="\t"), on=["id", "strain"], how="left")
runs["frac_deleted"] = runs.deleted_bases / runs.scorable_bases.replace(0, np.nan)
runs["allbase"] = np.select([runs.scorable_bases < 30, runs.frac_deleted >= 0.8, runs.frac_deleted >= 0.2],
                            ["not assessable (repeat / uncallable)", "confirmed deleted", "partly deleted"], "not deleted")
# A SNP-based call that overlaps an all-base deletion region of the same strain counts as at least partly deleted, also
# when the base count says "not deleted" or "not assessable": the count uses only the uniquely mappable bases inside the
# SNP interval, which can be few. The base-count result is kept in allbase_base_count.
runs = runs.reset_index(drop=True)
runs["allbase_base_count"] = runs.allbase
ov = np.zeros(len(runs), bool)
for (s, ch), sub in runs.groupby(["strain", "chrom"]):
    cc = calls[(calls.strain == s) & (calls.chrom == ch)]
    if cc.empty:
        continue
    a, b = cc.start.values, cc.end.values
    for i, r in zip(sub.index, sub.itertuples()):
        ov[i] = bool(((a < r.end) & (b > r.start - 1)).any())
runs["overlaps_allbase_call"] = ov
runs.loc[ov & runs.allbase.isin(["not deleted", "not assessable (repeat / uncallable)"]), "allbase"] = "partly deleted"
runs["flank_ok_either"] = (runs.dist_above < 3000) | (runs.dist_below < 3000)
runs["flank_ok_both"] = (runs.dist_above < 3000) & (runs.dist_below < 3000)
r5 = runs[runs.snp_count >= 5]

# ------------------------------------------------------------------ all-base regions: could / did the SNP method see them?
def count_snps(c, s, e):
    p = snp_pos[c]; return int(np.searchsorted(p, e, "right") - np.searchsorted(p, s + 1, "left"))
calls["snps_inside"] = [count_snps(c, s, e) for c, s, e in calls[["chrom", "start", "end"]].itertuples(index=False)]
hit = np.zeros(len(calls), bool); bs = np.full(len(calls), np.nan); be = np.full(len(calls), np.nan)
for (s, ch), sub in calls.groupby(["strain", "chrom"]):
    rr = r5[(r5.strain == s) & (r5.chrom == ch)]
    if rr.empty:
        continue
    a, b = rr.start.values - 1, rr.end.values
    for i, r in zip(sub.index, sub.itertuples()):
        o = (a < r.outer_end) & (b > r.outer_start)
        if o.any():
            hit[i] = True; bs[i] = a[o].min() - r.start; be[i] = b[o].max() - r.end
calls["found_by_snp_method"], calls["snp_start_minus_allbase_start"], calls["snp_end_minus_allbase_end"] = hit, bs, be
calls["snp_class"] = pd.cut(calls.snps_inside, [-1, 0, 4, 1e9], labels=["0 SNPs", "1-4 SNPs", ">=5 SNPs"])
calls["size_class"] = pd.cut(calls.deleted_bases, [0, 500, 1000, 5000, 20000, 1e9], labels=["100-500 bp", "0.5-1 kb", "1-5 kb", "5-20 kb", ">20 kb"])
E = E.merge(calls.groupby("event").agg(snps_inside=("snps_inside", "max"), found_by_snp_method=("found_by_snp_method", "any")),
            left_on="event", right_index=True)

# ------------------------------------------------------------------ Data S10 crosswalk
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")
p10 = pd.read_csv(C + "DataS10.per_strain_depth.tsv.gz", sep="\t")
S10_FAM = {"Grasshopper": "Grasshopper", "Pot2": "Pot2", "Fosbury": "Fosbury", "Pyret": "Pyret", "Inago2": "Inago2"}
rows = []
for r in s10.itertuples():
    sts = r.strains.split(",")
    x = p10[(p10.id == r.del_id) & p10.strain.isin(sts)]
    fr = x.deleted_bases / x.scorable_bases.replace(0, np.nan)
    ok = x.scorable_bases >= 30
    # as for the SNP-based calls: a listed strain whose all-base deletion overlaps the row counts as at least partly deleted
    xov = x.strain.isin(set(calls[(calls.chrom == r.chrom) & calls.strain.isin(sts) & (calls.start < r.end) & (calls.end > r.start - 1)].strain))
    full = (fr >= 0.8) & ok
    ndel, npart = int(full.sum()), int((~full & (((fr >= 0.2) & ok) | xov)).sum())
    nno, nun = int(((fr < 0.2) & ok & ~xov).sum()), int((~ok & ~xov).sum())
    bp, bn, bu = int(((fr >= 0.2) & (fr < 0.8) & ok).sum()), int(((fr < 0.2) & ok).sum()), int((~ok).sum())    # base count alone
    disp_base = ("not assessable (repeat / uncallable)" if bu == len(sts) else "supported" if ndel + bu == len(sts) else
                 "supported (partly deleted in some strains)" if bn == 0 else "not supported" if ndel + bp == 0 else "supported in some strains only")
    o = E[(E.chrom == r.chrom) & (E.outer_start < r.end) & (E.outer_end > r.start - 1)]
    oc = calls[(calls.chrom == r.chrom) & calls.strain.isin(sts) & (calls.outer_start < r.end) & (calls.outer_end > r.start - 1)]
    disp = ("not assessable (repeat / uncallable)" if nun == len(sts) else
            "supported" if ndel + nun == len(sts) else
            "supported (partly deleted in some strains)" if nno == 0 else
            "not supported" if ndel + npart == 0 else "supported in some strains only")
    allx = p10[(p10.id == r.del_id) & p10.strain.isin(cohort) & ~p10.strain.isin(sts)]
    extra = allx[(allx.deleted_bases / allx.scorable_bases.replace(0, np.nan) >= 0.8) & (allx.scorable_bases >= 30)].strain
    fams = sorted(set(",".join(o.te_family).split(",")) - {""})
    rows.append(dict(del_id=r.del_id, chrom=r.chrom, start=r.start, end=r.end, s10_te=r.te, s10_strains=r.strains,
                     n_listed=len(sts), n_deleted=ndel, n_partly=npart, n_not_deleted=nno, n_not_assessable=nun,
                     disposition=disp, disposition_base_count=disp_base, allbase_events=",".join(o.event_id),
                     allbase_te_class=";".join(sorted(set(o.te_class))), allbase_te_evidence=";".join(sorted(set(o.te_evidence))),
                     allbase_te_family=",".join(fams),
                     same_family_as_s10=(S10_FAM.get(r.te, r.te) in [x.replace("MAGGY_Fosbury", "Fosbury") for x in fams]) if r.te != "No_TE" else None,   # MAGGY and Fosbury are one family
                     junction_supported=bool(oc.junction.any()), extra_cohort_strains_deleted=",".join(extra),
                     allbase_start=int(o.start.min()) if len(o) else -1, allbase_end=int(o.end.max()) if len(o) else -1,
                     s10_min_size=r.min_size, s10_max_size=r.max_size))
X = pd.DataFrame(rows)
X["allbase_any_te"] = X.allbase_te_class.str.contains("TE at")

# ------------------------------------------------------------------ 100-kb window tracks
w = pd.read_csv(C + "windows100kb.bwa.uniq.tsv.gz", sep="\t")
w = w[w.strain.isin(cohort)]
W = w.groupby(["chrom", "start", "end"]).agg(scorable_bases=("scorable_bases", "first"),
                                             allbase_deleted_bp_sum=("deleted_bases_in_calls", "sum")).reset_index()
W["allbase_deleted_frac"] = W.allbase_deleted_bp_sum / (W.scorable_bases * len(cohort)).replace(0, np.nan)
W["n_snps"] = [count_snps(c, s, e) for c, s, e in W[["chrom", "start", "end"]].itertuples(index=False)]
W["snps_per_scorable_kb"] = W.n_snps / (W.scorable_bases / 1000).replace(0, np.nan)
W["repeat_frac"] = [n_in(rep_cs, c, s, e) / (e - s) for c, s, e in W[["chrom", "start", "end"]].itertuples(index=False)]
W["named_te_frac"] = [n_in(te_cs, c, s, e) / (e - s) for c, s, e in W[["chrom", "start", "end"]].itertuples(index=False)]
gmid = (genes.start + genes.end) // 2
W["n_genes"] = [int(((genes.chrom == c) & (gmid >= s) & (gmid < e)).sum()) for c, s, e in W[["chrom", "start", "end"]].itertuples(index=False)]
W["n_effectors"] = [int(((genes.chrom == c) & (gmid >= s) & (gmid < e) & genes.effector).sum()) for c, s, e in W[["chrom", "start", "end"]].itertuples(index=False)]
# SNP-based track: share of (SNP x strain) cells that fall inside a >=5-SNP missing run
snpdel = np.zeros(len(W)); key = {(c, s): i for i, (c, s) in enumerate(zip(W.chrom, W.start))}
for r in r5.itertuples():
    p = snp_pos[r.chrom]; q = p[(p >= r.start) & (p <= r.end)]
    for wv, n in zip(*np.unique((q // 100000) * 100000, return_counts=True)):
        snpdel[key[(r.chrom, int(wv))]] += n
W["snpbased_deleted_frac"] = snpdel / (W.n_snps * len(cohort)).replace(0, np.nan)
okw = (W.scorable_bases >= 20000)
tracks = ["allbase_deleted_frac", "snpbased_deleted_frac", "snps_per_scorable_kb", "repeat_frac", "named_te_frac"]
corr = W[okw][tracks].corr("spearman")
both0 = W[okw & W.snpbased_deleted_frac.notna()]
S["windows"] = dict(n_used=int(okw.sum()), n_with_snps=len(both0),
                    allbase_any=int((both0.allbase_deleted_frac > 0).sum()), snp_any=int((both0.snpbased_deleted_frac > 0).sum()),
                    both_any=int(((both0.allbase_deleted_frac > 0) & (both0.snpbased_deleted_frac > 0)).sum()),
                    allbase_only=int(((both0.allbase_deleted_frac > 0) & ~(both0.snpbased_deleted_frac > 0)).sum()),
                    snp_only=int((~(both0.allbase_deleted_frac > 0) & (both0.snpbased_deleted_frac > 0)).sum()),
                    windows_without_any_snp=int((okw & (W.n_snps == 0)).sum()),
                    allbase_deleted_bp_in_windows_without_snp=int(W[okw & (W.n_snps == 0)].allbase_deleted_bp_sum.sum()),
                    spearman=corr.round(3).to_dict())

# ------------------------------------------------------------------ is a reference TE at the breakpoints more frequent than expected?
# Null: breakpoints dropped at random scorable positions; a random position has a "TE at the breakpoint" if an
# annotated TE lies within PAD bp on the outer side - the same rule as for real edges (reference evidence only).
te_all = pd.concat([te[["chrom", "start", "end"]], interspersed[["chrom", "start", "end"]]])
te_m, _ = cover(te_all)
nearL, nearR = {}, {}
for c in CHROMS:
    cs = np.concatenate(([0], np.cumsum(te_m[c]))); pos = np.arange(sizes[c])
    nearL[c] = cs[pos] - cs[np.maximum(pos - PAD - 1, 0)] > 0
    nearR[c] = cs[np.minimum(pos + PAD + 2, sizes[c])] - cs[pos + 1] > 0
sc_pos = {c: np.flatnonzero(scor_m[c]) for c in CHROMS}
gl_L = np.concatenate([nearL[c][sc_pos[c]] for c in CHROMS]).mean(); gl_R = np.concatenate([nearR[c][sc_pos[c]] for c in CHROMS]).mean()
obs = E[~E.telomeric]
refL = (obs.L_te_ref != "") | obs.L_rm_ref; refR = (obs.R_te_ref != "") | obs.R_rm_ref
loc = []
for r in obs.itertuples():
    pair = []
    for edge, arr in ((r.start, nearL), (r.end - 1, nearR)):
        w0 = (edge // 100000) * 100000
        p = sc_pos[r.chrom]; p = p[(p >= w0) & (p < w0 + 100000)]
        pair.append(arr[r.chrom][p].mean() if len(p) else np.nan)
    loc.append(pair)
loc = np.array(loc)
S["te_enrichment_reference_annotation"] = dict(
    n_events_non_telomeric=len(obs),
    observed_te_at_either_breakpoint=round(float((refL | refR).mean()), 4),
    expected_random_genomewide=round(float(1 - (1 - gl_L) * (1 - gl_R)), 4),
    expected_random_same_100kb_windows=round(float(np.nanmean(1 - (1 - loc[:, 0]) * (1 - loc[:, 1]))), 4),
    observed_te_at_both_breakpoints=round(float((refL & refR).mean()), 4),
    expected_both_genomewide=round(float(gl_L * gl_R), 5),
    expected_both_same_100kb_windows=round(float(np.nanmean(loc[:, 0] * loc[:, 1])), 5))

# ------------------------------------------------------------------ effectors among deleted genes (hypergeometric tail)
allg = sorted(set(",".join(E.genes).split(",")) - {""})
k = int(genes[genes.gene.isin(allg)].effector.sum()); n = len(allg); K = int(genes.effector.sum()); N = len(genes)
pval = sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(k, min(K, n) + 1)) / math.comb(N, n)
S["genes"] = dict(genes_at_least_half_deleted=n, effectors=k, effector_share=round(k / max(n, 1), 4),
                  genome_effector_share=round(K / N, 4), hypergeometric_p=float(f"{pval:.3g}"))
for label, sub in (("events_with_TE_evidence", E[E.te_class != "no TE evidence"]), ("events_without_TE_evidence", E[E.te_class == "no TE evidence"])):
    g = sorted(set(",".join(sub.genes).split(",")) - {""}); kk = int(genes[genes.gene.isin(g)].effector.sum())
    S["genes"][label] = dict(genes=len(g), effectors=kk, effector_share=round(kk / max(len(g), 1), 4))

# ------------------------------------------------------------------ headline numbers
S["regions_strain_level"] = len(calls); S["events"] = len(E)
S["regions_per_strain"] = dict(median=float(calls.groupby("strain").size().median()), min=int(calls.groupby("strain").size().min()), max=int(calls.groupby("strain").size().max()))
S["deleted_scorable_bp_per_strain"] = dict(median=float(calls.groupby("strain").deleted_bases.sum().median()),
                                           min=int(calls.groupby("strain").deleted_bases.sum().min()), max=int(calls.groupby("strain").deleted_bases.sum().max()))
S["near_reference_fragments_per_strain"] = near.groupby("strain").size().to_dict()
S["events_private_to_one_strain"] = int((E.n_strains == 1).sum())
S["events_by_size"] = E.size_class.value_counts().sort_index().to_dict()
S["events_with_junction"] = int(E.junction.sum())
S["te_class_counts"] = E.te_class.value_counts().to_dict()
S["te_class_by_size"] = pd.crosstab(E.size_class, E.te_class).to_dict("index")
S["te_evidence_counts"] = E.te_evidence.value_counts().to_dict()
S["te_family_counts"] = E[E.te_family != ""].te_family.str.split(",").explode().value_counts().to_dict()
S["te_any_by_snp_detection_and_size"] = (E.assign(te=E.te_class != "no TE evidence").groupby(["size_class", "found_by_snp_method"], observed=True)
                                         .te.agg(["size", "mean"]).round(3).reset_index().to_dict("records"))
S["te_any_by_snp_detection"] = E.assign(te=E.te_class != "no TE evidence").groupby("found_by_snp_method").te.agg(["size", "mean"]).round(4).to_dict("index")
S["snp_method_recall"] = dict(
    regions=round(float(calls.found_by_snp_method.mean()), 4),
    deleted_bp=round(float(calls[calls.found_by_snp_method].deleted_bases.sum() / calls.deleted_bases.sum()), 4),
    events=round(float(E.found_by_snp_method.mean()), 4),
    by_snps_inside=calls.groupby("snp_class", observed=True).found_by_snp_method.agg(["size", "mean"]).round(4).to_dict("index"),
    by_size=calls.groupby("size_class", observed=True).found_by_snp_method.agg(["size", "mean"]).round(4).to_dict("index"))
S["snp_candidates_ge5"] = dict(n=len(r5), allbase=r5.allbase.value_counts().to_dict(),
                               allbase_base_count=r5.allbase_base_count.value_counts().to_dict())
f = r5[(r5.min_size > 100) & r5.flank_ok_either]
S["snp_candidates_ge5_min100_flank3kb"] = dict(n=len(f), allbase=f.allbase.value_counts().to_dict())
m = calls[calls.found_by_snp_method]
S["boundary_shift_bp_snp_vs_allbase"] = dict(
    start_abs_q25_q50_q75_q90=[float(v) for v in m.snp_start_minus_allbase_start.abs().quantile([.25, .5, .75, .9])],
    end_abs_q25_q50_q75_q90=[float(v) for v in m.snp_end_minus_allbase_end.abs().quantile([.25, .5, .75, .9])],
    share_snp_interval_shorter_at_start=round(float((m.snp_start_minus_allbase_start > 0).mean()), 3),
    share_snp_interval_shorter_at_end=round(float((m.snp_end_minus_allbase_end < 0).mean()), 3))
S["DataS10"] = dict(disposition=X.disposition.value_counts().to_dict(),
                    disposition_base_count=X.disposition_base_count.value_counts().to_dict(),
                    disposition_by_te=pd.crosstab(X.s10_te != "No_TE", X.disposition).rename(index={True: "TE", False: "No_TE"}).to_dict("index"),
                    junction_supported=int(X.junction_supported.sum()),
                    te_rows=int((X.s10_te != "No_TE").sum()),
                    te_rows_with_allbase_te=int(((X.s10_te != "No_TE") & X.allbase_any_te).sum()),
                    te_rows_same_family=int((X.same_family_as_s10 == True).sum()),
                    no_te_rows=int((X.s10_te == "No_TE").sum()),
                    no_te_rows_with_allbase_te=int(((X.s10_te == "No_TE") & X.allbase_any_te).sum()),
                    te_evidence_for_te_rows=X[X.s10_te != "No_TE"].allbase_te_evidence.value_counts().to_dict())

# reference-annotation TE share by size class (comparable across sensitivity runs)
_ref = (E.L_te_ref != "") | E.L_rm_ref | (E.R_te_ref != "") | E.R_rm_ref
S["ref_te_at_breakpoint_by_size"] = E.assign(te=_ref).groupby("size_class", observed=True).te.agg(["size", "mean"]).round(3).to_dict("index")
S["ref_te_at_breakpoint_events_ge1kb"] = round(float(_ref[E.deleted_bases >= 1000].mean()), 4)
S["ref_te_by_snp_detection"] = E.assign(te=_ref).groupby("found_by_snp_method").te.agg(["size", "mean"]).round(4).to_dict("index")

calls.drop(columns=["event"]).to_csv(OUT + f"allbase_calls_strain_level{SUFFIX}.tsv", sep="\t", index=False)
E.drop(columns=["event"]).to_csv(OUT + f"allbase_events{SUFFIX}.tsv", sep="\t", index=False)
if PRIMARY:
    frag.to_csv(R + f"allbase_fragments_strain_level{SUFFIX}.tsv", sep="\t", index=False)
    runs.to_csv(R + "snp_candidates_checked.tsv", sep="\t", index=False)
    X.to_csv(R + f"DataS10_crosswalk{SUFFIX}.tsv", sep="\t", index=False)
    W.to_csv(R + "windows100kb_tracks.tsv", sep="\t", index=False)
json.dump(S, open(OUT + f"summary{SUFFIX}.json", "w"), indent=1, default=str)
print(json.dumps(S, indent=1, default=str))
