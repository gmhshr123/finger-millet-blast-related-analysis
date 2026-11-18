#!/usr/bin/env python3
"""Follow-up checks on recurrent BUSCO missingness in the final RagTag cohort.

Builds on analyze_busco_missingness_randomness.py (reads its status matrix) and adds:
  1. E2 reference status and ODB12 descriptions for every marker.
  2. Influence of the six <50%-complete assemblies.
  3. Row-preserving permutation test after removing markers that are absent from the
     E2 reference itself or missing in >=90% of assemblies ("cohort-wide absent").
  4. Concentration of Missing calls in the most-recurrent markers vs the null.
  5. E2 chromosome / position / length associations with per-marker missing frequency.
  6. Ancestry check (Mantel test) on assemblies that are also tips of the SNP phylogeny.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT = Path(__file__).resolve().parents[1]
BASE = PROJECT / "results/busco_missingness_randomness"
OUT = BASE / "followup"
E2_TABLE = PROJECT / "results/busco_E2_fullgenome_ascomycota_odb12/E2_fullgenome_ascomycota_odb12/run_ascomycota_odb12/full_table.tsv"
SCAFFOLD_STATS = PROJECT / "results/manifests_and_qc/final_ragtag_fasta_scaffold_stats.tsv"
SUMMARY = PROJECT / "results/manifests_and_qc/busco_final_ragtag_summary.tsv"
TREE = PROJECT / "results/iqtree_rooted/MAFge003_missing30/finger_millet.MAFge003_missing30.rooted_CKF3589.treefile"
NEAR_UNIVERSAL = 0.90


# ---------- helpers ----------
def row_preserving_null(matrix: np.ndarray, permutations: int, seed: int):
    """Shuffle each sample's Missing calls among markers; return per-permutation column counts summary."""
    n, m = matrix.shape
    k = matrix.sum(axis=1).astype(int)
    rng = np.random.default_rng(seed)
    expected = k.sum() / m
    obs_counts = matrix.sum(axis=0)
    top = max(1, int(round(0.01 * m)))
    obs_T = float(np.square(obs_counts - expected).sum())
    obs_top_share = float(np.sort(obs_counts)[::-1][:top].sum() / max(k.sum(), 1))
    T = np.empty(permutations)
    top_share = np.empty(permutations)
    max_count = np.empty(permutations)
    null_counts_sum = np.zeros(m)
    for it in range(permutations):
        counts = np.zeros(m, dtype=np.int32)
        for ki in k:
            if ki:
                counts[rng.choice(m, size=int(ki), replace=False)] += 1
        T[it] = np.square(counts - expected).sum()
        s = np.sort(counts)[::-1]
        top_share[it] = s[:top].sum() / max(k.sum(), 1)
        max_count[it] = s[0]
        null_counts_sum += counts
    return {
        "samples": n,
        "markers": m,
        "total_missing_calls": int(k.sum()),
        "expected_missing_per_marker": round(expected, 4),
        "observed_T": round(obs_T, 2),
        "null_mean_T": round(float(T.mean()), 2),
        "null_99pct_T": round(float(np.quantile(T, 0.99)), 2),
        "T_ratio_obs_over_null": round(obs_T / T.mean(), 2),
        "empirical_p": (int((T >= obs_T).sum()) + 1) / (permutations + 1),
        "observed_max_samples_sharing_one_missing_marker": int(obs_counts.max()),
        "null_mean_max_sharing": round(float(max_count.mean()), 2),
        "null_99pct_max_sharing": float(np.quantile(max_count, 0.99)),
        "top1pct_markers": top,
        "observed_share_of_missing_calls_in_top1pct_markers": round(obs_top_share, 4),
        "null_mean_share_in_top1pct": round(float(top_share.mean()), 4),
        "permutations": permutations,
        "seed": seed,
    }


def parse_newick(text: str):
    """Minimal Newick parser -> dict tip -> list of (node_id) path with branch lengths."""
    text = text.strip().rstrip(";")
    pos = 0
    parent: dict[int, int] = {}
    blen: dict[int, float] = {}
    tips: dict[str, int] = {}
    counter = [0]

    def new_node():
        counter[0] += 1
        return counter[0]

    def read_label():
        nonlocal pos
        start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        return text[start:pos]

    def read_len():
        nonlocal pos
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",();":
                pos += 1
            return float(text[start:pos])
        return 0.0

    def node(par):
        nonlocal pos
        nid = new_node()
        if par is not None:
            parent[nid] = par
        if text[pos] == "(":
            pos += 1
            node(nid)
            while text[pos] == ",":
                pos += 1
                node(nid)
            assert text[pos] == ")"
            pos += 1
            read_label()
        else:
            tips[read_label()] = nid
        blen[nid] = read_len()
        return nid

    node(None)
    return parent, blen, tips


def patristic(parent, blen, tips, names):
    def path(n):
        d, out = 0.0, {}
        while True:
            out[n] = d
            if n not in parent:
                return out
            d += blen[n]
            n = parent[n]

    paths = {nm: path(tips[nm]) for nm in names}
    D = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        pa = paths[a]
        for j in range(i + 1, len(names)):
            pb = paths[names[j]]
            D[i, j] = D[j, i] = min(pa[x] + pb[x] for x in pa.keys() & pb.keys())
    return D


def mantel(A, B, permutations, seed):
    iu = np.triu_indices_from(A, 1)
    a = stats.rankdata(A[iu])
    r_obs = stats.spearmanr(A[iu], B[iu])[0]
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    hits = 0
    for _ in range(permutations):
        p = rng.permutation(n)
        r = np.corrcoef(a, stats.rankdata(B[p][:, p][iu]))[0, 1]
        hits += r >= r_obs
    return float(r_obs), (hits + 1) / (permutations + 1)


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261002)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    status = pd.read_csv(BASE / "busco_status_matrix.tsv", sep="\t", index_col=0)
    qc = pd.read_csv(BASE / "sample_qc_and_counts.tsv", sep="\t", index_col=0).loc[status.index]
    summ = pd.read_csv(SUMMARY, sep="\t", index_col=0).loc[status.index]
    markers = list(status.columns)
    M = (status.values == "Missing").astype(np.uint8)
    F = (status.values == "Fragmented")

    # E2 reference annotation (one row per marker; keep first non-Missing hit)
    e2 = pd.read_csv(E2_TABLE, sep="\t", comment=None, skiprows=2, dtype=str)
    e2.columns = [c.lstrip("# ").strip() for c in e2.columns]
    e2 = e2.drop_duplicates("Busco id").set_index("Busco id").reindex(markers)
    assert e2["Status"].notna().all(), "E2 table does not cover all markers"
    desc = pd.read_csv(OUT / "ascomycota_odb12_links_to_ODB12.txt", sep="\t", header=None, index_col=0)[1]
    anc_len = pd.read_csv(OUT / "ascomycota_odb12_ancestral_protein_lengths.tsv", sep="\t", header=None, index_col=0)[1]

    lt50 = qc.index[qc["complete_percent"] < 50].tolist()
    ge95 = qc.index[qc["complete_percent"] >= 95].tolist()
    cohorts = {
        "all_250": list(status.index),
        "all_minus_6_lt50pct": [s for s in status.index if s not in lt50],
        "shared_pre216_216": qc.index[qc["dataset_shared216"]].tolist(),
        "high_completeness_ge95pct": ge95,
    }
    idx = {s: i for i, s in enumerate(status.index)}

    # ---- marker annotation table ----
    ann = pd.DataFrame(index=markers)
    ann["odb12_description"] = desc.reindex(markers).values
    ann["ancestral_protein_length_aa"] = anc_len.reindex(markers).values
    ann["E2_reference_status"] = e2["Status"].values
    ann["E2_sequence"] = e2["Sequence"].values
    ann["E2_gene_start"] = pd.to_numeric(e2["Gene Start"]).values
    ann["E2_gene_end"] = pd.to_numeric(e2["Gene End"]).values
    ann["E2_reference_hit_length"] = pd.to_numeric(e2["Length"]).values
    for name, members in cohorts.items():
        rows = [idx[s] for s in members]
        ann[f"missing_{name}"] = M[rows].sum(axis=0)
        ann[f"missing_frac_{name}"] = (M[rows].mean(axis=0)).round(4)
    ann["fragmented_all_250"] = F.sum(axis=0)
    ann["missing_in_the_6_lt50pct"] = M[[idx[s] for s in lt50]].sum(axis=0)
    ann["E2_ragtag_assembly_status"] = status.loc["E2.ragtag"].values if "E2.ragtag" in status.index else ""
    ann["class"] = np.where(
        ann["E2_reference_status"] == "Missing", "absent_in_E2_reference",
        np.where(ann["missing_frac_all_250"] >= NEAR_UNIVERSAL, "missing_in_ge90pct_assemblies", "variable"),
    )
    ann.index.name = "busco_id"
    ann.sort_values(["missing_all_250", "missing_high_completeness_ge95pct"], ascending=False).to_csv(
        OUT / "marker_annotation_all.tsv", sep="\t")
    top = ann.sort_values(["missing_high_completeness_ge95pct", "missing_all_250"], ascending=False).head(60)
    top.to_csv(OUT / "top_recurrent_markers_annotated.tsv", sep="\t")

    # ---- low-quality influence ----
    low = qc.loc[lt50, ["complete_percent", "missing_count", "fragmented_count"]].copy()
    low["share_of_all_missing_calls"] = (low["missing_count"] / M.sum()).round(4)
    low.to_csv(OUT / "low_completeness_assemblies.tsv", sep="\t")

    # ---- permutation tests on full and reduced marker sets ----
    variable = ann.index[ann["class"] == "variable"].tolist()
    col = {b: j for j, b in enumerate(markers)}
    marker_sets = {"all_2826_markers": markers, "variable_markers_only": variable}
    perm_rows = []
    offset = 10
    for cname, members in cohorts.items():
        rows = [idx[s] for s in members]
        for mname, mlist in marker_sets.items():
            sub = M[np.ix_(rows, [col[b] for b in mlist])]
            res = row_preserving_null(sub, args.permutations, args.seed + offset)
            offset += 1
            perm_rows.append({"cohort": cname, "marker_set": mname, **res})
            print(cname, mname, res["T_ratio_obs_over_null"], res["empirical_p"], flush=True)
    perm = pd.DataFrame(perm_rows)
    perm.to_csv(OUT / "permutation_tests_by_cohort_and_marker_set.tsv", sep="\t", index=False)

    # ---- recurrence distribution among high-completeness assemblies (variable markers) ----
    hc_rows = [idx[s] for s in ge95]
    hc_counts = M[np.ix_(hc_rows, [col[b] for b in variable])].sum(axis=0)
    dist = pd.Series(hc_counts).value_counts().sort_index()
    rng = np.random.default_rng(args.seed + 99)
    k = M[np.ix_(hc_rows, [col[b] for b in variable])].sum(axis=1)
    null_hist = np.zeros(len(ge95) + 1)
    reps = 1000
    for _ in range(reps):
        c = np.zeros(len(variable), dtype=int)
        for ki in k:
            c[rng.choice(len(variable), ki, replace=False)] += 1
        null_hist += np.bincount(c, minlength=len(ge95) + 1)
    null_hist /= reps
    bins = [(0, 0), (1, 1), (2, 2), (3, 4), (5, 9), (10, 19), (20, 49), (50, 10**6)]
    recur = []
    for lo, hi in bins:
        recur.append({
            "assemblies_missing_marker": f"{lo}" if lo == hi else (f"{lo}+" if hi > 1000 else f"{lo}-{hi}"),
            "observed_markers": int(dist[(dist.index >= lo) & (dist.index <= hi)].sum()),
            "expected_markers_under_null": round(float(null_hist[lo:min(hi, len(ge95)) + 1].sum()), 1),
        })
    pd.DataFrame(recur).to_csv(OUT / "recurrence_histogram_ge95_variable_markers.tsv", sep="\t", index=False)

    # ---- position / length associations (high-completeness, variable markers with E2 coordinates) ----
    pos = ann.loc[variable].copy()
    pos = pos[pos["E2_sequence"].notna()]
    pos["lo"] = pos[["E2_gene_start", "E2_gene_end"]].min(axis=1)
    pos["hi"] = pos[["E2_gene_start", "E2_gene_end"]].max(axis=1)
    chrom_len = pos.groupby("E2_sequence")["hi"].max()
    pos["dist_to_nearest_end_bp"] = np.minimum(pos["lo"], chrom_len.reindex(pos["E2_sequence"]).values - pos["hi"])
    y = pos["missing_high_completeness_ge95pct"].values
    assoc = []
    main_chr = pos["E2_sequence"].value_counts()
    keep = main_chr.index[main_chr >= 20]
    grp = [y[pos["E2_sequence"].values == c] for c in keep]
    kw = stats.kruskal(*grp)
    by_chr = pos[pos["E2_sequence"].isin(keep)].groupby("E2_sequence")["missing_high_completeness_ge95pct"].agg(
        markers="size", mean_missing="mean", markers_missing_ge5=lambda v: int((v >= 5).sum()))
    by_chr.to_csv(OUT / "missing_by_E2_chromosome_ge95.tsv", sep="\t")
    assoc.append({"test": "Kruskal-Wallis: missing count ~ E2 chromosome", "n_markers": int(sum(len(g) for g in grp)),
                  "statistic": round(kw.statistic, 3), "p_value": kw.pvalue})
    for var in ["dist_to_nearest_end_bp", "ancestral_protein_length_aa", "E2_reference_hit_length"]:
        r = stats.spearmanr(pos[var], y)
        assoc.append({"test": f"Spearman: missing count ~ {var}", "n_markers": len(pos),
                      "statistic": round(r.statistic, 4), "p_value": r.pvalue})
    # neighbour similarity along E2 chromosomes vs shuffling counts within chromosome
    pos = pos.sort_values(["E2_sequence", "lo"])
    def neighbour_r(vals, seqs):
        same = seqs[1:] == seqs[:-1]
        return stats.spearmanr(vals[1:][same], vals[:-1][same]).statistic
    v = pos["missing_high_completeness_ge95pct"].values.astype(float)
    sq = pos["E2_sequence"].values
    r_obs = neighbour_r(v, sq)
    rng = np.random.default_rng(args.seed + 7)
    null_r = []
    for _ in range(2000):
        vv = v.copy()
        for c in np.unique(sq):
            m_ = sq == c
            vv[m_] = rng.permutation(vv[m_])
        null_r.append(neighbour_r(vv, sq))
    null_r = np.array(null_r)
    assoc.append({"test": "Neighbouring-marker Spearman r along E2 chromosomes (shuffle within chromosome, 2000x)",
                  "n_markers": len(pos), "statistic": round(float(r_obs), 4),
                  "p_value": (int((null_r >= r_obs).sum()) + 1) / 2001})
    # assembly-level (secondary)
    scaf = pd.read_csv(SCAFFOLD_STATS, sep="\t", index_col=0)
    samp = qc.loc[ge95].copy()
    samp["scaffold_n50_bp"] = scaf.reindex([s.replace(".ragtag", "") for s in samp.index])["scaffold_n50_bp"].values
    samp["scaffold_count"] = summ.loc[ge95, "scaffold_count"].values
    samp["percent_gaps"] = summ.loc[ge95, "percent_gaps"].values
    for var in ["scaffold_n50_bp", "scaffold_count", "percent_gaps"]:
        r = stats.spearmanr(samp[var], samp["missing_count"], nan_policy="omit")
        assoc.append({"test": f"Spearman across ge95 assemblies: missing count ~ {var}", "n_markers": len(samp),
                      "statistic": round(r.statistic, 4), "p_value": r.pvalue})
    pd.DataFrame(assoc).to_csv(OUT / "position_length_assembly_associations.tsv", sep="\t", index=False)

    # ---- ancestry check (Mantel) ----
    parent, blen, tips = parse_newick(TREE.read_text())
    in_tree = [s for s in cohorts["all_minus_6_lt50pct"] if s.replace(".ragtag", "") in tips]
    anc = {"tree": str(TREE.relative_to(PROJECT)), "assemblies_in_tree_excluding_lt50pct": len(in_tree),
           "of_which_ge95pct": len([s for s in in_tree if s in ge95])}
    if len(in_tree) >= 10:
        names = [s.replace(".ragtag", "") for s in in_tree]
        D = patristic(parent, blen, tips, names)
        sub = M[np.ix_([idx[s] for s in in_tree], [col[b] for b in variable])].astype(float)
        kk = sub.sum(axis=1)
        shared = sub @ sub.T
        excess = shared - np.outer(kk, kk) / sub.shape[1]
        # similarity should rise as phylogenetic distance falls -> test negative distance
        r, p = mantel(-D, excess, args.permutations, args.seed + 55)
        anc.update({"statistic": "Spearman Mantel r between excess shared Missing calls and (negative) patristic distance",
                    "r": round(r, 4), "one_sided_p": p, "markers": len(variable), "permutations": args.permutations})
    (OUT / "ancestry_mantel_test.json").write_text(json.dumps(anc, indent=2) + "\n")

    meta = {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
            "pandas": pd.__version__, "permutations": args.permutations, "seed": args.seed,
            "near_universal_threshold": NEAR_UNIVERSAL, "lt50_samples": lt50,
            "n_variable_markers": len(variable),
            "class_counts": ann["class"].value_counts().to_dict(),
            "command": " ".join(sys.argv)}
    (OUT / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta["class_counts"]), json.dumps(anc), flush=True)


if __name__ == "__main__":
    main()
