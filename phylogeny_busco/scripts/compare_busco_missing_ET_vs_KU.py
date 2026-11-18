#!/usr/bin/env python3
"""Compare BUSCO Missing calls between ET (E*, T*) and KU (K*, U*, UK*) assemblies (>=95% complete).

Group from sample-name prefix: E = Ethiopia, T = Tanzania, K = Kenya, U = Uganda, UK = (kept with KU;
sensitivity run without it). CKF* assemblies are excluded.
Tests:
  1. Per-assembly Missing count: Mann-Whitney U.
  2. Global marker-composition test: permute ET/KU labels among assemblies (10,000x). Statistic = sum over
     markers of the squared difference in Missing proportion between groups. Each assembly keeps its own
     Missing set, so this asks whether the group labels matter beyond chance.
     Run on raw Missing calls and on a version where each assembly's Missing row is divided by its total
     Missing count (composition only, removes overall-level differences).
  3. Per-marker Fisher exact test with BH FDR; Cochran-Mantel-Haenszel (stratified by Missing-count tercile)
     for markers missing in >=5 assemblies.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

PROJECT = Path(__file__).resolve().parents[1]
BASE = PROJECT / "results/busco_missingness_randomness"
OUT = BASE / "followup/ET_vs_KU"
OUT.mkdir(parents=True, exist_ok=True)
PERM, SEED = 10000, 20261002

status = pd.read_csv(BASE / "busco_status_matrix.tsv", sep="\t", index_col=0)
qc = pd.read_csv(BASE / "sample_qc_and_counts.tsv", sep="\t", index_col=0)
ann = pd.read_csv(BASE / "followup/marker_annotation_all.tsv", sep="\t", index_col=0)


def prefix(s):
    name = s.replace(".ragtag", "")
    for p in ("CKF", "UK", "E", "T", "K", "U"):
        if name.startswith(p):
            return p
    return "other"


qc["prefix"] = [prefix(s) for s in qc.index]
qc["group"] = qc["prefix"].map({"E": "ET", "T": "ET", "K": "KU", "U": "KU", "UK": "KU"})
good = qc[qc["complete_percent"] >= 95]
membership = good[["prefix", "group", "complete_percent", "missing_count"]]
membership.to_csv(OUT / "group_membership_ge95.tsv", sep="\t")


def global_perm(M, g, seed, normalise):
    X = M.astype(float)
    if normalise:
        X = X / X.sum(axis=1, keepdims=True)
    def stat(lab):
        return float(np.square(X[lab].mean(0) - X[~lab].mean(0)).sum())
    obs = stat(g)
    rng = np.random.default_rng(seed)
    null = np.array([stat(rng.permutation(g)) for _ in range(PERM)])
    return obs, null.mean(), (int((null >= obs).sum()) + 1) / (PERM + 1)


def run(include_uk: bool, tag: str):
    sub = good[good["group"].notna()]
    if not include_uk:
        sub = sub[sub["prefix"] != "UK"]
    M = (status.loc[sub.index] == "Missing").values
    g = (sub["group"] == "ET").values
    res = {"comparison": tag, "ET_n": int(g.sum()), "KU_n": int((~g).sum()),
           "ET_median_missing": float(np.median(M[g].sum(1))), "KU_median_missing": float(np.median(M[~g].sum(1))),
           "ET_median_complete_pct": float(sub.loc[g, "complete_percent"].median()),
           "KU_median_complete_pct": float(sub.loc[~g, "complete_percent"].median())}
    res["mannwhitney_p_missing_count"] = stats.mannwhitneyu(M[g].sum(1), M[~g].sum(1)).pvalue
    for norm, lab in [(False, "raw"), (True, "composition")]:
        o, nm, p = global_perm(M, g, SEED + (1 if norm else 0), norm)
        res[f"perm_{lab}_obs_over_null"] = round(o / nm, 3)
        res[f"perm_{lab}_p"] = p

    # per-marker
    a, b = M[g].sum(0), M[~g].sum(0)
    nE, nK = int(g.sum()), int((~g).sum())
    keep = (a + b) >= 1
    fisher_p = np.ones(M.shape[1])
    for j in np.where(keep)[0]:
        fisher_p[j] = stats.fisher_exact([[a[j], nE - a[j]], [b[j], nK - b[j]]]).pvalue
    q = np.ones_like(fisher_p)
    q[keep] = multipletests(fisher_p[keep], method="fdr_bh")[1]
    # CMH stratified by assembly Missing-count tercile
    k = M.sum(1)
    terc = pd.qcut(k, 3, labels=False, duplicates="drop")
    cmh_p = np.full(M.shape[1], np.nan)
    for j in np.where((a + b) >= 5)[0]:
        tables = []
        for t in np.unique(terc):
            s_ = terc == t
            ge, gk = g & s_, (~g) & s_
            tab = np.array([[M[ge, j].sum(), ge.sum() - M[ge, j].sum()],
                            [M[gk, j].sum(), gk.sum() - M[gk, j].sum()]], float)
            if tab.sum(1).min() > 0:
                tables.append(tab + 0.5)
        if tables:
            cmh_p[j] = StratifiedTable(tables).test_null_odds().pvalue
    tested = ~np.isnan(cmh_p)
    cmh_q = np.full_like(cmh_p, np.nan)
    cmh_q[tested] = multipletests(cmh_p[tested], method="fdr_bh")[1]
    mk = pd.DataFrame({
        "busco_id": status.columns, "description": ann.reindex(status.columns)["odb12_description"].values,
        "ET_missing": a, "ET_n": nE, "ET_pct": (100 * a / nE).round(1),
        "KU_missing": b, "KU_n": nK, "KU_pct": (100 * b / nK).round(1),
        "fisher_p": fisher_p, "fisher_q_BH": q, "cmh_tercile_p": cmh_p, "cmh_tercile_q_BH": cmh_q,
        "missing_pct_all180": ann.reindex(status.columns)["missing_frac_high_completeness_ge95pct"].values * 100,
    })
    mk = mk.sort_values("fisher_p")
    mk.to_csv(OUT / f"marker_ET_vs_KU_{tag}.tsv", sep="\t", index=False)
    res["markers_tested"] = int(keep.sum())
    res["markers_fisher_q_lt_0.05"] = int((q[keep] < 0.05).sum())
    res["markers_cmh_q_lt_0.05"] = int((cmh_q[tested] < 0.05).sum())
    return res, mk


rows = []
for inc, tag in [(True, "KU_with_UK"), (False, "KU_without_UK")]:
    r, mk = run(inc, tag)
    rows.append(r)
    if inc:
        mk_main = mk
pd.DataFrame(rows).to_csv(OUT / "ET_vs_KU_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).T.to_string())
print(membership.groupby(["group", "prefix"], dropna=False).size().to_string())
print(mk_main.head(25).to_string(index=False))
freq = mk_main[mk_main["missing_pct_all180"] >= 50].sort_values("missing_pct_all180", ascending=False)
print(freq[["busco_id", "description", "ET_pct", "KU_pct", "fisher_q_BH", "cmh_tercile_q_BH", "missing_pct_all180"]].to_string(index=False))
