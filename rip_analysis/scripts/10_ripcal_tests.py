#!/usr/bin/env python3
"""Significance tests for the RIPCaL counts, one row per group. Writes results/ripcal/ripcal_tests.tsv.

Two questions, each tested twice (RIPCaL's copy-level counts, and the event-level counts where a change shared by
several copies is counted once, because shared changes are not independent observations):

  direction  Are C-to-T changes (C>T + G>A) more common than T-to-C changes (T>C + A>G) than the base composition of the
             reference predicts? If a C is as likely to change as a T, the expected ratio is (C+G sites)/(A+T sites) of the
             reference, not 1. One-sided exact binomial test.
  site       Are the C-to-T changes spread over the four sites (C followed by A, T, C, G) in proportion to how often each
             site occurs in the reference? Chi-square goodness of fit, 3 degrees of freedom. This is the RIP-specific
             test: ordinary mutation can favour C-to-T, but it has no reason to prefer a C followed by T or A.

P values are adjusted over the groups by Benjamini-Hochberg. Only the standard library, numpy and pandas are needed.
"""
import math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
SITES = ["ct", "ca", "cc", "cg"]
PAIR = {"ct": ("CT->TT", "AG->AA"), "ca": ("CA->TA", "TG->TA"), "cc": ("CC->TC", "GG->GA"), "cg": ("CG->TG", "CG->CA")}
CONTEXT = {"CpT": "ct", "ApG": "ct", "CpA": "ca", "TpG": "ca", "CpC": "cc", "GpG": "cc", "CpG": "cg"}


def binom_upper(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p), summed in log space so that very small values do not underflow to 0."""
    if k <= 0:
        return 1.0
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p) + (n - i) * math.log1p(-p)
            for i in range(k, n + 1)]
    m = max(logs)
    return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in logs)) if m > -700 else 10 ** (m / math.log(10))


def chi2_sf_3df(x):
    """Upper tail of the chi-square distribution with 3 degrees of freedom (closed form)."""
    if x > 1400:
        return 0.0
    return math.erfc(math.sqrt(x / 2)) + math.sqrt(2 * x / math.pi) * math.exp(-x / 2)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def site_test(obs, sites):
    exp = np.asarray(sites, float) / sum(sites) * sum(obs)
    return float(((np.asarray(obs, float) - exp) ** 2 / exp).sum())


pc = pd.read_csv(RES / "ripcal" / "ripcal_per_copy.tsv", sep="\t")
pc = pc[pc.reference == "consensus"]
groups = pd.read_csv(RES / "ripcal" / "ripcal_vs_ours.tsv", sep="\t").group
fam = pd.read_csv(RES / "four_genomes" / "family_summary.tsv", sep="\t")
fam["group"] = fam.genome + "." + fam.family + "." + fam.subfamily
fam = fam.set_index("group")
ctx = pd.read_csv(RES / "four_genomes" / "context_table.tsv", sep="\t")
ctx["group"] = ctx.genome + "." + ctx.family + "." + ctx.subfamily
ctx["site"] = ctx.context.map(CONTEXT)
ctx = ctx.groupby(["group", "site"])[["consensus_sites", "events"]].sum()

rows = []
for g in groups:
    d = pc[pc.group == g]
    s = d.sum(numeric_only=True)
    m = d.iloc[0]
    gc, at = m["Model C Count"] + m["Model G Count"], m["Model A Count"] + m["Model T Count"]
    fwd, rev = s["C->T"] + s["G->A"], s["T->C"] + s["A->G"]
    f = fam.loc[g]
    e_fwd, e_rev = int(f.event_C_to_T + f.event_G_to_A), int(f.event_T_to_C + f.event_A_to_G)
    p0_e = f.consensus_CG_sites / (f.consensus_CG_sites + f.consensus_AT_sites)
    sites = [ctx.loc[(g, x), "consensus_sites"] for x in SITES]
    obs = [s[PAIR[x][0]] + s[PAIR[x][1]] for x in SITES]
    ev = [ctx.loc[(g, x), "events"] for x in SITES]
    r = dict(group=g, copies=len(d), reference_gc=gc / (gc + at), expected_ratio=gc / at,
             c_to_t=fwd, t_to_c=rev, ratio=fwd / rev, fold_over_expected=fwd / rev / (gc / at),
             p_direction=binom_upper(int(round(fwd)), int(round(fwd)) + int(round(rev)), gc / (gc + at)),
             event_c_to_t=e_fwd, event_t_to_c=e_rev, event_ratio=e_fwd / e_rev if e_rev else np.nan,
             event_fold_over_expected=(e_fwd / e_rev) / (p0_e / (1 - p0_e)) if e_rev else np.nan,
             p_direction_event=binom_upper(e_fwd, e_fwd + e_rev, p0_e),
             ti_tv=s["Transitions"] / s["Transversions"],
             chi2_site=site_test(obs, sites), chi2_site_event=site_test(ev, sites))
    for x, o, e, n in zip(SITES, obs, ev, sites):
        r[f"sites_{x}"] = n
        r[f"enrich_{x}"] = (o / sum(obs)) / (n / sum(sites))            # share of changes / share of sites; 1 = no preference
        r[f"event_enrich_{x}"] = (e / sum(ev)) / (n / sum(sites))
    r["p_site"], r["p_site_event"] = chi2_sf_3df(r["chi2_site"]), chi2_sf_3df(r["chi2_site_event"])
    rows.append(r)
T = pd.DataFrame(rows)
for c in ("p_direction", "p_direction_event", "p_site", "p_site_event"):
    T[c.replace("p_", "q_")] = bh(T[c])
T["excess_both_countings"] = (T.q_direction < 0.05) & (T.q_direction_event < 0.05)
T.to_csv(RES / "ripcal" / "ripcal_tests.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(T[["group", "copies", "expected_ratio", "ratio", "fold_over_expected", "q_direction", "event_ratio", "event_fold_over_expected",
         "q_direction_event", "ti_tv", "enrich_ct", "enrich_ca", "enrich_cc", "enrich_cg", "q_site", "q_site_event", "excess_both_countings"]]
      .sort_values("fold_over_expected", ascending=False).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
