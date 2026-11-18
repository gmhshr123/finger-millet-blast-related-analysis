#!/usr/bin/env python3
"""Summary table, workbook and figure for the RIP analysis (run from manuscripts_all/)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = "rip_analysis/results/four_genomes/"
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
f = pd.read_csv(R + "family_summary.tsv", sep="\t", keep_default_na=False, na_values=["", "nan"])
named = pd.read_csv(R + "consensus_vs_named.tsv", sep="\t", header=None, names="q s pident length qlen slen bits".split())
edta = pd.read_csv(R + "consensus_vs_edta.tsv", sep="\t", header=None, names="q s pident length qlen slen bits".split())


def identity(g, fam):
    q = f"{g}.{fam}.all_consensus"
    n = named[named.q == q]
    if len(n) and n.pident.iloc[0] >= 95 and n.length.iloc[0] >= 1500:
        name = n.s.iloc[0].split("#")[0].replace("MAGGY", "MAGGY/Fosbury")
        return f"{name} ({n.pident.iloc[0]:.1f}%)"
    e = edta[edta.q == q]
    lib = e.s.iloc[0].split("#")[0].replace("__TE_", ".TE").replace("_INT", "").replace("_LTR", "") if len(e) else ""
    if len(n) and n.length.iloc[0] >= 1500:
        return f"{n.s.iloc[0].split('#')[0]}-like ({n.pident.iloc[0]:.0f}%); library {lib}"
    return f"no published name; library {lib}"


f["identity"] = [identity(g, fam) for g, fam in zip(f.genome, f.family)]
# Pot2 sets hold two lineages. Identity of each subfamily consensus to GenBank Pot2 (Z33638), from BLASTN:
POT2 = {("E2", "sub1"): "Pot2 (97.2%)", ("MZ5", "sub1"): "Pot2 (97.2%)", ("B71", "sub2"): "Pot2 (98.5%)", ("PoP131", "sub1"): "Pot2 (99.9%)",
        ("B71", "sub1"): "second Pot2-like lineage (82% to Pot2)", ("PoP131", "sub2"): "second Pot2-like lineage (84% to Pot2)"}
f["identity"] = [POT2.get((g, s), i) if fam == "Pot2" and s != "all" else ("mixed: Pot2 and a second Pot2-like lineage" if fam == "Pot2" and g in ("B71", "PoP131") else i)
                 for g, fam, s, i in zip(f.genome, f.family, f.subfamily, f.identity)]
# rows to report: the subfamilies where a set was split, otherwise the whole set
split = f.groupby(["genome", "family"]).subfamily.transform(lambda s: (s != "all").any())
rep = f[(split & (f.subfamily != "all")) | (~split)].copy()
rep["group"] = np.where(rep.subfamily == "all", "whole set", rep.subfamily)
GEN = {"E2": "E2", "MZ5": "MZ5-1-6", "B71": "B71", "PoP131": "P131"}
rep["strain"] = rep.genome.map(GEN)
rep["order"] = rep.genome.map({"E2": 0, "MZ5": 1, "B71": 2, "PoP131": 3})
rep = rep.sort_values(["order", "family", "subfamily"])
cols = {"strain": "Strain", "family": "Family", "identity": "Identity of the consensus", "group": "Group", "n_copies": "Copies",
        "n_set": "Copies in the set", "consensus_length": "Consensus length (bp)", "median_divergence": "Median divergence from consensus",
        "event_substitutions": "Substitution events", "event_C_to_T": "C>T", "event_G_to_A": "G>A", "event_T_to_C": "T>C",
        "event_A_to_G": "A>G", "event_transversions": "Transversions", "event_rip_type_share": "RIP-type share (C>T + G>A)",
        "event_rip_vs_other_transitions": "RIP-type : other transitions", "event_ti_tv": "Transitions : transversions",
        "event_rip_per_CG_site": "RIP-type events per C/G site", "event_other_ti_per_AT_site": "Other transitions per A/T site",
        "event_rip_in_CpA_share": "Share of RIP-type events at CpA/TpG", "CpA_share_of_CG_sites": "CpA/TpG share of C/G sites",
        "event_rip_in_CpT_share": "Share of RIP-type events at CpT/ApG", "CpT_share_of_CG_sites": "CpT/ApG share of C/G sites",
        "singleton_substitutions": "Events seen in one copy only", "singleton_rip_type_share": "RIP-type share, one-copy events",
        "consensus_gc": "Consensus GC", "median_composite_index": "Median composite RIP index of copies",
        "copies_composite_above_0": "Copies with composite index > 0", "tsd_confirmed": "Flanking duplication confirmed"}
T = rep[list(cols)].rename(columns=cols)
T.to_csv("rip_analysis/results/RIP_summary.tsv", sep="\t", index=False, float_format="%.4f")
ctx = pd.read_csv(R + "context_table.tsv", sep="\t")
ctx = ctx.merge(rep[["genome", "family", "subfamily"]], on=["genome", "family", "subfamily"])
ctx = ctx[~ctx.context.str.contains("N")]
ctx["events_per_site"] = ctx.events / ctx.consensus_sites
copy = pd.read_csv(R + "copy_table.tsv", sep="\t").merge(rep[["genome", "family", "subfamily"]], on=["genome", "family", "subfamily"])
with pd.ExcelWriter("rip_analysis/results/RIP_summary.xlsx") as xw:
    T.to_excel(xw, sheet_name="Summary", index=False)
    ctx.to_excel(xw, sheet_name="Sequence context", index=False)
    copy.to_excel(xw, sheet_name="Per copy", index=False)
    f.drop(columns=["identity"]).to_excel(xw, sheet_name="All groups incl. unsplit", index=False)

short = T[["Strain", "Family", "Identity of the consensus", "Group", "Copies", "Median divergence from consensus", "Substitution events",
           "RIP-type share (C>T + G>A)", "RIP-type : other transitions", "Transitions : transversions",
           "Share of RIP-type events at CpA/TpG", "Share of RIP-type events at CpT/ApG", "Median composite RIP index of copies"]]
print(short.round(3).to_string(index=False))
print("\ncontext, events per available site (median over groups):")
print(ctx.groupby(["substitution", "context"]).events_per_site.median().round(3).to_string())
print("CpA share available vs observed, median:", T["CpA/TpG share of C/G sites"].median().round(3), T["Share of RIP-type events at CpA/TpG"].median().round(3),
      "| CpT:", T["CpT/ApG share of C/G sites"].median().round(3), T["Share of RIP-type events at CpT/ApG"].median().round(3))

# ---- figure: substitution spectrum per group (events)
lab = [f"{s}  {fam}  {g if g != 'whole set' else ''}  (n={n})" for s, fam, g, n in zip(T.Strain, T.Family, T.Group, T.Copies)]
tot = T["Substitution events"].values.astype(float)
parts = [("C>T", "#b8433f"), ("G>A", "#e08a5b"), ("T>C", "#2f6db0"), ("A>G", "#7fb0de"), ("Transversions", "#9aa1a9")]
fig, ax = plt.subplots(figsize=(10, 0.36 * len(T) + 1.6))
left = np.zeros(len(T)); y = np.arange(len(T))[::-1]
for col, color in parts:
    v = T[col].values / tot
    ax.barh(y, v, left=left, color=color, label=col, height=0.72, edgecolor="white", linewidth=0.6); left += v
ax.axvline(1 / 6, color="black", lw=1, ls="--"); ax.text(1 / 6, len(T) - 0.35, " expected RIP-type share\n if all changes were equally likely", fontsize=7.5, va="bottom")
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8); ax.set_xlim(0, 1); ax.set_xlabel("Share of substitution events against the family consensus")
for yi, n in zip(y, tot):
    ax.text(1.005, yi, f"{int(n)} events", fontsize=7.5, va="center")
ax.legend(ncol=5, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.09 - 0.6 / len(T)), frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Substitutions in full-length TE copies: C>T and G>A (RIP-type) against all other changes", fontsize=10, pad=22)
fig.tight_layout(); fig.savefig("rip_analysis/results/RIP_spectrum.png", dpi=200); fig.savefig("rip_analysis/results/RIP_spectrum.pdf")
