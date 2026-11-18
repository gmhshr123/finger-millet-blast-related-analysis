#!/usr/bin/env python3
"""Data S10 check, and genome distribution and size of two deletion groups (run from manuscripts_all/, after scripts 10, 16 and 29).

1. Data S10: is every listed deletion (each Data S10 row x each listed strain) also called by the all-base method?
   A listing is "found" when an all-base deletion region of the same strain overlaps the Data S10 interval. Otherwise
   the read-depth test of the interval (uniquely mappable bases, share with low depth; script 10) gives the reason.
2. Genome distribution of two groups of deletions, with the labels of the breakpoint browser:
     A  read evidence "simple deletion" and "flanks joined" (opposite flank, direct junction) at one or both ends
     B  any read evidence except "DNA still present" (not a deletion), with a TE at one or both ends
   Plotted per deletion event (events shared by strains are counted once; an event is in a group when at least one
   strain's region is). Events in both groups are drawn in both and left out of the tests.
3. Size of both groups: median, mean and histogram (event size = region length on E2).

Input:  deletion_allbase/inputs/DataS10_events.tsv
        deletion_allbase/results_edta/allbase_calls_strain_level_read_evidence.tsv, allbase_events_read_evidence.tsv,
        regions_for_browser.tsv, cluster_uniq/DataS10.per_strain_depth.tsv.gz, cluster_uniq/E2.chr.genome
Output: deletion_allbase/DataS10_check_and_deletion_tracks/
        Fig1_genome_distribution.png/.pdf, Fig2_distance_to_chromosome_end.png/.pdf, Fig3_size_distribution.png/.pdf
        Deletion_groups_tables.xlsx, summary.md
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats
from openpyxl.styles import Font

R = "deletion_allbase/results_edta/"
I = "deletion_allbase/inputs/"
OUT = "deletion_allbase/DataS10_check_and_deletion_tracks/"
os.makedirs(OUT, exist_ok=True)

# colours: categorical slots 1 and 2 of the dataviz reference palette; text stays in ink colours
CA, CB, INK, INK2, GRID, CHR = "#2a78d6", "#eb6834", "#1d1f23", "#5b6068", "#e3e1db", "#c9c6bd"
LA, LB = "Simple deletion, flanks joined", "TE at one or both ends"
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": CHR, "svg.fonttype": "none", "pdf.fonttype": 42})

sizes = dict(l.split() for l in open(R + "cluster_uniq/E2.chr.genome"))
sizes = {k: int(v) for k, v in sizes.items()}
CHROMS = sorted(sizes)

# ------------------------------------------------------------------ 1. Data S10
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")
calls = pd.read_csv(R + "allbase_calls_strain_level_read_evidence.tsv", sep="\t")
p10 = pd.read_csv(R + "cluster_uniq/DataS10.per_strain_depth.tsv.gz", sep="\t")
FOUND, EDGE = "found: overlaps an all-base deletion of this strain", "found at the edge: in the breakpoint zone of an all-base deletion of this strain"
UNTEST, SHORT = "not testable: fewer than 30 uniquely mappable bases", "partly low depth, but fewer than 100 low-depth bases (the minimum all-base call)"
NORMAL = "not found: read depth is normal (under 20% of the bases with low depth)"
rows = []
for r in s10.itertuples():
    for s in r.strains.split(","):
        c = calls[(calls.strain == s) & (calls.chrom == r.chrom)]
        core = c[(c.start < r.end) & (c.end > r.start - 1)]                     # Data S10 positions are 1-based
        zone = c[(c.outer_start < r.end) & (c.outer_end > r.start - 1)]          # region plus its breakpoint zones
        x = p10[(p10.id == r.del_id) & (p10.strain == s)].iloc[0]
        sc, de = int(x.scorable_bases), int(x.deleted_bases)
        frac = de / sc if sc else np.nan
        if len(core):
            res = FOUND
        elif len(zone):
            res = EDGE
        elif sc < 30:
            res = UNTEST
        elif frac >= 0.2:
            assert de < 100, (r.del_id, s)
            res = SHORT
        else:
            res = NORMAL
        rows.append({"Data S10 ID": r.del_id, "Strain": s, "Chromosome": r.chrom, "Start (Data S10)": r.start, "End (Data S10)": r.end,
                     "Length (bp)": r.end - r.start + 1, "Data S10 TE": r.te, "Result": res,
                     "All-base events": ",".join((core if len(core) else zone).event_id),
                     "Uniquely mappable bases": sc, "Low-depth bases": de, "Share with low depth": round(frac, 3) if sc else "",
                     "Relative depth (unique reads)": round(float(x.rel_depth_unique_scorable), 3) if pd.notna(x.rel_depth_unique_scorable) else ""})
L10 = pd.DataFrame(rows)
order = [FOUND, EDGE, UNTEST, SHORT, NORMAL]
lst_counts = L10.Result.value_counts().reindex(order, fill_value=0)


def per_deletion(g):
    f = g.Result.isin([FOUND, EDGE])
    t = ~g.Result.eq(UNTEST)
    return pd.Series({"Listed strains": len(g), "Found": int(f.sum()), "Not testable": int((~t).sum()),
                      "Testable, not found": int((t & ~f).sum()),
                      "Status": "found in all listed strains" if f.all() else
                                "found in all testable strains" if (f | ~t).all() and f.any() else
                                "not testable in any strain" if (~t).all() else
                                "found in some listed strains" if f.any() else "not found in any testable strain"})


D10 = L10.groupby("Data S10 ID", sort=False).apply(per_deletion, include_groups=False).reset_index()
D10 = s10[["del_id", "chrom", "start", "end", "te"]].rename(columns={"del_id": "Data S10 ID", "chrom": "Chromosome", "start": "Start",
                                                                   "end": "End", "te": "Data S10 TE"}).merge(D10, on="Data S10 ID")
del_counts = D10.Status.value_counts()

# ------------------------------------------------------------------ 2-3. the two groups
g = pd.read_csv(R + "regions_for_browser.tsv", sep="\t", keep_default_na=False)
DIR, TE = "opposite flank (direct deletion junction)", "transposable element (TE)"
g["A"] = (g["Read evidence"] == "simple deletion") & ((g["Left origin"] == DIR) | (g["Right origin"] == DIR))
g["B"] = (g["Read evidence"] != "DNA still present") & ((g["Left origin"] == TE) | (g["Right origin"] == TE))
n_b_dna = int(((g["Left origin"] == TE) | (g["Right origin"] == TE)).sum() - g.B.sum())
ev = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t").set_index("event_id")
E = g.groupby("Deletion event ID").agg(A=("A", "any"), B=("B", "any"), strains_A=("A", "sum"), strains_B=("B", "sum"))
E = E[E.A | E.B].join(ev[["chrom", "start", "end", "n_strains"]])
E["size"] = E.end - E.start
E["mid"] = (E.start + E.end) / 2
E["dist_end"] = np.minimum(E.mid, E.chrom.map(sizes) - E.mid)
E["group"] = np.select([E.A & E.B, E.A], ["both", "A"], "B")
EA, EB = E[E.A], E[E.B]
TA, TB = E[E.group == "A"], E[E.group == "B"]          # for the tests: events in one group only

# tests
chi = stats.chi2_contingency(pd.crosstab(pd.concat([TA.chrom, TB.chrom]), ["A"] * len(TA) + ["B"] * len(TB)))
mw_dist = stats.mannwhitneyu(TA.dist_end, TB.dist_end)
ks_dist = stats.ks_2samp(TA.dist_end, TB.dist_end)
mw_size = stats.mannwhitneyu(TA["size"], TB["size"])
RA, RB = g[g.A], g[g.B]

# ------------------------------------------------------------------ Fig 1: positions on the chromosomes
def fig1(name, note, true_size):
    """true_size False: marks at least 12 kb wide, centred on the event. True: each mark spans the event's real start
    and end; a hairline outline keeps events shorter than a pixel (about 4.6 kb at 200 dpi) visible as thin lines."""
    fig, ax = plt.subplots(figsize=(11, 6.2))
    for i, c in enumerate(CHROMS):
        y = -i
        ax.add_patch(Rectangle((0, y - 0.035), sizes[c] / 1e6, 0.07, color=CHR, lw=0))
        for sub, col, sgn in ((EA, CA, 1), (EB, CB, -1)):
            for r in sub[sub.chrom == c].itertuples():
                yy = y + 0.07 if sgn > 0 else y - 0.07 - 0.26
                if true_size:
                    ax.add_patch(Rectangle((r.start / 1e6, yy), r.size / 1e6, 0.26, facecolor=col, edgecolor=col, lw=0.4, alpha=.85))
                else:
                    w = max(r.size, 12000) / 1e6
                    ax.add_patch(Rectangle((r.mid / 1e6 - w / 2, yy), w, 0.26, color=col, lw=0, alpha=.85))
        na, nb = int((EA.chrom == c).sum()), int((EB.chrom == c).sum())
        ax.text(8.75, y + 0.2, f"{na}", ha="right", va="center", fontsize=9, color=INK2)
        ax.text(8.75, y - 0.2, f"{nb}", ha="right", va="center", fontsize=9, color=INK2)
        ax.text(-0.12, y, c, ha="right", va="center", fontsize=10, color=INK)
    ax.text(8.75, 0.62, "events", ha="right", va="center", fontsize=9, color=INK2)
    ax.set_xlim(0, 8.8); ax.set_ylim(-len(CHROMS) + 0.45, 0.75)
    ax.set_yticks([]); ax.set_xlabel("Position on the E2 chromosome (Mb)")
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=CA), Rectangle((0, 0), 1, 1, color=CB)],
              labels=[f"{LA} (above the line; {len(EA)} events)", f"{LB} (below the line; {len(EB)} events)"],
              loc="upper center", bbox_to_anchor=(0.45, 1.09), ncol=2, frameon=False, handlelength=1.2)
    fig.text(0.01, 0.005, note, fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(OUT + f"{name}.{ext}", dpi=200)
    plt.close(fig)


fig1("Fig1_genome_distribution", "Each mark is one deletion event.", False)
fig1("Fig1b_genome_distribution_true_size", "Each mark is one deletion event, drawn at its true length on E2. "
     "Events shorter than about 5 kb appear as thin lines of the minimum visible width.", True)

# ------------------------------------------------------------------ Fig 2: distance to the nearest chromosome end
bins = [0, 100e3, 500e3, 1e6, 2e6, 5e6]
blab = ["< 100 kb", "100-500 kb", "0.5-1 Mb", "1-2 Mb", "> 2 Mb"]
pa = pd.cut(EA.dist_end, bins, labels=blab).value_counts(normalize=True).reindex(blab) * 100
pb = pd.cut(EB.dist_end, bins, labels=blab).value_counts(normalize=True).reindex(blab) * 100
fig, ax = plt.subplots(figsize=(7.5, 4.2))
x = np.arange(len(blab)); w = 0.38
for off, p, col, lab, n in ((-w / 2 - 0.01, pa, CA, LA, len(EA)), (w / 2 + 0.01, pb, CB, LB, len(EB))):
    ax.bar(x + off, p.values, w, color=col, label=f"{lab} ({n} events)")
    for xi, v in zip(x + off, p.values):
        ax.text(xi, v + 0.8, f"{v:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK2)
ax.set_xticks(x, blab); ax.set_xlabel("Distance from the event to the nearest chromosome end"); ax.set_ylabel("Share of events (%)")
ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for sp in ("right", "top", "left"):
    ax.spines[sp].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT + f"Fig2_distance_to_chromosome_end.{ext}", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ Fig 3: sizes
lo, hi = 2, np.ceil(np.log10(E["size"].max()) * 4) / 4
sbins = 10 ** np.arange(lo, hi + 0.01, 0.25)
fig, axs = plt.subplots(2, 1, figsize=(7.5, 5.4), sharex=True, sharey=True)
for ax, sub, col, lab in ((axs[0], EA, CA, LA), (axs[1], EB, CB, LB)):
    h, _ = np.histogram(sub["size"], sbins)
    ax.bar(sbins[:-1], h / len(sub) * 100, width=np.diff(sbins) * 0.92, align="edge", color=col)
    med = sub["size"].median()
    ax.axvline(med, color=INK, lw=1.2, ls=(0, (4, 3)))
    ax.set_title(f"{lab}: {len(sub)} events, median {med / 1000:.1f} kb", loc="left", fontsize=10.5, color=INK)
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("right", "top", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0); ax.set_ylabel("Share of events (%)")
ymax = axs[0].get_ylim()[1]
for ax, sub in ((axs[0], EA), (axs[1], EB)):
    ax.text(sub["size"].median() * 1.1, ymax * 0.85, "median", fontsize=8.5, color=INK2)
axs[1].set_xscale("log")
axs[1].set_xticks([100, 1e3, 1e4, 1e5, 1e6], ["100 bp", "1 kb", "10 kb", "100 kb", "1 Mb"])
axs[1].set_xlim(sbins[0], sbins[-1])
axs[1].set_xlabel("Event size (region length on E2, log scale)")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT + f"Fig3_size_distribution.{ext}", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ tables
def sz(s):
    return dict(n=len(s), median=float(s.median()), mean=float(s.mean()), q1=float(s.quantile(.25)), q3=float(s.quantile(.75)),
                min=int(s.min()), max=int(s.max()))


size_rows = []
for unit, a, b in (("deletion events (counted once)", EA["size"], EB["size"]),
                   ("strain-level regions (as in the browser)", RA["Region length (bp)"], RB["Region length (bp)"])):
    for lab, s in ((LA, a), (LB, b)):
        d = sz(s)
        size_rows.append({"Unit": unit, "Group": lab, "N": d["n"], "Median (bp)": round(d["median"]), "Mean (bp)": round(d["mean"]),
                          "25th percentile (bp)": round(d["q1"]), "75th percentile (bp)": round(d["q3"]), "Smallest (bp)": d["min"], "Largest (bp)": d["max"]})
SIZE = pd.DataFrame(size_rows)
CH = pd.DataFrame({"Chromosome": CHROMS, "Length (Mb)": [round(sizes[c] / 1e6, 2) for c in CHROMS],
                   f"{LA}: events": [int((EA.chrom == c).sum()) for c in CHROMS],
                   f"{LB}: events": [int((EB.chrom == c).sum()) for c in CHROMS]})
EV = E.reset_index().rename(columns={"index": "Deletion event ID"})
EV = EV[["Deletion event ID", "group", "chrom", "start", "end", "size", "dist_end", "n_strains", "strains_A", "strains_B"]].rename(columns={
    "group": "Group (A, B or both)", "chrom": "Chromosome", "start": "Start (0-based)", "end": "End", "size": "Size (bp)",
    "dist_end": "Distance to nearest chromosome end (bp)", "n_strains": "Strains with this event",
    "strains_A": "Strains in group A", "strains_B": "Strains in group B"})
EV["Distance to nearest chromosome end (bp)"] = EV["Distance to nearest chromosome end (bp)"].round().astype(int)
TESTS = pd.DataFrame([
    {"Question": "Same distribution over the 7 chromosomes?", "Test": "chi-square test of independence (2 groups x 7 chromosomes)",
     "Statistic": round(chi.statistic, 2), "P value": float(f"{chi.pvalue:.3g}")},
    {"Question": "Same distance to the nearest chromosome end?", "Test": "Mann-Whitney U", "Statistic": float(mw_dist.statistic), "P value": float(f"{mw_dist.pvalue:.3g}")},
    {"Question": "Same distance to the nearest chromosome end?", "Test": "Kolmogorov-Smirnov (two samples)", "Statistic": round(ks_dist.statistic, 3), "P value": float(f"{ks_dist.pvalue:.3g}")},
    {"Question": "Same size?", "Test": "Mann-Whitney U", "Statistic": float(mw_size.statistic), "P value": float(f"{mw_size.pvalue:.3g}")}])
NOTES = pd.DataFrame({"Item": [
    "Group A", "Group B", "Unit", "Events in both groups", "Event position and size", "Data S10 listing", "Found", "Found at the edge",
    "Not testable", "Partly low depth", "Read depth normal"], "Definition": [
    "Read evidence 'simple deletion' and 'flanks joined' (direct deletion junction) at one or both ends, as in the breakpoint browser.",
    f"Any read evidence except 'DNA still present', with a TE joined to the flank at one or both ends. {n_b_dna} regions with 'DNA still present' were left out.",
    "Deletion events (regions in different strains that overlap by 80% or more). An event is in a group when the region of at least one strain is.",
    f"{int((E.group == 'both').sum())} events: different strains show different ends. Drawn in both groups; left out of the tests.",
    "Median start and end over the strains with the event, on E2. Distance to the chromosome end is measured from the event midpoint.",
    "One Data S10 deletion in one listed strain.",
    "An all-base deletion region of the same strain overlaps the Data S10 interval.",
    "The interval has no uniquely mappable bases and lies in the breakpoint zone (masked or not uniquely mappable sequence) next to an all-base deletion of the same strain.",
    "The interval has fewer than 30 uniquely mappable bases, so read depth cannot test it.",
    "20% or more of the uniquely mappable bases have low depth, but fewer than 100 bases, the minimum size of an all-base deletion call.",
    "Fewer than 20% of the uniquely mappable bases have low depth (depth of uniquely mapped reads at most 5% of the strain's median)."]})
with pd.ExcelWriter(OUT + "Deletion_groups_tables.xlsx", engine="openpyxl") as xw:
    for name, df in (("DataS10_by_deletion", D10), ("DataS10_by_strain", L10), ("Sizes", SIZE), ("Events_by_chromosome", CH),
                     ("Tests", TESTS), ("Events_A_and_B", EV), ("Definitions", NOTES)):
        df.to_excel(xw, sheet_name=name, index=False)
        ws = xw.sheets[name]
        for row in ws.iter_rows():
            for cell in row:
                cell.font = Font(name="Arial", size=10, bold=cell.row == 1)
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(10, max(len(str(c.value or "")) for c in col[:60]) + 2))
        ws.freeze_panes = "A2"

# ------------------------------------------------------------------ summary.md (all numbers computed here)
f = lambda v: f"{v:,.0f}"
pct = lambda a, b: f"{100 * a / b:.0f}%" if 100 * a / b >= 1 else f"{100 * a / b:.1f}%"
nl = len(L10)
sa, sb = SIZE.iloc[0], SIZE.iloc[1]
ra, rb = SIZE.iloc[2], SIZE.iloc[3]
md = f"""# Data S10 check and deletion groups

Produced by `deletion_allbase/scripts/35_deletion_groups.py`. Tables: `Deletion_groups_tables.xlsx`.

## 1. Are the Data S10 deletions found by the all-base analysis?

Data S10 lists {len(s10)} deletions in {nl} strain listings (one deletion in one strain). For each listing it was tested whether
an all-base deletion region of the same strain overlaps the Data S10 interval.

| Result | Listings |
|---|---|
""" + "\n".join(f"| {k} | {v} ({pct(v, nl)}) |" for k, v in lst_counts.items()) + f"""

Per Data S10 deletion:

| Status | Deletions |
|---|---|
""" + "\n".join(f"| {k} | {v} |" for k, v in del_counts.items()) + f"""

The listings with normal read depth are all short intervals ({L10[L10.Result == NORMAL]['Length (bp)'].min()}-{L10[L10.Result == NORMAL]['Length (bp)'].max()} bp):
{", ".join(sorted(set(L10[L10.Result == NORMAL]['Data S10 ID']), key=lambda s: int(s[3:])))}.

## 2. Genome distribution

Group A: {LA.lower()} ({len(EA)} events, {len(RA)} strain-level regions).
Group B: {LB} ({len(EB)} events, {len(RB)} strain-level regions).
{int((E.group == 'both').sum())} events are in both groups and are left out of the tests below.

- Chromosomes (Fig. 1): chi-square = {chi.statistic:.2f}, P = {chi.pvalue:.3g}.
- Distance to the nearest chromosome end (Fig. 2): median {TA.dist_end.median() / 1e3:,.0f} kb (A) and {TB.dist_end.median() / 1e3:,.0f} kb (B);
  share within 100 kb of an end: {pa.iloc[0]:.0f}% (A) and {pb.iloc[0]:.0f}% (B). Mann-Whitney P = {mw_dist.pvalue:.3g}; Kolmogorov-Smirnov P = {ks_dist.pvalue:.3g}.

## 3. Sizes (Fig. 3)

| Unit | Group | N | Median | Mean | 25th-75th percentile |
|---|---|---|---|---|---|
""" + "\n".join(f"| {r['Unit']} | {r['Group']} | {r['N']} | {f(r['Median (bp)'])} bp | {f(r['Mean (bp)'])} bp | {f(r['25th percentile (bp)'])}-{f(r['75th percentile (bp)'])} bp |"
                for _, r in SIZE.iterrows()) + f"""

Mann-Whitney test on event sizes: P = {mw_size.pvalue:.3g}.
"""
open(OUT + "summary.md", "w", encoding="utf-8").write(md)
print(md)
