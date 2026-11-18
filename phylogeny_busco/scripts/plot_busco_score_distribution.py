#!/usr/bin/env python3
"""Plot the distribution of BUSCO completeness for the 250 final RagTag assemblies."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
BASE = PROJECT / "results/busco_missingness_randomness"
OUT = BASE / "followup"
E2_TABLE = PROJECT / "results/busco_E2_fullgenome_ascomycota_odb12/E2_fullgenome_ascomycota_odb12/run_ascomycota_odb12/full_table.tsv"
N_MARKERS = 2826

BLUE, GRAY, INK, INK2, GRID = "#2a78d6", "#a3a29c", "#0b0b0b", "#52514e", "#e6e5e0"

status = pd.read_csv(BASE / "busco_status_matrix.tsv", sep="\t", index_col=0)
counts = pd.DataFrame({s: (status == s).sum(axis=1) for s in ["Complete", "Duplicated", "Fragmented", "Missing"]})
counts["complete_pct"] = 100 * (counts["Complete"] + counts["Duplicated"]) / N_MARKERS
counts["detected"] = counts["Complete"] + counts["Duplicated"] + counts["Fragmented"]
counts.to_csv(OUT / "per_assembly_busco_counts.tsv", sep="\t")

e2 = pd.read_csv(E2_TABLE, sep="\t", skiprows=2, dtype=str)
e2.columns = [c.lstrip("# ").strip() for c in e2.columns]
e2s = e2.drop_duplicates("Busco id")["Status"].value_counts()
e2_complete_pct = 100 * (e2s.get("Complete", 0) + e2s.get("Duplicated", 0)) / N_MARKERS
ceiling = 100 * (N_MARKERS - 31) / N_MARKERS  # if the 31 markers BUSCO misses in E2 are never found

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": GRAY,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.25]})
for ax in (a1, a2):
    ax.set_facecolor("#fcfcfb")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

x = counts["complete_pct"].values
bins = np.arange(0, 102.5, 2.5)
n, b, patches = a1.hist(x, bins=bins, color=BLUE, edgecolor="#fcfcfb", linewidth=1)
for p, left in zip(patches, b[:-1]):
    if left < 50:
        p.set_facecolor(GRAY)
a1.set_xlim(0, 100)
a1.set_xlabel("BUSCO complete (%)")
a1.set_ylabel("Assemblies")
a1.set_title("All 250 final RagTag assemblies", loc="left", color=INK, fontsize=11)
a1.annotate("6 assemblies <50%\n(E23, K19, K22, K5, T28, U47)", xy=(4, 5), xytext=(8, 60),
            color=INK2, fontsize=9, arrowprops=dict(arrowstyle="-", color=GRAY))

hi = x[x >= 80]
bins2 = np.arange(80, 100.01, 0.5)
a2.hist(hi, bins=bins2, color=BLUE, edgecolor="#fcfcfb", linewidth=1)
a2.set_xlim(80, 100)
a2.set_xlabel("BUSCO complete (%)")
a2.set_title(f"Zoom: the {len(hi)} assemblies at ≥80%", loc="left", color=INK, fontsize=11)
a2.set_ylim(0, a2.get_ylim()[1] * 1.25)
ymax = a2.get_ylim()[1]
box = dict(facecolor="#fcfcfb", edgecolor="none", pad=1.5)
a2.axvline(np.median(x), color=INK2, lw=1, ls=(0, (3, 3)))
a2.text(np.median(x) - 0.2, ymax * 0.97, f"cohort median\n{np.median(x):.1f}%", ha="right", va="top",
        fontsize=8.5, color=INK2, bbox=box)
a2.axvline(e2_complete_pct, color=INK2, lw=1, ls=(0, (1, 2)))
a2.text(e2_complete_pct - 0.2, ymax * 0.80,
        f"E2 reference genome {e2_complete_pct:.1f}%\n(practical ceiling: 31 markers\nare Missing even in E2)",
        ha="right", va="top", fontsize=8.5, color=INK2, bbox=box)

fig.suptitle("BUSCO completeness, ascomycota_odb12 (2,826 markers), BUSCO 6.0.0",
             x=0.01, ha="left", color=INK, fontsize=12)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"busco_completeness_distribution.{ext}", dpi=200, facecolor="white")

print(f"E2 reference: {dict(e2s)}; complete {e2_complete_pct:.2f}%")
print(counts[["complete_pct", "detected", "Missing"]].describe().round(2).to_string())
hq = counts[counts.complete_pct >= 50]
print("excluding 6 failed:", hq["detected"].describe().round(1).to_dict())
never = (status == "Missing").all(axis=0).sum()
print("markers never detected in any assembly:", never, "| detected in >=1:", N_MARKERS - never)
print("bins >=95:", (x >= 95).sum(), "90-95:", ((x >= 90) & (x < 95)).sum(), "50-90:", ((x >= 50) & (x < 90)).sum())
