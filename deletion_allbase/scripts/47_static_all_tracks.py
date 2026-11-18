#!/usr/bin/env python3
"""Static copy of interactive_all_tracks_version3.html (deletion tracks along E2), run from manuscripts_all/.

View: all chromosomes, 20 kb windows. The figure is drawn from the data
embedded in the interactive page with the page's own rules (build() and draw() in its script), so the two match:
  windows of W bp summed from the 10-kb bins, 300 kb of empty space between chromosomes;
  all-base deletions = deleted uniquely mappable bases / (uniquely mappable bases x 33 strains), 0 in windows with
  < 20% uniquely mappable sequence; SNP-based = SNP x strain cells in runs of >= 5 missing SNPs / (SNPs x 33);
  event counts = events overlapping the window; SNP density per kb of uniquely mappable sequence (0 in low windows);
  repeat fraction; effector genes by midpoint; bars scaled so the track maximum fills 96% of the height.
Output: static_figures_from_interactive_pages/all_tracks_all_chr_20kb.{png,pdf,svg}
Usage: 47_static_all_tracks.py [window_bp]   (default 20000)
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "deletion_allbase/DataS10_check_and_deletion_tracks/interactive_all_tracks_version3.html"
OUTD = "deletion_allbase/static_figures_from_interactive_pages/"
# --seven: TE track limited to the seven TE families selected in the manuscript. Reads the page written by
# 37_window_tracks_interactive.py --seven (all other tracks are identical to version3); output suffix _sevenTE
SEVEN_ONLY = "--seven" in sys.argv
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
if SEVEN_ONLY:
    SRC = "deletion_allbase/DataS10_check_and_deletion_tracks/Fig4_window_tracks_two_groups_sevenTE_interactive.html"
W = int(ARGS[0]) if ARGS else 20000
os.makedirs(OUTD, exist_ok=True)

html = open(SRC, encoding="utf-8").read()
start = html.index("const D = ") + len("const D = ")
D, _ = json.JSONDecoder().raw_decode(html[start:])

k, gap = W // D["bin"], round(300000 / W)
L, x, starts = [], 0, {}
for c in D["chroms"]:
    if L:
        x += gap
    starts[c] = x
    b = {key: np.asarray(v, float) for key, v in D["bins"][c].items()}
    evs = [e for e in D["events"] if e[1] == c]
    for i in range(-(-D["sizes"][c] // W)):
        s, e = i * W, min((i + 1) * W, D["sizes"][c])
        lo, hi = i * k, min((i + 1) * k, len(b["sc"]))
        sc, rp, sn, ef, dl, sd = (b[key][lo:hi].sum() for key in ("sc", "rp", "sn", "ef", "dl", "sd"))
        inw = [ev for ev in evs if ev[2] < e and ev[3] > s]
        low = sc < 0.2 * (e - s)
        L.append(dict(x=x, all=dl / (sc * D["n_strains"]) if not low and sc else 0,
                      snpd=sd / (sn * D["n_strains"]) if sn else 0,
                      A=sum(1 for ev in inw if ev[4]), B=sum(1 for ev in inw if ev[5]),
                      snp=sn / (sc / 1000) if not low and sc else 0, rep=rp / (e - s), eff=ef, low=float(low)))
        x += 1
N = x

TRACKS = [("all", "All-base\ndeletions", "#3b6fb6"), ("snpd", "SNP-based\ndeletions", "#c2571a"),
          ("A", "Simple deletion,\nflanks joined", "#2a78d6"), ("B", "Seven TEs at one\nor both ends" if SEVEN_ONLY else "TE at one or\nboth ends", "#eb6834"),
          ("snp", "SNP density", "#4a4f57"), ("rep", "Repeat fraction", "#7a4fa3"), ("eff", "Effector genes", "#2f8f6b"),
          ("low", "< 20% uniquely\nmappable", "#80868f")]


def ticks(mx, integer):
    raw = mx / 4
    p = 10 ** np.floor(np.log10(raw))
    step = next(m * p for m in (1, 2, 2.5, 5, 10) if m * p >= raw)
    if integer:
        step = max(1, np.ceil(step))
    return np.arange(0, mx + 1e-9, step)


heights = [1.0] * 7 + [0.28]
fig, axes = plt.subplots(len(TRACKS), 1, figsize=(16, 9.6), sharex=True,
                         gridspec_kw=dict(height_ratios=heights, hspace=0.18, left=0.135, right=0.99, top=0.985, bottom=0.05))
xs = np.array([d["x"] for d in L])
seps = [starts[c] - 0.5 * gap for c in D["chroms"][1:]]
for ax, (key, name, col) in zip(axes, TRACKS):
    v = np.array([d[key] for d in L])
    mx = v.max() or 1
    h = np.where(v > 0, np.maximum(v, 0.006 / 0.96 * mx), 0)   # the page draws at least 0.6% of the height
    ax.bar(xs + 0.5, h, width=1.0, color=col, linewidth=0, align="center")
    ax.set_ylim(0, mx / 0.96)
    for sx in seps:
        ax.axvline(sx, color="#5b6068", lw=1.0)
    ax.set_xlim(0, N)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if key == "low":
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
    else:
        ax.set_yticks(ticks(mx, key in ("A", "B", "eff")))
        ax.tick_params(axis="y", labelsize=8)
    ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=10, color="#464b53", labelpad=34)
    ax.tick_params(axis="x", length=0)
axes[-1].set_xticks([starts[c] + (-(-D["sizes"][c] // W)) / 2 for c in D["chroms"]])
axes[-1].set_xticklabels(D["chroms"], fontsize=10)

base = OUTD + f"all_tracks_all_chr_{W // 1000}kb" + ("_sevenTE" if SEVEN_ONLY else "")
for ext, kw in [("png", dict(dpi=300)), ("pdf", {}), ("svg", {})]:
    fig.savefig(f"{base}.{ext}", **kw)
print(len(L), "windows;", {key: round(max(d[key] for d in L), 4) for key, *_ in TRACKS}, "; written", base + ".{png,pdf,svg}")
