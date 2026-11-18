#!/usr/bin/env python3
"""Render a cladogram with right-aligned tips.

Layout:  tree | SampleID Host Origin   (tightly packed, coloured text only)

- Tree is compact horizontally.
- Sample ID and Host use host-coloured font.
- Origin uses origin-coloured font.
- No background colours, no pills, no bands.
- All bootstrap values shown.
- Pure white, no axes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from render_tree_publication import (
    ORIGIN_ABBREVIATIONS,
    ORIGIN_COLORS,
    ORIGIN_ORDER,
    clade_depths,
    display_host,
    load_manifest,
    origin_from_location,
    sample_display_name,
)


# ---------------------------------------------------------------------------
# Host colours
# ---------------------------------------------------------------------------

HOST_COLORS = {
    "Outgroup": "#000000",          # black
    "Eleusine coracana": "#0072B2",  # blue (Okabe-Ito)
    "Eleusine indica": "#009E73",    # bluish green
    "Eleusine": "#56B4E9",           # sky blue
    "Eleusine spp.": "#56B4E9",      # sky blue
    "Triticum aestivum": "#D55E00",  # vermillion
    "Oryza sativa": "#CC79A7",       # reddish purple
    "Digitaria sanguinalis": "#0072B2",  # blue (same genus as outgroup)
    "Digitaria spp.": "#332288",     # indigo (Tol)
    "Setaria italica": "#E69F00",    # orange
    "Setaria viridis": "#999933",    # olive (Tol)
    "Bromus spp.": "#882255",        # wine (Tol)
    "Festuca parvigluma": "#AA4499", # rose (Tol)
    "Lolium arundinaceum": "#AA4499",# rose
    "Lolium perenne": "#CC6677",     # dusty rose (Tol)
    "Eragrostis curvula": "#117733", # green (Tol)
    "Poa trivialis": "#44AA99",      # teal (Tol)
    "Panicum repens": "#AA3377",     # purple (Tol bright)
    "Pennisetum sp.": "#882255",     # wine
    "Echinochloa sp.": "#333333",    # dark grey
    "Leersia hexandra": "#D55E00",   # vermillion
    "Unknown": "#888888",            # grey
}

# Saturated origin colours for text
ORIGIN_TEXT_COLORS = {
    "Ethiopia":    "#B91C1C",
    "Uganda":      "#15803D",
    "Kenya":       "#1D4ED8",
    "Tanzania":    "#7E22CE",
    "Japan":       "#C2410C",
    "India":       "#BE185D",
    "Nepal":       "#6D28D9",
    "Philippines": "#0F766E",
    "Uruguay":     "#0369A1",
    "Brazil":      "#15803D",
    "Bolivia":     "#A16207",
    "China":       "#BE123C",
    "Bangladesh":  "#4338CA",
    "USA":         "#475569",
    "Paraguay":    "#92400E",
    "Mali":        "#78350F",
    "South Korea": "#0E7490",
    "Unknown":     "#6B7280",
}

HOST_ORDER = {name: idx for idx, name in enumerate(HOST_COLORS)}
OUTGROUPS = {"CKF3589", "DS0505", "W97-11"}

SUBFAMILY_COLORS = {
    "Chloridoideae": "#0072B2",  # Okabe-Ito blue
    "Panicoideae": "#009E73",    # Okabe-Ito bluish green
    "Pooideae": "#D55E00",       # Okabe-Ito vermillion
    "Oryzoideae": "#CC79A7",     # Okabe-Ito reddish purple
    "Unknown": "#666666",
}
OUTGROUP_COLOR = "#000000"

SUBFAMILY_ORDER = {name: idx for idx, name in enumerate(SUBFAMILY_COLORS)}


def host_group(row: dict[str, str], sample: str = "") -> str:
    if sample in OUTGROUPS:
        return "Outgroup"
    host = display_host(row).replace("\xa0", " ").strip()
    host = " ".join(host.split())
    if not host or host == "Unknown/VCF only":
        return "Unknown"
    if host == "Setaria vividis":
        return "Setaria viridis"
    if host == "Pennisetum sp. -":
        return "Pennisetum sp."
    if host.startswith("Pennisetum sp."):
        return "Pennisetum sp."
    if host == "Eleusine coracana subsp. africana":
        return "Eleusine coracana"
    return host


def subfamily_group(row: dict[str, str], sample: str = "") -> str:
    host = display_host(row).replace("\xa0", " ").strip()
    host = " ".join(host.split())
    if host == "Setaria vividis":
        host = "Setaria viridis"
    if host == "Pennisetum sp. -" or host.startswith("Pennisetum sp."):
        host = "Pennisetum sp."
    if host == "Eleusine coracana subsp. africana":
        host = "Eleusine coracana"
    genus = host.split()[0] if host else ""
    if genus in {"Eleusine", "Digitaria", "Eragrostis"}:
        return "Chloridoideae"
    if genus in {"Echinochloa", "Panicum", "Pennisetum", "Setaria"}:
        return "Panicoideae"
    if genus in {"Bromus", "Festuca", "Lolium", "Poa", "Triticum"}:
        return "Pooideae"
    if genus in {"Leersia", "Oryza"}:
        return "Oryzoideae"
    return "Unknown"


def short_host(full_host: str) -> str:
    parts = full_host.replace("\xa0", " ").strip().split()
    if not parts:
        return full_host
    if full_host in ("Unknown", "Outgroup", "Unknown/VCF only"):
        return full_host
    genus = parts[0]
    if len(parts) >= 2:
        species = parts[1]
        if species in ("sp.", "spp.", "sp"):
            return f"{genus} sp."
        return f"{genus[0]}. {species}"
    return genus


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def make_layout(tree, y_step: float):
    terminals = tree.get_terminals()
    y_pos = {clade: idx * y_step for idx, clade in enumerate(terminals)}
    def assign_y(clade):
        if clade in y_pos:
            return y_pos[clade]
        child_ys = [assign_y(child) for child in clade.clades]
        y_pos[clade] = sum(child_ys) / len(child_ys)
        return y_pos[clade]
    assign_y(tree.root)
    return y_pos


def support_label(clade) -> str:
    if clade.name and "/" in clade.name:
        _alrt, bootstrap = clade.name.split("/", 1)
        try:
            value = float(bootstrap)
        except ValueError:
            return ""
        return f"{value:g}"
    if clade.confidence is not None:
        return f"{float(clade.confidence):g}"
    if clade.name:
        try:
            value = float(clade.name)
        except ValueError:
            return ""
        return f"{value:g}"
    return ""


def bootstrap_value(clade) -> float | None:
    label = support_label(clade)
    if not label:
        return None
    try:
        return float(label)
    except ValueError:
        return None


def rotate_for_host_blocks(tree, manifest, color_mode="host_origin"):
    def leaf_key(sample):
        row = manifest.get(sample, {})
        group = subfamily_group(row, sample) if color_mode == "subfamily" else host_group(row, sample)
        origin = origin_from_location(row.get("location", ""))
        order = SUBFAMILY_ORDER if color_mode == "subfamily" else HOST_ORDER
        return (order.get(group, len(order)),
                ORIGIN_ORDER.get(origin, len(ORIGIN_ORDER)),
                sample_display_name(sample, manifest))

    def clade_key(clade):
        terms = clade.get_terminals()
        keys = [leaf_key(t.name) for t in terms]
        og = sum(1 for t in terms if t.name in OUTGROUPS)
        hc, oc = {}, {}
        for h, o, _ in keys:
            hc[h] = hc.get(h, 0) + 1
            oc[o] = oc.get(o, 0) + 1
        dh = min(hc, key=lambda k: (-hc[k], k))
        do = min(oc, key=lambda k: (-oc[k], k))
        return (0 if og else 1, -og/len(keys), dh, -hc[dh]/len(keys),
                do, -oc[do]/len(keys), min(k[2] for k in keys))

    def rotate(clade):
        for child in clade.clades:
            rotate(child)
        clade.clades.sort(key=clade_key)
    rotate(tree.root)


def rotate_eibtmp13_after_ta_group(tree):
    """Display-only rotation: place (Ta053i, Ta117, Ta37) before EiBTMP13_1."""
    eibtmp_tree_id = "ERR2612753"  # display name EiBTMP13_1 in the manifest
    ta_group = {"Ta053i", "Ta117", "Ta37"}
    names = {t.name for t in tree.get_terminals()}
    if eibtmp_tree_id not in names or not ta_group.issubset(names):
        return

    mrca = tree.common_ancestor(eibtmp_tree_id, *sorted(ta_group))

    def child_rank(child):
        child_names = {t.name for t in child.get_terminals()}
        if ta_group.issubset(child_names):
            return 0
        if eibtmp_tree_id in child_names:
            return 1
        return 2

    mrca.clades.sort(key=child_rank)


def draw_tree(ax, tree, depths, y_pos, support_font_size, branch_line_width):
    for clade in tree.find_clades(order="preorder"):
        x, y = depths[clade], y_pos[clade]
        for child in clade.clades:
            cx, cy = depths[child], y_pos[child]
            ax.plot([x, cx], [cy, cy], color="#253041", lw=branch_line_width, solid_capstyle="round")
            ax.plot([x, x], [y, cy], color="#253041", lw=branch_line_width, solid_capstyle="round")
            s = bootstrap_value(child)
            if s is not None:
                ax.text(x, cy - 0.20, f"{s:g}", fontsize=support_font_size,
                        color="#5b6472", ha="right", va="bottom")


# ---------------------------------------------------------------------------
# Main render
# ---------------------------------------------------------------------------

def render(
    tree_path,
    manifest_path,
    out_prefix,
    layout,
    formats,
    color_mode="host_origin",
    sample_font_size=12.0,
    host_font_size=12.0,
    origin_font_size=10.0,
    header_font_size=14.0,
    support_font_size=6.5,
    legend_font_size=13.0,
    legend_title_font_size=14.0,
    branch_line_width=1.05,
    legend_line_width=4.0,
    char_width=0.038,
    column_gap=0.15,
    figure_width=18.0,
    swap_eibtmp_ta_group=False,
):
    tree = Phylo.read(tree_path, "newick")
    manifest = load_manifest(manifest_path)
    rotate_for_host_blocks(tree, manifest, color_mode=color_mode)
    if swap_eibtmp_ta_group:
        rotate_eibtmp13_after_ta_group(tree)

    y_step = 1.24
    depths = clade_depths(tree, layout)

    # Compress tree horizontally
    raw_max = max(depths.values()) if depths else 1.0
    TREE_W = 6.0
    sc = TREE_W / raw_max if raw_max else 1.0
    depths = {k: v * sc for k, v in depths.items()}

    if layout == "cladogram":
        td = max(depths[c] for c in tree.get_terminals())
        for c in tree.get_terminals():
            depths[c] = td

    y_pos = make_layout(tree, y_step)
    terminals = tree.get_terminals()
    max_x = max(depths.values())

    # Collect text for width estimation
    rows_data = []
    for clade in terminals:
        sample = clade.name
        row = manifest.get(sample, {})
        origin = origin_from_location(row.get("location", ""))
        sname = sample_display_name(sample, manifest)
        hshort = short_host(display_host(row))
        rows_data.append((sample, row, origin, sname, hshort))

    max_slen = max(len(r[3]) for r in rows_data)
    max_hlen = max(len(r[4]) for r in rows_data)

    # Character width in data-units — calibrated for 18" figure
    # 18 in * 72 pt/in = 1296 pt.  Data range ~8 units → 1 du ≈ 162 pt.
    # 10pt font, avg char ≈ 5.2pt → 5.2/162 ≈ 0.032 du/char
    CW = char_width
    GAP = column_gap

    label_x     = max_x + 0.06
    host_x      = label_x + max_slen * CW + GAP
    origin_x    = host_x  + max_hlen * CW + GAP
    right_edge  = origin_x + 4 * CW + 0.3

    fig_w, fig_h = figure_width, max(20, len(terminals) * 0.43)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")
    ax.set_facecolor("white")

    draw_tree(
        ax,
        tree,
        depths,
        y_pos,
        support_font_size=support_font_size,
        branch_line_width=branch_line_width,
    )

    for clade, (sample, row, origin, sname, hshort) in zip(terminals, rows_data):
        hgroup = host_group(row, sample)
        if color_mode == "subfamily":
            color = OUTGROUP_COLOR if sample in OUTGROUPS else SUBFAMILY_COLORS.get(
                subfamily_group(row, sample), SUBFAMILY_COLORS["Unknown"]
            )
            hcol = color
            ocol = color
        else:
            hcol = HOST_COLORS.get(hgroup, HOST_COLORS["Unknown"])
            ocol = ORIGIN_TEXT_COLORS.get(origin, ORIGIN_TEXT_COLORS["Unknown"])
        origin_abbr = ORIGIN_ABBREVIATIONS.get(origin, origin)
        y = y_pos[clade]



        # Sample ID — host colour, large
        ax.text(label_x, y, sname, fontsize=sample_font_size, fontweight="bold",
                color=hcol, ha="left", va="center", clip_on=False)
        # Host — host colour, same style as Sample ID
        ax.text(host_x, y, hshort, fontsize=host_font_size, fontweight="bold",
                fontstyle="italic",
                color=hcol, ha="left", va="center", clip_on=False)
        # Origin — origin colour
        ax.text(origin_x, y, origin_abbr, fontsize=origin_font_size, fontweight="bold",
                color=ocol, ha="left", va="center", clip_on=False)

    # Headers — larger than data text
    hy = -y_step * 1.6
    ax.text(label_x, hy, "Sample ID", fontsize=header_font_size, fontweight="bold",
            color="#1e293b", ha="left", va="center")
    ax.text(host_x, hy, "Host", fontsize=header_font_size, fontweight="bold",
            color="#1e293b", ha="left", va="center")
    ax.text(origin_x, hy, "Origin", fontsize=header_font_size, fontweight="bold",
            color="#1e293b", ha="left", va="center")

    if color_mode == "subfamily":
        handles = [
            Line2D([0], [0], color=SUBFAMILY_COLORS[name], lw=legend_line_width, label=name)
            for name in ("Chloridoideae", "Panicoideae", "Pooideae", "Oryzoideae")
        ]
        handles.append(Line2D([0], [0], color=OUTGROUP_COLOR, lw=legend_line_width, label="Outgroup"))
        fig.legend(
            handles=handles,
            title="Host subfamily",
            loc="upper center",
            bbox_to_anchor=(0.5, 0.06),
            ncol=5,
            frameon=False,
            fontsize=legend_font_size,
            title_fontsize=legend_title_font_size,
            handlelength=2.2,
            columnspacing=1.3,
        )

    # Pure white
    ax.set_yticks([]); ax.set_xticks([]); ax.set_xlabel("")
    ax.set_xlim(-0.15, right_edge)
    ax.set_ylim(y_pos[terminals[-1]] + y_step * 0.85, hy - y_step * 0.5)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

    bottom_margin = 0.064 if color_mode == "subfamily" else 0.01
    fig.tight_layout(rect=(0, bottom_margin, 1, 0.998))
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(out_prefix.with_name(f"{out_prefix.name}.{fmt}"))
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--out-prefix", required=True, type=Path)
    p.add_argument("--layout", choices=["phylogram", "cladogram"], default="cladogram")
    p.add_argument("--formats", nargs="+", choices=["pdf", "svg"], default=["pdf"])
    p.add_argument("--color-mode", choices=["host_origin", "subfamily"], default="host_origin")
    p.add_argument("--sample-font-size", type=float, default=12.0)
    p.add_argument("--host-font-size", type=float, default=12.0)
    p.add_argument("--origin-font-size", type=float, default=10.0)
    p.add_argument("--header-font-size", type=float, default=14.0)
    p.add_argument("--support-font-size", type=float, default=6.5)
    p.add_argument("--legend-font-size", type=float, default=13.0)
    p.add_argument("--legend-title-font-size", type=float, default=14.0)
    p.add_argument("--branch-line-width", type=float, default=1.05)
    p.add_argument("--legend-line-width", type=float, default=4.0)
    p.add_argument("--char-width", type=float, default=0.038)
    p.add_argument("--column-gap", type=float, default=0.15)
    p.add_argument("--figure-width", type=float, default=18.0)
    p.add_argument("--swap-eibtmp-ta-group", action="store_true")
    a = p.parse_args()
    render(
        a.tree,
        a.manifest,
        a.out_prefix,
        a.layout,
        a.formats,
        color_mode=a.color_mode,
        sample_font_size=a.sample_font_size,
        host_font_size=a.host_font_size,
        origin_font_size=a.origin_font_size,
        header_font_size=a.header_font_size,
        support_font_size=a.support_font_size,
        legend_font_size=a.legend_font_size,
        legend_title_font_size=a.legend_title_font_size,
        branch_line_width=a.branch_line_width,
        legend_line_width=a.legend_line_width,
        char_width=a.char_width,
        column_gap=a.column_gap,
        figure_width=a.figure_width,
        swap_eibtmp_ta_group=a.swap_eibtmp_ta_group,
    )
    for fmt in a.formats:
        print(a.out_prefix.with_name(f"{a.out_prefix.name}.{fmt}"))

if __name__ == "__main__":
    main()
