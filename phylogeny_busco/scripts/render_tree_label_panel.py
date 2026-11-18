#!/usr/bin/env python3
"""Render a clean phylogeny with a spacious right-side label panel."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

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


HOST_COLORS = {
    "Outgroup": "#111827",
    "Eleusine coracana": "#00429D",
    "Eleusine indica": "#009E73",
    "Eleusine": "#56B4E9",
    "Eleusine spp.": "#56B4E9",
    "Triticum aestivum": "#D55E00",
    "Oryza sativa": "#CC79A7",
    "Digitaria sanguinalis": "#00BFC4",
    "Digitaria spp.": "#00A2C7",
    "Setaria italica": "#E69F00",
    "Setaria viridis": "#7CAE00",
    "Bromus spp.": "#A6761D",
    "Festuca parvigluma": "#8C564B",
    "Lolium arundinaceum": "#8C564B",
    "Lolium perenne": "#BC80BD",
    "Eragrostis curvula": "#7570B3",
    "Poa trivialis": "#1B9E77",
    "Panicum repens": "#66A61E",
    "Pennisetum sp.": "#E7298A",
    "Echinochloa sp.": "#525252",
    "Leersia hexandra": "#D95F02",
    "Unknown": "#6B7280",
}

HOST_ORDER = {name: idx for idx, name in enumerate(HOST_COLORS)}
OUTGROUPS = {"CKF3589", "DS0505", "W97-11"}


def host_group(row: dict[str, str], sample: str = "") -> str:
    """Group by host identity only, not geography or analysis section."""

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
    """Return the IQ-TREE ultrafast bootstrap support, when present.

    IQ-TREE writes SH-aLRT/UFBoot support as an internal node label such as
    "93.7/80". Bio.Phylo stores that as clade.name, not clade.confidence.
    """

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


def rotate_for_host_blocks(tree, manifest: dict[str, dict[str, str]]) -> None:
    """Rotate/sort sister clades to cluster host groups without changing topology."""

    def leaf_key(sample: str) -> tuple[int, int, str]:
        row = manifest.get(sample, {})
        host = host_group(row, sample)
        origin = origin_from_location(row.get("location", ""))
        return (
            HOST_ORDER.get(host, len(HOST_ORDER)),
            ORIGIN_ORDER.get(origin, len(ORIGIN_ORDER)),
            sample_display_name(sample, manifest),
        )

    def clade_key(clade):
        terminals = clade.get_terminals()
        keys = [leaf_key(term.name) for term in terminals]
        outgroup_count = sum(1 for term in terminals if term.name in OUTGROUPS)
        host_counts: dict[int, int] = {}
        origin_counts: dict[int, int] = {}
        for host, origin, _label in keys:
            host_counts[host] = host_counts.get(host, 0) + 1
            origin_counts[origin] = origin_counts.get(origin, 0) + 1
        dominant_host = min(host_counts, key=lambda key: (-host_counts[key], key))
        dominant_origin = min(origin_counts, key=lambda key: (-origin_counts[key], key))
        host_purity = -host_counts[dominant_host] / len(keys)
        origin_purity = -origin_counts[dominant_origin] / len(keys)
        first_label = min(key[2] for key in keys)
        outgroup_score = -outgroup_count / len(keys)
        return (0 if outgroup_count else 1, outgroup_score, dominant_host, host_purity, dominant_origin, origin_purity, first_label)

    def rotate(clade) -> None:
        for child in clade.clades:
            rotate(child)
        clade.clades.sort(key=clade_key)

    rotate(tree.root)


def draw_tree(ax, tree, depths, y_pos, support_font_size: float) -> None:
    for clade in tree.find_clades(order="preorder"):
        x = depths[clade]
        y = y_pos[clade]
        for child in clade.clades:
            child_x = depths[child]
            child_y = y_pos[child]
            ax.plot([x, child_x], [child_y, child_y], color="#253041", lw=1.05, solid_capstyle="round")
            ax.plot([x, x], [y, child_y], color="#253041", lw=1.05, solid_capstyle="round")
            support = bootstrap_value(child)
            if support is not None and support >= 70:
                ax.text(
                    x,
                    child_y - 0.20,
                    f"{support:g}",
                    fontsize=support_font_size,
                    color="#5b6472",
                    ha="right",
                    va="bottom",
                )


def contiguous_runs(terminals, y_pos, manifest):
    runs = []
    current_group = None
    start_y = None
    last_y = None
    for clade in terminals:
        row = manifest.get(clade.name, {})
        group = host_group(row, clade.name)
        y = y_pos[clade]
        if current_group is None:
            current_group = group
            start_y = y
        elif group != current_group:
            runs.append((current_group, start_y, last_y))
            current_group = group
            start_y = y
        last_y = y
    if current_group is not None:
        runs.append((current_group, start_y, last_y))
    return runs


def render(
    tree_path: Path,
    manifest_path: Path,
    out_prefix: Path,
    layout: str,
    formats: list[str],
) -> None:
    tree = Phylo.read(tree_path, "newick")
    manifest = load_manifest(manifest_path)
    rotate_for_host_blocks(tree, manifest)

    y_step = 1.24
    depths = clade_depths(tree, layout)
    if layout == "cladogram":
        terminal_depth = max(depths[clade] for clade in tree.get_terminals())
        for clade in tree.get_terminals():
            depths[clade] = terminal_depth
    y_pos = make_layout(tree, y_step)
    terminals = tree.get_terminals()
    max_x = max(depths.values())

    panel_x = max_x + max_x * 0.07
    origin_x = panel_x + max_x * 0.01
    host_x = panel_x + max_x * 0.045
    label_x = panel_x + max_x * 0.085
    origin_text_x = panel_x + max_x * 0.43
    host_text_x = panel_x + max_x * 0.55
    panel_right = max_x * 2.35
    fig_height = max(20, len(terminals) * 0.43)
    fig, ax = plt.subplots(figsize=(22, fig_height), facecolor="white")
    ax.set_facecolor("white")

    for group, y1, y2 in contiguous_runs(terminals, y_pos, manifest):
        color = HOST_COLORS.get(group, "#6B7280")
        band_y = y1 - y_step * 0.47
        band_h = (y2 - y1) + y_step * 0.94
        ax.add_patch(
            Rectangle(
                (panel_x - max_x * 0.015, band_y),
                panel_right - panel_x + max_x * 0.015,
                band_h,
                facecolor=color,
                edgecolor="none",
                alpha=0.16,
                zorder=0,
            )
        )
        ax.add_patch(
            Rectangle(
                (panel_x - max_x * 0.018, band_y),
                max_x * 0.006,
                band_h,
                facecolor=color,
                edgecolor="none",
                alpha=0.95,
                zorder=1,
            )
        )

    draw_tree(ax, tree, depths, y_pos, support_font_size=7.2)

    for clade in terminals:
        sample = clade.name
        row = manifest.get(sample, {})
        origin = origin_from_location(row.get("location", ""))
        origin_abbr = ORIGIN_ABBREVIATIONS.get(origin, origin)
        host = display_host(row)
        hgroup = host_group(row, sample)
        y = y_pos[clade]
        origin_color = ORIGIN_COLORS.get(origin, ORIGIN_COLORS["Unknown"])
        host_color = HOST_COLORS.get(hgroup, HOST_COLORS["Unknown"])
        leader_start = depths[clade] + max_x * 0.012
        leader_end = origin_x - max_x * 0.018
        if leader_end > leader_start:
            ax.plot(
                [leader_start, leader_end],
                [y, y],
                color="#d6dce5",
                lw=0.7,
                linestyle=(0, (1.0, 2.0)),
                zorder=1,
            )
        ax.scatter(depths[clade], y, s=92, color=host_color, edgecolor="white", linewidth=1.0, zorder=4)
        ax.scatter(origin_x, y, s=66, color=origin_color, edgecolor="white", linewidth=0.9, zorder=5)
        ax.scatter(host_x, y, s=66, marker="s", color=host_color, edgecolor="white", linewidth=0.9, zorder=5)
        ax.text(
            label_x,
            y,
            sample_display_name(sample, manifest),
            fontsize=11.5,
            color="#111827",
            ha="left",
            va="center",
            clip_on=False,
        )
        ax.text(
            origin_text_x,
            y,
            origin_abbr,
            fontsize=10.5,
            color="#3f4652",
            ha="left",
            va="center",
            clip_on=False,
        )
        ax.text(
            host_text_x,
            y,
            host,
            fontsize=10.5,
            color="#3f4652",
            ha="left",
            va="center",
            clip_on=False,
        )

    ax.axvline(panel_x - max_x * 0.035, color="#e5e7eb", lw=0.9, zorder=0)
    ax.set_yticks([])
    ax.set_xlabel("Topological depth" if layout == "cladogram" else "Substitutions per site", fontsize=10, color="#344054")
    ax.set_xlim(-max_x * 0.015, panel_right)
    ax.set_ylim(y_pos[terminals[-1]] + y_step * 0.85, -y_step * 0.85)
    ax.grid(axis="x", color="#eef2f7", lw=0.65)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#cbd5e1")
    ax.tick_params(axis="x", labelsize=9, colors="#475569")

    origins_in_tree = []
    hosts_in_tree = []
    for clade in terminals:
        row = manifest.get(clade.name, {})
        origin = origin_from_location(row.get("location", ""))
        group = host_group(row, clade.name)
        if origin not in origins_in_tree:
            origins_in_tree.append(origin)
        if group not in hosts_in_tree:
            hosts_in_tree.append(group)

    origin_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color=ORIGIN_COLORS.get(origin, "#9E9E9E"),
            label=ORIGIN_ABBREVIATIONS.get(origin, origin),
            markersize=8,
        )
        for origin in sorted(origins_in_tree)
    ]
    host_handles = [
        Line2D([0], [0], marker="s", linestyle="", color=HOST_COLORS.get(group, "#6B7280"), label=group, markersize=8)
        for group in hosts_in_tree
    ]
    fig.legend(
        handles=origin_handles,
        title="Origin",
        loc="lower left",
        bbox_to_anchor=(0.012, 0.028),
        ncol=5,
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        handletextpad=0.45,
        columnspacing=0.9,
    )
    fig.legend(
        handles=host_handles,
        title="Host/source group",
        loc="lower right",
        bbox_to_anchor=(0.988, 0.028),
        ncol=2,
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        handletextpad=0.45,
        columnspacing=0.9,
    )
    fig.tight_layout(rect=(0, 0.145, 1, 0.995))
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(out_prefix.with_name(f"{out_prefix.name}.{fmt}"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--layout", choices=["phylogram", "cladogram"], default="phylogram")
    parser.add_argument("--formats", nargs="+", choices=["pdf", "svg"], default=["pdf"])
    args = parser.parse_args()
    render(args.tree, args.manifest, args.out_prefix, args.layout, args.formats)
    for fmt in args.formats:
        print(args.out_prefix.with_name(f"{args.out_prefix.name}.{fmt}"))


if __name__ == "__main__":
    main()
