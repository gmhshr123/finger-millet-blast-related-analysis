#!/usr/bin/env python3
"""Create a polished rectangular phylogeny figure with origin annotations."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import re
import textwrap
from pathlib import Path

from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ORIGIN_COLORS = {
    "Ethiopia": "#E41A1C",
    "Uganda": "#4DAF4A",
    "Kenya": "#377EB8",
    "Tanzania": "#984EA3",
    "Japan": "#FF7F00",
    "India": "#F781BF",
    "Nepal": "#A65628",
    "Philippines": "#00BFC4",
    "Uruguay": "#1B9E77",
    "Brazil": "#66A61E",
    "Bolivia": "#D95F02",
    "China": "#E7298A",
    "Bangladesh": "#7570B3",
    "USA": "#525252",
    "Paraguay": "#A6761D",
    "Mali": "#8C564B",
    "South Korea": "#2B8CBE",
    "Unknown": "#9E9E9E",
}

ORIGIN_ABBREVIATIONS = {
    "Ethiopia": "ETH",
    "Uganda": "UGA",
    "Kenya": "KEN",
    "Tanzania": "TZA",
    "Japan": "JPN",
    "India": "IND",
    "Nepal": "NPL",
    "Philippines": "PHL",
    "Uruguay": "URY",
    "Brazil": "BRA",
    "Bolivia": "BOL",
    "China": "CHN",
    "Bangladesh": "BGD",
    "USA": "USA",
    "Paraguay": "PRY",
    "Mali": "MLI",
    "South Korea": "KOR",
    "Unknown": "UNK",
}

HOST_COLORS = {
    "Eastern African Eleusine coracana": "#00429D",
    "External Eleusine coracana": "#E69F00",
    "Eleusine indica": "#009E73",
    "Other Eleusine": "#CC79A7",
    "Other host": "#D55E00",
    "Unknown/VCF only": "#6B7280",
}

HOST_ORDER = {name: idx for idx, name in enumerate(HOST_COLORS)}
ORIGIN_ORDER = {name: idx for idx, name in enumerate(ORIGIN_COLORS)}


def load_manifest(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = row.get("vcf_sample", "")
            if sample and sample not in rows:
                rows[sample] = row
    return rows


def origin_from_location(location: str) -> str:
    if not location:
        return "Unknown"
    loc = location.replace("\xa0", " ").strip()
    if loc in {"Ethiopia", "Uganda", "Kenya", "Tanzania"}:
        return loc
    if "USA" in loc or loc.endswith(", PA") or loc.endswith(", MD"):
        return "USA"
    if "Brasil" in loc or "Brazil" in loc:
        return "Brazil"
    if "South Korea" in loc or loc == "Korea":
        return "South Korea"
    for country in (
        "China",
        "Japan",
        "India",
        "Nepal",
        "Philippines",
        "Uruguay",
        "Bolivia",
        "Bangladesh",
        "Paraguay",
        "Mali",
    ):
        if country in loc:
            return country
    if "," in loc:
        parts = [part.strip() for part in loc.split(",") if part.strip()]
        if parts:
            return parts[-1].replace("USA", "USA")
    if ":" in loc:
        return loc.split(":", 1)[0].strip()
    return loc.strip()


def host_group(row: dict[str, str]) -> str:
    host = row.get("host", "")
    if row.get("section") == "eastern_africa_selected":
        return "Eastern African Eleusine coracana"
    if "Eleusine coracana" in host:
        return "External Eleusine coracana"
    if "Eleusine indica" in host:
        return "Eleusine indica"
    if "Eleusine" in host:
        return "Other Eleusine"
    if host:
        return "Other host"
    return "Unknown/VCF only"


def display_host(row: dict[str, str]) -> str:
    host = re.sub(r"\s+", " ", row.get("host", "")).strip()
    host = host.replace("Setaria vividis", "Setaria viridis")
    host = host.replace("Pennisetum sp. -", "Pennisetum sp.")
    if not host:
        return host_group(row)
    return host


def sample_display_name(sample: str, manifest: dict[str, dict[str, str]]) -> str:
    row = manifest.get(sample, {})
    isolate = row.get("display_id") or row.get("isolate_id") or sample
    if row.get("section") == "unmapped_vcf_sample" or not row:
        return sample
    return isolate


def tip_label(sample: str, manifest: dict[str, dict[str, str]]) -> str:
    row = manifest.get(sample, {})
    isolate = sample_display_name(sample, manifest)
    origin = origin_from_location(row.get("location", ""))
    return f"{isolate} | {origin}"


def inline_tip_label(sample: str, manifest: dict[str, dict[str, str]]) -> str:
    row = manifest.get(sample, {})
    isolate = sample_display_name(sample, manifest)
    origin = origin_from_location(row.get("location", ""))
    host = display_host(row)
    return f"{isolate} | Origin: {origin} | Host: {host}"


def rotate_for_host(tree, manifest: dict[str, dict[str, str]]) -> None:
    """Rotate sister branches to place similar host groups near each other.

    This only changes child order at internal nodes, so the topology, branch
    lengths, and support values are unchanged.
    """

    def leaf_key(sample: str) -> tuple[int, int, str]:
        row = manifest.get(sample, {})
        host = host_group(row)
        origin = origin_from_location(row.get("location", ""))
        return (
            HOST_ORDER.get(host, len(HOST_ORDER)),
            ORIGIN_ORDER.get(origin, len(ORIGIN_ORDER)),
            sample_display_name(sample, manifest),
        )

    def clade_key(clade):
        terminals = clade.get_terminals()
        keys = [leaf_key(term.name) for term in terminals]
        host_counts = Counter(key[0] for key in keys)
        origin_counts = Counter(key[1] for key in keys)
        dominant_host = min(host_counts, key=lambda key: (-host_counts[key], key))
        dominant_origin = min(origin_counts, key=lambda key: (-origin_counts[key], key))
        host_purity = -host_counts[dominant_host] / len(keys)
        origin_purity = -origin_counts[dominant_origin] / len(keys)
        first_label = min(key[2] for key in keys)
        return (dominant_host, host_purity, dominant_origin, origin_purity, first_label)

    def rotate(clade) -> None:
        for child in clade.clades:
            rotate(child)
        clade.clades.sort(key=clade_key)

    rotate(tree.root)


def clade_depths(tree, layout: str):
    if layout == "cladogram":
        depths = tree.depths(unit_branch_lengths=True)
    else:
        depths = tree.depths()
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
    return depths


def make_layout(tree):
    terminals = tree.get_terminals()
    y_pos = {clade: idx for idx, clade in enumerate(terminals)}

    def assign_y(clade):
        if clade in y_pos:
            return y_pos[clade]
        child_ys = [assign_y(child) for child in clade.clades]
        y_pos[clade] = sum(child_ys) / len(child_ys)
        return y_pos[clade]

    assign_y(tree.root)
    return y_pos


def draw_tree(ax, tree, depths, y_pos, support_font_size: float):
    for clade in tree.find_clades(order="preorder"):
        x = depths[clade]
        y = y_pos[clade]
        for child in clade.clades:
            child_x = depths[child]
            child_y = y_pos[child]
            ax.plot([x, child_x], [child_y, child_y], color="#1f2937", lw=0.85, solid_capstyle="round")
            ax.plot([x, x], [y, child_y], color="#1f2937", lw=0.85, solid_capstyle="round")
            if child.confidence and child.confidence >= 90:
                ax.text(
                    x,
                    child_y - 0.18,
                    f"{child.confidence:g}",
                    fontsize=support_font_size,
                    color="#64748b",
                    ha="right",
                    va="bottom",
                )


def render(
    tree_path: Path,
    manifest_path: Path,
    out_prefix: Path,
    title: str,
    layout: str,
    formats: list[str],
    annotation_mode: str,
    rotate_hosts: bool,
    tip_font_size: float,
    support_font_size: float,
    row_height: float | None,
    show_header: bool,
) -> None:
    tree = Phylo.read(tree_path, "newick")
    manifest = load_manifest(manifest_path)
    if rotate_hosts:
        rotate_for_host(tree, manifest)
    depths = clade_depths(tree, layout)
    y_pos = make_layout(tree)
    terminals = tree.get_terminals()

    max_x = max(depths.values())
    label_x = max_x + max_x * 0.018
    if annotation_mode == "columns":
        origin_x = max_x + max_x * 0.37
        host_x = max_x + max_x * 0.42
        right_edge = max_x + max_x * 0.49
    else:
        marker_origin_x = max_x + max_x * 0.032
        marker_host_x = max_x + max_x * 0.052
        text_x = max_x + max_x * 0.075
        right_edge = max_x * 1.92

    if row_height is None:
        row_height = 0.235 if annotation_mode == "columns" else 0.39
    fig_height = max(18, len(terminals) * row_height)
    fig_width = 16 if annotation_mode == "columns" else 24
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor="white")
    ax.set_facecolor("white")
    draw_tree(ax, tree, depths, y_pos, support_font_size)

    strip_w = max_x * 0.014
    for clade in terminals:
        sample = clade.name
        row = manifest.get(sample, {})
        origin = origin_from_location(row.get("location", ""))
        hgroup = host_group(row)
        y = y_pos[clade]
        origin_color = ORIGIN_COLORS.get(origin, ORIGIN_COLORS["Unknown"])
        host_color = HOST_COLORS.get(hgroup, HOST_COLORS["Unknown/VCF only"])
        leader_start = depths[clade] + max_x * 0.012
        if annotation_mode == "inline":
            leader_end = marker_origin_x - max_x * 0.014
        else:
            leader_end = label_x - max_x * 0.012
        if leader_end > leader_start:
            ax.plot(
                [leader_start, leader_end],
                [y, y],
                color="#cbd5e1",
                lw=0.55,
                linestyle=(0, (1.2, 2.2)),
                zorder=1,
            )
        ax.scatter(
            depths[clade],
            y,
            s=88 if annotation_mode == "inline" else 34,
            color=host_color if annotation_mode == "inline" else origin_color,
            edgecolor="white",
            linewidth=0.95,
            zorder=5,
        )
        if annotation_mode == "columns":
            ax.text(label_x, y, tip_label(sample, manifest), va="center", ha="left", fontsize=tip_font_size, color="#111827")
            ax.add_patch(Rectangle((origin_x, y - 0.38), strip_w, 0.76, facecolor=origin_color, edgecolor="none"))
            ax.add_patch(Rectangle((host_x, y - 0.38), strip_w, 0.76, facecolor=host_color, edgecolor="none"))
        else:
            ax.scatter(
                marker_origin_x,
                y,
                s=62,
                marker="o",
                color=origin_color,
                edgecolor="white",
                linewidth=0.8,
                zorder=6,
            )
            ax.scatter(
                marker_host_x,
                y,
                s=62,
                marker="s",
                color=host_color,
                edgecolor="white",
                linewidth=0.8,
                zorder=6,
            )
            ax.text(
                text_x,
                y,
                inline_tip_label(sample, manifest),
                va="center",
                ha="left",
                fontsize=tip_font_size,
                color="#111827",
                clip_on=False,
            )

    if annotation_mode == "columns":
        ax.text(origin_x + strip_w / 2, -2.2, "Origin", ha="center", va="bottom", fontsize=8, rotation=90, color="#334155")
        ax.text(host_x + strip_w / 2, -2.2, "Host", ha="center", va="bottom", fontsize=8, rotation=90, color="#334155")

    top_limit = -4.8
    if show_header:
        ax.set_title(title, loc="left", fontsize=16, weight="bold", pad=18, color="#111827")
        subtitle = "IQ-TREE GTR+ASC consensus tree from BUSCO-gene SNPs; node labels show support where >=90."
        if rotate_hosts:
            subtitle += " Sister branches were rotated to group host colors where topology allows."
        ax.text(0, -3.8, subtitle, ha="left", va="bottom", fontsize=9, color="#475569")
    else:
        top_limit = -1.0

    xlabel = "Topological depth (cladogram)" if layout == "cladogram" else "Substitutions per site"
    ax.set_xlabel(xlabel, fontsize=9, color="#334155")
    ax.set_yticks([])
    ax.set_xlim(-max_x * 0.015, right_edge)
    ax.set_ylim(len(terminals) - 0.3, top_limit)
    ax.grid(axis="x", color="#e5e7eb", lw=0.55)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#cbd5e1")
    ax.tick_params(axis="x", labelsize=8, colors="#475569")

    origins_in_tree = []
    hosts_in_tree = []
    for clade in terminals:
        row = manifest.get(clade.name, {})
        origin = origin_from_location(row.get("location", ""))
        hgroup = host_group(row)
        if origin not in origins_in_tree:
            origins_in_tree.append(origin)
        if hgroup not in hosts_in_tree:
            hosts_in_tree.append(hgroup)

    origin_handles = [
        Line2D([0], [0], marker="s", linestyle="", color=ORIGIN_COLORS.get(origin, "#9ca3af"), label=origin, markersize=8)
        for origin in sorted(origins_in_tree)
    ]
    host_handles = [
        Line2D([0], [0], marker="s", linestyle="", color=HOST_COLORS.get(group, "#9ca3af"), label=group, markersize=8)
        for group in hosts_in_tree
    ]
    legend_font_size = 7 if annotation_mode == "columns" else 10
    legend_title_size = 8 if annotation_mode == "columns" else 11
    if annotation_mode == "inline":
        origin_legend = fig.legend(
            handles=origin_handles,
            title="Origin",
            loc="lower left",
            bbox_to_anchor=(0.012, 0.028),
            bbox_transform=fig.transFigure,
            ncol=5,
            frameon=False,
            fontsize=legend_font_size,
            title_fontsize=legend_title_size,
            handletextpad=0.4,
            columnspacing=0.9,
        )
        host_legend = fig.legend(
            handles=host_handles,
            title="Host/source group",
            loc="lower right",
            bbox_to_anchor=(0.988, 0.028),
            bbox_transform=fig.transFigure,
            ncol=2,
            frameon=False,
            fontsize=legend_font_size,
            title_fontsize=legend_title_size,
            handletextpad=0.4,
            columnspacing=0.9,
        )
    else:
        leg1 = ax.legend(
        handles=origin_handles,
        title="Origin",
        loc="lower left",
        bbox_to_anchor=(0.0, -0.115),
        ncol=5,
        frameon=False,
        fontsize=legend_font_size,
        title_fontsize=legend_title_size,
        handletextpad=0.4,
        columnspacing=0.9,
        )
        ax.add_artist(leg1)
        ax.legend(
            handles=host_handles,
            title="Host/source group",
            loc="lower right",
            bbox_to_anchor=(1.0, -0.115),
            ncol=2,
            frameon=False,
            fontsize=legend_font_size,
            title_fontsize=legend_title_size,
            handletextpad=0.4,
            columnspacing=0.9,
        )

    note = "\n".join(
        textwrap.wrap(
            "Tip labels use the local strain manifest; run IDs were replaced with isolate names where BioSample/SRA mappings were available.",
            120,
        )
    )
    fig.text(0.01, 0.006, note, fontsize=7, color="#64748b", ha="left", va="bottom")
    bottom_margin = 0.145 if annotation_mode == "inline" else 0.045
    fig.tight_layout(rect=(0, bottom_margin, 1, 0.985))
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        if fmt == "png":
            fig.savefig(out_prefix.with_suffix(".png"), dpi=320)
        elif fmt in {"pdf", "svg"}:
            fig.savefig(out_prefix.with_suffix(f".{fmt}"))
        else:
            raise ValueError(f"Unsupported format: {fmt}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--title", default="Magnaporthe BUSCO SNP Consensus Tree")
    parser.add_argument("--layout", choices=["phylogram", "cladogram"], default="phylogram")
    parser.add_argument("--formats", nargs="+", choices=["pdf", "png", "svg"], default=["pdf"])
    parser.add_argument("--annotation-mode", choices=["columns", "inline"], default="columns")
    parser.add_argument("--rotate-hosts", action="store_true")
    parser.add_argument("--tip-font-size", type=float, default=None)
    parser.add_argument("--support-font-size", type=float, default=5.2)
    parser.add_argument("--row-height", type=float, default=None)
    parser.add_argument("--no-header", action="store_true")
    args = parser.parse_args()
    tip_font_size = args.tip_font_size
    if tip_font_size is None:
        tip_font_size = 6.0 if args.annotation_mode == "columns" else 10.5
    render(
        args.tree,
        args.manifest,
        args.out_prefix,
        args.title,
        args.layout,
        args.formats,
        args.annotation_mode,
        args.rotate_hosts,
        tip_font_size,
        args.support_font_size,
        args.row_height,
        not args.no_header,
    )
    for fmt in args.formats:
        print(args.out_prefix.with_suffix(f".{fmt}"))


if __name__ == "__main__":
    main()
