#!/usr/bin/env python3
"""Parse BUSCO short_summary files into one TSV table."""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path


def parse_n50(text: str):
    match = re.search(r"^\s*([\d.]+)\s+([kM]?bp)\s+Scaffold N50", text, re.M)
    if not match:
        return "", ""
    value = float(match.group(1))
    unit = match.group(2)
    mult = {"bp": 1, "kbp": 1_000, "Mbp": 1_000_000}[unit]
    return f"{match.group(1)} {unit}", int(value * mult)


def parse_summary(path: str, lineage: str):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    run_dir = os.path.basename(os.path.dirname(path))
    sample = run_dir
    suffix = f"_{lineage}"
    if sample.endswith(suffix):
        sample = sample[: -len(suffix)]

    busco_match = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%",
        text,
    )
    if not busco_match:
        raise ValueError(f"Could not parse BUSCO percentages from {path}")

    def int_field(pattern: str):
        match = re.search(pattern, text, re.M)
        return int(match.group(1)) if match else ""

    n50_text, n50_bp = parse_n50(text)

    return {
        "sample": sample,
        "complete_percent": busco_match.group(1),
        "single_copy_percent": busco_match.group(2),
        "duplicated_percent": busco_match.group(3),
        "fragmented_percent": busco_match.group(4),
        "missing_percent": busco_match.group(5),
        "complete_count": int_field(r"^\s*(\d+)\s+Complete BUSCOs"),
        "single_copy_count": int_field(r"^\s*(\d+)\s+Complete and single-copy"),
        "duplicated_count": int_field(r"^\s*(\d+)\s+Complete and duplicated"),
        "fragmented_count": int_field(r"^\s*(\d+)\s+Fragmented"),
        "missing_count": int_field(r"^\s*(\d+)\s+Missing"),
        "scaffold_count": int_field(r"^\s*(\d+)\s+Number of scaffolds"),
        "contig_count": int_field(r"^\s*(\d+)\s+Number of contigs"),
        "total_length_bp": int_field(r"^\s*(\d+)\s+Total length"),
        "percent_gaps": (
            re.search(r"^\s*([\d.]+)%\s+Percent gaps", text, re.M).group(1)
            if re.search(r"^\s*([\d.]+)%\s+Percent gaps", text, re.M)
            else ""
        ),
        "scaffold_n50": n50_text,
        "scaffold_n50_bp": n50_bp,
        "summary_file": path,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--busco-outbase", required=True, help="BUSCO output base directory")
    parser.add_argument("--lineage", default="ascomycota_odb12")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    pattern = os.path.join(
        args.busco_outbase,
        "*",
        f"short_summary.specific.{args.lineage}.*.txt",
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"No BUSCO summary files found with pattern: {pattern}")

    rows = [parse_summary(path, args.lineage) for path in files]
    headers = list(rows[0].keys())

    Path(args.outfile).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outfile, "w", encoding="utf-8") as out:
        out.write("\t".join(headers) + "\n")
        for row in rows:
            out.write("\t".join(str(row[h]) for h in headers) + "\n")

    print(f"Parsed {len(rows)} BUSCO summaries -> {args.outfile}")


if __name__ == "__main__":
    main()
