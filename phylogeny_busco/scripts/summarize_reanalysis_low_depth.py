#!/usr/bin/env python3
"""List read-backed samples below a target mean-depth threshold."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def summarize_coverage(path: Path) -> tuple[float, float, int]:
    total = 0
    covered = 0
    reads = 0
    depth_sum = 0.0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length = int(row["endpos"]) - int(row["startpos"]) + 1
            total += length
            covered += int(row["covbases"])
            depth_sum += float(row["meandepth"]) * length
            reads += int(row["numreads"])
    if not total:
        return 0.0, 0.0, reads
    return covered / total, depth_sum / total, reads


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="/path/to/scratch/finger_millet_blast_phylogeny")
    parser.add_argument("--threshold", type=float, default=20.0)
    args = parser.parse_args()

    qc_dir = Path(args.base) / "reanalysis_deepvariant" / "qc" / "bam"
    rows = []
    for path in sorted(qc_dir.glob("*.coverage.tsv")):
        sample = path.name.removesuffix(".coverage.tsv")
        coverage, depth, reads = summarize_coverage(path)
        if depth < args.threshold:
            rows.append((depth, coverage, reads, sample))
    for depth, coverage, reads, sample in sorted(rows):
        print(f"{sample}\tcoverage={coverage:.4f}\tmean_depth={depth:.2f}\tnumreads={reads}")


if __name__ == "__main__":
    main()
