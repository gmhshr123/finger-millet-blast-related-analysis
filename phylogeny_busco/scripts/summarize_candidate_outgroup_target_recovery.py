#!/usr/bin/env python3
"""Summarize candidate outgroup target recovery from samtools coverage files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def summarize(path: Path) -> tuple[float, float, int]:
    total_len = 0
    covered = 0
    depth_sum = 0.0
    alignments = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length = int(row["endpos"]) - int(row["startpos"]) + 1
            total_len += length
            covered += int(row["covbases"])
            depth_sum += float(row["meandepth"]) * length
            alignments += int(row["numreads"])
    if total_len == 0:
        return 0.0, 0.0, alignments
    return covered / total_len, depth_sum / total_len, alignments


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="/path/to/scratch/finger_millet_blast_phylogeny")
    args = parser.parse_args()

    out_root = Path(args.base) / "reanalysis_assembly_alignments" / "candidate_outgroups"
    rows = []
    for coverage in out_root.glob("*/*/*.coverage.tsv"):
        sample = coverage.parents[1].name
        preset = coverage.parent.name
        cov_frac, depth, alignments = summarize(coverage)
        bam = coverage.with_suffix("").with_suffix(".target.sorted.bam")
        rows.append((sample, preset, cov_frac, depth, alignments, str(bam)))

    out = out_root / "candidate_outgroup_target_recovery.corrected.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "preset", "coverage_fraction", "mean_depth", "num_alignments", "bam"])
        writer.writerows(sorted(rows, key=lambda x: (x[1], -x[2], x[0])))

    print(out)
    with out.open("r", encoding="utf-8") as handle:
        print(handle.read(), end="")


if __name__ == "__main__":
    main()
