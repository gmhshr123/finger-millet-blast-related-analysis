#!/usr/bin/env python3
"""Audit mixed-input mapping outputs and logs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


BAD_PATTERNS = re.compile(
    r"(ERROR:|Error:|error:|failed|Failed|Cannot|No such file|Killed|Traceback|segmentation|exited with)",
    re.IGNORECASE,
)
BENIGN_PATTERNS = re.compile(
    r"(markdup_failed_falling_back|skip orientation|not enough pairs)",
    re.IGNORECASE,
)


def exists_nonempty(path: str) -> bool:
    return bool(path) and Path(path).is_file() and Path(path).stat().st_size > 0


def scan_log(path: Path, max_hits: int = 5) -> list[str]:
    hits: list[str] = []
    if not path.exists():
        return [f"MISSING_LOG\t{path}"]
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not BAD_PATTERNS.search(line):
                continue
            if BENIGN_PATTERNS.search(line):
                continue
            hits.append(f"{path.name}:{line_no}:{line.rstrip()[:250]}")
            if len(hits) >= max_hits:
                break
    return hits


def summarize_coverage(path: Path) -> tuple[int, float, float, int]:
    total_len = 0
    covered = 0
    depth_sum = 0.0
    reads = 0
    if not path.exists():
        return 0, 0.0, 0.0, 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            length = int(row["endpos"]) - int(row["startpos"]) + 1
            covbases = int(row["covbases"])
            meandepth = float(row["meandepth"])
            total_len += length
            covered += covbases
            depth_sum += meandepth * length
            reads += int(row["numreads"])
    coverage_fraction = covered / total_len if total_len else 0.0
    weighted_depth = depth_sum / total_len if total_len else 0.0
    return total_len, coverage_fraction, weighted_depth, reads


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="/path/to/scratch/finger_millet_blast_phylogeny")
    parser.add_argument("--job", default="45600964")
    args = parser.parse_args()

    base = Path(args.base)
    manifest = base / "reanalysis_inputs" / "manifests" / "mapping_manifest.tsv"
    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))

    counts: dict[str, int] = {}
    missing_outputs: list[str] = []
    suspicious_logs: list[str] = []
    read_qc: list[tuple[str, str, float, float, int]] = []
    asm_qc: list[tuple[str, float, float, int]] = []

    for row in rows:
        kind = row["input_kind"]
        counts[kind] = counts.get(kind, 0) + 1
        sample = row["sample"]
        task_id = row["id"]

        if kind in {"short_reads", "bam_extract_short_reads", "long_reads"}:
            bam = row["bam"]
            for label, path in [
                ("bam", bam),
                ("bai", f"{bam}.bai"),
                ("done", str(base / "reanalysis_deepvariant" / "bam" / f"{sample}.mapping.done")),
                ("flagstat", str(base / "reanalysis_deepvariant" / "qc" / "bam" / f"{sample}.flagstat.txt")),
                ("coverage", str(base / "reanalysis_deepvariant" / "qc" / "bam" / f"{sample}.coverage.tsv")),
            ]:
                if label == "done":
                    if not Path(path).exists():
                        missing_outputs.append(f"{sample}\t{kind}\tmissing_{label}\t{path}")
                elif not exists_nonempty(path):
                    missing_outputs.append(f"{sample}\t{kind}\tmissing_{label}\t{path}")
            total_len, cov_frac, depth, reads = summarize_coverage(
                base / "reanalysis_deepvariant" / "qc" / "bam" / f"{sample}.coverage.tsv"
            )
            if total_len:
                read_qc.append((sample, kind, cov_frac, depth, reads))
        elif kind == "assembly":
            asm_bam = base / "reanalysis_assembly_alignments" / "samples" / sample / f"{sample}.assembly_to_target.sorted.bam"
            for label, path in [
                ("asm_bam", asm_bam),
                ("asm_bai", Path(f"{asm_bam}.bai")),
                ("asm_paf", base / "reanalysis_assembly_alignments" / "samples" / sample / f"{sample}.assembly_to_target.paf"),
                ("asm_coverage", base / "reanalysis_assembly_alignments" / "samples" / sample / f"{sample}.coverage.tsv"),
            ]:
                if not path.is_file() or path.stat().st_size == 0:
                    missing_outputs.append(f"{sample}\t{kind}\tmissing_{label}\t{path}")
            total_len, cov_frac, depth, reads = summarize_coverage(
                base / "reanalysis_assembly_alignments" / "samples" / sample / f"{sample}.coverage.tsv"
            )
            if total_len:
                asm_qc.append((sample, cov_frac, depth, reads))

        log = base / "log" / f"fm_map_mixed_{args.job}_{task_id}.err"
        suspicious_logs.extend(scan_log(log))

    print("INPUT_KIND_COUNTS")
    for key in sorted(counts):
        print(f"{key}\t{counts[key]}")
    print("MISSING_OUTPUTS")
    if missing_outputs:
        print("\n".join(missing_outputs))
    else:
        print("none")
    print("SUSPICIOUS_LOG_LINES")
    if suspicious_logs:
        print("\n".join(suspicious_logs[:200]))
        if len(suspicious_logs) > 200:
            print(f"... {len(suspicious_logs) - 200} additional suspicious lines omitted")
    else:
        print("none")
    print("READ_QC_LOWEST_DEPTH")
    for sample, kind, cov_frac, depth, reads in sorted(read_qc, key=lambda x: x[3])[:15]:
        print(f"{sample}\t{kind}\tcoverage={cov_frac:.4f}\tmean_depth={depth:.2f}\tnumreads={reads}")
    print("READ_QC_DEPTH_COUNTS")
    for threshold in [1, 5, 10, 20]:
        n = sum(1 for _, _, _, depth, _ in read_qc if depth < threshold)
        print(f"mean_depth_lt_{threshold}\t{n}")
    print("READ_QC_COVERAGE_COUNTS")
    for threshold in [0.5, 0.8, 0.9]:
        n = sum(1 for _, _, cov_frac, _, _ in read_qc if cov_frac < threshold)
        print(f"coverage_lt_{threshold}\t{n}")
    print("ASSEMBLY_QC")
    for sample, cov_frac, depth, reads in sorted(asm_qc):
        print(f"{sample}\tcoverage={cov_frac:.4f}\tmean_depth={depth:.2f}\tnum_alignments={reads}")


if __name__ == "__main__":
    main()
