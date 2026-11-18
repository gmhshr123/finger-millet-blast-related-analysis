#!/usr/bin/env python3
"""Build final read-backed and assembly-backed manifests after QC decisions."""

from __future__ import annotations

import csv
from pathlib import Path


BASE = Path("/path/to/scratch/finger_millet_blast_phylogeny")
MAPPING_MANIFEST = BASE / "reanalysis_inputs" / "manifests" / "mapping_manifest.tsv"
OUTDIR = BASE / "reanalysis_inputs" / "manifests"
READ_OUT = OUTDIR / "final_read_backed_samples.tsv"
ASM_OUT = OUTDIR / "final_assembly_backed_samples.tsv"
LOW_DEPTH_OUT = OUTDIR / "final_low_depth_read_backed_samples.tsv"
ASM_SUMMARY = BASE / "reanalysis_assembly_alignments" / "final_assembly_backed" / "final_assembly_backed_target_recovery.tsv"

READ_DEPTH_MIN = 10.0
ASSEMBLY_BACKED = {"B51", "G22", "PH42", "chw", "ho", "p29", "pg1213", "W97-11", "DS0505"}


def coverage_summary(path: Path) -> tuple[float, float, int]:
    total = covered = reads = 0
    depth_sum = 0.0
    if not path.exists():
        return 0.0, 0.0, 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length = int(row["endpos"]) - int(row["startpos"]) + 1
            total += length
            covered += int(row["covbases"])
            reads += int(row["numreads"])
            depth_sum += float(row["meandepth"]) * length
    if total == 0:
        return 0.0, 0.0, reads
    return covered / total, depth_sum / total, reads


def read_assembly_summary() -> dict[str, dict[str, str]]:
    summaries: dict[str, dict[str, str]] = {}
    if not ASM_SUMMARY.exists():
        return summaries
    with ASM_SUMMARY.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            summaries[row["sample"]] = row
    return summaries


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    mapping_rows = list(csv.DictReader(MAPPING_MANIFEST.open(), delimiter="\t"))
    by_sample = {row["sample"]: row for row in mapping_rows}
    asm_summaries = read_assembly_summary()

    read_rows = []
    low_depth_rows = []
    for row in mapping_rows:
        sample = row["sample"]
        if sample in ASSEMBLY_BACKED:
            continue
        if row["input_kind"] not in {"short_reads", "long_reads", "bam_extract_short_reads"}:
            continue
        cov, depth, reads = coverage_summary(BASE / "reanalysis_deepvariant" / "qc" / "bam" / f"{sample}.coverage.tsv")
        output = {
            "sample": sample,
            "input_kind": row["input_kind"],
            "mapper": row["mapper"],
            "deepvariant_model": row["deepvariant_model"],
            "bam": row["bam"],
            "gvcf": row["gvcf"],
            "vcf": row["vcf"],
            "coverage_fraction": f"{cov:.6f}",
            "mean_depth": f"{depth:.4f}",
            "numreads": str(reads),
        }
        output["decision"] = "include_read_backed"
        output["depth_flag"] = "review_low_depth_lt_10x" if depth < READ_DEPTH_MIN else "pass_depth_ge_10x"
        read_rows.append(output)
        if depth < READ_DEPTH_MIN:
            low_depth_rows.append(output)

    asm_rows = []
    for sample in sorted(ASSEMBLY_BACKED):
        summary = asm_summaries.get(sample, {})
        original = by_sample.get(sample, {})
        asm_rows.append(
            {
                "sample": sample,
                "role": summary.get("role", "outgroup_assembly" if sample in {"W97-11", "DS0505"} else "ingroup_assembly"),
                "preset": summary.get("preset", "asm20" if sample in {"W97-11", "DS0505"} else "asm5"),
                "assembly_fasta": summary.get("assembly_fasta", original.get("assembly_fasta", "")),
                "bam": summary.get("bam", ""),
                "coverage_fraction": summary.get("coverage_fraction", ""),
                "mean_depth": summary.get("mean_depth", ""),
                "num_alignments": summary.get("num_alignments", ""),
                "decision": "include_assembly_backed",
            }
        )

    read_fields = [
        "sample",
        "input_kind",
        "mapper",
        "deepvariant_model",
        "bam",
        "gvcf",
        "vcf",
        "coverage_fraction",
        "mean_depth",
        "numreads",
        "depth_flag",
        "decision",
    ]
    asm_fields = [
        "sample",
        "role",
        "preset",
        "assembly_fasta",
        "bam",
        "coverage_fraction",
        "mean_depth",
        "num_alignments",
        "decision",
    ]

    for path, fields, rows in [
        (READ_OUT, read_fields, sorted(read_rows, key=lambda r: r["sample"])),
        (LOW_DEPTH_OUT, read_fields, sorted(low_depth_rows, key=lambda r: float(r["mean_depth"]))),
        (ASM_OUT, asm_fields, asm_rows),
    ]:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    print(f"read_backed_include_now\t{len(read_rows)}\t{READ_OUT}")
    print(f"read_backed_low_depth_flagged\t{len(low_depth_rows)}\t{LOW_DEPTH_OUT}")
    print(f"assembly_backed\t{len(asm_rows)}\t{ASM_OUT}")


if __name__ == "__main__":
    main()
