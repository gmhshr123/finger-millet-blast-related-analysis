#!/usr/bin/env python3
"""Build a per-sample mapping manifest from staged reanalysis inputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FASTQ_SUFFIXES = (".fq", ".fastq", ".fq.gz", ".fastq.gz")
ASSEMBLY_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")


def read_key_value(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            values[key] = value
    return values


def is_fastq(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(FASTQ_SUFFIXES) and not name.endswith(".done")


def is_assembly(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(ASSEMBLY_SUFFIXES) and "genomic" in name and not name.endswith(".done")


def pair_fastqs(paths: list[Path]) -> tuple[list[tuple[Path, Path]], list[Path]]:
    by_name = {p.name: p for p in paths}
    used: set[Path] = set()
    pairs: list[tuple[Path, Path]] = []

    tokens = [
        ("_1.fastq.gz", "_2.fastq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
        ("_1.fastq", "_2.fastq"),
        ("_1.fq", "_2.fq"),
        (".1.fastq.gz", ".2.fastq.gz"),
        (".1.fq.gz", ".2.fq.gz"),
        (".1.fastq", ".2.fastq"),
        (".1.fq", ".2.fq"),
    ]
    for path in sorted(paths):
        if path in used:
            continue
        mate = None
        for left, right in tokens:
            if path.name.endswith(left):
                mate_name = path.name[: -len(left)] + right
                mate = by_name.get(mate_name)
                break
        if mate is not None:
            used.add(path)
            used.add(mate)
            pairs.append((path, mate))

    singles = [p for p in sorted(paths) if p not in used]
    return pairs, singles


def encode_paths(paths: list[Path]) -> str:
    return ";".join(str(path) for path in paths)


def encode_pairs(pairs: list[tuple[Path, Path]]) -> str:
    return ";".join(f"{left},{right}" for left, right in pairs)


def build_manifest(base: Path) -> list[dict[str, str]]:
    samples_dir = base / "reanalysis_inputs" / "samples"
    bam_dir = base / "reanalysis_deepvariant" / "bam"
    gvcf_dir = base / "reanalysis_deepvariant" / "gvcf"
    vcf_dir = base / "reanalysis_deepvariant" / "vcf"
    rows: list[dict[str, str]] = []

    for row_id, sample_dir in enumerate(sorted(p for p in samples_dir.iterdir() if p.is_dir()), start=1):
        sample = sample_dir.name
        short_pairs: list[tuple[Path, Path]] = []
        short_singles: list[Path] = []
        long_reads: list[Path] = []
        notes: list[str] = []

        root_fastqs = [p for p in sample_dir.iterdir() if p.is_file() and is_fastq(p)]
        root_pairs, root_singles = pair_fastqs(root_fastqs)
        if root_pairs or root_singles:
            short_pairs.extend(root_pairs)
            short_singles.extend(root_singles)
            notes.append("root_fastq_assumed_short_read")

        public_reads = sample_dir / "public_reads"
        if public_reads.exists():
            for run_dir in sorted(p for p in public_reads.iterdir() if p.is_dir()):
                info = read_key_value(run_dir / "run_info.tsv")
                platform = info.get("platform", "").upper()
                run_fastqs = [p for p in run_dir.iterdir() if p.is_file() and is_fastq(p)]
                pairs, singles = pair_fastqs(run_fastqs)
                if "ILLUMINA" in platform:
                    short_pairs.extend(pairs)
                    short_singles.extend(singles)
                elif "OXFORD_NANOPORE" in platform or "PACBIO" in platform:
                    long_reads.extend(run_fastqs)
                else:
                    short_pairs.extend(pairs)
                    short_singles.extend(singles)
                    if run_fastqs:
                        notes.append(f"{run_dir.name}_unknown_platform_treated_as_short")

        source_bams = sorted(sample_dir.glob("*.bam"))
        assemblies = sorted(p for p in (sample_dir / "assembly").glob("*") if p.is_file() and is_assembly(p))

        if short_pairs or short_singles:
            input_kind = "short_reads"
            mapper = "bwa_mem2"
            deepvariant_model = "WGS"
            if long_reads:
                notes.append("long_reads_present_not_primary")
            if source_bams:
                notes.append("source_bam_present_not_primary")
        elif long_reads:
            input_kind = "long_reads"
            mapper = "minimap2_map_ont"
            deepvariant_model = "ONT_R104"
        elif source_bams:
            input_kind = "bam_extract_short_reads"
            mapper = "samtools_fastq_then_bwa_mem2"
            deepvariant_model = "WGS"
        elif assemblies:
            input_kind = "assembly"
            mapper = "minimap2_asm5"
            deepvariant_model = ""
        else:
            input_kind = "missing"
            mapper = "none"
            deepvariant_model = ""

        bam_name = (
            f"{sample}.target.markdup.bam"
            if input_kind in {"short_reads", "bam_extract_short_reads", "long_reads"}
            else ""
        )

        rows.append(
            {
                "id": str(row_id),
                "sample": sample,
                "input_kind": input_kind,
                "mapper": mapper,
                "deepvariant_model": deepvariant_model,
                "short_pairs": encode_pairs(short_pairs),
                "short_singles": encode_paths(short_singles),
                "long_reads": encode_paths(long_reads),
                "source_bams": encode_paths(source_bams),
                "assembly_fasta": str(assemblies[0]) if assemblies else "",
                "bam": str(bam_dir / bam_name) if bam_name else "",
                "gvcf": str(gvcf_dir / f"{sample}.dv.g.vcf.gz") if deepvariant_model else "",
                "vcf": str(vcf_dir / f"{sample}.dv.vcf.gz") if deepvariant_model else "",
                "notes": ";".join(notes),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="/path/to/scratch/finger_millet_blast_phylogeny")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    base = Path(args.base)
    out = Path(args.out) if args.out else base / "reanalysis_inputs" / "manifests" / "mapping_manifest.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = build_manifest(base)
    fieldnames = [
        "id",
        "sample",
        "input_kind",
        "mapper",
        "deepvariant_model",
        "short_pairs",
        "short_singles",
        "long_reads",
        "source_bams",
        "assembly_fasta",
        "bam",
        "gvcf",
        "vcf",
        "notes",
    ]
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["input_kind"]] = counts.get(row["input_kind"], 0) + 1
    for key in sorted(counts):
        print(f"{key}\t{counts[key]}")
    print(f"manifest\t{out}")


if __name__ == "__main__":
    main()
