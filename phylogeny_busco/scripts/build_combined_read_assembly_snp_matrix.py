#!/usr/bin/env python3
"""Build a combined read-backed plus assembly-backed SNP matrix.

Read-backed samples come from a filtered biallelic SNP VCF. Assembly-backed
samples are projected onto the same SNP coordinates from assembly-to-reference
BAMs; non-covered or conflicting assembly bases are encoded as N.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
from collections import Counter
from pathlib import Path


IUPAC = {
    frozenset(["A", "C"]): "M",
    frozenset(["A", "G"]): "R",
    frozenset(["A", "T"]): "W",
    frozenset(["C", "G"]): "S",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "T"]): "K",
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open(encoding="utf-8", errors="replace", newline="")


def gt_to_base(gt: str, ref: str, alt: str) -> str:
    if gt in {".", "./.", ".|."} or gt.startswith("./") or gt.startswith(".|") or "." in gt:
        return "N"
    alleles = re.split(r"[|/]", gt)
    bases = []
    for allele in alleles:
        if allele == "0":
            bases.append(ref)
        elif allele == "1":
            bases.append(alt)
        else:
            return "N"
    unique = frozenset(bases)
    if len(unique) == 1:
        return next(iter(unique))
    return IUPAC.get(unique, "N")


def read_vcf(vcf: Path):
    samples: list[str] = []
    sites: list[dict[str, str | int]] = []
    seqs: dict[str, list[str]] = {}
    with open_text(vcf) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                seqs = {sample: [] for sample in samples}
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            chrom, pos_s, _id, ref, alt = fields[:5]
            if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"}:
                continue
            if "," in alt:
                continue
            pos = int(pos_s)
            sites.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "id": f"{chrom}:{pos}:{ref}>{alt}"})
            for sample, value in zip(samples, fields[9 : 9 + len(samples)]):
                gt = value.split(":", 1)[0]
                seqs[sample].append(gt_to_base(gt, ref, alt))
    if not samples:
        raise RuntimeError(f"No #CHROM header found in {vcf}")
    return samples, sites, seqs


def load_assembly_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        rows = [{k.strip().strip("\ufeff"): (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]
    rows = [row for row in rows if row.get("decision", "").startswith("include") and row.get("bam")]
    role_order = {"ingroup_assembly": 0, "outgroup_assembly": 1}
    return sorted(rows, key=lambda row: (role_order.get(row.get("role", ""), 9), row.get("sample", "")))


def write_sites_bed(sites: list[dict[str, str | int]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for site in sites:
            pos = int(site["pos"])
            handle.write(f"{site['chrom']}\t{pos - 1}\t{pos}\n")


def parse_pileup_bases(bases: str, ref: str) -> tuple[str, int, str]:
    counts: Counter[str] = Counter()
    skipped = 0
    i = 0
    while i < len(bases):
        char = bases[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            i += 1
            start = i
            while i < len(bases) and bases[i].isdigit():
                i += 1
            length = int(bases[start:i] or "0")
            i += length
            continue
        if char in ".,":
            counts[ref] += 1
        elif char in "ACGTacgt":
            counts[char.upper()] += 1
        elif char in "*<>Nn":
            skipped += 1
        else:
            skipped += 1
        i += 1

    informative = sum(counts.values())
    observed = [base for base, count in counts.items() if count > 0]
    if informative == 0:
        return "N", informative, "no_base"
    if skipped:
        return "N", informative, "deletion_or_skip"
    if len(observed) == 1:
        return observed[0], informative, "called"
    return "N", informative, "conflict"


def assembly_sequence(sample: str, bam: Path, ref: Path, bed: Path, sites: list[dict[str, str | int]]):
    command = [
        "samtools",
        "mpileup",
        "-aa",
        "-A",
        "-B",
        "-Q",
        "0",
        "-q",
        "0",
        "-f",
        str(ref),
        "-l",
        str(bed),
        str(bam),
    ]
    by_site: dict[tuple[str, int], tuple[str, int, str]] = {}
    proc = subprocess.run(command, check=True, text=True, capture_output=True)
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos_s, ref_base, depth_s, pile_bases = fields[:5]
        base, informative_depth, status = parse_pileup_bases(pile_bases, ref_base.upper())
        if int(depth_s) == 0:
            base, status = "N", "no_coverage"
        by_site[(chrom, int(pos_s))] = (base, informative_depth, status)

    seq: list[str] = []
    status_counts: Counter[str] = Counter()
    informative_depth_sum = 0
    informative_depth_n = 0
    for site in sites:
        key = (str(site["chrom"]), int(site["pos"]))
        base, informative_depth, status = by_site.get(key, ("N", 0, "not_reported"))
        seq.append(base)
        status_counts[status] += 1
        if informative_depth:
            informative_depth_sum += informative_depth
            informative_depth_n += 1
    mean_depth = informative_depth_sum / informative_depth_n if informative_depth_n else 0
    return seq, status_counts, mean_depth


def sample_qc(name: str, seq: list[str]) -> dict[str, str | int | float]:
    counts = Counter(seq)
    called = sum(counts[base] for base in "ACGT")
    missing = counts["N"]
    ambiguity = sum(counts[base] for base in "MRWSYK")
    total = len(seq)
    return {
        "sample": name,
        "sites": total,
        "called_acgt": called,
        "ambiguity_iupac": ambiguity,
        "missing_N": missing,
        "missing_pct": f"{(100 * missing / total) if total else 0:.4f}",
    }


def write_fasta(path: Path, seqs: dict[str, list[str]], sample_order: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for sample in sample_order:
            sequence = "".join(seqs[sample])
            handle.write(f">{sample}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def write_phylip(path: Path, seqs: dict[str, list[str]], sample_order: list[str]) -> None:
    length = len(seqs[sample_order[0]]) if sample_order else 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"{len(sample_order)} {length}\n")
        for sample in sample_order:
            handle.write(f"{sample.replace(' ', '_')} {''.join(seqs[sample])}\n")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--assembly-manifest", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--prefix", default="finger_millet.read90_assembly9")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    read_samples, sites, seqs = read_vcf(args.vcf)
    if not sites:
        raise RuntimeError(f"No biallelic SNP sites loaded from {args.vcf}")

    bed = args.outdir / f"{args.prefix}.sites.bed"
    write_sites_bed(sites, bed)

    assembly_rows = load_assembly_manifest(args.assembly_manifest)
    assembly_qc_rows = []
    for row in assembly_rows:
        sample = row["sample"]
        bam = Path(row["bam"])
        if not bam.exists():
            raise FileNotFoundError(bam)
        seq, status_counts, mean_depth = assembly_sequence(sample, bam, args.reference, bed, sites)
        seqs[sample] = seq
        qc = sample_qc(sample, seq)
        qc.update(
            {
                "role": row.get("role", ""),
                "preset": row.get("preset", ""),
                "mean_informative_pileup_depth": f"{mean_depth:.4f}",
                "called_status": status_counts.get("called", 0),
                "no_coverage_status": status_counts.get("no_coverage", 0),
                "conflict_status": status_counts.get("conflict", 0),
                "deletion_or_skip_status": status_counts.get("deletion_or_skip", 0),
            }
        )
        assembly_qc_rows.append(qc)

    sample_order = read_samples + [row["sample"] for row in assembly_rows]
    read_qc_rows = [dict(sample_qc(sample, seqs[sample]), role="read_backed", preset="", mean_informative_pileup_depth="") for sample in read_samples]

    fasta = args.outdir / f"{args.prefix}.fasta"
    phy = args.outdir / f"{args.prefix}.phy"
    site_table = args.outdir / f"{args.prefix}.sites.tsv"
    sample_qc_table = args.outdir / f"{args.prefix}.sample_qc.tsv"
    summary = args.outdir / f"{args.prefix}.summary.tsv"

    write_fasta(fasta, seqs, sample_order)
    write_phylip(phy, seqs, sample_order)
    write_tsv(site_table, sites)
    write_tsv(sample_qc_table, read_qc_rows + assembly_qc_rows)
    write_tsv(
        summary,
        [
            {"metric": "read_samples", "value": len(read_samples)},
            {"metric": "assembly_samples", "value": len(assembly_rows)},
            {"metric": "total_samples", "value": len(sample_order)},
            {"metric": "sites", "value": len(sites)},
            {"metric": "fasta", "value": str(fasta)},
            {"metric": "phylip", "value": str(phy)},
            {"metric": "sample_qc", "value": str(sample_qc_table)},
        ],
    )
    print(f"Wrote {phy}")
    print(f"samples={len(sample_order)} sites={len(sites)}")


if __name__ == "__main__":
    main()
