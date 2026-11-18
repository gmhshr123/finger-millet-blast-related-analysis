#!/usr/bin/env python3
"""Compute simple scaffold-count statistics for FASTA assemblies."""

from __future__ import annotations

import argparse
from pathlib import Path


def fasta_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    current = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    return lengths


def sample_name(path: Path, suffix: str) -> str:
    name = path.name
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    return path.stem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--glob", default="*.fasta")
    parser.add_argument("--sample-suffix", default="")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    paths = sorted(input_dir.glob(args.glob))
    if not paths:
        raise SystemExit(f"No FASTA files found: {input_dir}/{args.glob}")

    out = Path(args.outfile)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "sample",
                    "source_file",
                    "scaffold_count",
                    "scaffolds_gt_1000bp",
                    "scaffolds_ge_1000bp",
                    "total_length_bp",
                ]
            )
            + "\n"
        )
        for path in paths:
            lengths = fasta_lengths(path)
            handle.write(
                "\t".join(
                    [
                        sample_name(path, args.sample_suffix),
                        str(path),
                        str(len(lengths)),
                        str(sum(length > 1000 for length in lengths)),
                        str(sum(length >= 1000 for length in lengths)),
                        str(sum(lengths)),
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
