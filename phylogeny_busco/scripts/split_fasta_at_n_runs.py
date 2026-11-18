#!/usr/bin/env python3
"""Split FASTA records at N runs for BUSCO validation."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def write_wrapped(handle, seq: str, width: int = 80) -> None:
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def split_record(name: str, seq: str, min_len: int):
    start = 0
    part = 1
    for match in re.finditer(r"N+", seq.upper()):
        end = match.start()
        if end - start >= min_len:
            yield f"{name}_part{part:06d}_{start + 1}_{end}", seq[start:end]
            part += 1
        start = match.end()
    if len(seq) - start >= min_len:
        yield f"{name}_part{part:06d}_{start + 1}_{len(seq)}", seq[start:]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-len", type=int, default=200)
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_records = 0
    total_parts = 0
    total_bases = 0

    with open(args.input, "r", encoding="utf-8") as inp, open(
        out_path, "w", encoding="utf-8"
    ) as out:
        name = None
        chunks: list[str] = []

        def flush() -> None:
            nonlocal total_records, total_parts, total_bases
            if name is None:
                return
            total_records += 1
            seq = "".join(chunks)
            for part_name, part_seq in split_record(name, seq, args.min_len):
                out.write(f">{part_name}\n")
                write_wrapped(out, part_seq)
                total_parts += 1
                total_bases += len(part_seq)

        for raw_line in inp:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        flush()

    print(
        f"records={total_records}\tparts={total_parts}\t"
        f"non_n_bases_written={total_bases}\tmin_len={args.min_len}"
    )


if __name__ == "__main__":
    main()
