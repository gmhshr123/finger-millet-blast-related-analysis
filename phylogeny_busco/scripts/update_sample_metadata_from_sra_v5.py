#!/usr/bin/env python3
"""Update local tree sample manifests from SRA_list_v5_KMDChecked.xlsx."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


SPECIAL_BY_VCF_SAMPLE = {
    "CKF3589": "DsCKF3589",
    "DS0505": "DsDS0505",
    "W97-11": "DsW97-11",
}

DISPLAY_OVERRIDES = {
    "Ei_EI9411_F1": "EiEI9411",
    "SiMG03": "SiMG08",
}


def clean(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    return " ".join(text.replace("\xa0", " ").strip().split())


def load_updates(xlsx_path: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    df = pd.read_excel(xlsx_path, dtype=str).fillna("")
    df = df[df["Isolate ID"].map(clean).ne("")]
    df = df[~df["Isolate ID"].map(clean).str.startswith("Added by")]

    by_isolate: dict[str, dict[str, str]] = {}
    by_accession: dict[str, dict[str, str]] = {}
    for _, raw in df.iterrows():
        row = {
            "isolate_id": clean(raw["Isolate ID"]),
            "id_on_tree": clean(raw["ID on tree"]),
            "host": clean(raw["Host"]),
            "location": clean(raw["Location"]),
            "accession": clean(raw["SRA Accession No"]),
            "reference": clean(raw["Reference"]),
            "notes": clean(raw["Notes"]),
        }
        by_isolate[row["isolate_id"]] = row
        if row["accession"]:
            by_accession[row["accession"]] = row
    return by_isolate, by_accession


def find_update(row: dict[str, str], by_isolate: dict[str, dict[str, str]], by_accession: dict[str, dict[str, str]]):
    sample = row.get("vcf_sample", "")
    if sample in SPECIAL_BY_VCF_SAMPLE:
        return by_isolate.get(SPECIAL_BY_VCF_SAMPLE[sample])

    for key in ("accession", "sra"):
        value = clean(row.get(key, ""))
        if value and value in by_accession:
            return by_accession[value]

    isolate = clean(row.get("isolate_id", ""))
    if isolate in by_isolate:
        return by_isolate[isolate]
    return None


def display_id_for(row: dict[str, str], update: dict[str, str]) -> str:
    sample = row.get("vcf_sample", "")
    if sample in DISPLAY_OVERRIDES:
        return DISPLAY_OVERRIDES[sample]
    if update.get("id_on_tree"):
        return update["id_on_tree"]
    return update["isolate_id"]


def update_manifest(path: Path, by_isolate: dict[str, dict[str, str]], by_accession: dict[str, dict[str, str]]):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    if "display_id" not in fieldnames:
        fieldnames.append("display_id")

    changes: list[dict[str, str]] = []
    for row in rows:
        before = {key: row.get(key, "") for key in ("isolate_id", "host", "location", "accession", "display_id")}
        update = find_update(row, by_isolate, by_accession)
        if not update:
            row.setdefault("display_id", row.get("isolate_id", ""))
            continue

        row["isolate_id"] = update["isolate_id"]
        row["host"] = update["host"]
        row["location"] = update["location"]
        if update["accession"]:
            row["accession"] = update["accession"]
        row["display_id"] = display_id_for(row, update)

        after = {key: row.get(key, "") for key in ("isolate_id", "host", "location", "accession", "display_id")}
        if before != after:
            changes.append(
                {
                    "manifest": str(path),
                    "vcf_sample": row.get("vcf_sample", ""),
                    "before_isolate_id": before["isolate_id"],
                    "after_isolate_id": after["isolate_id"],
                    "before_host": before["host"],
                    "after_host": after["host"],
                    "before_location": before["location"],
                    "after_location": after["location"],
                    "before_display_id": before["display_id"],
                    "after_display_id": after["display_id"],
                }
            )

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return changes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, default=Path("SRA_list_v5_KMDChecked.xlsx"))
    parser.add_argument("--report", type=Path, default=Path("01_from_scratch_reanalysis/results/manifests_and_qc/sra_v5_metadata_update_report.tsv"))
    parser.add_argument("manifests", nargs="+", type=Path)
    args = parser.parse_args()

    by_isolate, by_accession = load_updates(args.xlsx)
    all_changes: list[dict[str, str]] = []
    for manifest in args.manifests:
        all_changes.extend(update_manifest(manifest, by_isolate, by_accession))

    args.report.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "manifest",
        "vcf_sample",
        "before_isolate_id",
        "after_isolate_id",
        "before_host",
        "after_host",
        "before_location",
        "after_location",
        "before_display_id",
        "after_display_id",
    ]
    with args.report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_changes)
    print(f"Updated {len(args.manifests)} manifests; wrote {len(all_changes)} changed rows to {args.report}")


if __name__ == "__main__":
    main()
