#!/usr/bin/env python3
"""Resolve missing phylogeny samples to public read or assembly downloads.

The input is the local staging plan, where rows marked input_method=missing
carry the sample/isolate names and public accessions from SRA_list_v4.
Outputs are TSV manifests that can be used on the cluster.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


ENA_FIELDS = [
    "run_accession",
    "sample_accession",
    "study_accession",
    "instrument_platform",
    "instrument_model",
    "library_layout",
    "fastq_ftp",
    "submitted_ftp",
    "submitted_format",
    "read_count",
    "base_count",
    "scientific_name",
]


def fetch_text(url: str, tries: int = 3) -> str:
    last_exc: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "sample-resolver/1.0"})
            with urllib.request.urlopen(req, timeout=60) as handle:
                return handle.read().decode("utf-8", errors="replace")
        except Exception as exc:  # pragma: no cover - operational retries
            last_exc = exc
            if attempt < tries - 1:
                time.sleep(2 + attempt)
    raise RuntimeError(f"failed to fetch {url}: {last_exc}")


def ena_read_runs(accession: str) -> list[dict[str, str]]:
    query = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "download": "true",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(query)
    text = fetch_text(url)
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) <= 1:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def ncbi_runinfo(accession: str) -> list[dict[str, str]]:
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?" + urllib.parse.urlencode(
        {"acc": accession}
    )
    text = fetch_text(url)
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) <= 1 or lines[0].startswith("Error"):
        return []
    return list(csv.DictReader(lines))


def ncbi_assembly_summaries(accession: str) -> list[dict[str, str]]:
    term = f"{accession}[All Fields]"
    esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "assembly", "term": term, "retmode": "json"}
    )
    search = json.loads(fetch_text(esearch))
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "assembly", "id": ",".join(ids), "retmode": "json"}
    )
    summary = json.loads(fetch_text(esummary))
    result = summary.get("result", {})
    rows = []
    for uid in result.get("uids", []):
        item = result.get(uid, {})
        ftp_path = item.get("ftppath_refseq") or item.get("ftppath_genbank") or ""
        if ftp_path:
            prefix = ftp_path.rsplit("/", 1)[-1]
            genomic_fna = f"{ftp_path}/{prefix}_genomic.fna.gz"
        else:
            genomic_fna = ""
        rows.append(
            {
                "assembly_uid": uid,
                "assembly_accession": item.get("assemblyaccession", ""),
                "assembly_name": item.get("assemblyname", ""),
                "organism": item.get("organism", ""),
                "ftp_path": ftp_path,
                "genomic_fna_ftp": genomic_fna,
            }
        )
    return rows


def http_url(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    if value.startswith(("http://", "https://", "ftp://")):
        return value
    return "https://" + value


def read_missing(plan: Path) -> list[dict[str, str]]:
    with plan.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row for row in rows if row.get("input_method") == "missing"]


def main() -> int:
    repo = Path.cwd()
    plan = repo / "results" / "staging_plan_96_plus_CKF3589.tsv"
    out_reads = repo / "results" / "missing_public_read_downloads.tsv"
    out_asm = repo / "results" / "missing_public_assembly_downloads.tsv"
    out_status = repo / "results" / "missing_public_resolution_status.tsv"

    read_rows: list[dict[str, str]] = []
    asm_rows: list[dict[str, str]] = []
    status_rows: list[dict[str, str]] = []

    for sample in read_missing(plan):
        accession = sample.get("accession", "").strip()
        sample_name = sample.get("sample", "").strip()
        isolate_id = sample.get("isolate_id", "").strip()
        status = {
            "sample": sample_name,
            "isolate_id": isolate_id,
            "accession": accession,
            "ena_read_rows": "0",
            "ncbi_runinfo_rows": "0",
            "assembly_rows": "0",
            "resolution": "unresolved",
            "notes": "",
        }

        try:
            ena_rows = ena_read_runs(accession) if accession else []
        except Exception as exc:
            ena_rows = []
            status["notes"] += f"ENA error: {exc}; "
        status["ena_read_rows"] = str(len(ena_rows))

        if ena_rows:
            status["resolution"] = "read_fastq"
            for row in ena_rows:
                fastq_urls = [http_url(x) for x in row.get("fastq_ftp", "").split(";") if x.strip()]
                submitted_urls = [http_url(x) for x in row.get("submitted_ftp", "").split(";") if x.strip()]
                urls = fastq_urls if fastq_urls else submitted_urls
                read_rows.append(
                    {
                        "sample": sample_name,
                        "isolate_id": isolate_id,
                        "source_accession": accession,
                        "run_accession": row.get("run_accession", ""),
                        "sample_accession": row.get("sample_accession", ""),
                        "study_accession": row.get("study_accession", ""),
                        "platform": row.get("instrument_platform", ""),
                        "instrument": row.get("instrument_model", ""),
                        "layout": row.get("library_layout", ""),
                        "url_1": urls[0] if len(urls) > 0 else "",
                        "url_2": urls[1] if len(urls) > 1 else "",
                        "url_3": urls[2] if len(urls) > 2 else "",
                        "read_count": row.get("read_count", ""),
                        "base_count": row.get("base_count", ""),
                    }
                )
        else:
            try:
                runinfo_rows = ncbi_runinfo(accession) if accession else []
            except Exception as exc:
                runinfo_rows = []
                status["notes"] += f"NCBI runinfo error: {exc}; "
            status["ncbi_runinfo_rows"] = str(len(runinfo_rows))
            if runinfo_rows:
                status["resolution"] = "read_run_no_direct_fastq_url"
                status["notes"] += "NCBI runinfo found run metadata but no direct URLs in runinfo; "

        try:
            assemblies = ncbi_assembly_summaries(accession) if accession else []
        except Exception as exc:
            assemblies = []
            status["notes"] += f"NCBI assembly error: {exc}; "
        status["assembly_rows"] = str(len(assemblies))

        for assembly in assemblies:
            asm_rows.append(
                {
                    "sample": sample_name,
                    "isolate_id": isolate_id,
                    "source_accession": accession,
                    **assembly,
                    "genomic_fna_url": http_url(assembly.get("genomic_fna_ftp", "")),
                }
            )
        if status["resolution"] == "unresolved" and assemblies:
            status["resolution"] = "assembly"

        status_rows.append(status)

    def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(
        out_reads,
        read_rows,
        [
            "sample",
            "isolate_id",
            "source_accession",
            "run_accession",
            "sample_accession",
            "study_accession",
            "platform",
            "instrument",
            "layout",
            "url_1",
            "url_2",
            "url_3",
            "read_count",
            "base_count",
        ],
    )
    write_tsv(
        out_asm,
        asm_rows,
        [
            "sample",
            "isolate_id",
            "source_accession",
            "assembly_uid",
            "assembly_accession",
            "assembly_name",
            "organism",
            "ftp_path",
            "genomic_fna_ftp",
            "genomic_fna_url",
        ],
    )
    write_tsv(
        out_status,
        status_rows,
        [
            "sample",
            "isolate_id",
            "accession",
            "ena_read_rows",
            "ncbi_runinfo_rows",
            "assembly_rows",
            "resolution",
            "notes",
        ],
    )

    print(f"wrote {out_reads} ({len(read_rows)} read/run rows)")
    print(f"wrote {out_asm} ({len(asm_rows)} assembly rows)")
    print(f"wrote {out_status} ({len(status_rows)} sample statuses)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
