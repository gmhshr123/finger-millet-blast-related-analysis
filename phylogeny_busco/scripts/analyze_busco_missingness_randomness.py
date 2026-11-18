#!/usr/bin/env python3
"""Test recurrence of missing BUSCO groups in the final RagTag cohort."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import statistics
import subprocess
import tarfile
from collections import Counter
from pathlib import Path, PurePosixPath

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = PROJECT / "results/manifests_and_qc/busco_final_ragtag_summary.tsv"
DEFAULT_SHARED = PROJECT / "results/manifests_and_qc/busco_pre216_vs_final_ragtag_comparison.tsv"
DEFAULT_E2_TABLE = PROJECT / "results/busco_E2_fullgenome_ascomycota_odb12/E2_fullgenome_ascomycota_odb12/run_ascomycota_odb12/full_table.tsv"
DEFAULT_OUT = PROJECT / "results/busco_missingness_randomness"
REMOTE_OUT = "/path/to/scratch/finger_millet_blast_phylogeny/results/busco_final_ragtag_ascomycota_odb12"
LINEAGE_SUFFIX = "_ascomycota_odb12"
EXPECTED_BUSCOS = 2826
STATES = ("Complete", "Duplicated", "Fragmented", "Missing")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_full_tables(ssh_host: str, remote_root: str) -> dict[str, str]:
    remote_command = (
        f"find {remote_root} -type f -path '*/run_ascomycota_odb12/full_table.tsv' -print0 "
        "| tar --null -T - -czf -"
    )
    result = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=30", ssh_host, remote_command],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    tables: dict[str, str] = {}
    with tarfile.open(fileobj=io.BytesIO(result.stdout), mode="r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith("/full_table.tsv"):
                continue
            parts = PurePosixPath(member.name).parts
            run_dir = parts[-3]
            if not run_dir.endswith(LINEAGE_SUFFIX):
                continue
            sample = run_dir[: -len(LINEAGE_SUFFIX)]
            stream = archive.extractfile(member)
            if stream is None:
                continue
            tables[sample] = stream.read().decode("utf-8")
    return tables


def parse_full_table(text: str, sample: str) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    statuses: dict[str, str] = {}
    header: list[str] | None = None
    annotations: dict[str, dict[str, str]] = {}
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    for fields in reader:
        if not fields:
            continue
        if fields[0].lstrip("# ").strip().lower() == "busco id":
            header = [field.lstrip("# ").strip() for field in fields]
            continue
        if fields[0].startswith("#") or header is None:
            continue
        row = dict(zip(header, fields))
        busco_id = row.get("Busco id", "").strip()
        status = row.get("Status", "").strip()
        if not busco_id or status not in STATES:
            raise ValueError(f"Malformed BUSCO row in {sample}: {fields[:3]}")
        previous = statuses.get(busco_id)
        if previous is not None and previous != status:
            raise ValueError(f"Conflicting statuses for BUSCO ID {busco_id} in {sample}: {previous}, {status}")
        statuses[busco_id] = status
        if sample.startswith("E2") and status != "Missing" and busco_id not in annotations:
            annotations[busco_id] = row
    if len(statuses) != EXPECTED_BUSCOS:
        raise ValueError(f"{sample}: expected {EXPECTED_BUSCOS} BUSCO rows, found {len(statuses)}")
    return statuses, annotations


def poisson_binomial_logpmf(probabilities: np.ndarray) -> np.ndarray:
    log_pmf = np.array([0.0])
    for probability in probabilities:
        next_log_pmf = np.full(len(log_pmf) + 1, -np.inf, dtype=float)
        if probability < 1.0:
            next_log_pmf[:-1] = log_pmf + math.log1p(-probability)
        if probability > 0.0:
            next_log_pmf[1:] = np.logaddexp(
                next_log_pmf[1:], log_pmf + math.log(probability)
            )
        log_pmf = next_log_pmf
    return log_pmf


def benjamini_hochberg_log(log_pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(log_pvalues)
    ranked = log_pvalues[order]
    adjusted = ranked + math.log(len(ranked)) - np.log(np.arange(1, len(ranked) + 1))
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 0.0)
    return result


def probability_for_tsv(log_probability: float) -> float:
    return math.exp(max(log_probability, math.log(np.finfo(float).tiny)))


def marker_statistics(
    sample_ids: list[str],
    matrix: np.ndarray,
    busco_ids: list[str],
    annotations: dict[str, dict[str, str]],
    e2_statuses: dict[str, str],
    dataset: str,
) -> list[dict[str, object]]:
    row_missing = matrix.sum(axis=1).astype(int)
    probabilities = row_missing / matrix.shape[1]
    log_pmf = poisson_binomial_logpmf(probabilities)
    log_upper = np.array([np.logaddexp.reduce(log_pmf[index:]) for index in range(len(log_pmf))])
    log_lower = np.array([np.logaddexp.reduce(log_pmf[: index + 1]) for index in range(len(log_pmf))])
    observed = matrix.sum(axis=0).astype(int)
    expected = float(row_missing.sum()) / matrix.shape[1]
    log_p_upper = np.array([log_upper[min(int(value), len(log_upper) - 1)] for value in observed])
    log_p_lower = np.array([log_lower[min(int(value), len(log_lower) - 1)] for value in observed])
    log_p_two = np.minimum(0.0, math.log(2.0) + np.minimum(log_p_upper, log_p_lower))
    log_q_upper = benjamini_hochberg_log(log_p_upper)
    log_q_two = benjamini_hochberg_log(log_p_two)
    rows: list[dict[str, object]] = []
    for index, busco_id in enumerate(busco_ids):
        annotation = annotations.get(busco_id, {})
        start = annotation.get("Gene Start", "")
        end = annotation.get("Gene End", "")
        rows.append(
            {
                "dataset": dataset,
                "busco_id": busco_id,
                "missing_samples": int(observed[index]),
                "sample_count": len(sample_ids),
                "expected_missing_samples": round(expected, 5),
                "fold_vs_expected": round(float(observed[index]) / expected, 4) if expected else "NA",
                "poisson_binomial_p_upper": probability_for_tsv(log_p_upper[index]),
                "poisson_binomial_p_upper_log10": float(log_p_upper[index] / math.log(10)),
                "poisson_binomial_q_upper_BH": probability_for_tsv(log_q_upper[index]),
                "poisson_binomial_q_upper_BH_log10": float(log_q_upper[index] / math.log(10)),
                "poisson_binomial_p_two_sided": probability_for_tsv(log_p_two[index]),
                "poisson_binomial_p_two_sided_log10": float(log_p_two[index] / math.log(10)),
                "poisson_binomial_q_two_sided_BH": probability_for_tsv(log_q_two[index]),
                "poisson_binomial_q_two_sided_BH_log10": float(log_q_two[index] / math.log(10)),
                "E2_status": e2_statuses.get(busco_id, "Not found"),
                "E2_sequence": annotation.get("Sequence", ""),
                "E2_gene_start": start,
                "E2_gene_end": end,
                "E2_strand": annotation.get("Strand", ""),
                "E2_reference_hit_length": annotation.get("Length", ""),
                "description": annotation.get("Description", ""),
            }
        )
    rows.sort(key=lambda row: (-int(row["missing_samples"]), str(row["busco_id"])))
    return rows


def permutation_test(matrix: np.ndarray, permutations: int, seed: int) -> dict[str, object]:
    n_samples, n_buscos = matrix.shape
    row_missing = matrix.sum(axis=1).astype(int)
    observed_counts = matrix.sum(axis=0).astype(int)
    expected_count = float(row_missing.sum()) / n_buscos
    observed_stat = float(np.square(observed_counts - expected_count).sum())
    rng = np.random.default_rng(seed)
    simulated = np.empty(permutations, dtype=float)
    for iteration in range(permutations):
        counts = np.zeros(n_buscos, dtype=np.int32)
        for k in row_missing:
            if k:
                counts[rng.choice(n_buscos, size=int(k), replace=False)] += 1
        simulated[iteration] = np.square(counts - expected_count).sum()
    exceed = int(np.count_nonzero(simulated >= observed_stat))
    null_var = float(np.var(simulated, ddof=1)) if permutations > 1 else 0.0
    return {
        "sample_count": n_samples,
        "busco_count": n_buscos,
        "total_missing_calls": int(row_missing.sum()),
        "mean_missing_per_busco_null": expected_count,
        "observed_dispersion_T": observed_stat,
        "null_mean_T": float(np.mean(simulated)),
        "observed_to_null_mean_ratio": observed_stat / float(np.mean(simulated)),
        "null_sd_T": math.sqrt(null_var),
        "null_95pct_T": float(np.quantile(simulated, 0.95)),
        "null_99pct_T": float(np.quantile(simulated, 0.99)),
        "empirical_p_upper": (exceed + 1) / (permutations + 1),
        "permutations": permutations,
        "seed": seed,
        "observed_buscos_missing_zero": int(np.count_nonzero(observed_counts == 0)),
        "observed_buscos_missing_once": int(np.count_nonzero(observed_counts == 1)),
        "observed_buscos_missing_ge_5": int(np.count_nonzero(observed_counts >= 5)),
        "observed_buscos_missing_ge_10": int(np.count_nonzero(observed_counts >= 10)),
        "observed_max_samples_missing_same_busco": int(observed_counts.max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ssh-host", default="user@cluster")
    parser.add_argument("--remote-output", default=REMOTE_OUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--shared-comparison", type=Path, default=DEFAULT_SHARED)
    parser.add_argument("--e2-full-table", type=Path, default=DEFAULT_E2_TABLE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261002)
    args = parser.parse_args()
    if args.permutations < 999:
        parser.error("Use at least 999 permutations; 10,000 is recommended.")

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_rows = read_tsv(args.summary)
    expected = {row["sample"]: row for row in summary_rows}
    shared_rows = read_tsv(args.shared_comparison)
    shared_ids = {
        row["sample"] + ".ragtag"
        for row in shared_rows
        if row.get("in_pre_216") == "yes" and row.get("in_final_ragtag") == "yes"
    }
    print(f"Fetching BUSCO full tables from {args.ssh_host}:{args.remote_output}", flush=True)
    remote_tables = get_full_tables(args.ssh_host, args.remote_output)
    table_ids = {sample for sample in remote_tables if sample in expected}
    missing_tables = sorted(set(expected) - table_ids)
    unexpected_tables = sorted(set(remote_tables) - set(expected))
    if missing_tables or len(table_ids) != len(expected):
        raise RuntimeError(
            f"Found {len(table_ids)}/{len(expected)} expected full tables; missing={missing_tables[:10]}"
        )

    parsed: dict[str, dict[str, str]] = {}
    e2_annotations: dict[str, dict[str, str]] = {}
    busco_ids: list[str] | None = None
    validation_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for sample in sorted(expected):
        statuses, _ = parse_full_table(remote_tables[sample], sample)
        ids = sorted(statuses)
        if busco_ids is None:
            busco_ids = ids
        elif ids != busco_ids:
            raise ValueError(f"BUSCO ID set differs for {sample}")
        parsed[sample] = statuses
        counts = Counter(statuses.values())
        summary = expected[sample]
        complete = counts["Complete"] + counts["Duplicated"]
        summary_complete = int(summary["complete_count"])
        summary_missing = int(summary["missing_count"])
        valid = complete == summary_complete and counts["Missing"] == summary_missing
        validation_rows.append(
            {
                "sample": sample,
                "valid_counts_vs_short_summary": valid,
                "full_table_complete": complete,
                "summary_complete": summary_complete,
                "full_table_missing": counts["Missing"],
                "summary_missing": summary_missing,
                "full_table_fragmented": counts["Fragmented"],
                "summary_fragmented": int(summary["fragmented_count"]),
                "full_table_single_copy": counts["Complete"],
                "full_table_duplicated": counts["Duplicated"],
            }
        )
        if not valid or counts["Fragmented"] != int(summary["fragmented_count"]):
            raise ValueError(f"Full table counts do not reproduce short summary for {sample}")
        sample_rows.append(
            {
                "sample": sample,
                "dataset_shared216": sample in shared_ids,
                "complete_percent": summary["complete_percent"],
                "missing_percent": summary["missing_percent"],
                "missing_count": counts["Missing"],
                "complete_count": counts["Complete"] + counts["Duplicated"],
                "fragmented_count": counts["Fragmented"],
                "duplicated_count": counts["Duplicated"],
            }
        )

    e2_text = args.e2_full_table.read_text(encoding="utf-8-sig")
    e2_statuses, e2_annotations = parse_full_table(e2_text, "E2_reference")
    assert busco_ids is not None
    samples = sorted(parsed)
    matrix = np.array(
        [[parsed[sample][busco_id] == "Missing" for busco_id in busco_ids] for sample in samples],
        dtype=np.uint8,
    )
    write_tsv(args.outdir / "sample_qc_and_counts.tsv", list(sample_rows[0]), sample_rows)
    write_tsv(args.outdir / "full_table_validation.tsv", list(validation_rows[0]), validation_rows)
    status_path = args.outdir / "busco_status_matrix.tsv"
    with status_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", *busco_ids])
        for sample in samples:
            writer.writerow([sample, *(parsed[sample][busco_id] for busco_id in busco_ids)])

    datasets: list[tuple[str, list[str]]] = [
        ("all_250", samples),
        ("shared_pre216_216", [sample for sample in samples if sample in shared_ids]),
        (
            "high_completeness_ge95pct",
            [sample for sample in samples if float(expected[sample]["complete_percent"]) >= 95.0],
        ),
    ]
    rng_seed_offsets = {"all_250": 0, "shared_pre216_216": 1, "high_completeness_ge95pct": 2}
    global_rows: list[dict[str, object]] = []
    marker_rows: list[dict[str, object]] = []
    print(f"Validated {len(samples)} full tables. Running {args.permutations:,} permutations per cohort.", flush=True)
    for name, subset in datasets:
        if not subset:
            raise ValueError(f"Empty sensitivity cohort: {name}")
        indices = [samples.index(sample) for sample in subset]
        subset_matrix = matrix[indices]
        result = permutation_test(subset_matrix, args.permutations, args.seed + rng_seed_offsets[name])
        result["dataset"] = name
        global_rows.append(result)
        marker_rows.extend(
            marker_statistics(subset, subset_matrix, busco_ids, e2_annotations, e2_statuses, name)
        )
        print(f"{name}: observed T={result['observed_dispersion_T']:.2f}; empirical p={result['empirical_p_upper']:.6g}", flush=True)

    write_tsv(args.outdir / "cohort_randomization_tests.tsv", list(global_rows[0]), global_rows)
    marker_fields = list(marker_rows[0])
    write_tsv(args.outdir / "marker_missingness_results.tsv", marker_fields, marker_rows)
    write_tsv(
        args.outdir / "top_recurrent_markers_all250.tsv",
        marker_fields,
        [row for row in marker_rows if row["dataset"] == "all_250"][:100],
    )

    status_counts = Counter(status for sample in samples for status in parsed[sample].values())
    high_n = len(datasets[2][1])
    shared_n = len(datasets[1][1])
    top = [row for row in marker_rows if row["dataset"] == "all_250"][:10]
    all250_by_id = {row["busco_id"]: row for row in marker_rows if row["dataset"] == "all_250"}
    high_by_id = {
        row["busco_id"]: row
        for row in marker_rows
        if row["dataset"] == "high_completeness_ge95pct"
    }
    recurrent_e2_detected = [
        row for row in all250_by_id.values() if row["E2_status"] != "Missing"
    ]
    recurrent_e2_detected.sort(key=lambda row: (-int(row["missing_samples"]), str(row["busco_id"])))
    all_missing_ids = [
        busco_ids[index]
        for index, missing_count in enumerate(matrix.sum(axis=0))
        if int(missing_count) == len(samples)
    ]
    all_missing_in_e2 = [busco_id for busco_id in all_missing_ids if e2_statuses.get(busco_id) == "Missing"]
    e2_missing_recurrence = [
        int(matrix[:, index].sum())
        for index, busco_id in enumerate(busco_ids)
        if e2_statuses.get(busco_id) == "Missing"
    ]
    report = [
        "# BUSCO Missingness Randomness Analysis",
        "",
        "Analysis of the 250 final RagTag assemblies with BUSCO 6.0.0 and `ascomycota_odb12`.",
        "",
        "## Data validation",
        "",
        f"- Parsed and validated {len(samples)} of {len(expected)} expected per-sample full tables; unexpected remote tables not in the summary: {len(unexpected_tables)}.",
        f"- Each sample had {EXPECTED_BUSCOS} unique BUSCO IDs, and all full-table Complete, Duplicated, Fragmented, and Missing counts agreed with its short summary.",
        f"- Total calls across the cohort: Complete {status_counts['Complete']:,}; Duplicated {status_counts['Duplicated']:,}; Fragmented {status_counts['Fragmented']:,}; Missing {status_counts['Missing']:,}.",
            f"- BUSCOs missing from all 250 RagTag assemblies: {len(all_missing_ids)}; among these, {len(all_missing_in_e2)} are also Missing in the E2 full-genome BUSCO table.",
        "",
        "## Randomization test",
        "",
        "For each assembly, the test held its observed number of Missing BUSCOs fixed and randomly reassigned those calls among the 2,826 BUSCO groups. The statistic is the sum of squared deviations in missing-call counts across groups. The empirical upper-tail p-value compares the observed statistic with the randomized distribution.",
        "",
        "| Cohort | Samples | Observed dispersion | Randomized mean | Observed / null mean | Randomized 99th percentile | Empirical p |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in global_rows:
        report.append(
            f"| {row['dataset']} | {row['sample_count']} | {row['observed_dispersion_T']:.2f} | {row['null_mean_T']:.2f} | {row['observed_to_null_mean_ratio']:.1f}x | {row['null_99pct_T']:.2f} | {row['empirical_p_upper']:.6g} |"
        )
    report.extend(
        [
            "",
            f"Permutations per cohort: {args.permutations:,}. Seed: {args.seed}; cohort offsets 0, 1, and 2 were used in the order shown above.",
            "",
            "## Most recurrent markers",
            "",
            "| BUSCO | Missing assemblies | E2 BUSCO status | E2 sequence | E2 coordinates | Description |",
            "|---|---:|---|---|---|---|",
        ]
    )
    for row in top:
        coords = f"{row['E2_gene_start']}-{row['E2_gene_end']}" if row["E2_gene_start"] else "NA"
        report.append(
            f"| {row['busco_id']} | {row['missing_samples']} | {row['E2_status']} | {row['E2_sequence'] or 'NA'} | {coords} | {row['description'] or 'NA'} |"
        )
    report.extend(
        [
            "",
            "Marker-level upper-tail and two-sided Poisson-binomial p-values and BH-adjusted values are in `marker_missingness_results.tsv`. These marker-level p-values treat sample rows as independent Bernoulli draws with their sample-specific missing fraction; interpret them cautiously because strains may share ancestry. The cohort-wide row-preserving randomization test is the primary result.",
            "",
            "A Missing BUSCO call does not prove biological gene deletion. Recurrence indicates marker-level non-random detection failure under this null; use the E2 coordinates, descriptions, sample identities, and assembly evidence to investigate causes.",
            "",
            f"The E2 full-genome table marks {len(e2_missing_recurrence)} of {EXPECTED_BUSCOS} groups as Missing. In the RagTag cohort, those E2-Missing markers recur in {min(e2_missing_recurrence)} to {max(e2_missing_recurrence)} assemblies (median {statistics.median(e2_missing_recurrence):g}); {sum(count >= 100 for count in e2_missing_recurrence)} occur as Missing in at least 100 assemblies. This overlap suggests a substantial shared reference/lineage detectability component. It does not by itself show whether the cause is true shared gene absence or systematic BUSCO detection failure.",
            "",
            "BUSCO IDs missing in every RagTag assembly: " + (", ".join(all_missing_ids) if all_missing_ids else "none") + ".",
            "Among these, IDs also Missing in E2: " + (", ".join(all_missing_in_e2) if all_missing_in_e2 else "none") + ".",
            "",
            "Some frequently missing groups are Complete in E2, so excess recurrence also affects groups the reference BUSCO run detected:",
            "",
            "| BUSCO | Missing in all 250 | Missing in >=95% cohort | E2 sequence | E2 coordinates | Description |",
            "|---|---:|---:|---|---|---|",
        ]
    )
    for row in recurrent_e2_detected[:5]:
        high_row = high_by_id[row["busco_id"]]
        coords = f"{row['E2_gene_start']}-{row['E2_gene_end']}" if row["E2_gene_start"] else "NA"
        report.append(
            f"| {row['busco_id']} | {row['missing_samples']} | {high_row['missing_samples']} | {row['E2_sequence'] or 'NA'} | {coords} | {row['description'] or 'NA'} |"
        )
    report.extend(
        [
            "",
            "Very small marker-level p-values are represented by the smallest positive floating-point value in the numeric p/q columns; use the accompanying `_log10` columns for their actual scale.",
            "",
            "## Files",
            "",
            "- `cohort_randomization_tests.tsv`: primary and sensitivity tests.",
            "- `marker_missingness_results.tsv`: all markers in all three cohorts, with E2 annotations.",
            "- `top_recurrent_markers_all250.tsv`: top 100 recurrent markers in the full cohort.",
            "- `busco_status_matrix.tsv`: validated categorical BUSCO states for each sample and marker.",
            "- `sample_qc_and_counts.tsv`: sample scores and sensitivity cohort membership.",
            "- `full_table_validation.tsv`: per-sample agreement check against BUSCO short summaries.",
            "",
            f"Shared pre-216 subset: {shared_n} assemblies. High-completeness subset (>=95% complete): {high_n} assemblies.",
        ]
    )
    (args.outdir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    metadata = {
        "ssh_host": args.ssh_host,
        "remote_output": args.remote_output,
        "summary": str(args.summary),
        "shared_comparison": str(args.shared_comparison),
        "e2_full_table": str(args.e2_full_table),
        "outdir": str(args.outdir),
        "samples": len(samples),
        "busco_groups": EXPECTED_BUSCOS,
        "permutations_per_cohort": args.permutations,
        "seed": args.seed,
        "unexpected_remote_tables": unexpected_tables,
    }
    (args.outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote report and result tables to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
