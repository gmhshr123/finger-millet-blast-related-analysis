from __future__ import annotations

import csv
import time
import urllib.parse
import urllib.request
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill


HERE = Path(__file__).resolve().parent
WORKBOOK = HERE / "Fosbury_FullLengthElements_LTR_identity_with_TSD.xlsx"
OUT_TSV = HERE / "Fosbury_LTR_pairwise_indel_analysis.tsv"
OUT_XLSX = HERE / "Fosbury_LTR_pairwise_indel_analysis.xlsx"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_fasta(text: str) -> str:
    return "".join(
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()


def fetch_sequence(seqid: str, start: int, end: int) -> str:
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": seqid,
            "seq_start": start,
            "seq_stop": end,
            "strand": 1,
            "rettype": "fasta",
            "retmode": "text",
            "tool": "fosbury_ltr_indel_analysis",
        }
    )
    url = f"{NCBI_EFETCH}?{params}"
    last_error = None
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                sequence = parse_fasta(response.read().decode("utf-8"))
            if not sequence:
                raise RuntimeError(f"empty sequence returned for {seqid}:{start}-{end}")
            time.sleep(0.35)
            return sequence
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(attempt)
    raise RuntimeError(f"failed to fetch {seqid}:{start}-{end}: {last_error}")


def global_edit_alignment(a: str, b: str):
    """Needleman-Wunsch alignment minimizing edit distance."""
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = i
        ptr[i][0] = "D"
    for j in range(1, m + 1):
        dp[0][j] = j
        ptr[0][j] = "I"

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            sub_cost = 0 if ai == b[j - 1] else 1
            choices = [
                (dp[i - 1][j - 1] + sub_cost, "M" if sub_cost == 0 else "S"),
                (dp[i - 1][j] + 1, "D"),
                (dp[i][j - 1] + 1, "I"),
            ]
            # Prefer matches/substitutions over gaps when edit cost ties.
            best_cost, best_op = min(choices, key=lambda x: (x[0], 0 if x[1] in {"M", "S"} else 1))
            dp[i][j] = best_cost
            ptr[i][j] = best_op

    ops = []
    i, j = n, m
    while i > 0 or j > 0:
        op = ptr[i][j]
        if op in {"M", "S"}:
            ops.append(op)
            i -= 1
            j -= 1
        elif op == "D":
            ops.append("D")
            i -= 1
        elif op == "I":
            ops.append("I")
            j -= 1
        else:
            raise RuntimeError(f"unexpected traceback op {op!r} at {i},{j}")
    ops.reverse()
    return ops


def summarize_ops(ops):
    matches = ops.count("M")
    substitutions = ops.count("S")
    insertions = ops.count("I")
    deletions = ops.count("D")
    indel_runs = []
    current = None
    length = 0
    for op in ops + ["X"]:
        if op in {"I", "D"} and op == current:
            length += 1
        else:
            if current in {"I", "D"}:
                indel_runs.append((current, length))
            current = op if op in {"I", "D"} else None
            length = 1 if op in {"I", "D"} else 0
    indel_bases = insertions + deletions
    return {
        "alignment_length": len(ops),
        "matches": matches,
        "substitutions": substitutions,
        "insertion_bases_vs_5ltr": insertions,
        "deletion_bases_vs_5ltr": deletions,
        "indel_bases_total": indel_bases,
        "indel_events": len(indel_runs),
        "max_indel_event_bp": max([run_len for _, run_len in indel_runs], default=0),
        "indel_run_summary": ";".join(f"{kind}{run_len}" for kind, run_len in indel_runs) or "none",
        "alignment_identity_pct": matches / len(ops) if ops else 0,
        "indel_only_differences": substitutions == 0 and indel_bases > 0,
        "all_indel_events_single_bp": bool(indel_runs) and all(run_len == 1 for _, run_len in indel_runs),
        "single_1bp_indel_event_only": substitutions == 0 and len(indel_runs) == 1 and indel_bases == 1,
        "perfect_alignment": substitutions == 0 and indel_bases == 0,
    }


def row_group(identity_pct: float) -> str:
    if identity_pct == 1:
        return "100%"
    if identity_pct >= 0.99:
        return ">99_to_<100%"
    return "<99%"


def main() -> None:
    workbook = load_workbook(WORKBOOK, data_only=True)
    sheet = workbook["LTR_identity"]
    headers = [cell.value for cell in sheet[1]]
    idx = {name: i for i, name in enumerate(headers)}
    output_rows = []

    for sheet_row, row in enumerate(sheet.iter_rows(min_row=2, values_only=True), start=2):
        seqid = row[idx["seqid"]]
        five_start = int(row[idx["five_ltr_start"]])
        five_end = int(row[idx["five_ltr_end"]])
        three_start = int(row[idx["three_ltr_start"]])
        three_end = int(row[idx["three_ltr_end"]])

        # Motif-aware corrections already established in the TSD update.
        if sheet_row == 24:
            three_end += 1
        elif sheet_row == 71:
            three_end -= 1

        five_ltr = fetch_sequence(seqid, five_start, five_end)
        three_ltr = fetch_sequence(seqid, three_start, three_end)
        ops = global_edit_alignment(five_ltr, three_ltr)
        summary = summarize_ops(ops)
        identity_pct = float(row[idx["identity_pct"]])
        direct_identical = int(row[idx["identical_bases"]])
        direct_compare_len = int(row[idx["compare_len_bp"]])

        output_rows.append(
            {
                "source_sheet_row": sheet_row,
                "seqid": seqid,
                "five_ltr_coords": f"{five_start}-{five_end}",
                "three_ltr_coords_used": f"{three_start}-{three_end}",
                "direct_identical_bases": direct_identical,
                "direct_compare_len_bp": direct_compare_len,
                "direct_identity_pct": f"{identity_pct:.6f}",
                "identity_group": row_group(identity_pct),
                "five_ltr_len_bp": len(five_ltr),
                "three_ltr_len_bp": len(three_ltr),
                **summary,
            }
        )

        if (sheet_row - 1) % 10 == 0 or sheet_row == sheet.max_row:
            print(f"Aligned {sheet_row - 1}/{sheet.max_row - 1}")

    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    out_wb = Workbook()
    out_ws = out_wb.active
    out_ws.title = "LTR_indel_analysis"
    out_ws.append(list(output_rows[0]))
    for row in output_rows:
        out_ws.append([row[key] for key in output_rows[0]])

    fill = PatternFill("solid", fgColor="1F4E79")
    font = Font(color="FFFFFF", bold=True)
    for cell in out_ws[1]:
        cell.fill = fill
        cell.font = font
        cell.alignment = Alignment(wrap_text=True, vertical="center")
    out_ws.freeze_panes = "A2"
    out_ws.auto_filter.ref = out_ws.dimensions
    for column in out_ws.columns:
        letter = column[0].column_letter
        out_ws.column_dimensions[letter].width = min(max(max(len(str(cell.value or "")) for cell in column) + 2, 10), 24)
    out_wb.save(OUT_XLSX)

    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_XLSX}")


if __name__ == "__main__":
    main()
