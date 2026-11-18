from __future__ import annotations

import csv
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill


HERE = Path(__file__).resolve().parent
INPUT = HERE / "Fosbury_FullLengthElements_LTR_identity_with_TSD.xlsx"
ANALYSIS = HERE / "Fosbury_LTR_pairwise_indel_analysis.tsv"
OUTPUT = HERE / "Fosbury_FullLengthElements_LTR_identity_with_TSD_and_alignment.xlsx"


def main() -> None:
    workbook = load_workbook(INPUT, data_only=True)
    sheet = workbook["LTR_identity"]

    with ANALYSIS.open(encoding="utf-8", newline="") as handle:
        analysis_rows = {
            int(row["source_sheet_row"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }

    # Keep this compact: original direct identity remains; these columns show
    # the indel-aware alignment result and what kinds of differences were found.
    new_columns = [
        "aligned_identity_pct",
        "substitutions",
        "indel_bases",
        "indel_events",
        "indel_summary",
    ]

    # Remove old copies if rerun.
    for header in new_columns:
        current_headers = {
            sheet.cell(row=1, column=col).value: col
            for col in range(1, sheet.max_column + 1)
        }
        if header in current_headers:
            sheet.delete_cols(current_headers[header])

    first_new_col = sheet.max_column + 1
    for offset, header in enumerate(new_columns):
        sheet.cell(row=1, column=first_new_col + offset).value = header

    for sheet_row in range(2, sheet.max_row + 1):
        analysis = analysis_rows[sheet_row]
        values = [
            float(analysis["alignment_identity_pct"]),
            int(analysis["substitutions"]),
            int(analysis["indel_bases_total"]),
            int(analysis["indel_events"]),
            analysis["indel_run_summary"],
        ]
        for offset, value in enumerate(values):
            sheet.cell(row=sheet_row, column=first_new_col + offset).value = value

    header_fill = PatternFill("solid", fgColor="1F4E79")
    header_font = Font(color="FFFFFF", bold=True)
    for cell in sheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(wrap_text=True, vertical="center")

    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions

    for col in range(first_new_col, first_new_col + len(new_columns)):
        letter = sheet.cell(row=1, column=col).column_letter
        sheet.column_dimensions[letter].width = 16 if col != first_new_col + 4 else 22
    for cell in sheet[sheet.cell(row=1, column=first_new_col).column_letter][1:]:
        cell.number_format = "0.00%"

    workbook.save(OUTPUT)
    load_workbook(OUTPUT, data_only=True).close()
    print(OUTPUT)


if __name__ == "__main__":
    main()
