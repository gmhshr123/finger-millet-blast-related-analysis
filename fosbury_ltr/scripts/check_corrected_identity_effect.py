from pathlib import Path
import csv

from openpyxl import load_workbook


base = Path(__file__).resolve().parent
wb = load_workbook(base / "Fosbury_FullLengthElements_LTR_identity.xlsx", data_only=True)
ws = wb["LTR_identity"]
headers = [cell.value for cell in ws[1]]
idx = {name: i for i, name in enumerate(headers)}
rows = {i: row for i, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2)}

with (base / "Fosbury_LTR_lt99_TSD_recheck_with_LTR_motifs.tsv").open(encoding="utf-8") as handle:
    for record in csv.DictReader(handle, delimiter="\t"):
        sheet_row = int(record["source_sheet_row"])
        row = rows[sheet_row]
        print(
            sheet_row,
            row[idx["identical_bases"]],
            row[idx["compare_len_bp"]],
            f"{float(row[idx['identity_pct']]) * 100:.2f}%",
            "start_shift",
            record["start_shift_bp"],
            "end_shift",
            record["end_shift_bp"],
        )
