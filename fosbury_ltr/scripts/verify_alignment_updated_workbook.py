from openpyxl import load_workbook
from pathlib import Path


path = Path(__file__).resolve().parent / "Fosbury_FullLengthElements_LTR_identity_with_TSD_and_alignment.xlsx"
wb = load_workbook(path, data_only=True)
ws = wb["LTR_identity"]
headers = [cell.value for cell in ws[1]]
idx = {name: i for i, name in enumerate(headers)}

low = []
for n, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2):
    if row[idx["identity_pct"]] < 0.99:
        low.append(
            (
                n,
                round(row[idx["identity_pct"]] * 100, 2),
                round(row[idx["aligned_identity_pct"]] * 100, 2),
                row[idx["substitutions"]],
                row[idx["indel_bases"]],
                row[idx["indel_summary"]],
            )
        )
print(low)
print("rows", ws.max_row - 1, "cols", ws.max_column)
