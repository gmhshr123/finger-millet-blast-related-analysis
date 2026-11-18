from openpyxl import load_workbook
from pathlib import Path


path = Path(__file__).resolve().parent / "Fosbury_FullLengthElements_LTR_identity_with_TSD.xlsx"
wb = load_workbook(path, data_only=True)
ws = wb["LTR_identity"]
headers = [cell.value for cell in ws[1]]
idx = {name: i for i, name in enumerate(headers)}

for n, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2):
    if row[idx["tsd_match"]] != "yes":
        print(
            n,
            {
                key: row[idx[key]]
                for key in [
                    "seqid",
                    "element_start",
                    "element_end",
                    "five_ltr_start",
                    "three_ltr_end",
                    "identical_bases",
                    "identity_pct",
                    "left_5bp_tsd",
                    "right_5bp_tsd",
                    "tsd_match",
                ]
            },
        )
