from pathlib import Path
import csv

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill


here = Path(__file__).resolve().parent
source = here / "Fosbury_LTR_lt99_TSD_check.tsv"
output = here / "Fosbury_LTR_lt99_TSD_check.xlsx"

with source.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

workbook = Workbook()
sheet = workbook.active
sheet.title = "TSD_check"

headers = list(rows[0])
sheet.append(headers)
for row in rows:
    sheet.append([row[header] for header in headers])

header_fill = PatternFill("solid", fgColor="1F4E79")
header_font = Font(color="FFFFFF", bold=True)
for cell in sheet[1]:
    cell.fill = header_fill
    cell.font = header_font
    cell.alignment = Alignment(wrap_text=True, vertical="center")

sheet.freeze_panes = "A2"
sheet.auto_filter.ref = sheet.dimensions
for column, width in {
    "A": 16,
    "B": 14,
    "C": 10,
    "D": 18,
    "E": 14,
    "F": 14,
    "G": 18,
    "H": 18,
    "I": 18,
    "J": 16,
    "K": 14,
    "L": 14,
    "M": 12,
}.items():
    sheet.column_dimensions[column].width = width

workbook.save(output)
print(output)
