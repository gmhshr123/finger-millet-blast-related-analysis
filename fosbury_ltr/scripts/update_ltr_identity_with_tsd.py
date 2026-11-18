from __future__ import annotations

import time
import urllib.parse
import urllib.request
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill


HERE = Path(__file__).resolve().parent
INPUT = HERE / "Fosbury_FullLengthElements_LTR_identity.xlsx"
OUTPUT = HERE / "Fosbury_FullLengthElements_LTR_identity_with_TSD.xlsx"
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
            "tool": "fosbury_tsd_update",
        }
    )
    url = f"{NCBI_EFETCH}?{params}"
    last_error = None
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                text = response.read().decode("utf-8")
            seq = parse_fasta(text)
            if not seq:
                raise RuntimeError(f"empty sequence returned for {seqid}:{start}-{end}")
            time.sleep(0.35)
            return seq
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(attempt)
    raise RuntimeError(f"failed to fetch {seqid}:{start}-{end}: {last_error}")


def slice_seq(window_seq: str, window_start: int, start: int, end: int) -> str | None:
    a = start - window_start
    b = end - window_start + 1
    if a < 0 or b > len(window_seq) or a >= b:
        return None
    return window_seq[a:b]


def count_identical(left: str, right: str) -> tuple[int, int, float]:
    compare_len = min(len(left), len(right))
    identical = sum(1 for a, b in zip(left[:compare_len], right[:compare_len]) if a == b)
    identity = identical / compare_len if compare_len else 0
    return compare_len, identical, identity


def best_tsd_candidate(window_seq: str, window_start: int, five_ltr_start: int, three_ltr_end: int):
    best = None
    for start_shift in range(-5, 6):
        candidate_start = five_ltr_start + start_shift
        left_tsd = slice_seq(window_seq, window_start, candidate_start - 5, candidate_start - 1)
        first_two = slice_seq(window_seq, window_start, candidate_start, candidate_start + 1)
        if left_tsd is None or first_two is None:
            continue
        for end_shift in range(-5, 6):
            candidate_end = three_ltr_end + end_shift
            last_two = slice_seq(window_seq, window_start, candidate_end - 1, candidate_end)
            right_tsd = slice_seq(window_seq, window_start, candidate_end + 1, candidate_end + 5)
            if last_two is None or right_tsd is None:
                continue
            tsd_match = left_tsd == right_tsd
            motif_match = first_two == "TG" and last_two == "CA"
            score = (
                1 if tsd_match else 0,
                1 if motif_match else 0,
                -(abs(start_shift) + abs(end_shift)),
                1 if start_shift == 0 and end_shift == 0 else 0,
            )
            candidate = {
                "left_tsd": left_tsd,
                "right_tsd": right_tsd,
                "tsd_match": tsd_match,
                "motif_match": motif_match,
                "start_shift": start_shift,
                "end_shift": end_shift,
                "five_ltr_start": candidate_start,
                "three_ltr_end": candidate_end,
                "score": score,
            }
            if best is None or candidate["score"] > best["score"]:
                best = candidate
    if best is None:
        raise RuntimeError("No TSD candidate found")
    return best


def main() -> None:
    workbook = load_workbook(INPUT, data_only=True)
    sheet = workbook["LTR_identity"]
    headers = [cell.value for cell in sheet[1]]
    idx = {name: i + 1 for i, name in enumerate(headers)}

    # Remove old TSD columns if the script has been run before.
    existing = {sheet.cell(row=1, column=col).value: col for col in range(1, sheet.max_column + 1)}
    for header in ["left_5bp_tsd", "right_5bp_tsd", "tsd_match"]:
        if header in existing:
            sheet.delete_cols(existing[header])
            existing = {sheet.cell(row=1, column=col).value: col for col in range(1, sheet.max_column + 1)}

    left_col = sheet.max_column + 1
    right_col = sheet.max_column + 2
    match_col = sheet.max_column + 3
    sheet.cell(row=1, column=left_col).value = "left_5bp_tsd"
    sheet.cell(row=1, column=right_col).value = "right_5bp_tsd"
    sheet.cell(row=1, column=match_col).value = "tsd_match"

    changed_identity_rows = []
    shifted_rows = []
    for row_number in range(2, sheet.max_row + 1):
        seqid = sheet.cell(row=row_number, column=idx["seqid"]).value
        five_start = int(sheet.cell(row=row_number, column=idx["five_ltr_start"]).value)
        five_end = int(sheet.cell(row=row_number, column=idx["five_ltr_end"]).value)
        three_start = int(sheet.cell(row=row_number, column=idx["three_ltr_start"]).value)
        three_end = int(sheet.cell(row=row_number, column=idx["three_ltr_end"]).value)

        window_start = max(1, five_start - 12)
        window_end = three_end + 12
        window_seq = fetch_sequence(seqid, window_start, window_end)
        best = best_tsd_candidate(window_seq, window_start, five_start, three_end)

        sheet.cell(row=row_number, column=left_col).value = best["left_tsd"]
        sheet.cell(row=row_number, column=right_col).value = best["right_tsd"]
        sheet.cell(row=row_number, column=match_col).value = "yes" if best["tsd_match"] else "no"

        if best["start_shift"] != 0 or best["end_shift"] != 0:
            shifted_rows.append((row_number, best["start_shift"], best["end_shift"]))
            five_ltr = fetch_sequence(seqid, best["five_ltr_start"], five_end)
            three_ltr = fetch_sequence(seqid, three_start, best["three_ltr_end"])
            compare_len, identical, identity = count_identical(five_ltr, three_ltr)

            old_identical = sheet.cell(row=row_number, column=idx["identical_bases"]).value
            old_identity = sheet.cell(row=row_number, column=idx["identity_pct"]).value
            if old_identical != identical or abs(float(old_identity) - identity) > 1e-12:
                changed_identity_rows.append((row_number, old_identical, identical, old_identity, identity))

            sheet.cell(row=row_number, column=idx["5_ltr_len_bp"]).value = len(five_ltr)
            sheet.cell(row=row_number, column=idx["3_ltr_len_bp"]).value = len(three_ltr)
            sheet.cell(row=row_number, column=idx["compare_len_bp"]).value = compare_len
            sheet.cell(row=row_number, column=idx["identical_bases"]).value = identical
            sheet.cell(row=row_number, column=idx["identity_pct"]).value = identity

        if (row_number - 1) % 10 == 0 or row_number == sheet.max_row:
            print(f"Processed {row_number - 1}/{sheet.max_row - 1}")

    header_fill = PatternFill("solid", fgColor="1F4E79")
    header_font = Font(color="FFFFFF", bold=True)
    for cell in sheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(wrap_text=True, vertical="center")
    for column in [left_col, right_col, match_col]:
        letter = sheet.cell(row=1, column=column).column_letter
        sheet.column_dimensions[letter].width = 14
    for cell in sheet["O"][1:]:
        cell.number_format = "0.00%"

    workbook.save(OUTPUT)
    load_workbook(OUTPUT, data_only=True).close()

    print(f"Wrote {OUTPUT}")
    print(f"Rows with boundary shifts: {shifted_rows}")
    print(f"Rows with changed identity values: {changed_identity_rows}")


if __name__ == "__main__":
    main()
