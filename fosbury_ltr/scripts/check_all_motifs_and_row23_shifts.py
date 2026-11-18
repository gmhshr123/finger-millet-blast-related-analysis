from __future__ import annotations

import time
import urllib.parse
import urllib.request
from pathlib import Path

from openpyxl import load_workbook


HERE = Path(__file__).resolve().parent
WORKBOOK = HERE / "Fosbury_FullLengthElements_LTR_identity.xlsx"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_fasta(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">")).upper()


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
            "tool": "fosbury_motif_check",
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
                raise RuntimeError("empty sequence")
            time.sleep(0.35)
            return seq
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(attempt)
    raise RuntimeError(f"failed to fetch {seqid}:{start}-{end}: {last_error}")


def subseq(window: str, window_start: int, start: int, end: int) -> str:
    a = start - window_start
    b = end - window_start + 1
    if a < 0 or b > len(window) or a >= b:
        return ""
    return window[a:b]


def matching_positions(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x == y)


def main() -> None:
    workbook = load_workbook(WORKBOOK, data_only=True)
    sheet = workbook["LTR_identity"]
    headers = [cell.value for cell in sheet[1]]
    idx = {name: i for i, name in enumerate(headers)}

    motif_failures = []
    current_tsd_failures = []
    rows = list(sheet.iter_rows(min_row=2, values_only=True))
    for sheet_row, row in enumerate(rows, start=2):
        seqid = row[idx["seqid"]]
        five_start = int(row[idx["five_ltr_start"]])
        three_end = int(row[idx["three_ltr_end"]])
        window_start = max(1, five_start - 8)
        window_end = three_end + 8
        window = fetch_sequence(seqid, window_start, window_end)

        first_two = subseq(window, window_start, five_start, five_start + 1)
        last_two = subseq(window, window_start, three_end - 1, three_end)
        left_tsd = subseq(window, window_start, five_start - 5, five_start - 1)
        right_tsd = subseq(window, window_start, three_end + 1, three_end + 5)

        if first_two != "TG" or last_two != "CA":
            motif_failures.append((sheet_row, seqid, first_two, last_two))
        if left_tsd != right_tsd:
            current_tsd_failures.append((sheet_row, seqid, left_tsd, right_tsd, matching_positions(left_tsd, right_tsd)))

        if (sheet_row - 1) % 20 == 0 or sheet_row == len(rows) + 1:
            print(f"Checked {sheet_row - 1}/{len(rows)}")

    row23 = rows[21]  # Excel row 23
    seqid = row23[idx["seqid"]]
    five_start = int(row23[idx["five_ltr_start"]])
    three_end = int(row23[idx["three_ltr_end"]])
    window_start = max(1, five_start - 15)
    window_end = three_end + 15
    window = fetch_sequence(seqid, window_start, window_end)

    row23_candidates = []
    for start_shift in range(-5, 6):
        candidate_start = five_start + start_shift
        left_tsd = subseq(window, window_start, candidate_start - 5, candidate_start - 1)
        first_two = subseq(window, window_start, candidate_start, candidate_start + 1)
        for end_shift in range(-5, 6):
            candidate_end = three_end + end_shift
            right_tsd = subseq(window, window_start, candidate_end + 1, candidate_end + 5)
            last_two = subseq(window, window_start, candidate_end - 1, candidate_end)
            row23_candidates.append(
                {
                    "start_shift": start_shift,
                    "end_shift": end_shift,
                    "left_tsd": left_tsd,
                    "right_tsd": right_tsd,
                    "matches": matching_positions(left_tsd, right_tsd),
                    "tsd_match": left_tsd == right_tsd,
                    "first_two": first_two,
                    "last_two": last_two,
                    "motif_match": first_two == "TG" and last_two == "CA",
                }
            )

    exact = [c for c in row23_candidates if c["tsd_match"]]
    best = sorted(
        row23_candidates,
        key=lambda c: (c["matches"], c["motif_match"], -(abs(c["start_shift"]) + abs(c["end_shift"]))),
        reverse=True,
    )[:10]

    print("motif_failures", motif_failures)
    print("current_tsd_failures", current_tsd_failures)
    print("row23_exact_tsd_candidates", exact)
    print("row23_best_candidates", best)


if __name__ == "__main__":
    main()
