from __future__ import annotations

import csv
import textwrap
import time
import urllib.parse
import urllib.request
from pathlib import Path

from openpyxl import load_workbook


HERE = Path(__file__).resolve().parent
WORKBOOK = HERE / "Fosbury_FullLengthElements_LTR_identity.xlsx"
FULL_FASTA = HERE / "Fosbury_LTR_lt99_full_length_sequences.fasta"
FLANK_FASTA = HERE / "Fosbury_LTR_lt99_full_length_with_5bp_flanks.fasta"
TSD_TSV = HERE / "Fosbury_LTR_lt99_TSD_check.tsv"
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
            "tool": "fosbury_tsd_check",
        }
    )
    url = f"{NCBI_EFETCH}?{params}"
    last_error = None
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                text = response.read().decode("utf-8")
            sequence = parse_fasta(text)
            if not sequence:
                raise RuntimeError(f"empty sequence returned for {seqid}:{start}-{end}")
            time.sleep(0.35)
            return sequence
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(attempt)
    raise RuntimeError(f"failed to fetch {seqid}:{start}-{end}: {last_error}")


def wrap_fasta(sequence: str) -> str:
    return "\n".join(textwrap.wrap(sequence, width=70))


def main() -> None:
    workbook = load_workbook(WORKBOOK, data_only=True)
    sheet = workbook["LTR_identity"]
    headers = [cell.value for cell in sheet[1]]
    idx = {name: i for i, name in enumerate(headers)}

    selected = []
    for excel_row, row in enumerate(sheet.iter_rows(min_row=2, values_only=True), start=2):
        identity_pct = row[idx["identity_pct"]]
        if identity_pct is not None and float(identity_pct) < 0.99:
            selected.append((excel_row, row))

    full_fasta_records = []
    flank_fasta_records = []
    tsd_rows = []

    for record_number, (excel_row, row) in enumerate(selected, start=1):
        genome = row[idx["genome"]]
        seqid = row[idx["seqid"]]
        element_start = int(row[idx["element_start"]])
        element_end = int(row[idx["element_end"]])
        identical_bases = int(row[idx["identical_bases"]])
        compare_len = int(row[idx["compare_len_bp"]])
        identity_pct = float(row[idx["identity_pct"]])

        element_sequence = fetch_sequence(seqid, element_start, element_end)
        flank_start = max(1, element_start - 5)
        flank_end = element_end + 5
        flanked_sequence = fetch_sequence(seqid, flank_start, flank_end)

        left_tsd = flanked_sequence[: element_start - flank_start]
        full_from_flank = flanked_sequence[element_start - flank_start : element_start - flank_start + len(element_sequence)]
        right_tsd = flanked_sequence[element_start - flank_start + len(element_sequence) :]
        tsd_match = left_tsd == right_tsd and len(left_tsd) == 5 and len(right_tsd) == 5

        if full_from_flank != element_sequence:
            raise RuntimeError(f"flanked sequence did not contain expected element for row {excel_row}")

        record_id = (
            f"Fosbury_lt99_{record_number:02d}|sheet_row={excel_row}|"
            f"{genome}|{seqid}:{element_start}-{element_end}|"
            f"LTR_identity={identity_pct:.2%}|identical={identical_bases}/{compare_len}|"
            f"TSD_left={left_tsd}|TSD_right={right_tsd}|TSD_match={tsd_match}"
        )
        full_fasta_records.append(f">{record_id}\n{wrap_fasta(element_sequence)}\n")
        flank_fasta_records.append(
            f">{record_id}|with_5bp_flanks={flank_start}-{flank_end}\n{wrap_fasta(flanked_sequence)}\n"
        )
        tsd_rows.append(
            {
                "record_id": f"Fosbury_lt99_{record_number:02d}",
                "source_sheet_row": excel_row,
                "genome": genome,
                "seqid": seqid,
                "element_start": element_start,
                "element_end": element_end,
                "element_length_bp": len(element_sequence),
                "ltr_identical_bases": identical_bases,
                "ltr_compare_len_bp": compare_len,
                "ltr_identity_pct": f"{identity_pct:.6f}",
                "left_5bp_flank": left_tsd,
                "right_5bp_flank": right_tsd,
                "has_5bp_tsd": "yes" if tsd_match else "no",
            }
        )
        print(f"Fetched {record_number}/{len(selected)}: row {excel_row} {seqid}:{element_start}-{element_end}")

    FULL_FASTA.write_text("\n".join(full_fasta_records), encoding="utf-8")
    FLANK_FASTA.write_text("\n".join(flank_fasta_records), encoding="utf-8")
    with TSD_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(tsd_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(tsd_rows)

    print(f"Selected elements: {len(selected)}")
    print(f"Wrote {FULL_FASTA}")
    print(f"Wrote {FLANK_FASTA}")
    print(f"Wrote {TSD_TSV}")


if __name__ == "__main__":
    main()
