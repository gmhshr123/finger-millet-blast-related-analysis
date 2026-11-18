from pathlib import Path
import csv
import re

from openpyxl import load_workbook, Workbook
from openpyxl.styles import Alignment, Font, PatternFill


here = Path(__file__).resolve().parent
workbook_path = here / "Fosbury_FullLengthElements_LTR_identity.xlsx"
flank_fasta = here / "Fosbury_LTR_lt99_full_length_with_5bp_flanks.fasta"
out_tsv = here / "Fosbury_LTR_lt99_TSD_recheck_with_LTR_motifs.tsv"
out_xlsx = here / "Fosbury_LTR_lt99_TSD_recheck_with_LTR_motifs.xlsx"
corrected_fasta = here / "Fosbury_LTR_lt99_corrected_element_sequences_TG_to_CA.fasta"
corrected_with_tsd_fasta = here / "Fosbury_LTR_lt99_corrected_element_sequences_with_5bp_TSDs.fasta"


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def get_subseq(seq, window_start, coord_start, coord_end):
    a = coord_start - window_start
    b = coord_end - window_start + 1
    if a < 0 or b > len(seq) or a >= b:
        return None
    return seq[a:b]


def wrap_fasta(sequence):
    return "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70))


def main():
    wb = load_workbook(workbook_path, data_only=True)
    ws = wb["LTR_identity"]
    headers = [cell.value for cell in ws[1]]
    idx = {name: i for i, name in enumerate(headers)}

    rows = []
    for excel_row, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2):
        if row[idx["identity_pct"]] is not None and float(row[idx["identity_pct"]]) < 0.99:
            rows.append((excel_row, row))

    fasta_records = read_fasta(flank_fasta)
    by_record_id = {}
    for header, seq in fasta_records.items():
        record_id = header.split("|", 1)[0]
        match = re.search(r"with_5bp_flanks=(\d+)-(\d+)", header)
        if not match:
            raise ValueError(f"Could not parse flanked coordinates from {header}")
        by_record_id[record_id] = (int(match.group(1)), int(match.group(2)), seq)

    output_rows = []
    corrected_records = []
    corrected_with_tsd_records = []
    for n, (excel_row, row) in enumerate(rows, start=1):
        record_id = f"Fosbury_lt99_{n:02d}"
        window_start, window_end, seq = by_record_id[record_id]

        element_start_with_tsd = int(row[idx["element_start"]])
        element_end_with_tsd = int(row[idx["element_end"]])
        annotated_ltr_start = int(row[idx["five_ltr_start"]])
        annotated_ltr_end = int(row[idx["three_ltr_end"]])

        best = None
        candidates = []
        for start_shift in range(-5, 6):
            candidate_ltr_start = annotated_ltr_start + start_shift
            left_tsd = get_subseq(seq, window_start, candidate_ltr_start - 5, candidate_ltr_start - 1)
            first_two = get_subseq(seq, window_start, candidate_ltr_start, candidate_ltr_start + 1)
            if left_tsd is None or first_two is None:
                continue
            for end_shift in range(-5, 6):
                candidate_ltr_end = annotated_ltr_end + end_shift
                last_two = get_subseq(seq, window_start, candidate_ltr_end - 1, candidate_ltr_end)
                right_tsd = get_subseq(seq, window_start, candidate_ltr_end + 1, candidate_ltr_end + 5)
                if last_two is None or right_tsd is None:
                    continue
                has_tsd = left_tsd == right_tsd
                has_motif = first_two == "TG" and last_two == "CA"
                score = (1 if has_tsd else 0, 1 if has_motif else 0, -(abs(start_shift) + abs(end_shift)))
                candidate = {
                    "candidate_ltr_start": candidate_ltr_start,
                    "candidate_ltr_end": candidate_ltr_end,
                    "start_shift": start_shift,
                    "end_shift": end_shift,
                    "left_5bp_tsd": left_tsd,
                    "right_5bp_tsd": right_tsd,
                    "has_5bp_tsd": "yes" if has_tsd else "no",
                    "first_two_bases": first_two,
                    "last_two_bases": last_two,
                    "has_TG_CA_motif": "yes" if has_motif else "no",
                    "score": score,
                }
                candidates.append(candidate)
                if best is None or candidate["score"] > best["score"]:
                    best = candidate

        if best is None:
            raise RuntimeError(f"No candidates found for {record_id}")

        output_rows.append(
            {
                "record_id": record_id,
                "source_sheet_row": excel_row,
                "seqid": row[idx["seqid"]],
                "element_start_with_tsd": element_start_with_tsd,
                "element_end_with_tsd": element_end_with_tsd,
                "annotated_5_ltr_start": annotated_ltr_start,
                "annotated_3_ltr_end": annotated_ltr_end,
                "best_5_ltr_start": best["candidate_ltr_start"],
                "best_3_ltr_end": best["candidate_ltr_end"],
                "start_shift_bp": best["start_shift"],
                "end_shift_bp": best["end_shift"],
                "left_5bp_tsd": best["left_5bp_tsd"],
                "right_5bp_tsd": best["right_5bp_tsd"],
                "has_5bp_tsd": best["has_5bp_tsd"],
                "first_two_bases": best["first_two_bases"],
                "last_two_bases": best["last_two_bases"],
                "has_TG_CA_motif": best["has_TG_CA_motif"],
                "ltr_identity_pct": f"{float(row[idx['identity_pct']]):.6f}",
                "ltr_identical_bases": int(row[idx["identical_bases"]]),
                "ltr_compare_len_bp": int(row[idx["compare_len_bp"]]),
            }
        )
        element_sequence = get_subseq(
            seq,
            window_start,
            best["candidate_ltr_start"],
            best["candidate_ltr_end"],
        )
        element_with_tsd_sequence = get_subseq(
            seq,
            window_start,
            best["candidate_ltr_start"] - 5,
            best["candidate_ltr_end"] + 5,
        )
        fasta_header = (
            f"{record_id}|sheet_row={excel_row}|{row[idx['seqid']]}:"
            f"{best['candidate_ltr_start']}-{best['candidate_ltr_end']}|"
            f"LTR_identity={float(row[idx['identity_pct']]):.2%}|"
            f"identical={int(row[idx['identical_bases']])}/{int(row[idx['compare_len_bp']])}|"
            f"TSD={best['left_5bp_tsd']}|start={best['first_two_bases']}|end={best['last_two_bases']}"
        )
        corrected_records.append(f">{fasta_header}\n{wrap_fasta(element_sequence)}\n")
        corrected_with_tsd_records.append(
            f">{fasta_header}|with_TSDs={best['candidate_ltr_start'] - 5}-{best['candidate_ltr_end'] + 5}\n"
            f"{wrap_fasta(element_with_tsd_sequence)}\n"
        )

    with out_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    corrected_fasta.write_text("\n".join(corrected_records), encoding="utf-8")
    corrected_with_tsd_fasta.write_text("\n".join(corrected_with_tsd_records), encoding="utf-8")

    out_wb = Workbook()
    out_ws = out_wb.active
    out_ws.title = "TSD_recheck"
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
    out_wb.save(out_xlsx)

    print(f"Wrote {out_tsv}")
    print(f"Wrote {out_xlsx}")
    print(f"Wrote {corrected_fasta}")
    print(f"Wrote {corrected_with_tsd_fasta}")
    print(f"Rows with best TSD and TG...CA motif: {sum(r['has_5bp_tsd'] == 'yes' and r['has_TG_CA_motif'] == 'yes' for r in output_rows)}/{len(output_rows)}")
    for row in output_rows:
        print(row["record_id"], row["left_5bp_tsd"], row["right_5bp_tsd"], row["has_5bp_tsd"], row["first_two_bases"], row["last_two_bases"], row["has_TG_CA_motif"], row["start_shift_bp"], row["end_shift_bp"])


if __name__ == "__main__":
    main()
