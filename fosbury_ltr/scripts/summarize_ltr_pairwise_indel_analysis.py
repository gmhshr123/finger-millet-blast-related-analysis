from __future__ import annotations

import collections
import csv
from pathlib import Path


path = Path(__file__).resolve().parent / "Fosbury_LTR_pairwise_indel_analysis.tsv"
with path.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

groups = collections.defaultdict(list)
for row in rows:
    groups[row["identity_group"]].append(row)

for group in ["<99%", ">99_to_<100%", "100%"]:
    records = groups[group]
    print(group, len(records))
    print("  indel_only:", sum(row["indel_only_differences"] == "True" for row in records))
    print("  all_indel_events_single_bp:", sum(row["all_indel_events_single_bp"] == "True" for row in records))
    print("  single_1bp_indel_only:", sum(row["single_1bp_indel_event_only"] == "True" for row in records))
    print("  has_substitutions:", sum(int(row["substitutions"]) > 0 for row in records))
    print("  perfect_alignment:", sum(row["perfect_alignment"] == "True" for row in records))
    nonperfect = [
        (
            row["source_sheet_row"],
            row["direct_identity_pct"],
            row["substitutions"],
            row["indel_run_summary"],
            row["indel_bases_total"],
            row["alignment_identity_pct"],
        )
        for row in records
        if row["perfect_alignment"] != "True"
    ]
    print("  nonperfect_rows:", nonperfect)
