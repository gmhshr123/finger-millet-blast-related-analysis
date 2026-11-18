#!/usr/bin/env python3
"""Supplementary files for the all-base deletion analysis (run from manuscripts_all/ after script 26).

Writes deletion_allbase/supplement/:
  Data_S10_updated_all-base_deletions.xlsx   all tables, one sheet each, with a README sheet
  Figure_S4_updated_100kb_tracks.png/.pdf     33 strains, all-base and SNP-based deletions with SNP density and repeats
  Figure_S4_all_226_strains.png/.pdf          all 226 strains
  Supplement_legends.md                       legends for the figures and the workbook

All values are copied from the result tables; nothing is recomputed here except simple per-strain sums.
Coordinates follow the depth tracks: start is 0-based and end is exclusive (BED convention).
"""
import json, os, shutil
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

R = os.environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")
R226 = "deletion_allbase/results/"
I = "deletion_allbase/inputs/"
OUT = "deletion_allbase/supplement/"
os.makedirs(OUT, exist_ok=True)
S = json.load(open(R + "summary.json"))
cohort = [l.strip() for l in open(I + "cohort33.txt")]

LABEL = {"simple deletion": "simple deletion", "other sequence at edge": "other sequence at the edge",
         "DNA still present": "DNA still present", "depth only": "depth only"}
fam = lambda s: s.str.replace("MAGGY_Fosbury", "MAGGY/Fosbury", regex=False)
yn = lambda s: s.map({True: "yes", False: "no"})

# ------------------------------------------------------------------ events
e = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t")
events = pd.DataFrame({
    "event_id": e.event_id, "chromosome": e.chrom, "start": e.start, "end": e.end,
    "breakpoint_zone_left_start": e.outer_start, "breakpoint_zone_right_end": e.outer_end,
    "deleted_bases": e.deleted_bases, "size_class": e.size_class,
    "n_strains": e.n_strains, "strains": e.strains,
    "read_evidence": e.read_class.map(LABEL),
    "confirmed_deletion": yn(e.read_class.isin(["simple deletion", "other sequence at edge"])),
    "poorly_mappable": yn(e.poorly_mappable.astype(bool)),
    "TE_at_breakpoint": e.te_class, "TE_evidence": e.te_evidence, "TE_family": fam(e.te_family.fillna("")),
    "TE_family_left": fam(e.L_family.fillna("")), "TE_family_right": fam(e.R_family.fillna("")),
    "SNPs_inside": e.snps_inside, "detected_by_SNP_method": yn(e.found_by_snp_method.astype(bool)),
    "n_genes": e.n_genes, "genes": e.genes.fillna(""), "effector_genes": e.effector_genes.fillna("")})

# ------------------------------------------------------------------ strain-level regions
c = pd.read_csv(R + "allbase_calls_strain_level_read_evidence.tsv", sep="\t")
regions = pd.DataFrame({
    "strain": c.strain, "event_id": c.event_id, "chromosome": c.chrom, "start": c.start, "end": c.end,
    "breakpoint_zone_left_start": c.outer_start, "breakpoint_zone_right_end": c.outer_end,
    "deleted_bases": c.deleted_bases, "size_class": c.size_class,
    "strain_median_depth": c.strain_median_depth,
    "read_evidence": c.read_class.map(LABEL),
    "junction_Delly_or_split_reads": yn(c.junction.astype(bool)), "discordant_pairs": yn(c.pairs_deletion.astype(bool)),
    "concordant_pairs": yn(c.pairs_present.astype(bool)),
    "soft_clipped_reads_match_TE": yn(c.edge_te_reads.astype(bool)),
    "TE_family_left_reads": fam(c.L_te_family_reads.fillna("")), "TE_family_right_reads": fam(c.R_te_family_reads.fillna("")),
    "mate_maps_elsewhere": yn(c.edge_mate_far.astype(bool)), "mate_unmapped": yn(c.edge_mate_unmapped.astype(bool)),
    "poorly_mappable": yn(c.poorly_mappable.astype(bool)),
    "SNPs_inside": c.snps_inside, "detected_by_SNP_method": yn(c.found_by_snp_method.astype(bool))})

# ------------------------------------------------------------------ per strain
strains = (c.groupby("strain").agg(median_depth=("strain_median_depth", "first"), median_insert_size=("lib_median", "median"),
                                   deleted_regions=("start", "size"), deleted_Mb=("deleted_bases", lambda v: round(v.sum() / 1e6, 3)),
                                   confirmed_regions=("read_class", lambda v: int(v.isin(["simple deletion", "other sequence at edge"]).sum())))
           .reindex(cohort).reset_index())

# ------------------------------------------------------------------ Data S10 rows checked
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")
x = pd.read_csv(R + "DataS10_crosswalk_read_evidence.tsv", sep="\t")
m = s10.merge(x, on="del_id", suffixes=("", "_x"))
DISP = {"supported": "supported", "supported (partly deleted in some strains)": "supported (partly deleted in some strains)",
        "not assessable (repeat / uncallable)": "cannot be tested (repeat or sequence that is not uniquely mappable)", "not supported": "not supported"}
s10_out = pd.DataFrame({
    "DataS10_deletion_id": m.del_id, "chromosome": m.chrom, "DataS10_start": m.start, "DataS10_end": m.end,
    "DataS10_strains": m.strains, "DataS10_TE": m.te, "DataS10_SNP_count": m.snp_count,
    "DataS10_min_size": m.min_size, "DataS10_max_size": m.max_size,
    "strains_listed": m.n_listed, "strains_deleted": m.n_deleted, "strains_partly_deleted": m.n_partly,
    "strains_not_deleted": m.n_not_deleted, "strains_not_testable": m.n_not_assessable,
    "read_depth_result": m.disposition.map(DISP),
    "all-base_events": m.allbase_events.fillna(""), "all-base_read_evidence": m.read_class.replace({"other sequence at edge": "other sequence at the edge"}),
    "all-base_TE_at_breakpoint": m.allbase_te_class.fillna(""), "all-base_TE_family": fam(m.allbase_te_family.fillna("")),
    "same_TE_family_as_DataS10": m.same_family_as_s10.map({True: "yes", False: "no"}).fillna("")})

# ------------------------------------------------------------------ SNP-based calls checked (runs of 5 or more missing SNPs)
r = pd.read_csv(R + "snp_candidates_checked.tsv", sep="\t")
r = r[(r.snp_count >= 5) & r.strain.isin(cohort)]
assert len(r) == S["snp_candidates_ge5"]["n"]
snp_out = pd.DataFrame({
    "strain": r.strain, "chromosome": r.chrom, "start_first_missing_SNP": r.start, "end_last_missing_SNP": r.end,
    "missing_SNPs": r.snp_count, "min_size": r.min_size, "max_size": r.max_size,
    "uniquely_mappable_bases": r.scorable_bases, "deleted_bases": r.deleted_bases,
    "fraction_deleted": r.frac_deleted.round(3), "depth_uniquely_mapped_relative_to_median": r.rel_depth_unique_scorable,
    "depth_all_reads_relative_to_median": r.rel_depth_allreads,
    "read_depth_result": r.allbase.replace({"not assessable (repeat / uncallable)": "cannot be tested (repeat or sequence that is not uniquely mappable)"})})

# ------------------------------------------------------------------ Data S9 checked
s9 = pd.read_csv(R + "DataS9_crosswalk.tsv", sep="\t")
s9_out = pd.DataFrame({
    "DataS9_example": s9.example, "chromosome": s9.chrom, "DataS9_proximal_breakpoint": s9.s9_proximal, "DataS9_distal_breakpoint": s9.s9_distal,
    "strain": s9.strain, "DataS9_TE": s9.s9_te, "result": s9.result,
    "all-base_events": s9.events.fillna(""), "all-base_left_breakpoint_zone": s9.allbase_start_zone.fillna(""),
    "all-base_right_breakpoint_zone": s9.allbase_end_zone.fillna(""),
    "proximal_breakpoint_agrees": s9.proximal_edge_agrees.map({True: "yes", False: "no"}).fillna(""),
    "distal_breakpoint_agrees": s9.distal_edge_agrees.map({True: "yes", False: "no"}).fillna(""),
    "all-base_TE_family": fam(s9.allbase_te_families.fillna("")),
    "same_TE_family_as_DataS9": s9.same_family_as_s9.map({True: "yes", False: "no"}).fillna(""),
    "all-base_read_evidence": s9.read_class.fillna("").str.replace("other sequence at edge", "other sequence at the edge", regex=False)})

# ------------------------------------------------------------------ figure source data
w = pd.read_csv(R + "windows100kb_tracks.tsv", sep="\t")
win = pd.DataFrame({
    "chromosome": w.chrom, "window_start": w.start, "window_end": w.end, "uniquely_mappable_bases": w.scorable_bases,
    "all-base_fraction_of_uniquely_mappable_bases_deleted_mean_of_33_strains": w.allbase_deleted_frac.round(5),
    "SNP-based_fraction_of_SNP_calls_in_deletions": w.snpbased_deleted_frac.round(5),
    "SNPs": w.n_snps, "SNPs_per_uniquely_mappable_kb": w.snps_per_scorable_kb.round(3), "repeat_fraction": w.repeat_frac.round(4),
    "genes": w.n_genes, "effector_genes": w.n_effectors})
p226 = pd.read_csv(R226 + "figure_226_strains.per_strain.tsv", sep="\t")
p226["population"] = p226.population.replace({"neither (<90%)": "neither (less than 90% membership)"})
p226 = p226.rename(columns={"deleted_Mb": "deleted_Mb_read_depth_only"}).sort_values(["population", "strain"])

# ------------------------------------------------------------------ README
cls = e.read_class.value_counts()
readme = [
    ("All-base deletion analysis: supplementary tables", ""),
    ("", ""),
    ("Scope", "Deletions relative to the E2 reference in the 33 Ethiopian strains of Data S10, called from read depth at every base (all-base method) and tested with split reads and paired-end reads. The original Data S9 and Data S10 are not changed; their rows are listed here with the result of the all-base test."),
    ("Counts", f"{S['regions_strain_level']:,} strain-level regions; {S['events']} events; {int(cls.get('simple deletion', 0) + cls.get('other sequence at edge', 0))} events confirmed as deletions by split reads or paired-end reads; "
               f"{int(cls.get('DNA still present', 0))} present in the strain (not deletions); {int(cls.get('depth only', 0))} supported by read depth only."),
    ("Coordinates", "E2 assembly. start is 0-based and end is exclusive (BED convention). start and end are the first and last deleted uniquely mappable base; the true breakpoints lie between breakpoint_zone_left_start and start, and between end and breakpoint_zone_right_end."),
    ("", ""),
    ("Sheet", "Content"),
    ("Events", "One row per deletion event. Regions of different strains with at least 80% reciprocal overlap form one event."),
    ("Regions_per_strain", "One row per deleted region in one strain, with the read evidence found in that strain."),
    ("Strains", "Read depth, insert size and number of deleted regions for each of the 33 strains."),
    ("DataS10_checked", "The 124 SNP-based events of Data S10 with the result of the read-depth test and the overlapping all-base events."),
    ("SNP_calls_checked", f"The {S['snp_candidates_ge5']['n']:,} strain-level SNP-based calls (runs of five or more consecutive SNPs with a missing genotype call) with the result of the read-depth test."),
    ("DataS9_checked", "The manually inspected deletions of Data S9 that lie in the 33 strains, compared with the all-base calls."),
    ("Windows_100kb", "Values per 100-kb window; source data of the updated Figure S4 (33 strains)."),
    ("Strains_226", "Deleted sequence per strain for all 226 strains; source data of the 226-strain figure. Read depth only: not merged and not tested with split reads or paired-end reads."),
    ("Parameters", "All thresholds used."),
    ("", ""),
    ("Column", "Meaning"),
    ("deleted_bases", "Number of uniquely mappable bases in the region with a depth of uniquely mapped reads at or below the deletion cut-off. Uniquely mappable bases exclude annotated repeats, assembly gaps and bases where E2's own reads have less than 25% of the E2 median depth."),
    ("read_evidence", "simple deletion: discordant pairs or split reads connect the two flanks. other sequence at the edge: reads at a breakpoint continue into a TE or other sequence (soft-clipped part matches a TE, or mate maps elsewhere or is unmapped). DNA still present: concordant pairs map on opposite sides of the region; not a deletion. depth only: none of these."),
    ("confirmed_deletion", "yes for 'simple deletion' and 'other sequence at the edge'."),
    ("poorly_mappable", "yes when the depth of uniquely mapped E2 reads in the region is below 60% of the E2 median, or below 50% of the median in one of the four strains nearly identical to E2 (E5, E6, E19, E36)."),
    ("TE_at_breakpoint, TE_evidence", "A TE annotated in E2 within 20 bp of a breakpoint zone (reference annotation), or at least three soft-clipped reads at a breakpoint that match the TE library (reads)."),
    ("TE_family", "Name from the combined EDTA TE library (275 entries from B71, E2, MZ5-1-6 and P131). An entry has a published name when it matches that TE at 80% identity or more; other entries keep their EDTA label. MAGGY and Fosbury are one library entry."),
    ("detected_by_SNP_method", "yes when a run of five or more missing SNPs in the same strain overlaps the region (for an event: in at least one strain)."),
    ("genes, effector_genes", "E2 genes with at least half of their length inside the event; effectors as in Data S2."),
    ("read_depth_result", "confirmed deleted / supported: at least 80% of the uniquely mappable bases of the interval are deleted. partly deleted: 20% to 80%, or an all-base deletion of the same strain overlaps the interval. not deleted / not supported: below 20% and no overlapping all-base deletion. cannot be tested: fewer than 30 uniquely mappable bases and no overlapping all-base deletion."),
]
params = pd.DataFrame([
    ("Mapping quality for a uniquely mapped read", "20 or higher"),
    ("Depth of uniquely mapped E2 reads for a base to be uniquely mappable", "at least 25% of the E2 median"),
    ("Depth of uniquely mapped reads for a base to be called deleted", f"at most {S['depth_threshold_frac'] * 100:g}% of the strain's median"),
    ("Minimum deleted bases per call", "100"),
    ("Consecutive uniquely mappable bases that are present, ending a call", "25"),
    ("Distance below which adjacent calls of one strain are merged", f"{S['merge_gap_bp']} bp"),
    ("Reciprocal overlap for regions of different strains to form one event", "at least 80%"),
    ("Supporting reads for a junction, discordant pairs, or a TE at a breakpoint", "at least 3"),
    ("Reads at a breakpoint whose mate maps elsewhere or is unmapped", "at least 3, and at least 10% of the reads there"),
    ("Length of the soft-clipped part of a read", "at least 25 nt"),
    ("BLAST e-value for aligning a soft-clipped sequence to the TE library", "1e-5"),
    ("Distance of an annotated TE from a breakpoint", "20 bp or less"),
    ("Consecutive SNPs with a missing genotype call for an SNP-based call", "5 or more"),
    ("Uniquely mappable bases deleted for an SNP-based call to be confirmed", "at least 80%"),
    ("Fraction of a gene inside an event for the gene to count as deleted", "at least 50%"),
    ("Software", "BWA-MEM 0.7.18; SAMtools 1.21; BEDTools 2.31.1; Delly 1.3.3; RepeatMasker 4.1.9; BLAST+ 2.16.0"),
], columns=["Parameter", "Value"])

# ------------------------------------------------------------------ write workbook
xlsx = OUT + "Data_S10_updated_all-base_deletions.xlsx"
sheets = [("Events", events), ("Regions_per_strain", regions), ("Strains", strains), ("DataS10_checked", s10_out),
          ("SNP_calls_checked", snp_out), ("DataS9_checked", s9_out), ("Windows_100kb", win), ("Strains_226", p226), ("Parameters", params)]
with pd.ExcelWriter(xlsx, engine="openpyxl") as wr:
    pd.DataFrame(readme).to_excel(wr, sheet_name="README", index=False, header=False)
    for name, df in sheets:
        df.to_excel(wr, sheet_name=name, index=False)
wb = load_workbook(xlsx)
arial, bold = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True)
head = PatternFill("solid", start_color="DDE5EE")
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for cell in row:
            cell.font = arial
    if ws.title == "README":
        ws.column_dimensions["A"].width = 34; ws.column_dimensions["B"].width = 130
        ws["A1"].font = Font(name="Arial", size=12, bold=True)
        for row in ws.iter_rows(min_row=2):
            row[0].font = bold
            row[1].alignment = Alignment(wrap_text=True, vertical="top"); row[0].alignment = Alignment(vertical="top")
            if row[0].value in ("Sheet", "Column"):
                row[0].fill = head; row[1].fill = head; row[1].font = bold
        continue
    for cell in ws[1]:
        cell.font = bold; cell.fill = head; cell.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"; ws.auto_filter.ref = ws.dimensions
    for i, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), 1):
        width = max(len(str(cl.value)) if cl.value is not None else 0 for cl in col)
        ws.column_dimensions[get_column_letter(i)].width = max(9, min(width + 2, 45))
wb.save(xlsx)

# ------------------------------------------------------------------ figures
for src, dst in ((R + "figure_windows100kb", "Figure_S4_updated_100kb_tracks"), (R226 + "figure_226_strains", "Figure_S4_all_226_strains")):
    for ext in (".png", ".pdf"):
        shutil.copyfile(src + ext, OUT + dst + ext)

# ------------------------------------------------------------------ legends
W = S["windows"]; g = S["genes"]
te_all = 1 - S["te_class_counts"]["no TE evidence"] / S["events"]
p = pd.read_csv(R226 + "figure_226_strains.per_strain.tsv", sep="\t").groupby("population").deleted_Mb.agg(["size", "median"])
legends = f"""# Supplementary files for the all-base deletion analysis

The original Data S9 and Data S10 are not changed.

## Figure_S4_updated_100kb_tracks (replaces Figure S4)

**Proposed legend.** Distribution of deletions in 33 Ethiopian MoE strains along the seven chromosomes of the E2 reference, in consecutive windows of 100 kb. Each bar is one window. Tracks from top to bottom: (1) fraction of the uniquely mappable bases deleted according to read depth at every base (all-base method; mean of the 33 strains); (2) fraction of SNP genotype calls that lie in SNP-based deletions (runs of five or more consecutive SNPs with a missing genotype call, among the 198,454 SNPs); (3) number of SNPs per kb of uniquely mappable sequence; (4) fraction of the window annotated as repeat; (5) number of effector genes. Grey: windows with fewer than 20 kb of uniquely mappable bases. Across the {W['n_with_snps']} windows with SNPs, the deleted fractions from the two methods are correlated (Spearman {W['spearman']['allbase_deleted_frac']['snpbased_deleted_frac']:.2f}); {W['both_any']} of the {W['snp_any']} windows with SNP-based deletions also have all-base deletions.

Source data: sheet `Windows_100kb`.

## Figure_S4_all_226_strains (new; all 226 strains)

**Proposed legend.** Deletions relative to the E2 reference in all 226 strains, called from read depth at every base. Each row is one strain. Strains are grouped by population: ET (n = {int(p.loc['ET', 'size'])}) at the top, then the strains with less than 90% membership in either population (n = {int(p.loc['neither (<90%)', 'size'])}), then KU (n = {int(p.loc['KU', 'size'])}). Within a group, strains are sorted by their total deleted sequence. The E2 genome is divided into consecutive windows of 100 kb; each column is one window, in order along the chromosomes, and white vertical lines separate the chromosomes. Colour: the fraction of the bases in a window that are deleted in a strain. A base is deleted when its depth of uniquely mapped reads is at most {S['depth_threshold_frac'] * 100:g}% of the strain's median depth. Only uniquely mappable, non-repeat bases are counted: bases outside annotated repeats to which the reads of the reference strain E2 map uniquely. Grey: windows with fewer than 20 kb of uniquely mappable, non-repeat sequence. Bars on the right: total deleted sequence of each strain in Mb (median {p.loc['ET', 'median']:.2f} Mb in ET and {p.loc['KU', 'median']:.2f} Mb in KU strains).

Limits to state with this figure: the calls are from read depth only (not merged, not tested with split reads or paired-end reads); sequence absent in KU strains may have been gained by E2 and the ET strains; strains with low sequencing coverage have more calls.

Source data: sheet `Strains_226`.

## Data_S10_updated_all-base_deletions.xlsx

**Proposed legend.** Deletions in 33 Ethiopian MoE strains relative to the E2 reference, called from read depth at every base and tested with split reads and paired-end reads. Sheet `Events` lists the {S['events']} deletion events with their strains, read evidence, TEs at the breakpoints, genes and effector genes. Sheet `Regions_per_strain` lists the {S['regions_strain_level']:,} strain-level regions. Sheets `DataS10_checked`, `SNP_calls_checked` and `DataS9_checked` give the result of the all-base test for the 124 SNP-based events, for all {S['snp_candidates_ge5']['n']:,} strain-level SNP-based calls, and for the manually inspected deletions. Column definitions and all thresholds are in the sheets `README` and `Parameters`.

## Values in the manuscript text that the all-base analysis changes

| Manuscript | All-base analysis |
|---|---|
| 124 deletions in 33 strains (SNP-based) | {S['events']} events ({S['regions_strain_level']:,} strain-level regions); the SNP-based method detects those with five or more SNPs |
| about 45% of deletions TE-associated (56 of 124) | a TE at a breakpoint in {te_all * 100:.0f}% of all events |
| 194 genes in the deletions, 36 effectors (18.6%) | {g['genes_at_least_half_deleted']} genes, {g['effectors']} effectors ({g['effector_share'] * 100:.1f}%); {g['genome_effector_share'] * 100:.1f}% of all genes are effectors (p = {g['hypergeometric_p']:.2g}) |
| Figure S4 legend: "similar profiles in A and B" | Spearman {W['spearman']['allbase_deleted_frac']['snpbased_deleted_frac']:.2f} per 100-kb window; {W['both_any']} of {W['snp_any']} windows with SNP-based deletions also have all-base deletions |

The methods section also needs a paragraph on the all-base method (report, section 3). The manuscript text and the supplementary Word file were not edited.
"""
open(OUT + "Supplement_legends.md", "w", encoding="utf-8").write(legends)
print("written to", OUT)
for name, df in sheets:
    print(f"  {name}: {len(df)} rows, {df.shape[1]} columns")
