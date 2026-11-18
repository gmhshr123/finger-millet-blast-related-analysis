#!/usr/bin/env python3
"""Workbook of breakpoint sequences and TEs for every deleted region in the 33 strains (run from manuscripts_all/).

Input: deletion_allbase/results_edta/cluster_uniq/breakpoint_sequences.tsv (script 28, downloaded from the cluster)
and the result tables of scripts 10 and 16.
Output: deletion_allbase/supplement/Deletion_breakpoints_and_TEs_33_strains.xlsx
  Breakpoint_sequences   the main table: two rows per breakpoint (one row for a region without read-defined breakpoints),
                         E2 sequence above the strain sequence in one column, aligned base by base, and whether the SNP-based method found
                         the deletion; SNP-based calls with no all-base region are added as rows
  Read_me                definitions
  Summary                counts (computed here from the same rows)
and deletion_allbase/results_edta/breakpoints_long.tsv, regions_for_browser.tsv (input of script 31, the HTML browser).

Positions in this workbook are 1-based.
"""
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

R = os.environ.get("DEL_RESULTS", "deletion_allbase/results_edta/")
I = "deletion_allbase/inputs/"
OUT = "deletion_allbase/supplement/Deletion_breakpoints_and_TEs_33_strains.xlsx"
cohort = [l.strip() for l in open(I + "cohort33.txt")]

BLACK, RED, PURPLE, BLUE = "000000", "C00000", "7030A0", "0070C0"
mono = lambda color, bold=False: InlineFont(rFont="Courier New", sz=10, color=color, b=bold)
ARIAL, BOLD = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True)
HEAD = PatternFill("solid", start_color="DDE5EE")
TOP = Alignment(wrap_text=True, vertical="top")

TYPE = {"LTR/Gypsy": "LTR retrotransposon (Gypsy)", "LTR/Copia": "LTR retrotransposon (Copia)", "LTR/unknown": "LTR retrotransposon",
        "DNA/DTA": "DNA transposon (hAT)", "DNA/DTC": "DNA transposon (CACTA)", "DNA/DTH": "DNA transposon (PIF/Harbinger)",
        "DNA/DTM": "DNA transposon (Mutator)", "DNA/DTT": "DNA transposon (Tc1/Mariner)", "DNA/Helitron": "DNA transposon (Helitron)",
        "DNA/unknown": "DNA transposon", "MITE/DTA": "MITE (hAT)", "MITE/DTC": "MITE (CACTA)", "MITE/DTH": "MITE (PIF/Harbinger)",
        "MITE/DTM": "MITE (Mutator)", "LINE/CRE": "LINE (CRE)", "LINE/I-Jockey": "LINE (I/Jockey)", "LINE/Tad1": "LINE (Tad1)", "unknown": "unknown"}
PART = {"LTR": "long terminal repeat (LTR)", "INT": "internal region"}
nm = pd.read_csv(R + "cluster_uniq/EDTA_name_map.tsv", sep="\t")
te_type = {k: TYPE.get(v, v) for k, v in zip(nm.lib_id, nm.te_class)}
fam = lambda x: "" if pd.isna(x) else str(x).replace("MAGGY_Fosbury", "MAGGY/Fosbury").replace(",", ", ")
LABEL = {"simple deletion": "simple deletion", "other sequence at edge": "other sequence at the edge",
         "DNA still present": "DNA still present", "depth only": "depth only"}
CLASS_DEF = {
    "simple deletion": "Discordant read pairs or split reads join the two flanks directly. The region is deleted and nothing else is at the breakpoints.",
    "other sequence at the edge": "At one or both breakpoints, the soft-clipped reads continue into other sequence (a TE or other DNA) instead of the opposite flank. "
                                  "The region is missing at this position, and other sequence is joined to the flank.",
    "DNA still present": "Read pairs span the region at the normal insert size (concordant pairs). The low read depth is not caused by a deletion; this region is not a deletion.",
    "depth only": "No informative reads at the breakpoints. The deletion is supported by low read depth alone.",
}
# origin of the clipped sequence: (label in the workbook, definition)
ORIGIN = {
    "TE": ("transposable element (TE)", "Matches an entry of the TE library (BLASTN, e-value 1e-5)."),
    "the other flank of this deletion": ("opposite flank (direct deletion junction)",
                                         "Aligns in E2 at the other breakpoint of the same deletion: the two flanks are joined directly. "
                                         "When the opposite flank itself lies in a TE, the clipped sequence also matches the TE library; "
                                         "such breakpoints are still direct junctions, and 'TE match' notes the TE."),
    "E2 sequence from inside the deleted region": ("sequence from within the deleted interval",
                                                   "Aligns in E2 inside the deleted region itself. Part of the region is retained elsewhere: a more complex rearrangement than one deletion."),
    "other E2 sequence": ("sequence from elsewhere in E2 (1 to 4 copies)", "Aligns at another place in E2, at 1 to 4 sites, and is not a library TE."),
    "repeated E2 sequence (5 or more copies)": ("repetitive E2 sequence (5 or more copies)",
                                                "Aligns at 5 or more sites in E2 and is not a library TE: an unclassified repeat."),
    "sequence not found in E2": ("no match in E2 or the TE library",
                                 "15 bases or longer, but no BLASTN match in E2 or the TE library (e-value 1e-5). Sequence absent from E2, or a consensus with read errors."),
    "too short to identify": ("too short to classify (< 15 bases)", "Fewer than 15 bases; not searched."),
    "adapter": ("sequencing adapter (artefact)", "Contains the Illumina adapter (AGATCGGAAGAGC). The read continued past the end of a short DNA fragment into the adapter; "
                                                 "this is not sequence of the strain, and the breakpoint may be an artefact."),
}

b = pd.read_csv(R + "cluster_uniq/breakpoint_sequences.tsv", sep="\t", keep_default_na=False)
c = pd.read_csv(R + "allbase_calls_strain_level_read_evidence.tsv", sep="\t")
e = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t").set_index("event_id")
key = ["strain", "chrom", "start", "end"]
m = b.merge(c[key + ["deleted_bases", "read_class"]], on=key, how="left", validate="1:1")
assert len(m) == len(c) and m.read_class.notna().all()
m["order"] = m.strain.map({s: i for i, s in enumerate(cohort)})
m = m.sort_values(["order", "chrom", "start"]).reset_index(drop=True)


def follows(row, side):
    """what the soft-clipped consensus is: TE, the other flank, other E2 sequence, or sequence not found in E2"""
    cons = row[f"{side}_clip_consensus"]
    if not cons:
        return ""
    if "AGATCGGAAGAGC" in cons or "GCTCTTCCGATCT" in cons:   # Illumina adapter: read-through of a short fragment
        return "adapter"
    hit = row[f"{side}_E2_hit"]
    if hit:
        ch, se = hit.rsplit(":", 1); hs, he = map(int, se.split("-"))
        other = int(row["R_bp"] if side == "L" else row["L_bp"])
        lo, hi = (row["end"] - 300, row["outer_end"] + 300) if side == "L" else (row["outer_start"] - 300, row["start"] + 300)
        # checked before the TE match: when the opposite flank lies in a TE, the junction is still direct
        if ch == row["chrom"] and ((other >= 0 and hs - 300 <= other <= he + 300) or (hs < hi and he > lo)):
            return "the other flank of this deletion"
    if row[f"{side}_te_family"]:
        return "TE"
    if hit:
        if ch == row["chrom"] and hs < row["outer_end"] and he > row["outer_start"]:
            return "E2 sequence from inside the deleted region"
        return "other E2 sequence" if int(row[f"{side}_E2_n_hits"]) < 5 else "repeated E2 sequence (5 or more copies)"
    return "sequence not found in E2" if len(cons) >= 15 else "too short to identify"


GREY = "A6A6A6"
W = 40                      # bases shown on each side of the breakpoint


def line(parts):
    """one cell with one monospace line; parts = [(text, colour), ...]"""
    return CellRichText(*[TextBlock(mono(col, col in (RED, PURPLE, BLUE)), t) for t, col in parts if t])


def alignment(side, before, after, clip, col):
    """two cells, E2 and the strain, in adjacent columns of the same row: both are 40 bases | 40 bases,
    so '|' and each flank base are at the same character position in the two cells"""
    if side == "left":       # flank | deleted ; flank | clipped
        return (line([(before, BLACK), ("|", GREY), (after, RED)]),
                line([(before, BLACK), ("|", GREY), (clip[:W], col)]))
    c = clip[-W:]            # deleted | flank ; clipped | flank
    return (line([(before, RED), ("|", GREY), (after, BLACK)]),
            line([(" " * (W - len(c)) + c, col), ("|", GREY), (after, BLACK)]))   # padding joins the clip: Excel refuses a block of spaces only


# ------------------------------------------------------------------ SNP-based calls (five or more consecutive SNPs with a missing genotype)
snp = pd.read_csv(R + "snp_candidates_checked.tsv", sep="\t")
snp = snp[(snp.snp_count >= 5) & snp.strain.isin(cohort)].copy()
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")


def s10_id(x):
    hit = s10[(s10.chrom == x.chrom) & (s10.start <= x.end) & (s10.end >= x.start)]
    ids = [d for d, st in zip(hit.del_id, hit.strains) if x.strain in st.split(",")]
    return f"Data S10 {', '.join(ids)}" if ids else "not in Data S10"


snp["s10"] = [s10_id(x) for x in snp.itertuples()]
snp["used"] = False


def snp_text(x):
    return f"{x.s10}; {int(x.snp_count)} missing SNPs, {int(x.start):,}-{int(x.end):,}"


m = m.merge(c[key + ["snps_inside", "found_by_snp_method"]], on=key, how="left", validate="1:1")


def snp_for_region(r):
    o = snp[(snp.strain == r["strain"]) & (snp.chrom == r["chrom"]) & (snp.start - 1 < r["outer_end"]) & (snp.end > r["outer_start"])]
    snp.loc[o.index, "used"] = True
    assert bool(len(o)) == bool(r["found_by_snp_method"])
    if len(o):
        return "Yes: " + " | ".join(snp_text(x) for x in o.itertuples())
    n = int(r["snps_inside"])
    if n == 0:
        return "No: no SNPs in this region"
    if n < 5:
        return f"No: only {n} SNP{'s' if n > 1 else ''} in this region (5 needed)"
    return f"No: {n} SNPs in this region, but fewer than 5 consecutive SNPs with a missing genotype"


regions, bps, rows = [], [], []     # rows: (sort key, list of cells) for the main sheet
for r in m.to_dict("records"):
    ev = e.loc[r["event_id"]]
    region = [r["strain"], r["chrom"], r["event_id"], int(ev.n_strains), int(r["start"]) + 1, int(r["end"]), int(r["end"]) - int(r["start"]),
              LABEL[r["read_class"]], snp_for_region(r)]
    order = (r["order"], r["chrom"], int(r["start"]))
    side_info = {}
    for side, name in (("L", "left"), ("R", "right")):
        bp = int(r[f"{side}_bp"])
        if bp < 0:
            side_info[side] = None
            continue
        what = follows(r, side)
        cons = r[f"{side}_clip_consensus"]
        before, after = r[f"{side}_E2_before"], r[f"{side}_E2_after"]
        hit = r[f"{side}_E2_hit"]
        if hit:
            ch, se = hit.rsplit(":", 1); hs, he = se.split("-")
            hit = f"{ch}:{int(hs) + 1}-{he} ({r[f'{side}_E2_hit_strand']})"
        lib = r[f"{side}_te_entry"].split("#")[0] if r[f"{side}_te_entry"] else ""
        is_te = what == "TE"
        info = dict(strain=r["strain"], event=r["event_id"], chrom=r["chrom"], side=name, pos=bp if side == "L" else bp + 1,
                    reads=int(r[f"{side}_n_reads"]), clip=cons, origin=ORIGIN[what][0] if what else "",
                    family=fam(r[f"{side}_te_family"]) if is_te else "", tetype=te_type.get(lib, "") if is_te else "",
                    entry=lib if is_te else "", part=PART.get(lib.rsplit("_", 1)[-1], "") if is_te else "",
                    te_pos=f"{r[f'{side}_te_start']}-{r[f'{side}_te_end']} of {r[f'{side}_te_length']}" if is_te else "",
                    strand=r[f"{side}_te_strand"] if is_te else "",
                    ident=float(r[f"{side}_te_identity"]) if is_te and r[f"{side}_te_identity"] != "" else "",
                    e2hit=hit, e2n=int(r[f"{side}_E2_n_hits"]) if hit else "", e2ann=fam(ev[f"{side}_te_ref"]),
                    read_class=LABEL[r["read_class"]], e2_left=before, e2_right=after)
        te_desc = ""
        if is_te:
            te_desc = f"{info['tetype']}{', ' + info['part'] if info['part'] else ''}; bases {info['te_pos']} of {lib}, {info['strand']} strand, {info['ident']:g}% identity"
        info["flank_te"] = ""
        if what == "the other flank of this deletion" and r[f"{side}_te_family"]:   # direct junction; the opposite flank lies in a TE
            info["flank_te"] = (f"The opposite flank lies in a TE: {fam(r[f'{side}_te_family'])} ({te_type.get(lib, '')}), "
                                f"{float(r[f'{side}_te_identity']):g}% identity to {lib}")
            te_desc = info["flank_te"]
        bps.append(info); side_info[side] = info
        rows.append((order + (side,), region + [name, info["pos"], info["reads"],
                     *alignment(name, before, after, cons, PURPLE if is_te else BLUE),
                     info["origin"], info["family"], te_desc, info["e2ann"], cons]))
    if not side_info["L"] and not side_info["R"]:
        ann = "; ".join(f"{s}: {fam(ev[f'{k}_te_ref'])}" for k, s in (("L", "left"), ("R", "right")) if fam(ev[f"{k}_te_ref"]))
        rows.append((order + ("",), region + ["not defined by reads", "", "", "", "", "", "", "", ann, ""]))
    L, Rr = side_info["L"], side_info["R"]
    te_reads = "; ".join(f"{x['side']}: {x['family']}" for x in (L, Rr) if x and x["family"])
    regions.append([r["strain"], r["event_id"], int(ev.n_strains), r["chrom"], int(r["start"]) + 1, int(r["end"]), int(r["end"]) - int(r["start"]),
                    int(r["deleted_bases"]), LABEL[r["read_class"]], L["pos"] if L else "", Rr["pos"] if Rr else "",
                    Rr["pos"] - L["pos"] - 1 if L and Rr and Rr["pos"] > L["pos"] + 1 else "",
                    L["origin"] if L else "", Rr["origin"] if Rr else "", te_reads, fam(ev["L_te_ref"]), fam(ev["R_te_ref"]), region[8]])

# SNP-based calls that overlap no all-base region of the same strain: added as rows of their own
ALLBASE = {"not deleted": "not called by read depth: {p} of the {n} uniquely mappable bases have low depth",
           "partly deleted": "not called by read depth: {p} of the {n} uniquely mappable bases have low depth",
           "confirmed deleted": "not called by read depth as a region: {p} of the {n} uniquely mappable bases have low depth",
           "not assessable (repeat / uncallable)": "not tested by read depth: only {n} uniquely mappable bases (30 needed)"}
only = snp[~snp.used]
order_of = {s: i for i, s in enumerate(cohort)}
snp_only = []
for x in only.itertuples():
    frac = "" if pd.isna(x.frac_deleted) else f"{x.frac_deleted:.0%}"
    snp_only.append([x.strain, x.chrom, int(x.start), int(x.end), int(x.snp_count), x.s10, ALLBASE[x.allbase].format(p=frac, n=int(x.scorable_bases))])
    rows.append(((order_of[x.strain], x.chrom, int(x.start) - 1, ""),
                 [x.strain, x.chrom, "", "", int(x.start), int(x.end), int(x.end) - int(x.start) + 1, "SNP-based call only",
                  f"Yes (SNP-based call only): {snp_text(x)}; {ALLBASE[x.allbase].format(p=frac, n=int(x.scorable_bases))}",
                  "", "", "", "", "", "", "", "", "", ""]))
rows.sort(key=lambda t: (t[0][0], t[0][1], t[0][2], {"left": 0, "right": 1}.get(t[0][3], 0)))

wb = Workbook()

# ------------------------------------------------------------------ Breakpoint_sequences (main sheet)
ws = wb.active; ws.title = "Breakpoint_sequences"
cols = [("Strain", 7), ("Chromosome", 9), ("Deletion event ID", 9), ("Strains with this deletion", 8),
        ("Region start (E2)", 10), ("Region end (E2)", 10), ("Region length (bp)", 9), ("Read evidence", 25),
        ("Found by the SNP-based method?", 40),
        ("Breakpoint", 11), ("Breakpoint position (E2)", 12), ("Soft-clipped reads", 8),
        ("Sequence of", 9),
        ("Sequence at the breakpoint: E2 reference (upper row) and strain from reads (lower row), 40 bases on each side; | = breakpoint", 94),
        ("Sequence joined to the flank in the strain", 24), ("TE family", 15), ("TE match", 44),
        ("TE annotated in E2 at this breakpoint", 20), ("Full clipped sequence (from reads)", 40)]
ws.append([n for n, _ in cols])
for i, (_, width) in enumerate(cols, 1):
    cell = ws.cell(row=1, column=i); cell.font = BOLD; cell.alignment = TOP; cell.fill = HEAD
    ws.column_dimensions[get_column_letter(i)].width = width
ws.row_dimensions[1].height = 70
LABCOL, SEQ, CLIP, WRAP = 13, 14, 19, {9, 17, 18}
GREYTXT = "FF808080"
PAIR_END = Border(bottom=Side(style="thin", color="BFBFBF"))
# a breakpoint with sequences takes two rows: E2 above, the strain below, in the same column so the bases line up.
# The lower row repeats the other columns (grey text), so filtering and sorting keep the two rows together.
out = []                    # (cells, lower row of a pair?, last row of the breakpoint?)
for _, cells in rows:
    e2seq, stseq = cells[12], cells[13]
    if isinstance(e2seq, CellRichText):
        out.append((cells[:12] + ["E2", e2seq] + cells[14:], False, False))
        out.append((cells[:12] + [cells[0], stseq] + cells[14:], True, True))
    else:
        out.append((cells[:12] + ["", ""] + cells[14:], False, True))
for k, (cells, lower, last) in enumerate(out, 2):
    ws.append(cells)
    paired = isinstance(ws.cell(row=k, column=SEQ).value, CellRichText)
    for cell in ws[k]:
        if cell.column == SEQ:      # one line, never wrapped, left-aligned, so E2 and the strain line up base by base
            cell.alignment = Alignment(horizontal="left", wrap_text=False, vertical="center")
        else:
            colour = GREYTXT if lower and cell.column != LABCOL else None
            cell.font = Font(name="Courier New", size=9, color=colour) if cell.column == CLIP else Font(name="Arial", size=10, color=colour)
            cell.alignment = Alignment(wrap_text=cell.column in WRAP and not paired, vertical="center",   # paired rows: one line each
                                       horizontal="center" if cell.column == LABCOL else None)
        if last:
            cell.border = PAIR_END
    if paired:
        ws.row_dimensions[k].height = 15
ws.freeze_panes = "C2"; ws.auto_filter.ref = ws.dimensions

# ------------------------------------------------------------------ Summary (values computed from the rows above)
sm = wb.create_sheet("Summary")
bp_df = pd.DataFrame(bps)
reg_df = pd.DataFrame([{"snp": cells[8], "cls": cells[7]} for _, cells in rows if cells[7] != "SNP-based call only"])
reg_df = reg_df.loc[~pd.Series([(cells[0], cells[2]) for _, cells in rows if cells[7] != "SNP-based call only"]).duplicated().values]
assert len(reg_df) == len(regions)
found = reg_df.snp.str.startswith("Yes")
tables = [
    ("Deleted regions (all-base method) and the SNP-based method", ["", "Regions"], [
        ["Found by the SNP-based method", int(found.sum())],
        ["Not found: no SNPs in the region", int(reg_df.snp.str.startswith("No: no SNPs").sum())],
        ["Not found: 1 to 4 SNPs in the region", int(reg_df.snp.str.startswith("No: only").sum())],
        ["Not found: 5 or more SNPs, but fewer than 5 consecutive missing genotypes", int(reg_df.snp.str.contains("but fewer than 5").sum())],
        ["All deleted regions", len(reg_df)]]),
    ("SNP-based calls (5 or more consecutive SNPs with a missing genotype in one strain)", ["", "Calls"], [
        ["Overlap a deleted region of the same strain", int(snp.used.sum())],
        ["No deleted region: low depth over 20% to 80% of the bases", int((~snp.used & (snp.allbase == "partly deleted")).sum())],
        ["No deleted region: low depth over less than 20% of the bases", int((~snp.used & (snp.allbase == "not deleted")).sum())],
        ["No deleted region: low depth over 80% or more of the bases", int((~snp.used & (snp.allbase == "confirmed deleted")).sum())],
        ["No deleted region: fewer than 30 uniquely mappable bases, not testable", int((~snp.used & snp.allbase.str.startswith("not assessable")).sum())],
        ["All SNP-based calls", len(snp)]]),
    ("Read evidence of the deleted regions", ["", "Regions"], [[k, int((reg_df.cls == k).sum())] for k in CLASS_DEF] + [["All deleted regions", len(reg_df)]]),
    ("Sequence joined to the flank at the read-defined breakpoints", ["", "Left", "Right", "Total"],
     [[lab, int(((bp_df.origin == lab) & (bp_df.side == "left")).sum()), int(((bp_df.origin == lab) & (bp_df.side == "right")).sum()),
       int((bp_df.origin == lab).sum())] for lab, _ in ORIGIN.values()] + [["All breakpoints", int((bp_df.side == "left").sum()), int((bp_df.side == "right").sum()), len(bp_df)]]),
    ("TE families joined to the flank", ["TE family", "Breakpoints", "Deletion events", "Strains", "TE classification"],
     [[f_, len(g), g.event.nunique(), g.strain.nunique(), ", ".join(sorted(set(g.tetype)))]
      for f_, g in sorted(bp_df[bp_df.family != ""].groupby("family"), key=lambda t: -len(t[1]))]),
]
for title, head, body in tables:
    sm.append([title]); sm.cell(row=sm.max_row, column=1).font = Font(name="Arial", size=11, bold=True)
    sm.append(head)
    for cell in sm[sm.max_row]:
        cell.font = BOLD; cell.fill = HEAD
    for row in body:
        sm.append(row)
        for cell in sm[sm.max_row]:
            cell.font = BOLD if str(row[0]).startswith("All ") else ARIAL
    sm.append([])
sm.column_dimensions["A"].width = 70
for col in "BCD":
    sm.column_dimensions[col].width = 12
sm.column_dimensions["E"].width = 40

# ------------------------------------------------------------------ Read_me
rd = wb.create_sheet("Read_me")
lines = [
    ("title", "Deletion breakpoints and transposable elements (TEs) in the 33 Ethiopian strains", ""),
    ("", "", ""),
    ("", "What the main sheet shows", "Sheet Breakpoint_sequences has two rows per breakpoint of each region deleted in a strain relative to the E2 reference: "
                                      "the E2 sequence and, below it, the strain sequence. The lower row repeats the other columns in grey, so that "
                                      "filtering and sorting keep the two rows together; to count breakpoints, count only the rows with 'Sequence of' = E2. "
                                      "A region whose breakpoints are not defined by reads has one row. SNP-based calls that overlap no deleted region are added as rows "
                                      "('SNP-based call only'). Rows are sorted by strain, chromosome and position. Positions are on the E2 assembly, 1-based."),
    ("", "Deleted region", "Uniquely mappable bases with a depth of uniquely mapped reads of at most 5% of the strain's median depth (all-base method). "
                           "Regions in different strains that overlap by 80% or more share one deletion event ID."),
    ("", "Breakpoint", "Left: last retained base before the deletion. Right: first retained base after it. Defined where soft-clipped reads end at the same base "
                       "(reads whose remaining part, at least 10 bases, does not align to E2 at that place). Breakpoints with 1 or 2 reads are weakly supported."),
    ("", "Sequence at the breakpoint", "Two rows per breakpoint, aligned base by base in the same column. Upper row ('Sequence of' = E2): the E2 reference. "
                                       "Lower row ('Sequence of' = the strain): the strain, from its reads. "
                                       "Black = flank, the same in E2 and the strain. Red = E2 sequence deleted in the strain. "
                                       "Purple = sequence joined to the flank that matches a TE. Blue = other sequence joined to the flank. "
                                       "| = the breakpoint. 40 bases are shown on each side; the full clipped sequence is in the last column."),
    ("", "Clipped sequence", "Consensus of the soft-clipped parts of the reads at the breakpoint (majority base at each position, at most 120 bases). "
                             "It can contain errors where few reads cover a position."),
    ("", "Found by the SNP-based method?", "The SNP-based method calls a deletion where five or more consecutive SNPs have a missing genotype in a strain. "
                                           "'Yes' gives the Data S10 deletion ID, the number of missing SNPs, and the positions of the first and last missing SNP. "
                                           "'No' gives the reason: no SNPs, fewer than 5 SNPs, or SNPs with genotype calls. "
                                           "'SNP-based call only' rows give the share of the uniquely mappable bases between the first and last missing SNP that have low depth."),
    ("", "", ""),
    ("head", "Read evidence", ""),
] + [("", k, v) for k, v in CLASS_DEF.items()] + [
    ("", "", ""),
    ("head", "Sequence joined to the flank in the strain", ""),
] + [("", lab, d) for lab, d in ORIGIN.values()] + [
    ("", "", ""),
    ("head", "TE columns", ""),
    ("", "TE family", "Best BLASTN match (e-value 1e-5) in the combined EDTA library (275 entries from B71, E2, MZ5-1-6 and P131). "
                      "The published name when the entry matches that TE at 80% identity or more; otherwise the EDTA label."),
    ("", "TE match", "Classification of the entry; for LTR retrotransposons, whether the long terminal repeat (LTR) or the internal region matched; "
                     "the bases of the entry matched, strand, and identity."),
    ("", "TE annotated in E2 at this breakpoint", "TE families in the E2 annotation (RepeatMasker, EDTA library) within 20 bp of the breakpoint. From the reference, not from reads. "
                                                  "When E2 has the same TE as the strain, reads align through it and are not clipped, so the TE is seen only here."),
]
for style, a_, b_ in lines:
    rd.append([a_, b_])
    i = rd.max_row
    A, B = rd.cell(row=i, column=1), rd.cell(row=i, column=2)
    A.alignment = B.alignment = TOP
    A.font, B.font = BOLD, ARIAL
    if style == "title":
        A.font = Font(name="Arial", size=13, bold=True)
    elif style == "head":
        A.font = Font(name="Arial", size=11, bold=True); A.fill = B.fill = HEAD
rd.column_dimensions["A"].width = 36; rd.column_dimensions["B"].width = 130

os.makedirs(os.path.dirname(OUT), exist_ok=True)
wb.save(OUT)

# tables for the HTML browser (script 31)
pd.DataFrame([{k: v for k, v in x.items()} for x in bps]).to_csv(R + "breakpoints_long.tsv", sep="\t", index=False)
rnames = ["Strain", "Deletion event ID", "Strains with this event", "Chromosome", "Region start", "Region end", "Region length (bp)", "Deleted bases",
          "Read evidence", "Left breakpoint", "Right breakpoint", "Deletion length between the breakpoints (bp)", "Left origin", "Right origin",
          "TE at a breakpoint, from reads", "TE annotated in E2 left", "TE annotated in E2 right", "Found by the SNP-based method?"]
pd.DataFrame(regions, columns=rnames).to_csv(R + "regions_for_browser.tsv", sep="\t", index=False)
pd.DataFrame(snp_only, columns=["strain", "chrom", "first_missing_snp", "last_missing_snp", "n_missing_snps", "s10", "read_depth"]).to_csv(
    R + "snp_only_for_browser.tsv", sep="\t", index=False)
print("written", OUT, "| rows", len(rows), "| regions", len(regions), "| breakpoints", len(bps), "| SNP calls", len(snp), "overlapping", int(snp.used.sum()),
      "| SNP-only rows", len(only), only.allbase.value_counts().to_dict(), "| regions found by SNP", int(found.sum()))
