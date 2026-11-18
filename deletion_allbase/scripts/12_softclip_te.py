#!/usr/bin/env python3
"""Split-read evidence at the edges of every depth-based deletion call in one strain.

BWA-MEM soft-clips a read where it stops matching the reference. At a deletion breakpoint
the clipped part of the read is the sequence that sits next to the breakpoint in the strain:
  - if it is the other side of the deletion, the read's supplementary alignment (SA tag)
    lands at the other edge  -> a simple junction;
  - if it is a transposable element, the clipped part matches the TE library, even when the
    reference has no TE copy at that place (a TE that is only present in the strain).

For the left edge, reads are used whose alignment ends inside the left breakpoint zone
[outer_start, start] (clip on the read's right side); for the right edge, reads whose
alignment starts inside [end, outer_end] (clip on the read's left side).
"""
import argparse, collections, re, subprocess, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--strain", required=True)
ap.add_argument("--calls", default="calls_uniq/depth_calls.bwa.uniq.primary.tsv")
ap.add_argument("--bam-dir", default="bam_bwa")
ap.add_argument("--te-db", default="te_lib/named_TEs")
ap.add_argument("--name-map", default=None, help="table lib_id / label (24_te_name_map.py); default: name before #")
ap.add_argument("--out-dir", default="softclip")
ap.add_argument("--pad", type=int, default=30)
ap.add_argument("--min-clip", type=int, default=25)
a = ap.parse_args()
s = a.strain
calls = pd.read_csv(a.calls, sep="\t")
calls = calls[calls.strain == s].reset_index(drop=True)
calls["id"] = calls.chrom + ":" + calls.start.astype(str) + "-" + calls.end.astype(str)

zones = []   # chrom, zstart, zend, call index, side
for i, r in calls.iterrows():
    zones.append((r.chrom, max(0, r.outer_start - a.pad), r.start + a.pad, i, "L"))
    zones.append((r.chrom, max(0, r.end - a.pad), r.outer_end + a.pad, i, "R"))
bed = f"{a.out_dir}/{s}.zones.bed"
with open(bed, "w") as fh:
    for c, zs, ze, i, side in sorted(zones):
        fh.write(f"{c}\t{zs}\t{ze}\n")
by_chrom = collections.defaultdict(list)
for z in zones:
    by_chrom[z[0]].append(z)

cig = re.compile(r"(\d+)([MIDNSHP=X])")
recs = []            # (call index, side, clip position, read id, clipped sequence, SA chrom, SA pos, SA mapq)
p = subprocess.Popen(["samtools", "view", "-M", "-L", bed, "-F", "0x900", f"{a.bam_dir}/{s}.bwa.markdup.bam"],
                     stdout=subprocess.PIPE, text=True)
for n, line in enumerate(p.stdout):
    f = line.split("\t", 11)
    ops = cig.findall(f[5])
    if not ops or "S" not in f[5]:
        continue
    pos = int(f[3]) - 1
    aend = pos + sum(int(l) for l, o in ops if o in "MDN=X")
    sa = re.search(r"SA:Z:([^,]+),(\d+),[+-],[^,]+,(\d+)", f[11]) if len(f) > 11 else None
    sa = (sa.group(1), int(sa.group(2)), int(sa.group(3))) if sa else ("", -1, -1)
    seq = f[9]
    cand = []
    if ops[-1][1] == "S" and int(ops[-1][0]) >= a.min_clip:
        cand.append(("L", aend, seq[-int(ops[-1][0]):]))       # read runs out of reference at the left edge
    if ops[0][1] == "S" and int(ops[0][0]) >= a.min_clip:
        cand.append(("R", pos, seq[:int(ops[0][0])]))          # read starts at the right edge
    for side, bp, clip in cand:
        for c, zs, ze, i, zside in by_chrom[f[2]]:
            if zside == side and zs <= bp <= ze:
                recs.append((i, side, bp, f"r{n}", clip, *sa))
p.wait()

fa = f"{a.out_dir}/{s}.clips.fa"
with open(fa, "w") as fh:
    for k, r in enumerate(recs):
        fh.write(f">{k}\n{r[4]}\n")
hits = {}
labels = dict(pd.read_csv(a.name_map, sep="\t")[["lib_id", "label"]].values) if a.name_map else None
if recs:
    out = subprocess.run(["blastn", "-task", "blastn", "-query", fa, "-db", a.te_db, "-evalue", "1e-5",
                          "-outfmt", "6 qseqid sseqid pident length", "-max_target_seqs", "5", "-num_threads", "2"],
                         capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        q, sj, pid, ln = l.split("\t")
        fam = labels.get(sj.split("#")[0]) if labels else re.sub(r"(_\d+)?#.*", "", sj)
        hits.setdefault(int(q), fam or sj.split("#")[0])        # best hit only (first line per query)

rows = []
D = pd.DataFrame([(k, *r[:4], len(r[4]), *r[5:], hits.get(k, "")) for k, r in enumerate(recs)],
                 columns=["k", "call", "side", "bp", "read", "cliplen", "sa_chrom", "sa_pos", "sa_mapq", "te"])
for i, r in calls.iterrows():
    row = dict(strain=s, id=r.id)
    for side in "LR":
        d = D[(D.call == i) & (D.side == side)]
        row[f"{side}_clipped_reads"] = len(d)
        if len(d) == 0:
            row.update({f"{side}_bp": -1, f"{side}_bp_reads": 0, f"{side}_te_reads": 0, f"{side}_te": "", f"{side}_sa_other_edge": 0,
                        f"{side}_sa_elsewhere": 0})
            continue
        bp, nbp = d.bp.value_counts().index[0], int(d.bp.value_counts().iloc[0])
        d = d[(d.bp - bp).abs() <= 3]                                  # reads agreeing on the modal breakpoint
        te = d[d.te != ""].te.value_counts()
        # does the clipped part map at the other edge of this call?
        ozs, oze = (r.end - 200, r.outer_end + 200) if side == "L" else (r.outer_start - 200, r.start + 200)
        other = (d.sa_chrom == r.chrom) & (d.sa_pos >= ozs) & (d.sa_pos <= oze)
        row.update({f"{side}_bp": int(bp), f"{side}_bp_reads": len(d), f"{side}_te_reads": int(te.sum()),
                    f"{side}_te": ",".join(f"{k}:{v}" for k, v in te.items()), f"{side}_sa_other_edge": int(other.sum()),
                    f"{side}_sa_elsewhere": int(((d.sa_chrom != "") & ~other).sum())})
    rows.append(row)
pd.DataFrame(rows).to_csv(f"{a.out_dir}/{s}.softclip_edges.tsv", sep="\t", index=False)
print(s, len(calls), "calls;", len(recs), "clipped reads at edges;", len(hits), "clips matching a TE", file=sys.stderr)
