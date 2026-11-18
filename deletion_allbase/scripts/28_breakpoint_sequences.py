#!/usr/bin/env python3
"""Sequences at the breakpoints of every strain-level deleted region (run on the cluster, one strain per call).

For each region and each side with a breakpoint defined by soft-clipped reads (L_bp / R_bp from script 12):
  - the E2 sequence on both sides of the breakpoint (outside the deletion, and inside it);
  - the consensus of the soft-clipped parts of the reads that end at the breakpoint. This is the sequence
    that follows the flank in the strain (the other flank, a TE, or other sequence);
  - the best match of that consensus in the TE library (family, entry, position in the entry, strand);
  - where that consensus aligns in E2 (best hit and number of hits).

Coordinates are 0-based. L_bp is the first base after the left flank (first deleted base); R_bp is the
first base of the right flank (first base after the deletion). The deletion defined by reads is [L_bp, R_bp).
"""
import argparse, collections, re, subprocess, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--strain", required=True)
ap.add_argument("--regions", default="inputs/regions_breakpoints.tsv")
ap.add_argument("--bam-dir", default="bam_bwa")
ap.add_argument("--ref", default="ref/E2.fasta")
ap.add_argument("--te-db", default="te_lib/edta/EDTA")
ap.add_argument("--name-map", default="te_lib/edta/name_map.tsv")
ap.add_argument("--e2-db", default="tmp/E2_blastdb/E2")
ap.add_argument("--out-dir", default="bpseq")
ap.add_argument("--flank", type=int, default=40)
ap.add_argument("--min-clip", type=int, default=10)
ap.add_argument("--max-len", type=int, default=120)
a = ap.parse_args()
s = a.strain

# ---- reference
ref, name, buf = {}, None, []
for line in open(a.ref):
    if line.startswith(">"):
        if name:
            ref[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
ref[name] = "".join(buf)

reg = pd.read_csv(a.regions, sep="\t")
reg = reg[reg.strain == s].reset_index(drop=True)
bam = f"{a.bam_dir}/{s}.bwa.markdup.bam"
cig = re.compile(r"(\d+)([MIDNSHP=X])")


def clipped_reads(chrom, bp, side):
    """soft-clipped parts of the reads whose alignment ends (L) or starts (R) exactly at bp"""
    out = subprocess.run(["samtools", "view", "-F", "0x900", bam, f"{chrom}:{max(1, bp - 2)}-{bp + 3}"],
                         capture_output=True, text=True, check=True).stdout
    seqs = []
    for line in out.splitlines():
        f = line.split("\t", 10)
        ops = cig.findall(f[5])
        if not ops:
            continue
        pos = int(f[3]) - 1
        aend = pos + sum(int(l) for l, o in ops if o in "MDN=X")
        if side == "L" and ops[-1][1] == "S" and int(ops[-1][0]) >= a.min_clip and aend == bp:
            seqs.append(f[9][-int(ops[-1][0]):])
        if side == "R" and ops[0][1] == "S" and int(ops[0][0]) >= a.min_clip and pos == bp:
            seqs.append(f[9][:int(ops[0][0])])
    return seqs


def consensus(seqs, side):
    """majority base per column; L: columns counted from the breakpoint to the right, R: to the left"""
    if not seqs:
        return ""
    if side == "R":
        seqs = [x[::-1] for x in seqs]
    need = 2 if len(seqs) >= 2 else 1
    cons = []
    for i in range(min(a.max_len, max(len(x) for x in seqs))):
        col = collections.Counter(x[i] for x in seqs if len(x) > i)
        if sum(col.values()) < need:
            break
        cons.append(col.most_common(1)[0][0])
    c = "".join(cons)
    return c[::-1] if side == "R" else c


rows, fa = [], []
for r in reg.itertuples():
    row = dict(strain=s, event_id=r.event_id, chrom=r.chrom, start=r.start, end=r.end, outer_start=r.outer_start, outer_end=r.outer_end)
    seq = ref[r.chrom]
    for side, bp in (("L", int(r.L_bp)), ("R", int(r.R_bp))):
        row[f"{side}_bp"] = bp
        if bp < 0:
            row.update({f"{side}_n_reads": 0, f"{side}_E2_before": "", f"{side}_E2_after": "", f"{side}_clip_consensus": ""})
            continue
        reads = clipped_reads(r.chrom, bp, side)
        cons = consensus(reads, side)
        row.update({f"{side}_n_reads": len(reads), f"{side}_E2_before": seq[max(0, bp - a.flank):bp], f"{side}_E2_after": seq[bp:bp + a.flank],
                    f"{side}_clip_consensus": cons})
        if len(cons) >= 15:
            fa.append((f"{len(rows)}_{side}", cons))
    rows.append(row)

qfa = f"{a.out_dir}/{s}.clip_consensus.fa"
with open(qfa, "w") as fh:
    for k, c in fa:
        fh.write(f">{k}\n{c}\n")
labels = dict(pd.read_csv(a.name_map, sep="\t")[["lib_id", "label"]].values)
te, e2 = {}, {}
if fa:
    fmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
    out = subprocess.run(["blastn", "-task", "blastn", "-query", qfa, "-db", a.te_db, "-evalue", "1e-5", "-outfmt", fmt,
                          "-max_target_seqs", "5", "-num_threads", "2"], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        f = l.split("\t")
        if f[0] not in te:                     # best hit only (first line per query)
            lib = f[1].split("#")[0]
            te[f[0]] = dict(te_family=labels.get(lib) or lib, te_entry=f[1], te_identity=float(f[2]), te_aligned_bp=int(f[3]),
                            te_start=min(int(f[6]), int(f[7])), te_end=max(int(f[6]), int(f[7])), te_length=int(f[11]),
                            te_strand="+" if int(f[7]) >= int(f[6]) else "-")
    out = subprocess.run(["blastn", "-task", "blastn", "-query", qfa, "-db", a.e2_db, "-evalue", "1e-5", "-outfmt", fmt,
                          "-max_target_seqs", "50", "-max_hsps", "1", "-num_threads", "2"], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        f = l.split("\t")
        d = e2.setdefault(f[0], dict(n=0))
        d["n"] += 1
        if d["n"] == 1:
            d.update(E2_hit=f"{f[1]}:{min(int(f[6]), int(f[7])) - 1}-{max(int(f[6]), int(f[7]))}", E2_hit_strand="+" if int(f[7]) >= int(f[6]) else "-",
                     E2_hit_identity=float(f[2]), E2_hit_aligned_bp=int(f[3]))
for i, row in enumerate(rows):
    for side in "LR":
        t = te.get(f"{i}_{side}", {}); h = e2.get(f"{i}_{side}", {})
        for k in ("te_family", "te_entry", "te_identity", "te_aligned_bp", "te_start", "te_end", "te_length", "te_strand"):
            row[f"{side}_{k}"] = t.get(k, "")
        row[f"{side}_E2_hit"] = h.get("E2_hit", ""); row[f"{side}_E2_hit_strand"] = h.get("E2_hit_strand", "")
        row[f"{side}_E2_hit_identity"] = h.get("E2_hit_identity", ""); row[f"{side}_E2_hit_aligned_bp"] = h.get("E2_hit_aligned_bp", "")
        row[f"{side}_E2_n_hits"] = h.get("n", 0) if row[f"{side}_clip_consensus"] else ""
pd.DataFrame(rows).to_csv(f"{a.out_dir}/{s}.breakpoint_sequences.tsv", sep="\t", index=False)
print(s, len(rows), "regions;", len(fa), "breakpoints with a consensus;", len(te), "match a TE", file=sys.stderr)
