#!/usr/bin/env python3
"""Read-pair evidence for every depth-based deletion call in one strain.

A stretch with no uniquely mapped reads can mean two different things:
  (1) the DNA is missing in the strain (a real deletion), or
  (2) DNA of about the same length is still there, but it is too different from E2 for the
      reads to map (diverged sequence), or it has been replaced by something else.
Read depth cannot tell these apart. Read pairs can. The two reads of a pair come from the
two ends of one DNA fragment, and the fragment length of a library is known.
  - If the DNA is missing, pairs that straddle the gap look too far apart on E2 by exactly
    the deletion size, and there are about as many of them as at any ordinary position.
  - If same-length DNA is still there, only unusually long fragments can straddle the gap,
    so there are few such pairs and they do NOT look too far apart.
  - If something long has been inserted instead, no pair straddles the gap at all.

For every call the two explanations are compared with a log-likelihood ratio (LLR) built from
the strain's own fragment-length distribution. Only pairs in which neither read is soft-clipped
on the side facing the interval are used, so that the expected numbers are simple to work out:
  LLR = sum_i [ log f(t_i - d) - log f(t_i) ] - lambda * (Z1 - Z0)
where t_i are the E2 distances of the straddling pairs, d is the deletion size (taken as the
best value between the inner and the outer interval), f is the fragment-length distribution,
lambda the number of read pairs starting per base, and Z0 / Z1 the number of ways a fragment
can straddle the gap under "nothing missing" / "d bases missing". LLR > 0 favours a deletion.

Other columns: reads with a gap in their own alignment (CIGAR D) over the interval, reads
mapped inside the interval and how many mismatches they carry, and reads next to the
interval whose mate is unmapped or maps far away.
"""
import argparse, re, subprocess, sys
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--strain", required=True)
ap.add_argument("--calls", default="calls_uniq/depth_calls.bwa.uniq.primary.tsv")
ap.add_argument("--bam-dir", default="bam_bwa")
ap.add_argument("--out", required=True)
ap.add_argument("--flank", type=int, default=1000)
ap.add_argument("--slack", type=int, default=20)
a = ap.parse_args()
s = a.strain
bam = f"{a.bam_dir}/{s}.bwa.markdup.bam"
calls = pd.read_csv(a.calls, sep="\t")
calls = calls[calls.strain == s].reset_index(drop=True)
MAXI = 100_000
EPS = 0.01                      # share of pairs allowed to be chimeric / wrongly placed
cig = re.compile(r"(\d+)([MIDNSHP=X])")
CHROMS = [f"Chr0{i}" for i in range(1, 8)]


def view(regions, flags_excl="0xD04"):
    """Primary, mapped, non-duplicate reads overlapping the regions (each read once)."""
    p = subprocess.Popen(["samtools", "view", "-F", flags_excl, bam, *regions], stdout=subprocess.PIPE, text=True)
    seen = set()
    for line in p.stdout:
        f = line.split("\t")
        key = (f[0], int(f[1]) & 0xC0)
        if key in seen:
            continue
        seen.add(key)
        tags = "\t".join(f[11:])
        nm = re.search(r"NM:i:(\d+)", tags); mq = re.search(r"MQ:i:(\d+)", tags)
        mc = re.search(r"MC:Z:(\d+)S", tags)                     # mate's alignment starts with a soft clip
        ops = cig.findall(f[5])
        pos = int(f[3]) - 1
        alen = sum(int(l) for l, o in ops if o in "MDN=X")
        maxd = max([int(l) for l, o in ops if o == "D"], default=0)
        # "facing" soft clip: the read (forward) or its mate (reverse) stops matching on the side that faces the other read
        clip = (int(ops[-1][0]) if ops and ops[-1][1] == "S" else 0) >= 10 or (mc is not None and int(mc.group(1)) >= 10)
        yield (int(f[1]), pos, pos + alen, int(f[4]), f[6], int(f[7]) - 1, int(f[8]),
               int(nm.group(1)) if nm else 0, int(mq.group(1)) if mq else 60, maxd, alen, clip)
    p.wait()


# ---------------------------------------------------------------- the strain's library: fragment lengths and pair density
sizes = {}
for l in subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True, check=True).stdout.splitlines():
    c, n = l.split("\t")[:2]
    sizes[c] = int(n)
tl, binc, rlens = [], [], []
for c in CHROMS:
    mid = sizes[c] // 2
    lo = mid - 500_000
    cnt = np.zeros(100, int)
    for fl, pos, aend, mapq, rnext, pnext, tlen, nm, mq, maxd, alen, clip in view([f"{c}:{lo + 1}-{lo + 1_000_000}"], "0xD0C"):
        if mapq >= 20 and mq >= 20 and not clip and not fl & 0x10 and fl & 0x20 and rnext == "=" and 0 < tlen and lo <= pos < lo + 1_000_000:
            if tlen <= MAXI:
                tl.append(tlen)
            cnt[(pos - lo) // 10_000] += 1
            rlens.append(alen)
    binc.extend(cnt)
tl = np.array(tl)
lam = float(np.median(binc)) / 10_000                   # straddling-capable pairs starting per base
RL = int(np.median(rlens))
f = np.bincount(tl, minlength=MAXI + 1).astype(float)
k = np.ones(11) / 11                                    # light smoothing
f = np.convolve(f, k, mode="same")
f = f / f.sum()
logf = np.log((1 - EPS) * f + EPS / MAXI)               # a stray chimeric pair cannot decide the outcome on its own
F1 = np.cumsum(f[::-1])[::-1]                           # P(X >= k)
G = np.cumsum(F1[::-1])[::-1]                           # sum_x f(x) * (x - k + 1) for x >= k
lib_med, lib_q = float(np.median(tl)), np.percentile(tl, [5, 95])
print(f"{s}: pairs sampled {tl.size}; fragment median {lib_med:.0f} (5-95%: {lib_q[0]:.0f}-{lib_q[1]:.0f}); "
      f"read length {RL}; pairs per base {lam:.3f}", file=sys.stderr)


def Z(L):
    L = np.asarray(L)
    return G[np.clip(L, 0, MAXI)]


def lf(x):
    x = np.asarray(x)
    return np.where((x >= 0) & (x <= MAXI), logf[np.clip(x, 0, MAXI)], np.log(EPS / MAXI))


rows = []
for r in calls.itertuples():
    A, B = r.start + a.slack, r.end - a.slack
    if B <= A:
        A = B = (r.start + r.end) // 2
    span, outer = r.end - r.start, r.outer_end - r.outer_start
    left = f"{r.chrom}:{max(1, r.outer_start - a.flank + 1)}-{r.start + a.slack}"
    right = f"{r.chrom}:{max(1, r.end - a.slack + 1)}-{min(sizes[r.chrom], r.outer_end + a.flank)}"
    regions = [left, right] if span > 3000 else [f"{r.chrom}:{max(1, r.outer_start - a.flank + 1)}-{min(sizes[r.chrom], r.outer_end + a.flank)}"]
    t = []; n_dcig = 0; ins = []; anch = {"L": [0, 0, 0], "R": [0, 0, 0]}
    for fl, pos, aend, mapq, rnext, pnext, tlen, nm, mq, maxd, alen, clip in view(regions):
        fwd = not fl & 0x10
        mate_unmapped = bool(fl & 0x8)
        mate_far = (not mate_unmapped) and (rnext != "=" or abs(tlen) > outer + 20_000)
        if mapq >= 20:
            # pairs straddling the interval: forward read ends left of it, its mate starts right of it
            if (fwd and not clip and not mate_unmapped and fl & 0x20 and rnext == "=" and tlen > 0 and mq >= 20
                    and pos >= r.outer_start - a.flank and aend <= A and pnext >= B):
                t.append(tlen)
            if maxd >= 20 and pos < r.outer_end and aend > r.outer_start:
                n_dcig += 1
            # reads just outside the interval that point towards it: where is the mate?
            if fwd and r.outer_start - 400 <= pos and aend <= A:
                anch["L"][0] += 1; anch["L"][1] += mate_unmapped; anch["L"][2] += mate_far
            if not fwd and pos >= B and aend <= r.outer_end + 400:
                anch["R"][0] += 1; anch["R"][1] += mate_unmapped; anch["R"][2] += mate_far
        if span <= 3000 and pos >= r.start and aend <= r.end:
            ins.append((mapq, nm, alen))
    t = np.array(t, dtype=int)
    L0 = (B - A) + 2 * RL
    # deletion sizes that are possible: at least the deleted bases scored, at most the outer interval
    low = min(span, int(r.deleted_bases)) if hasattr(r, "deleted_bases") else span
    d_lo, d_hi = max(1, low - 3 * a.slack), outer + 3 * a.slack
    step = max(1, (d_hi - d_lo) // 400)
    ds = np.arange(d_lo, d_hi + 1, step)
    z0 = float(Z(L0))
    z1 = Z(np.maximum(L0 - ds, 2 * RL))
    base = lf(t).sum() if t.size else 0.0
    ll = np.array([lf(t - d).sum() if t.size else 0.0 for d in ds]) - base - lam * (z1 - z0)
    j = int(np.argmax(ll))
    ins = np.array(ins).reshape(-1, 3)
    row = dict(strain=s, id=f"{r.chrom}:{r.start}-{r.end}", span=span, outer_span=outer,
               n_straddling_pairs=int(t.size), median_pair_distance=float(np.median(t)) if t.size else np.nan,
               lib_median=lib_med, expected_pairs_if_present=round(lam * z0, 2), expected_pairs_if_deleted=round(lam * float(z1[j]), 2),
               best_deletion_size=int(ds[j]), llr_deletion_vs_present=round(float(ll[j]), 2),
               reads_with_alignment_gap=n_dcig,
               reads_inside=len(ins), reads_inside_mapq20=int((ins[:, 0] >= 20).sum()) if len(ins) else 0,
               mismatches_per_100bp_inside=round(100 * ins[:, 1].sum() / max(ins[:, 2].sum(), 1), 3) if len(ins) else np.nan)
    for side in "LR":
        n, nu, nf = anch[side]
        row.update({f"{side}_anchor_reads": n, f"{side}_mate_unmapped": int(nu), f"{side}_mate_far": int(nf)})
    rows.append(row)
pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
print(s, len(rows), "intervals written to", a.out, file=sys.stderr)
