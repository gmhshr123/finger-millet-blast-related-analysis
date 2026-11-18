#!/usr/bin/env python3
"""All-base deletion calling from per-base read depth.

Every base of the seven E2 chromosomes gets one of three states in each strain:
  present  - scorable base with unique-read depth above the threshold
  deleted  - scorable base with unique-read depth at or below the threshold
  unknown  - base that is not scorable

A base is "callable" if the reference strain's own reads (E2 mapped to E2) cover it with
uniquely mapped reads (MAPQ>=20) and it is not an assembly gap. If E2's own reads cannot
be placed there, short reads cannot tell presence from absence in any strain.
With --exclude-bed (annotated repeats) the scorable set is "callable and not repeat".
Depth inside repeats depends on which repeat copies a strain carries elsewhere in its
genome, so it is not used as evidence; repeats are bridged instead.

A deletion call is a stretch of deleted scorable bases. It may run through unknown bases
(repeats, gaps) and through short blips of present bases (< --min-present-run scorable
bases in a row, e.g. a few stray reads), but it is ended by a longer run of present bases.
start/end are the first and last deleted scorable base. outer_start/outer_end are the
nearest flanking runs of present scorable bases, so the true breakpoints lie in
[outer_start, start] and [end, outer_end].
"""
import argparse, os, sys
from multiprocessing import Pool
import numpy as np
import pandas as pd

CHROMS = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07"]
THR_FRACS = [0.0, 0.05, 0.10, 0.20]      # depth threshold as fraction of the strain's median
PRIMARY_THR = 0.05         # primary cut-off
MIN_LENS = [50, 100, 500, 1000]          # minimum number of deleted scorable bases
PRIMARY_MINLEN = 100
WIN = 100_000


def load_depth(path, sizes):
    """bedGraph (run-length) -> dict chrom -> uint16 depth array."""
    df = pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e", "d"],
                     dtype={"c": "category", "s": np.int32, "e": np.int32, "d": np.int32})
    out = {}
    for c in CHROMS:
        sub = df[df.c == c]
        arr = np.repeat(np.minimum(sub.d.values, 65535).astype(np.uint16), (sub.e - sub.s).values)
        assert arr.size == sizes[c], (path, c, arr.size, sizes[c])
        out[c] = arr
    return out


def runs(mask):
    """start, end (half-open) of True runs in a boolean array."""
    d = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def median_depth(depth, mask):
    return float(np.median(np.concatenate([depth[c][mask[c]] for c in CHROMS])))


def load_mask(bed, sizes, value=True):
    m = {c: np.full(sizes[c], not value, bool) for c in CHROMS}
    for l in open(bed):
        f = l.split("\t")
        if f[0] in m:
            m[f[0]][int(f[1]):int(f[2])] = value
    return m


def write_bed(mask, path):
    with open(path, "w") as fh:
        for c in CHROMS:
            st, en = runs(mask[c])
            for s, e in zip(st, en):
                fh.write(f"{c}\t{s}\t{e}\n")


def segments(deleted_c, min_present_run):
    """Deletion segments in compressed (scorable-only) coordinates.
    Returns arrays seg_start, seg_end of stretches not interrupted by a long present run."""
    present = ~deleted_c
    ps, pe = runs(present)
    breaker = present.copy()
    breaker[present] = np.repeat((pe - ps) >= min_present_run, pe - ps)   # short present blips do not break
    return runs(~breaker)


def call_sample(args):
    s, a = args
    sizes, scor = a["sizes"], a["scorable"]
    q20 = load_depth(f"{a['depth']}/{s}.{a['aligner']}.q20.bga.gz", sizes)
    q0 = load_depth(f"{a['depth']}/{s}.{a['aligner']}.q0.bga.gz", sizes)
    med = median_depth(q20, scor)
    calls, wins = [], []
    for thr in THR_FRACS:
        cut = thr * med
        for c in CHROMS:
            idx = np.flatnonzero(scor[c])                       # genome position of each scorable base
            dele = q20[c][idx] <= cut
            st, en = segments(dele, a["min_present_run"])
            csum = np.concatenate(([0], np.cumsum(dele, dtype=np.int64)))
            keep = (csum[en] - csum[st]) >= min(MIN_LENS)
            incall = np.zeros(idx.size, bool)
            for cs, ce in zip(st[keep], en[keep]):
                d = np.flatnonzero(dele[cs:ce])
                i_s, i_e = cs + d[0], cs + d[-1] + 1            # trimmed to deleted bases (compressed)
                ndel = int(csum[i_e] - csum[i_s])
                g_s, g_e = int(idx[i_s]), int(idx[i_e - 1]) + 1
                o_s = int(idx[i_s - 1]) + 1 if i_s > 0 else 0   # just after the last present base before
                o_e = int(idx[i_e]) if i_e < idx.size else sizes[c]
                calls.append((s, thr, c, g_s, g_e, o_s, o_e, g_e - g_s, ndel, i_e - i_s,
                              float(q0[c][g_s:g_e].mean()), float(q20[c][idx[i_s:i_e]].mean()), med))
                if ndel >= PRIMARY_MINLEN:
                    incall[i_s:i_e] = dele[i_s:i_e]
            if thr == PRIMARY_THR:
                # 100-kb window track: deleted scorable bases inside calls >= PRIMARY_MINLEN
                w = idx // WIN
                nwin = (sizes[c] + WIN - 1) // WIN
                n_sc = np.bincount(w, minlength=nwin); n_low = np.bincount(w, weights=dele, minlength=nwin)
                n_in = np.bincount(w, weights=incall, minlength=nwin)
                dsum = np.bincount(w, weights=q20[c][idx], minlength=nwin)
                for k in range(nwin):
                    wins.append((s, c, k * WIN, min((k + 1) * WIN, sizes[c]), int(n_sc[k]), int(n_low[k]), int(n_in[k]),
                                 dsum[k] / n_sc[k] if n_sc[k] else np.nan, med))
    cols = ["strain", "thr_frac", "chrom", "start", "end", "outer_start", "outer_end", "span",
            "deleted_bases", "scorable_bases", "mean_depth_allreads", "mean_depth_unique_scorable", "strain_median_depth"]
    wcols = ["strain", "chrom", "start", "end", "scorable_bases", "lowdepth_bases", "deleted_bases_in_calls",
             "mean_unique_depth_scorable", "strain_median_depth"]
    return pd.DataFrame(calls, columns=cols), pd.DataFrame(wins, columns=wcols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", default="depth")
    ap.add_argument("--ref", default="ref")
    ap.add_argument("--samples", default="meta/copy38.txt")
    ap.add_argument("--aligner", default="bwa", choices=["bwa", "bowtie2"])
    ap.add_argument("--self-sample", default="E2")
    ap.add_argument("--self-frac", type=float, default=0.25,
                    help="base is callable if E2's own unique depth >= this fraction of E2's median")
    ap.add_argument("--exclude-bed", nargs="*", default=[], help="BED files of bases not used as evidence (repeats)")
    ap.add_argument("--min-present-run", type=int, default=1,
                    help="a call is ended only by this many present scorable bases in a row")
    ap.add_argument("--tag", default=None, help="label for output files (default: aligner)")
    ap.add_argument("--out", default="calls")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tag = a.tag or a.aligner

    sizes = {l.split()[0]: int(l.split()[1]) for l in open(f"{a.ref}/E2.chr.genome")}

    # ---- callability from E2 self-mapping and assembly gaps
    self_q20 = load_depth(f"{a.depth}/{a.self_sample}.{a.aligner}.q20.bga.gz", sizes)
    nogap = load_mask(f"{a.ref}/E2.gaps.bed", sizes, value=False)
    self_med = median_depth(self_q20, nogap)
    callable_ = {c: nogap[c] & (self_q20[c] >= a.self_frac * self_med) for c in CHROMS}
    scorable = {c: callable_[c].copy() for c in CHROMS}
    for bed in a.exclude_bed:
        ex = load_mask(bed, sizes, value=True)
        for c in CHROMS:
            scorable[c] &= ~ex[c]
    write_bed(callable_, f"{a.out}/callable.{tag}.bed")
    write_bed(scorable, f"{a.out}/scorable.{tag}.bed")
    tot = sum(sizes.values())
    ncal = sum(int(v.sum()) for v in callable_.values()); nsc = sum(int(v.sum()) for v in scorable.values())
    with open(f"{a.out}/callable.{tag}.summary.txt", "w") as fh:
        fh.write(f"self_sample\t{a.self_sample}\nself_median_unique_depth\t{self_med}\nself_frac\t{a.self_frac}\n"
                 f"genome_bases_chr\t{tot}\ncallable_bases\t{ncal}\ncallable_fraction\t{ncal / tot:.5f}\n"
                 f"exclude_bed\t{','.join(a.exclude_bed)}\nscorable_bases\t{nsc}\nscorable_fraction\t{nsc / tot:.5f}\n"
                 f"min_present_run\t{a.min_present_run}\n")
    print(open(f"{a.out}/callable.{tag}.summary.txt").read(), flush=True)

    samples = [l.strip() for l in open(a.samples) if l.strip()]
    job = dict(sizes=sizes, scorable=scorable, depth=a.depth, aligner=a.aligner, min_present_run=a.min_present_run)
    with Pool(a.threads) as pool:
        res = pool.map(call_sample, [(s, job) for s in samples], chunksize=1)
    calls = pd.concat([r[0] for r in res if len(r[0])]); wins = pd.concat([r[1] for r in res])
    calls.to_csv(f"{a.out}/depth_calls.{tag}.all_thresholds.tsv.gz", sep="\t", index=False)
    wins.to_csv(f"{a.out}/windows100kb.{tag}.tsv.gz", sep="\t", index=False)
    prim = calls[(calls.thr_frac == PRIMARY_THR) & (calls.deleted_bases >= PRIMARY_MINLEN)]
    prim.to_csv(f"{a.out}/depth_calls.{tag}.primary.tsv", sep="\t", index=False)
    sens = (calls.assign(k=1).merge(pd.DataFrame({"min_len": MIN_LENS, "k": 1}), on="k")
            .query("deleted_bases >= min_len")
            .groupby(["thr_frac", "min_len"]).agg(n_calls=("span", "size"), deleted_bp=("deleted_bases", "sum"),
                                                  n_strains=("strain", "nunique")).reset_index())
    sens.to_csv(f"{a.out}/sensitivity.{tag}.tsv", sep="\t", index=False)
    print(sens.to_string(), flush=True)
    print(prim.groupby("strain").agg(n=("span", "size"), bp=("deleted_bases", "sum")).to_string())


if __name__ == "__main__":
    sys.exit(main())
