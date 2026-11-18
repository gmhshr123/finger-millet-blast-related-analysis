#!/usr/bin/env python3
"""Genome-wide RIP profile with the definitions and default settings of The RIPper (van Wyk et al. 2019, PeerJ 7:e7447).

Settings taken from the tool's own documentation (https://theripper.hawk.rocks, "Background"):
  sliding window 1000 bp, step 500 bp
  product index    TpA / ApT                       RIP-affected if > 1.1
  substrate index  (CpA + TpG) / (ApC + GpT)       RIP-affected if <= 0.75
  composite index  product - substrate             RIP-affected if > 0
  a window is RIP-positive when all three hold
  LRAR (large RIP-affected region): 7 or more consecutive RIP-positive windows (4000 bp or more)
This is a re-implementation of those definitions, not the RIPper program itself.

Usage: 05_ripper_genomewide.py --genome NAME=fasta [--genome ...] --out outdir [--bed NAME=te.bed]
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genome", action="append", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--window", type=int, default=1000)
ap.add_argument("--step", type=int, default=500)
ap.add_argument("--chain", type=int, default=7)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
CODE = np.full(256, 4, np.int8)
for i, b in enumerate("ACGT"):
    CODE[ord(b)] = i; CODE[ord(b.lower())] = i
DI = {"TA": 3 * 5 + 0, "AT": 0 * 5 + 3, "CA": 1 * 5 + 0, "TG": 3 * 5 + 2, "AC": 0 * 5 + 1, "GT": 2 * 5 + 3}


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    yield name, "".join(buf)


summ, lrars = [], []
for spec in a.genome:
    g, path = spec.split("=", 1)
    W = []
    for chrom, seq in read_fasta(path):
        x = CODE[np.frombuffer(seq.encode(), np.uint8)]
        L = x.size
        if L < a.window:
            continue
        di = x[:-1].astype(np.int16) * 5 + x[1:]                     # dinucleotide starting at each base
        cs = {k: np.concatenate(([0], np.cumsum(di == v))) for k, v in DI.items()}
        gc = np.concatenate(([0], np.cumsum((x == 1) | (x == 2)))); acgt = np.concatenate(([0], np.cumsum(x < 4)))
        st = np.arange(0, L - a.window + 1, a.step); en = st + a.window
        n = {k: v[en - 1] - v[st] for k, v in cs.items()}             # dinucleotides fully inside the window
        with np.errstate(divide="ignore", invalid="ignore"):
            prod = n["TA"] / n["AT"]; sub = (n["CA"] + n["TG"]) / (n["AC"] + n["GT"])
        comp = prod - sub
        pos = (prod > 1.1) & (sub <= 0.75) & (comp > 0)
        pos &= np.isfinite(prod) & np.isfinite(sub)
        d = pd.DataFrame({"genome": g, "chrom": chrom, "start": st, "end": en, "product": prod, "substrate": sub, "composite": comp,
                          "gc": (gc[en] - gc[st]) / np.maximum(acgt[en] - acgt[st], 1), "acgt": acgt[en] - acgt[st], "rip": pos})
        W.append(d)
        # LRARs: runs of >= chain consecutive RIP-positive windows
        e = np.diff(np.concatenate(([0], pos.astype(np.int8), [0])))
        for i, j in zip(np.flatnonzero(e == 1), np.flatnonzero(e == -1)):
            if j - i >= a.chain:
                r = d.iloc[i:j]
                lrars.append((g, chrom, int(st[i]), int(en[j - 1]), j - i, int(en[j - 1] - st[i]), r["product"].mean(), r.substrate.mean(),
                              r.composite.mean(), r.gc.mean()))
    W = pd.concat(W, ignore_index=True)
    W.to_csv(f"{a.out}/{g}.windows.tsv.gz", sep="\t", index=False, float_format="%.4f")
    ok = W[W.acgt >= 0.9 * a.window]                                   # windows that are mostly real sequence, not N
    lg = [x for x in lrars if x[0] == g]
    # bases covered by at least one RIP-positive window
    covered = 0
    for c, d in W[W.rip].groupby("chrom"):
        s, e = d.start.values, d.end.values
        covered += int((e - s).sum() - np.maximum(e[:-1] - s[1:], 0).sum())
    total_bp = int(W.groupby("chrom").end.max().sum())
    summ.append(dict(genome=g, sequences=W.chrom.nunique(), genome_bp_in_windows=total_bp, windows=len(W),
                     rip_positive_windows=int(W.rip.sum()), rip_positive_share_of_windows=W.rip.mean(),
                     bp_in_rip_positive_windows=covered, rip_share_of_genome_bp=covered / total_bp,
                     windows_composite_above_0=float((ok.composite > 0).mean()), windows_product_above_1_1=float((ok["product"] > 1.1).mean()),
                     windows_substrate_le_0_75=float((ok.substrate <= 0.75).mean()),
                     median_product=ok["product"].median(), median_substrate=ok.substrate.median(), median_composite=ok.composite.median(),
                     genome_gc=float((ok.gc * ok.acgt).sum() / ok.acgt.sum()), gc_of_rip_windows=float(W[W.rip].gc.mean()) if W.rip.any() else np.nan,
                     lrar_count=len(lg), lrar_total_bp=sum(x[5] for x in lg), lrar_mean_bp=float(np.mean([x[5] for x in lg])) if lg else 0,
                     lrar_max_bp=max([x[5] for x in lg]) if lg else 0, lrar_mean_gc=float(np.mean([x[9] for x in lg])) if lg else np.nan))
S = pd.DataFrame(summ)
S.to_csv(f"{a.out}/ripper_genome_profile.tsv", sep="\t", index=False, float_format="%.5f")
pd.DataFrame(lrars, columns=["genome", "chrom", "start", "end", "windows", "length", "mean_product", "mean_substrate", "mean_composite",
                             "mean_gc"]).to_csv(f"{a.out}/ripper_LRARs.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print(S.set_index("genome").T.to_string())
