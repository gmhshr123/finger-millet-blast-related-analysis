#!/usr/bin/env python3
"""Count all 16 dinucleotides in sliding windows along each genome.

The RIPper's indices use only six dinucleotides (TpA, ApT, CpA, TpG, ApC, GpT). With all 16 counts per window, any other
index can be calculated afterwards, for example one built on CpT -> TpT instead of CpA -> TpA (see 12_alternative_indices.py).

Sequence is converted to uppercase first (the NCBI assemblies are soft-masked). Windows are placed as in The RIPper:
full-length windows only, start 0, step --step. A dinucleotide is counted in a window when both bases are inside it.
Overlapping occurrences are counted (TTT holds two TpT).

Usage: 11_dinucleotide_windows.py --genome NAME=fasta [--genome ...] --out outdir [--window 1000 --step 500] [--ripcal-scan]
  --ripcal-scan  also write ripcal_scan_summary.tsv: the share of windows that pass RIPCaL's default index-scan thresholds
                 (TpA/ApT >= 0.89, or (CpA+TpG)/(ApC+GpT) <= 1.03). Meant for --window 300 --step 150, RIPCaL's defaults.
"""
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genome", action="append", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--window", type=int, default=1000)
ap.add_argument("--step", type=int, default=500)
ap.add_argument("--ripcal-scan", action="store_true")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
CODE = np.full(256, 4, np.int8)
for i, b in enumerate("ACGT"):
    CODE[ord(b)] = i; CODE[ord(b.lower())] = i
DI = [x + y for x in "ACGT" for y in "ACGT"]


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


scan = []
for spec in a.genome:
    g, path = spec.split("=", 1)
    W = []
    for chrom, seq in read_fasta(path):
        raw = np.frombuffer(seq.encode(), np.uint8)
        x = CODE[raw]
        L = x.size
        if L < a.window:
            continue
        di = x[:-1].astype(np.int16) * 5 + x[1:]                     # dinucleotide starting at each base
        st = np.arange(0, L - a.window + 1, a.step); en = st + a.window
        d = {"genome": g, "chrom": chrom, "start": st, "end": en}
        acgt = np.concatenate(([0], np.cumsum(x < 4))); gc = np.concatenate(([0], np.cumsum((x == 1) | (x == 2))))
        d["acgt"] = acgt[en] - acgt[st]; d["gc"] = gc[en] - gc[st]
        low = np.concatenate(([0], np.cumsum(raw >= 97)))            # lowercase = soft-masked repeat in the NCBI assemblies
        d["soft_masked"] = low[en] - low[st]
        for k in DI:
            cs = np.concatenate(([0], np.cumsum(di == "ACGT".index(k[0]) * 5 + "ACGT".index(k[1]))))
            d[k] = cs[en - 1] - cs[st]                                # dinucleotides fully inside the window
        W.append(pd.DataFrame(d))
    W = pd.concat(W, ignore_index=True)
    W.to_csv(f"{a.out}/{g}.dinuc_w{a.window}.tsv.gz", sep="\t", index=False)
    print(g, len(W), "windows", flush=True)
    if a.ripcal_scan:
        with np.errstate(divide="ignore", invalid="ignore"):
            p = W.TA / W.AT; s = (W.CA + W.TG) / (W.AC + W.GT)
        f1 = (p >= 0.89) & (W.AT > 0); f2 = (s <= 1.03) & ((W.AC + W.GT) > 0)
        scan.append(dict(genome=g, window=a.window, step=a.step, windows=len(W), pass_TpA_ApT=int(f1.sum()), pass_CpA_TpG=int(f2.sum()),
                         pass_either=int((f1 | f2).sum()), pass_both=int((f1 & f2).sum()), share_either=float((f1 | f2).mean()),
                         share_both=float((f1 & f2).mean())))
if scan:
    pd.DataFrame(scan).to_csv(f"{a.out}/ripcal_scan_summary.tsv", sep="\t", index=False, float_format="%.5f")
    print(pd.DataFrame(scan).to_string())
