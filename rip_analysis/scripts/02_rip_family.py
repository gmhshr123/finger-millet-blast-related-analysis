#!/usr/bin/env python3
"""RIP-type mutation analysis of full-length TE copies, one family in one genome at a time.

Steps for each (genome, family) set of full-length copies:
  1. cut the copies out of the genome assembly (reverse-complemented where the table gives a minus strand)
  2. align them with MAFFT (--adjustdirection, so copies on the other strand are turned round)
  3. split the set into subfamilies. Some sets mix lineages that separated long ago (for example two
     Pot2-like elements that differ at 17% of positions). A consensus of such a mixture is not an
     ancestral element, and the differences between lineages would be counted as mutations. Copies are
     grouped by average-linkage clustering on the TRANSVERSION distance (purine <-> pyrimidine changes
     only) and the tree is cut at --split-distance. RIP only makes transitions, so a copy that was hit
     hard by RIP stays with its relatives. Every group with at least --min-sub copies is analysed as
     its own subfamily; the unsplit set is analysed too and reported as subfamily "all".
  4. majority-rule consensus = model of the intact ancestral element. A column in which most copies
     have a gap is treated as an insertion in a few copies and is not part of the consensus; a column
     with a tie for the most common base is left undetermined (N) and not used.
  5. every difference between a copy and the consensus is one substitution, written consensus>copy.
     RIP converts C to T; on the other strand this shows as G to A. So:
        RIP-type transitions   C>T and G>A
        other transitions      T>C and A>G
        transversions          the eight remaining changes
     Each substitution is counted per copy ("copy level") and once per alignment column and derived
     base ("event level"), because a change shared by several copies may have happened only once.
     "singleton" counts use only the changes seen in exactly one copy.
  6. sequence context: for C>T the base after the C in the consensus is recorded, for G>A the base
     before the G (CpA on one strand is TpG on the other, CpT is ApG, and so on).
  7. composition-based RIP indices of each copy, whole sequence and sliding windows (as in The RIPper):
        product index    TpA / ApT                    (RIP raises it; RIPper cut-off > 1.1)
        substrate index  (CpA + TpG) / (ApC + GpT)    (RIP lowers it; RIPper cut-off < 0.75)
        composite index  product - substrate          (RIPper cut-off > 0)

Usage: 02_rip_family.py --copies full_length_copies.tsv --genome E2=ref/E2.fasta [--genome B71=...] --out outdir
"""
import argparse, collections, os, subprocess, sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ap = argparse.ArgumentParser()
ap.add_argument("--copies", required=True)
ap.add_argument("--genome", action="append", required=True, help="NAME=path/to/assembly.fasta (repeat for each genome)")
ap.add_argument("--out", required=True)
ap.add_argument("--min-copies", type=int, default=9, help="smallest (genome, family) set that is aligned")
ap.add_argument("--min-sub", type=int, default=5, help="smallest subfamily that is analysed")
ap.add_argument("--split-distance", type=float, default=0.01, help="transversion distance at which the tree is cut")
ap.add_argument("--window", type=int, default=500)
ap.add_argument("--step", type=int, default=100)
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
for d in ("alignments", "consensus", "copies"):
    os.makedirs(f"{a.out}/{d}", exist_ok=True)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "ACGT"


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    seqs[name] = "".join(buf)
    return seqs


def rip_indices(seq):
    c = collections.Counter(seq[i:i + 2] for i in range(len(seq) - 1))
    prod = c["TA"] / c["AT"] if c["AT"] else float("nan")
    sub = (c["CA"] + c["TG"]) / (c["AC"] + c["GT"]) if (c["AC"] + c["GT"]) else float("nan")
    gc = (seq.count("G") + seq.count("C")) / max(1, sum(seq.count(b) for b in BASES))
    return prod, sub, prod - sub, gc


def transversion_distance(aln, names):
    """pairwise share of compared positions where one copy has a purine and the other a pyrimidine"""
    X = np.array([np.frombuffer(aln[k].encode(), dtype="S1") for k in names])
    pur = np.isin(X, [b"A", b"G"]).astype(np.float32); pyr = np.isin(X, [b"C", b"T"]).astype(np.float32)
    valid = pur + pyr
    mism = pur @ pyr.T + pyr @ pur.T
    comp = valid @ valid.T
    D = np.where(comp > 0, mism / np.maximum(comp, 1), 0.0)
    np.fill_diagonal(D, 0.0)
    return (D + D.T) / 2


fam_rows, copy_rows, event_rows, win_rows, ctx_rows, member_rows = [], [], [], [], [], []


def analyse(gname, fam, subfam, names, aln, te_class, extra):
    """steps 4-7 for one group of aligned copies"""
    tag = f"{gname}.{fam}.{subfam}"
    n = len(names); L = len(aln[names[0]])
    cons = []
    for j in range(L):
        col = collections.Counter(aln[k][j] for k in names)
        if col["-"] > n / 2:
            cons.append("-"); continue
        top = sorted(((col[b], b) for b in BASES), reverse=True)
        cons.append(top[0][1] if top[0][0] > top[1][0] and top[0][0] > 0 else "N")
    kept = [j for j in range(L) if cons[j] != "-"]
    cseq = "".join(cons[j] for j in kept)
    with open(f"{a.out}/consensus/{tag}.consensus.fa", "w") as fh:
        fh.write(f">{tag}_consensus n_copies={n} length={len(cseq)}\n{cseq}\n")
    nxt = {kept[i]: (cons[kept[i + 1]] if i + 1 < len(kept) else "N") for i in range(len(kept))}
    prv = {kept[i]: (cons[kept[i - 1]] if i > 0 else "N") for i in range(len(kept))}
    site = collections.Counter(cseq)
    avail = collections.Counter()
    for j in kept:
        if cons[j] == "C":
            avail["C>T:Cp" + nxt[j]] += 1
        elif cons[j] == "G":
            avail["G>A:" + prv[j] + "pG"] += 1
    events = collections.defaultdict(int); tot = collections.Counter(); mine = []
    for k in names:
        s = aln[k]; c = collections.Counter(); compared = 0
        for j in kept:
            x, y = cons[j], s[j]
            if x not in BASES or y not in BASES:
                continue
            compared += 1
            if x != y:
                t = f"{x}>{y}"
                c[t] += 1
                ctx = ("Cp" + nxt[j]) if t == "C>T" else ((prv[j] + "pG") if t == "G>A" else "")
                events[(j, t, ctx)] += 1
        raw = s.replace("-", "")
        p, sb, comp, gc = rip_indices(raw)
        nsub = sum(c.values()); rip = c["C>T"] + c["G>A"]; oti = c["T>C"] + c["A>G"]
        row = dict(genome=gname, family=fam, subfamily=subfam, copy_id=k, length=len(raw), compared_sites=compared,
                   substitutions=nsub, divergence=nsub / max(compared, 1), C_to_T=c["C>T"], G_to_A=c["G>A"],
                   T_to_C=c["T>C"], A_to_G=c["A>G"], transversions=nsub - rip - oti,
                   rip_type_share=rip / nsub if nsub else float("nan"), gc=gc,
                   product_index=p, substrate_index=sb, composite_index=comp)
        mine.append(row); copy_rows.append(row); tot.update(c)
        if subfam != "all" or extra.get("n_subfamilies", 0) == 0:       # windows once per copy
            for w in range(0, max(1, len(raw) - a.window + 1), a.step):
                seg = raw[w:w + a.window]
                wp, ws, wc, wg = rip_indices(seg)
                win_rows.append((gname, fam, subfam, k, w, w + len(seg), wp, ws, wc, wg))
    ev = collections.Counter(); ev_ctx = collections.Counter(); cp_ctx = collections.Counter(); singles = collections.Counter()
    for (j, t, ctx), m in events.items():
        ev[t] += 1
        if m == 1:
            singles[t] += 1
        if ctx:
            ev_ctx[t + ":" + ctx] += 1; cp_ctx[t + ":" + ctx] += m
        event_rows.append((gname, fam, subfam, j, t, ctx, m))
    for key in sorted(set(avail) | set(ev_ctx)):
        ctx_rows.append((gname, fam, subfam, key.split(":")[0], key.split(":")[1], avail[key], ev_ctx[key], cp_ctx[key]))

    def spectrum(c, pre):
        nsub = sum(c.values()); rip = c["C>T"] + c["G>A"]; oti = c["T>C"] + c["A>G"]; tv = nsub - rip - oti
        return {pre + "substitutions": nsub, pre + "C_to_T": c["C>T"], pre + "G_to_A": c["G>A"], pre + "T_to_C": c["T>C"],
                pre + "A_to_G": c["A>G"], pre + "transversions": tv,
                pre + "ti_tv": (rip + oti) / tv if tv else float("nan"),
                pre + "rip_type_share": rip / nsub if nsub else float("nan"),
                pre + "rip_vs_other_transitions": rip / oti if oti else float("nan")}
    cg, at = site["C"] + site["G"], site["A"] + site["T"]
    pc, sc, cc, gcc = rip_indices(cseq.replace("N", ""))
    row = dict(genome=gname, family=fam, subfamily=subfam, te_class=te_class, n_copies=n, **extra,
               consensus_length=len(cseq), consensus_undetermined=cseq.count("N"), consensus_gc=gcc,
               consensus_CG_sites=cg, consensus_AT_sites=at, consensus_product_index=pc,
               consensus_substrate_index=sc, consensus_composite_index=cc)
    row.update(spectrum(tot, "copy_")); row.update(spectrum(ev, "event_")); row.update(spectrum(singles, "singleton_"))
    row["event_rip_per_CG_site"] = (ev["C>T"] + ev["G>A"]) / cg if cg else float("nan")
    row["event_other_ti_per_AT_site"] = (ev["T>C"] + ev["A>G"]) / at if at else float("nan")
    cpa = ev_ctx["C>T:CpA"] + ev_ctx["G>A:TpG"]; cpt = ev_ctx["C>T:CpT"] + ev_ctx["G>A:ApG"]; allrip = ev["C>T"] + ev["G>A"]
    row["event_rip_in_CpA_share"] = cpa / allrip if allrip else float("nan")
    row["CpA_share_of_CG_sites"] = (avail["C>T:CpA"] + avail["G>A:TpG"]) / cg if cg else float("nan")
    row["event_rip_in_CpT_share"] = cpt / allrip if allrip else float("nan")
    row["CpT_share_of_CG_sites"] = (avail["C>T:CpT"] + avail["G>A:ApG"]) / cg if cg else float("nan")
    dv = pd.Series([r["divergence"] for r in mine]); ci = pd.Series([r["composite_index"] for r in mine])
    row["median_divergence"] = float(dv.median()); row["max_divergence"] = float(dv.max())
    row["median_composite_index"] = float(ci.median()); row["copies_composite_above_0"] = int((ci > 0).sum())
    fam_rows.append(row)
    print(f"{tag}: {n} copies, consensus {len(cseq)} bp, median divergence {dv.median():.4f}, {sum(ev.values())} events, "
          f"{row['event_rip_type_share']:.2f} RIP-type", file=sys.stderr, flush=True)


copies = pd.read_csv(a.copies, sep="\t", keep_default_na=False)
copies = copies[copies.selected.astype(str).str.lower() == "true"]
for gname, gpath in dict(g.split("=", 1) for g in a.genome).items():
    G = read_fasta(gpath)
    for fam, sub in copies[copies.genome == gname].groupby("family"):
        if len(sub) < a.min_copies:
            continue
        tag = f"{gname}.{fam}"
        # ---- 1. cut out the copies; the listed target-site duplication should flank each LTR element
        seqs, tsd_ok, tsd_n = {}, 0, 0
        for r in sub.itertuples():
            chrom = G[r.chrom]
            assert r.end <= len(chrom), (r.copy_id, len(chrom))
            s = chrom[r.start - 1:r.end].upper()
            seqs[r.copy_id] = s.translate(COMP)[::-1] if r.strand == "-" else s
            if r.tsd:
                tsd_n += 1; k = len(r.tsd)
                left, right = chrom[r.start - 1 - k:r.start - 1].upper(), chrom[r.end:r.end + k].upper()
                tsd_ok += (left == right and left in {r.tsd.upper(), r.tsd.upper().translate(COMP)[::-1]})
        fa = f"{a.out}/copies/{tag}.fa"
        with open(fa, "w") as fh:
            for k, v in seqs.items():
                fh.write(f">{k}\n{v}\n")
        # ---- 2. align
        aln_path = f"{a.out}/alignments/{tag}.aln.fa"
        with open(aln_path, "w") as fh:
            subprocess.run(["mafft", "--auto", "--adjustdirection", "--thread", str(a.threads), "--quiet", fa], stdout=fh, check=True)
        raw_aln = read_fasta(aln_path)
        aln = {(k[3:] if k.startswith("_R_") else k): v.upper() for k, v in raw_aln.items()}
        names = list(aln)
        # ---- 3. subfamilies from transversion distance
        D = transversion_distance(aln, names)
        cl = fcluster(linkage(squareform(D, checks=False), "average"), a.split_distance, "distance")
        sizes = collections.Counter(cl)
        order = [c for c, m in sizes.most_common() if m >= a.min_sub]
        label = {c: f"sub{i + 1}" for i, c in enumerate(order)}
        for k, c in zip(names, cl):
            member_rows.append((gname, fam, k, label.get(c, "unassigned")))
        extra = dict(n_set=len(names), copies_turned_round=sum(k.startswith("_R_") for k in raw_aln),
                     tsd_confirmed=f"{tsd_ok}/{tsd_n}" if tsd_n else "", alignment_length=len(aln[names[0]]),
                     n_subfamilies=len(order), copies_unassigned=sum(m for c, m in sizes.items() if c not in label),
                     max_transversion_distance=float(D.max()))
        te_class = sub.te_class.iloc[0]
        analyse(gname, fam, "all", names, aln, te_class, extra)
        if not (len(order) == 1 and sizes[order[0]] == len(names)):       # a split changes something
            for c in order:
                analyse(gname, fam, label[c], [k for k, x in zip(names, cl) if x == c], aln, te_class, extra)

pd.DataFrame(fam_rows).to_csv(f"{a.out}/family_summary.tsv", sep="\t", index=False, float_format="%.4f")
pd.DataFrame(copy_rows).to_csv(f"{a.out}/copy_table.tsv", sep="\t", index=False, float_format="%.4f")
pd.DataFrame(member_rows, columns=["genome", "family", "copy_id", "subfamily"]).to_csv(f"{a.out}/subfamily_membership.tsv", sep="\t", index=False)
pd.DataFrame(event_rows, columns=["genome", "family", "subfamily", "alignment_column", "substitution", "context", "n_copies"]).to_csv(
    f"{a.out}/substitution_events.tsv.gz", sep="\t", index=False)
pd.DataFrame(ctx_rows, columns=["genome", "family", "subfamily", "substitution", "context", "consensus_sites", "events", "copy_level"]).to_csv(
    f"{a.out}/context_table.tsv", sep="\t", index=False)
pd.DataFrame(win_rows, columns=["genome", "family", "subfamily", "copy_id", "start", "end", "product_index", "substrate_index",
                                "composite_index", "gc"]).to_csv(f"{a.out}/sliding_windows.tsv.gz", sep="\t", index=False, float_format="%.4f")
