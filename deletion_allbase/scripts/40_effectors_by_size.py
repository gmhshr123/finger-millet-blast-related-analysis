#!/usr/bin/env python3
"""Size of the 191 one-TE-family deletions and the effector share of their deleted genes by deletion size and by the
number of strains sharing the deletion; overlap of their deleted genes with the Data S10 genes. Run from manuscripts_all/.

Same rules as 38_whole_gene_stats.py (copied here so that its outputs are not rewritten): event sets from
regions_for_browser.tsv, whole deleted genes (100% of the uniquely mappable bases), effectors from Data S2.
Event size = end - start of the event on E2 (allbase_events_read_evidence.tsv), as in Fig. 3 of script 35.
A gene deleted by several events is placed once, in the bin of its largest (or most widely shared) event, so the bins add up to the total.
Genes with an mRNA under 200 bp are left out of all counts; the genes kept are those of Data S2, which holds only genes of 200 bp or more.
101 genes have reversed exon coordinates in E2.gff3.gz and a false 4-bp span in E2.genes.bed; their span is taken here from all their exon coordinates.
Output: printed, deletion_allbase/supplement_figures_and_NoteS5/effectors_by_size_200bp.xlsx
"""
import collections, gzip, os, re, sys
from math import comb
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from depth_tracks import deletion_extent

R, I = "deletion_allbase/results_edta/", "deletion_allbase/inputs/"
OUT = "deletion_allbase/supplement_figures_and_NoteS5/effectors_by_size_200bp.xlsx"
cohort = [l.strip() for l in open(I + "cohort33.txt")]
TE = "transposable element (TE)"

# ---- the 191 events (as in script 38)
g = pd.read_csv(R + "regions_for_browser.tsv", sep="\t", keep_default_na=False)
B = g[(g["Left origin"] == TE) | (g["Right origin"] == TE)].copy()


def fams(s, side):
    out = set()
    for part in s.split(";"):
        part = part.strip()
        if part and part.split(":", 1)[0].strip() == side:
            out |= {x.strip() for x in part.split(":", 1)[1].split(",")}
    return out


B["L"] = B["TE at a breakpoint, from reads"].map(lambda s: fams(s, "left"))
B["R"] = B["TE at a breakpoint, from reads"].map(lambda s: fams(s, "right"))
B["dif"] = [bool(l and r and not (l & r)) for l, r in zip(B.L, B.R)]
ev = B.groupby("Deletion event ID").agg(L=("L", lambda x: set().union(*x)), R=("R", lambda x: set().union(*x)), dif=("dif", "any"))
ev["nfam"] = [len(l | r) for l, r in zip(ev.L, ev.R)]
IDS = set(ev[~ev.dif & (ev.nfam == 1)].index)
print("events:", len(IDS))

E = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t").set_index("event_id").loc[sorted(IDS)]
E["size"] = E.end - E.start

# ---- genes, effectors, deleted bases (as in script 38)
genes, ALLSEQ = [], set()
for l in open(R + "cluster_uniq/E2.genes.bed", "rb").read().decode().split("\n"):
    f = l.replace("\r", "").split("\t")
    if len(f) >= 4:
        ALLSEQ.add(f[3])
    if len(f) >= 4 and f[0].startswith("Chr"):
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))
G = pd.DataFrame(genes, columns=["chrom", "a", "b", "gene"]).drop_duplicates("gene")
ex = collections.defaultdict(list)
for l in gzip.open(I + "E2.gff3.gz", "rt"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 9 and f[2] == "exon":
        par = re.search(r"Parent=([^;,]+)", f[8]).group(1)
        if par.endswith("-mRNA-1"):
            ex[par[:-7]] += [int(f[3]), int(f[4])]
bad = G.b - G.a <= 4
G.loc[bad, "a"] = [min(ex[x]) - 1 for x in G.gene[bad]]
G.loc[bad, "b"] = [max(ex[x]) for x in G.gene[bad]]
print("gene spans corrected:", int(bad.sum()))
s2 = pd.read_csv(I + "DataS2_gene_sections.tsv", sep="\t")
S2G = set(s2.gene[s2.gene.astype(str).str.startswith("EcMO")])
SHORT = ALLSEQ - S2G
print("genes in the annotation:", len(G), "; not in Data S2 (under 200 bp), left out:", len(set(G.gene) & SHORT))
G = G[G.gene.isin(S2G)]
EFF = set(s2[s2.section.str.startswith("Effector", na=False)].gene) & set(G.gene)
# expected share over all sequences of the assembly, the 6.6% of the manuscript
S2EFF = set(s2[s2.section.str.startswith("Effector", na=False)].gene) & ALLSEQ
K, N = len(S2EFF), len(S2G)
print(f"expected effector share: {K}/{N} = {K / N * 100:.2f}%; seven chromosomes {len(EFF)}/{len(G)}")
sizes = {l.split()[0]: int(l.split()[1]) for l in open(R + "cluster_uniq/E2.chr.genome")}
scor = {c: np.zeros(n, bool) for c, n in sizes.items()}
for l in open(R + "cluster_uniq/scorable.bwa.uniq.bed"):
    c, a, b = l.split("\t")[:3]; scor[c][int(a):int(b)] = True
frag = pd.read_csv(R + "allbase_fragments_strain_level.tsv", sep="\t")
dele = {}
for (st, c), x in frag.groupby(["strain", "chrom"]):
    m = np.zeros(sizes[c], bool)
    for a, b in zip(x.start.values, x.end.values):
        m[a:b] = True
    dele[(st, c)] = m
Gc = {c: x for c, x in G.groupby("chrom")}
calls = pd.read_csv(R + "allbase_calls_strain_level.tsv", sep="\t")

EXT = {}                                           # event -> extent of each strain's region (bp), see whole_pairs


def whole_pairs(ids, skip=frozenset()):
    """(event, gene): whole gene deleted by this event in at least one strain; skip: (strain, event) regions left out"""
    pairs = set()
    sub = calls[calls.event_id.isin(ids) & calls.strain.isin(cohort)]
    sub = sub[[(s, e) not in skip for s, e in zip(sub.strain, sub.event_id)]]
    for (st, c), x in sub.groupby(["strain", "chrom"]):
        ext = [(s, e, eid) for s, e, eid, _, _ in deletion_extent(st, c, 0, sizes[c]) if eid in set(x.event_id)]
        for s, e, eid in ext:
            EXT.setdefault(eid, []).append(e - s)
        m = np.zeros(sizes[c], bool)
        for s, e, _ in ext:
            m[s:e] = True
        m &= dele[(st, c)] | ~scor[c]
        gx = Gc.get(c)
        for s, e, eid in ext:
            for a, b, gene in zip(gx.a, gx.b, gx.gene):
                if b <= s or a >= e:
                    continue
                sc = scor[c][a:b]; n = int(sc.sum())
                share = (m[a:b] & sc).sum() / n if n else (1.0 if (a >= s and b <= e) else 0.0)
                if share >= 1 - 1e-9:
                    pairs.add((eid, gene))
    return pairs


P = pd.DataFrame(sorted(whole_pairs(IDS)), columns=["event", "gene"])
# size between the read-defined breakpoints where they exist (the extent used for the genes), median over the strains;
# 'span' (first to last deleted uniquely mappable base) misses deleted repeats at the edges
E["span"] = E["size"]
E["size"] = [int(np.median(EXT[i])) for i in E.index]
P["eff"] = P.gene.isin(EFF)
P = P.join(E[["size", "span", "n_strains"]], on="event")
allg = set(P.gene)
print(f"unique whole-deleted genes {len(allg)}, effectors {len(allg & EFF)}")


def hyper(k, n):
    return sum(comb(K, i) * comb(N - K, n - i) for i in range(k, min(n, K) + 1)) / comb(N, n) if n else float("nan")


def table(col, bins, labels, name):
    E["bin"] = pd.cut(E[col], bins, labels=labels, right=False)
    P["bin"] = pd.cut(P[col], bins, labels=labels, right=False)
    one = P.loc[P.groupby("gene")[col].idxmax()]             # each gene once, in the bin of its largest value of col
    rows = []
    for lab in labels:
        e, p = E[E.bin == lab], P[P.bin == lab]
        gs = set(one[one.bin == lab].gene); k = len(gs & EFF)
        rows.append({name: lab, "Deletions": len(e), "Deletions carrying a whole gene": int(p.event.nunique()),
                     "Unique deleted genes": len(gs), "Effectors": k,
                     "Effector share (%)": round(k / len(gs) * 100, 1) if gs else None,
                     "Enrichment P (vs 6.6%)": float(f"{hyper(k, len(gs)):.2g}") if gs else None})
    return pd.DataFrame(rows)


SZ = table("size", [0, 1_000, 5_000, 20_000, 50_000, 10**9], ["< 1 kb", "1-5 kb", "5-20 kb", "20-50 kb", "> 50 kb"], "Deletion size")
SH = table("n_strains", [1, 2, 6, 16, 100], ["1 strain", "2-5 strains", "6-15 strains", "16-33 strains"], "Strains sharing the deletion")

# effector share against size, per deletion that carries genes (Spearman)
from scipy.stats import spearmanr, mannwhitneyu
pe = P.groupby("event").agg(size=("size", "first"), genes=("gene", "nunique"), effs=("eff", "sum"))
pe["share"] = pe.effs / pe.genes
rho = spearmanr(pe["size"], pe.share)
# compare sizes of deletions with and without an effector
mw = mannwhitneyu(pe[pe.effs > 0]["size"], pe[pe.effs == 0]["size"])

# overlap with the Data S10 genes
s10 = pd.read_csv(I + "DataS10_events.tsv", sep="\t")
S10G = set(",".join(s10.genes.dropna()).split(",")) - SHORT
OV = pd.DataFrame([
    {"Genes": "Deleted in the 191 deletions, also listed in Data S10", "Genes ": len(allg & S10G), "Effectors": len(allg & S10G & EFF)},
    {"Genes": "Deleted in the 191 deletions, not in Data S10", "Genes ": len(allg - S10G), "Effectors": len((allg - S10G) & EFF)},
    {"Genes": "Listed in Data S10 (all 124 deletions)", "Genes ": len(S10G), "Effectors": len(S10G & EFF)},
    {"Genes": "Listed in Data S10, not deleted in the 191 deletions", "Genes ": len(S10G - allg), "Effectors": len((S10G - allg) & EFF)}])
OV["Effector share (%)"] = (OV.Effectors / OV["Genes "] * 100).round(1)

# unique deleted genes, all deletions vs TE-associated, for the SNP-based (Data S10) and the all-base deletions
S10TE = set(",".join(s10[s10.te != "No_TE"].genes.dropna()).split(",")) - SHORT
ALLEV = set(calls.event_id.dropna())
SET245 = set(ev.index)
gen = {k: {gn for _, gn in whole_pairs(v)} for k, v in (("all", ALLEV), ("245", SET245))}
U = []
for lab, gs, n in [("SNP-based, all 124 Data S10 deletions (genes listed)", S10G, len(s10)),
                   (f"SNP-based, {int((s10.te != 'No_TE').sum())} TE-associated Data S10 deletions (genes listed)", S10TE, int((s10.te != "No_TE").sum())),
                   (f"All-base, all {len(ALLEV)} deletions (whole genes)", gen["all"], len(ALLEV)),
                   ("All-base, 245 deletions with a TE at a breakpoint (whole genes)", gen["245"], len(SET245)),
                   ("All-base, 191 deletions with one TE family (whole genes)", allg, len(IDS))]:
    k = len(gs & EFF)
    U.append({"Deletion set": lab, "Deletions": n, "Unique deleted genes": len(gs), "Effectors": k,
              "Effector share (%)": round(k / len(gs) * 100, 1), "Enrichment P (vs 6.6%)": float(f"{hyper(k, len(gs)):.2g}")})
U = pd.DataFrame(U)
X = pd.DataFrame([{"Data S10 TE-associated genes": len(S10TE), "effectors": len(S10TE & EFF),
                   "also deleted in the 191": len(S10TE & allg), "effectors also deleted in the 191": len(S10TE & allg & EFF),
                   "also deleted in the 245": len(S10TE & gen["245"]), "effectors also deleted in the 245": len(S10TE & gen["245"] & EFF),
                   "also deleted in any all-base deletion": len(S10TE & gen["all"]), "effectors, any": len(S10TE & gen["all"] & EFF)}])

SIZE = pd.DataFrame([{"Unit": u, "N": len(E), "Mean (bp)": round(E[k].mean()), "Median (bp)": int(E[k].median()),
                      "25th percentile": int(E[k].quantile(.25)), "75th percentile": int(E[k].quantile(.75)),
                      "Min": int(E[k].min()), "Max": int(E[k].max())}
                     for u, k in (("191 deletions, between the read-defined breakpoints (median over strains)", "size"),
                                  ("191 deletions, first to last deleted uniquely mappable base (event span, as in Fig. 3)", "span"))])
TESTS = pd.DataFrame([{"Test": "Spearman: effector share of a deletion's genes vs deletion size (deletions carrying a gene)",
                       "N": len(pe), "Statistic": round(rho.statistic, 3), "P": float(f"{rho.pvalue:.2g}")},
                      {"Test": "Mann-Whitney: size of deletions with vs without an effector (deletions carrying a gene)",
                       "N": len(pe), "Statistic": float(mw.statistic), "P": float(f"{mw.pvalue:.2g}"),
                       "Median size with effector": int(pe[pe.effs > 0]["size"].median()),
                       "Median size without effector": int(pe[pe.effs == 0]["size"].median())}])
# comparison table with the SNP-based deletions: real deletions only ('DNA still present' regions left out), seven TEs of the manuscript
# (Family13 = MZ5.TE00000001, Family1 = MZ5.TE00000027: tentative, from the Data S10 rows)
SEVEN = {"Pot2": "Pot2", "MAGGY/Fosbury": "Fosbury", "Grasshopper": "Grh", "Pyret": "Pyret", "MZ5.TE00000001": "Family13",
         "MZ5.TE00000027": "Family1", "Inago2": "Inago2"}
bl = pd.read_csv(R + "breakpoints_long.tsv", sep="\t", keep_default_na=False)
dna = g[g["Read evidence"] == "DNA still present"]
SKIP = frozenset(zip(dna.Strain, dna["Deletion event ID"]))
real = g[g["Read evidence"] != "DNA still present"]
te = bl[(bl.family != "") & np.array([(s, e) not in SKIP for s, e in zip(bl.strain, bl.event)])]
te7 = te[te.family.isin(SEVEN)]
nreg = real.groupby("Deletion event ID").Strain.nunique()
st = te[te.reads >= 3].groupby(["event", "family"]).strain.nunique().reset_index(name="k")
st = st[st.k / st.event.map(nreg) >= 0.5]
CMP = {"Any TE, any support": set(te.event), "Seven TEs, any support": set(te7.event),
       "Any TE, >=3 reads in >=50% of strains (Data S10 rule)": set(st.event),
       "Seven TEs, >=3 reads in >=50% of strains (Data S10 rule)": set(st[st.family.isin(SEVEN)].event),
       "One TE family (191), real deletions": IDS - {e for s, e in SKIP if not set(real[real["Deletion event ID"] == e].Strain) & set(te[te.event == e].strain)}}
C = []
for lab, ids in CMP.items():
    pp = whole_pairs(ids, SKIP); gs = {x for _, x in pp}; k = len(gs & EFF)
    C.append({"Deletion set": lab, "Deletions": len(ids), "Carrying >=1 whole gene": len({e for e, _ in pp}),
              "Unique deleted genes": len(gs), "Effectors": k, "Effector share (%)": round(k / len(gs) * 100, 1)})
C = pd.DataFrame(C)
F7 = pd.DataFrame([{"TE": SEVEN[f], "Deletions, any support": te7[te7.family == f].event.nunique(),
                    "Deletions, Data S10 rule": st[st.family == f].event.nunique()} for f in SEVEN])
# Note S5 of the supplement: size groups and sizes for the 187 seven-TE deletions (each gene once, in its largest deletion)
EV = pd.read_csv(R + "allbase_events_read_evidence.tsv", sep="\t").set_index("event_id")
EXT.clear()
pp7 = whole_pairs(CMP["Seven TEs, any support"], SKIP)
S7 = pd.DataFrame({"size": {k: int(np.median(v)) for k, v in EXT.items()}})
S7["span"] = (EV.end - EV.start).reindex(S7.index)
P7 = pd.DataFrame(sorted(pp7), columns=["event", "gene"]).join(S7, on="event")
P7["bin"] = pd.cut(P7["size"], [0, 1_000, 5_000, 20_000, 50_000, 10**9], labels=["< 1 kb", "1-5 kb", "5-20 kb", "20-50 kb", "> 50 kb"], right=False)
one7 = P7.loc[P7.groupby("gene")["size"].idxmax()]
SZ7 = one7.groupby("bin", observed=False).agg(genes=("gene", "nunique"), effectors=("gene", lambda s: len(set(s) & EFF)))
SZ7["share_%"] = (SZ7.effectors / SZ7.genes * 100).round(1)
SZ7["Enrichment P (vs 6.6%)"] = [float(f"{hyper(k, n):.2g}") for k, n in zip(SZ7.effectors, SZ7.genes)]
SZ7["deletions"] = pd.cut(S7["size"], [0, 1_000, 5_000, 20_000, 50_000, 10**9], labels=SZ7.index, right=False).value_counts().reindex(SZ7.index)
g7 = S7.groupby(pd.cut(S7["size"], [0, 1_000, 5_000, 20_000, 50_000, 10**9], labels=SZ7.index, right=False), observed=False)["size"]
SZ7["Mean size (kb)"] = (g7.mean() / 1000).round(1)
SZ7["Median size (kb)"] = (g7.median() / 1000).round(1)
SIZE7 = pd.DataFrame([{"Set": lab, "N": len(x), "Mean between breakpoints": round(x["size"].mean()), "Median between breakpoints": int(x["size"].median()),
                       "Mean depth span": round(x.span.mean()), "Median depth span": int(x.span.median())}
                      for lab, x in (("all 187", S7), ("under 50 kb", S7[S7["size"] < 50_000]), ("1-50 kb", S7[(S7["size"] >= 1_000) & (S7["size"] < 50_000)]))])
nf = te7.groupby("event").family.nunique()
sf = te7.groupby(["strain", "event"]).family.nunique()
print("seven-TE events with 2+ families:", int((nf > 1).sum()), "| of these with two families in one strain:", len({e for (_, e), v in sf.items() if v > 1}), "| family counts sum:", int(te7.groupby("family").event.nunique().sum()))
print(SZ7.to_string()); print(SIZE7.to_string(index=False)); print()
print("events with at least one real deletion region:", real["Deletion event ID"].nunique(), "of", g["Deletion event ID"].nunique())
for t in (C, F7):
    print(t.to_string(index=False)); print()
for t in (SIZE, SZ, SH, TESTS, OV, U, X.T):
    print(t.to_string(index=X.T is t)); print()
with pd.ExcelWriter(OUT) as w:
    for name, t in (("Size", SIZE), ("By_size", SZ), ("By_sharing", SH), ("Tests", TESTS), ("DataS10_overlap", OV),
                    ("Unique_genes_by_set", U), ("DataS10_TE_genes", X), ("Comparison_table", C), ("Seven_TEs", F7),
                    ("Size_groups_187", SZ7.reset_index()), ("Sizes_187", SIZE7)):
        t.to_excel(w, sheet_name=name, index=False)
    P.drop(columns="bin").to_excel(w, sheet_name="Event_gene_pairs", index=False)
print("written", OUT)
