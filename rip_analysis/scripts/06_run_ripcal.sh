#!/bin/bash
# Run RIPCaL 2.0 (Hane & Oliver 2008; alignment-based mode) on every group of aligned full-length TE copies,
# once against RIPCaL's degenerate consensus and once against the copy with the highest G+C content.
# Run from /path/to/scratch/finger_millet_blast_phylogeny/rip_analysis.
#
# Installation notes
#   source   ripcal/ripcal2_install.zip from https://sourceforge.net/projects/ripcal/ (RIPCAL_2.0)
#   Perl     conda environment ripcal/env: perl-bioperl-core, perl-bio-procedural, perl-gd, perl-math-round
#   Tk       RIPCaL loads the Tk graphical toolkit even in command-line mode. Tk is not installed; ripcal/stub
#            holds empty Tk modules so that the program starts. The RIPCaL script itself is unchanged.
set -euo pipefail
B=$PWD/ripcal
PERL="$B/env/bin/perl -I$B/stub"
RES=${RES:-results_4genomes}
module load SciPy-bundle/2024.05-gfbf-2024a

# 1. one alignment file per group (whole set, and each subfamily), without columns that are gaps in every copy
python - "$RES" "$B/groups" <<'EOF'
import os, sys
import pandas as pd
res, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
mem = pd.read_csv(f"{res}/subfamily_membership.tsv", sep="\t")
fam = pd.read_csv(f"{res}/family_summary.tsv", sep="\t")
def read_fasta(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; n = n[3:] if n.startswith("_R_") else n; d[n] = []
        else:
            d[n].append(l.strip().upper())
    return {k: "".join(v) for k, v in d.items()}
for r in fam.itertuples():
    aln = read_fasta(f"{res}/alignments/{r.genome}.{r.family}.aln.fa")
    m = mem[(mem.genome == r.genome) & (mem.family == r.family)]
    ids = list(m.copy_id) if r.subfamily == "all" else list(m[m.subfamily == r.subfamily].copy_id)
    keep = [j for j in range(len(aln[ids[0]])) if any(aln[k][j] != "-" for k in ids)]
    with open(f"{out}/{r.genome}.{r.family}.{r.subfamily}.fa", "w") as fh:
        for i, k in enumerate(ids):
            fh.write(f">{k.replace(':', '_').replace('..', '-')}\n{''.join(aln[k][j] for j in keep)}\n")
print(len(fam), "groups written")
EOF

# 2. RIPCaL, alignment-based, two reference models. Output names carry a time stamp, so each run has its own folder.
for f in "$B"/groups/*.fa; do
  g=$(basename "$f" .fa)
  for model in consensus gc; do
    d="$B/runs/$g/$model"; mkdir -p "$d"; cp "$f" "$d/input.fa"
    ( cd "$d" && rm -f ./*_ripalign.txt ./*.gif && $PERL "$B/src/ripcal2_install/perl/ripcal" -c -seq input.fa -type align -model $model > ripcal.log 2>&1 ) || echo "FAILED $g $model"
    t=$(ls "$d"/*_ripalign.txt 2>/dev/null | head -1)
    [ -n "$t" ] && cp "$t" "$d/ripalign.txt" && $B/env/bin/perl "$B/src/ripcal2_install/perl/ripcal_summarise" "$d/ripalign.txt" "$g" 1 > "$d/summary.txt" 2>/dev/null
  done
done

# 3. collect
mkdir -p results_ripcal
python - "$B/runs" results_ripcal <<'EOF'
import glob, os, sys
import pandas as pd
runs, out = sys.argv[1:3]
per_copy, summ = [], []
for d in sorted(glob.glob(f"{runs}/*/*")):
    g, model = d.split("/")[-2:]
    if not os.path.exists(f"{d}/ripalign.txt"):
        continue
    t = pd.read_csv(f"{d}/ripalign.txt", sep="\t"); t.columns = [c.strip() for c in t.columns]
    t.insert(0, "reference", model); t.insert(0, "group", g); per_copy.append(t)
    s = pd.read_csv(f"{d}/summary.txt", sep="\t"); s = s.loc[:, [c for c in s.columns if not c.startswith("Unnamed")]]
    s.columns = [c.split(" ", 1)[1] if c[:2].isdigit() else c for c in s.columns]
    s.insert(0, "reference", model); s.insert(0, "group", g); summ.append(s)
pd.concat(per_copy).to_csv(f"{out}/ripcal_per_copy.tsv", sep="\t", index=False)
pd.concat(summ).drop(columns=["model"], errors="ignore").to_csv(f"{out}/ripcal_summary.tsv", sep="\t", index=False)
print(len(summ), "RIPCaL runs collected")
EOF
