#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
MANIFEST="${MANIFEST:-$PROJECT/reanalysis_inputs/manifests/final_read_backed_samples.tsv}"
OUT="${OUT:-$PROJECT/reanalysis_deepvariant/gvcf/final_read_backed_gvcf_list.txt}"

mkdir -p "$(dirname "$OUT")"

awk '
  BEGIN { FS = "\t" }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "gvcf") {
        gvcf_col = i
      }
    }
    if (!gvcf_col) {
      print "ERROR: manifest is missing required gvcf column" > "/dev/stderr"
      exit 1
    }
    next
  }
  $gvcf_col != "" { print $gvcf_col }
' "$MANIFEST" > "$OUT.tmp"
missing=0
while read -r gvcf; do
  if [[ ! -s "$gvcf" ]]; then
    echo "MISSING: $gvcf" >&2
    missing=$((missing + 1))
  fi
done < "$OUT.tmp"

if [[ "$missing" -ne 0 ]]; then
  echo "ERROR: $missing gVCFs are missing; not writing final list" >&2
  exit 1
fi

mv "$OUT.tmp" "$OUT"
echo "Wrote $OUT with $(wc -l < "$OUT") gVCFs"
