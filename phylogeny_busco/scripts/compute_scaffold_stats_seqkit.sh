#!/bin/bash
set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
cd "$PROJECT"

module load SeqKit/2.9.0 >/dev/null 2>&1 || true

write_stats() {
  local input_dir="$1"
  local find_pattern="$2"
  local suffix="$3"
  local outfile="$4"

  mkdir -p "$(dirname "$outfile")"
  printf "sample\tsource_file\tscaffold_count\tscaffolds_gt_1000bp\tscaffolds_ge_1000bp\ttotal_length_bp\tscaffold_n50_bp\n" > "$outfile"

  find "$input_dir" -maxdepth 1 -type f -name "$find_pattern" -print0 |
    sort -z |
    while IFS= read -r -d '' fasta; do
    [[ -s "$fasta" ]] || continue
    local sample
    sample="$(basename "$fasta")"
    sample="${sample%"$suffix"}"
    seqkit fx2tab -n -l "$fasta" |
      awk -v sample="$sample" -v fasta="$fasta" '
        BEGIN { FS="\t"; count=0; gt1000=0; ge1000=0; total=0 }
        {
          len=$NF + 0
          count++
          if (len > 1000) gt1000++
          if (len >= 1000) ge1000++
          total += len
        }
        END {
          print sample "\t" fasta "\t" count "\t" gt1000 "\t" ge1000 "\t" total
        }
      ' > "${outfile}.tmp_row"
    local n50
    n50="$(seqkit stats -T -a "$fasta" | awk -F '\t' 'NR == 2 { print $13 }')"
    printf "%s\t%s\n" "$(cat "${outfile}.tmp_row")" "$n50" >> "$outfile"
    rm -f "${outfile}.tmp_row"
  done
}

write_stats \
  scaf \
  "*.scafSeq" \
  ".scafSeq" \
  results/manifests_and_qc/pre216_scaf_fasta_scaffold_stats.tsv

write_stats \
  final \
  "*.fasta" \
  ".ragtag.fasta" \
  results/manifests_and_qc/final_ragtag_fasta_scaffold_stats.tsv

wc -l \
  results/manifests_and_qc/pre216_scaf_fasta_scaffold_stats.tsv \
  results/manifests_and_qc/final_ragtag_fasta_scaffold_stats.tsv
