#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/scratch/finger_millet_blast_phylogeny}"
REF="$BASE/reanalysis_inputs/reference/E2.single_copy.fasta"
OUTROOT="$BASE/reanalysis_assembly_alignments/final_assembly_backed"
SUMMARY="$OUTROOT/final_assembly_backed_target_recovery.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

set +u
source ~/.bashrc || true
set -u
module load minimap2/2.29-GCCcore-13.3.0 || true
module load SAMtools/1.21-GCC-13.3.0 || true

mkdir -p "$OUTROOT"
printf "sample\trole\tpreset\tassembly_fasta\tcoverage_fraction\tmean_depth\tnum_alignments\tbam\n" > "$SUMMARY"

assembly_path() {
  local sample="$1"
  find "$BASE/reanalysis_inputs/samples/$sample/assembly" -maxdepth 1 -type f -name "*_genomic.fna.gz" | head -1
}

summarize_coverage() {
  local sample="$1"
  local role="$2"
  local preset="$3"
  local asm="$4"
  local coverage="$5"
  local bam="$6"
  awk -v sample="$sample" -v role="$role" -v preset="$preset" -v asm="$asm" -v bam="$bam" '
    BEGIN { FS = OFS = "\t" }
    NR > 1 {
      len = $3 - $2 + 1
      total_len += len
      covered += $5
      depth_sum += $7 * len
      reads += $4
    }
    END {
      cov = total_len ? covered / total_len : 0
      depth = total_len ? depth_sum / total_len : 0
      printf "%s\t%s\t%s\t%s\t%.6f\t%.4f\t%d\t%s\n", sample, role, preset, asm, cov, depth, reads, bam
    }' "$coverage" >> "$SUMMARY"
}

run_one() {
  local sample="$1"
  local role="$2"
  local preset="$3"
  local asm
  local out
  local bam
  local coverage

  asm="$(assembly_path "$sample")"
  [[ -s "$asm" ]] || { echo "ERROR: missing assembly for $sample" >&2; exit 1; }
  out="$OUTROOT/$sample"
  mkdir -p "$out"
  bam="$out/${sample}.${preset}.target.sorted.bam"
  coverage="$out/${sample}.${preset}.coverage.tsv"

  if [[ ! -s "$bam" ]]; then
    minimap2 -t "$THREADS" -ax "$preset" "$REF" "$asm" |
      samtools sort -@ "$THREADS" -o "$bam" -
    samtools index -@ "$THREADS" "$bam"
  elif [[ ! -s "$bam.bai" ]]; then
    samtools index -@ "$THREADS" "$bam"
  fi
  samtools coverage "$bam" > "$coverage"
  summarize_coverage "$sample" "$role" "$preset" "$asm" "$coverage" "$bam"
}

for sample in B51 G22 PH42 chw ho p29 pg1213; do
  run_one "$sample" "ingroup_assembly" "asm5"
done

for sample in W97-11 DS0505; do
  run_one "$sample" "outgroup_assembly" "asm20"
done

cat "$SUMMARY"
