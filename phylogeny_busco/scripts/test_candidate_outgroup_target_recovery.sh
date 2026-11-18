#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/scratch/finger_millet_blast_phylogeny}"
REF="$BASE/reanalysis_inputs/reference/E2.single_copy.fasta"
OUT="$BASE/reanalysis_assembly_alignments/candidate_outgroups"
SUMMARY="$OUT/candidate_outgroup_target_recovery.tsv"

set +u
source ~/.bashrc || true
set -u
module load minimap2/2.29-GCCcore-13.3.0 || true
module load SAMtools/1.21-GCC-13.3.0 || true

mkdir -p "$OUT"
printf "sample\tpreset\tcoverage_fraction\tmean_depth\tnum_alignments\tbam\n" > "$SUMMARY"

for sample in W97-11 DsLIZ VO107 DS0505 DS9461 W98-20; do
  asm="$(find "$BASE/reanalysis_inputs/samples/$sample/assembly" -maxdepth 1 -name "*_genomic.fna.gz" | head -1)"
  [[ -s "$asm" ]] || { echo "missing assembly for $sample" >&2; continue; }
  for preset in asm5 asm10 asm20; do
    sample_out="$OUT/$sample/$preset"
    mkdir -p "$sample_out"
    bam="$sample_out/${sample}.${preset}.target.sorted.bam"
    coverage="$sample_out/${sample}.${preset}.coverage.tsv"
    minimap2 -t "${SLURM_CPUS_PER_TASK:-8}" -ax "$preset" "$REF" "$asm" |
      samtools sort -@ "${SLURM_CPUS_PER_TASK:-8}" -o "$bam" -
    samtools index -@ "${SLURM_CPUS_PER_TASK:-8}" "$bam"
    samtools coverage "$bam" > "$coverage"
    awk -v sample="$sample" -v preset="$preset" -v bam="$bam" '
      BEGIN { FS = OFS = "\t" }
      NR > 1 {
        len = $3 - $2 + 1
        total_len += len
        covered += $6
        depth_sum += $7 * len
        reads += $4
      }
      END {
        cov = total_len ? covered / total_len : 0
        depth = total_len ? depth_sum / total_len : 0
        printf "%s\t%s\t%.6f\t%.4f\t%d\t%s\n", sample, preset, cov, depth, reads, bam
      }' "$coverage" >> "$SUMMARY"
  done
done

sort -k2,2 -k3,3nr "$SUMMARY" > "$SUMMARY.sorted"
mv "$SUMMARY.sorted" "$SUMMARY"
cat "$SUMMARY"
