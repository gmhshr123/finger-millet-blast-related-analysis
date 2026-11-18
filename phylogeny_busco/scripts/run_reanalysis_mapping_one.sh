#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/scratch/finger_millet_blast_phylogeny}"
MANIFEST="${MANIFEST:-$BASE/reanalysis_inputs/manifests/mapping_manifest.tsv}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-${1:-}}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
REF="${REF:-$BASE/reanalysis_inputs/reference/E2.single_copy.fasta}"
PROJECT="${PROJECT:-$BASE/reanalysis_deepvariant}"
ASM_OUTROOT="${ASM_OUTROOT:-$BASE/reanalysis_assembly_alignments}"

if [[ -z "$TASK_ID" ]]; then
  echo "Need SLURM_ARRAY_TASK_ID or task id argument" >&2
  exit 2
fi

mkdir -p "$PROJECT"/{bam,qc/bam,tmp,log} "$ASM_OUTROOT"/{samples,log}

line="$(awk -v task="$TASK_ID" 'BEGIN { FS = "\t"; OFS = "\034" } NR > 1 && $1 == task { for (i = 1; i <= 14; i++) printf "%s%s", $i, (i < 14 ? OFS : ORS); exit }' "$MANIFEST")"
if [[ -z "$line" ]]; then
  echo "No manifest row for task $TASK_ID" >&2
  exit 2
fi

IFS=$'\034' read -r row_id sample input_kind mapper deepvariant_model short_pairs short_singles long_reads source_bams assembly_fasta manifest_bam manifest_gvcf manifest_vcf notes <<< "$line"
sample_out="$PROJECT/qc/bam/$sample"
tmpbase="${SLURM_TMPDIR:-$PROJECT/tmp}"
tmp="$tmpbase/${sample}_map_${SLURM_JOB_ID:-manual}_${TASK_ID}"
mkdir -p "$sample_out" "$tmp"

done_marker="$PROJECT/bam/$sample.mapping.done"
if [[ -f "$done_marker" ]]; then
  echo "Already completed: $sample"
  exit 0
fi

log() {
  printf "%s\t%s\t%s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$sample" "$*" >&2
}

fastq_records() {
  local file="$1"
  if [[ ! -s "$file" ]]; then
    echo 0
    return 0
  fi
  if gzip -t "$file" >/dev/null 2>&1; then
    gzip -cd "$file" 2>/dev/null | awk 'END { print int(NR / 4) }'
  else
    awk 'END { print int(NR / 4) }' "$file"
  fi
}

sort_bam() {
  local out="$1"
  samtools sort -@ "$THREADS" -T "$tmp/${sample}.coord" -o "$out" -
}

markdup_final() {
  local final_bam="$1"
  shift
  local bams=("$@")
  local merged="$tmp/${sample}.merged.coord.bam"
  local markdup_stats="$PROJECT/qc/bam/${sample}.markdup.stats.txt"
  if [[ "${#bams[@]}" -eq 1 ]]; then
    mv "${bams[0]}" "$merged"
  else
    samtools merge -@ "$THREADS" -f "$merged" "${bams[@]}"
    rm -f "${bams[@]}"
  fi
  if samtools markdup -@ "$THREADS" -s "$merged" "$final_bam" 2> "$markdup_stats"; then
    :
  else
    rm -f "$final_bam"
    printf "markdup_failed_falling_back_to_coordinate_sorted_bam\n" >> "$markdup_stats"
    mv "$merged" "$final_bam"
    return 0
  fi
  rm -f "$merged"
}

qc_bam() {
  local bam="$1"
  samtools index -@ "$THREADS" "$bam"
  samtools flagstat -@ "$THREADS" "$bam" > "$PROJECT/qc/bam/$sample.flagstat.txt"
  samtools stats -@ "$THREADS" "$bam" > "$PROJECT/qc/bam/$sample.samtools_stats.txt"
  samtools coverage "$bam" > "$PROJECT/qc/bam/$sample.coverage.tsv"
  samtools idxstats "$bam" > "$PROJECT/qc/bam/$sample.idxstats.txt"
}

map_short_fastqs() {
  local bams=()
  local unit=0
  local rg out bam_name bam_fix bam_sorted

  if [[ -n "${short_pairs:-}" ]]; then
    IFS=';' read -r -a pairs <<< "$short_pairs"
    for pair in "${pairs[@]}"; do
      [[ -z "$pair" ]] && continue
      IFS=',' read -r r1 r2 <<< "$pair"
      unit=$((unit + 1))
      rg="@RG\tID:${sample}_${unit}\tSM:${sample}\tPL:ILLUMINA"
      bam_name="$tmp/${sample}.pair_${unit}.name.bam"
      bam_fix="$tmp/${sample}.pair_${unit}.fixmate.bam"
      bam_sorted="$tmp/${sample}.pair_${unit}.coord.bam"
      log "bwa-mem2 pair $unit"
      bwa-mem2 mem -t "$THREADS" -R "$rg" "$REF" "$r1" "$r2" |
        samtools view -@ "$THREADS" -u -F 4 - |
        samtools sort -n -@ "$THREADS" -T "$tmp/${sample}.pair_${unit}.name" -o "$bam_name" -
      samtools fixmate -m -@ "$THREADS" "$bam_name" "$bam_fix"
      samtools sort -@ "$THREADS" -T "$tmp/${sample}.pair_${unit}.coord" -o "$bam_sorted" "$bam_fix"
      rm -f "$bam_name" "$bam_fix"
      bams+=("$bam_sorted")
    done
  fi

  if [[ -n "${short_singles:-}" ]]; then
    IFS=';' read -r -a singles <<< "$short_singles"
    for r1 in "${singles[@]}"; do
      [[ -z "$r1" ]] && continue
      unit=$((unit + 1))
      rg="@RG\tID:${sample}_${unit}\tSM:${sample}\tPL:ILLUMINA"
      out="$tmp/${sample}.short_single_${unit}.bam"
      log "bwa-mem2 single $unit"
      bwa-mem2 mem -t "$THREADS" -R "$rg" "$REF" "$r1" |
        samtools view -@ "$THREADS" -u -F 4 - |
        sort_bam "$out"
      bams+=("$out")
    done
  fi

  if [[ "${#bams[@]}" -eq 0 ]]; then
    log "no short-read files found"
    exit 1
  fi

  final="$manifest_bam"
  markdup_final "$final" "${bams[@]}"
  qc_bam "$final"
}

map_bam_extracted_reads() {
  local bams=()
  local unit=0
  local r1_count r2_count single_count singleton_count
  IFS=';' read -r -a source_array <<< "$source_bams"
  for source_bam in "${source_array[@]}"; do
    [[ -z "$source_bam" ]] && continue
    unit=$((unit + 1))
    prefix="$tmp/${sample}.bam_${unit}"
    r1="$prefix.R1.fq.gz"
    r2="$prefix.R2.fq.gz"
    single="$prefix.single.fq.gz"
    singleton="$prefix.singleton.fq.gz"
    log "extract reads from source BAM $unit"
    samtools collate -@ "$THREADS" -u -O "$source_bam" |
      samtools fastq -@ "$THREADS" -1 "$r1" -2 "$r2" -0 "$single" -s "$singleton" -n -
    rg="@RG\tID:${sample}_bam_${unit}\tSM:${sample}\tPL:ILLUMINA"

    r1_count="$(fastq_records "$r1")"
    r2_count="$(fastq_records "$r2")"
    single_count="$(fastq_records "$single")"
    singleton_count="$(fastq_records "$singleton")"
    log "extracted FASTQ records pair=(${r1_count},${r2_count}) single=${single_count} singleton=${singleton_count}"

    if [[ "$r1_count" -gt 0 && "$r2_count" -gt 0 ]]; then
      out="$tmp/${sample}.bam_extract_${unit}.paired.bam"
      bam_name="$tmp/${sample}.bam_extract_${unit}.name.bam"
      bam_fix="$tmp/${sample}.bam_extract_${unit}.fixmate.bam"
      log "bwa-mem2 extracted BAM pair $unit"
      bwa-mem2 mem -t "$THREADS" -R "$rg" "$REF" "$r1" "$r2" |
        samtools view -@ "$THREADS" -u -F 4 - |
        samtools sort -n -@ "$THREADS" -T "$tmp/${sample}.bam_extract_${unit}.name" -o "$bam_name" -
      samtools fixmate -m -@ "$THREADS" "$bam_name" "$bam_fix"
      samtools sort -@ "$THREADS" -T "$tmp/${sample}.bam_extract_${unit}.coord" -o "$out" "$bam_fix"
      rm -f "$bam_name" "$bam_fix"
      bams+=("$out")
    fi
    if [[ "$single_count" -gt 0 ]]; then
      out="$tmp/${sample}.bam_extract_${unit}.single.bam"
      log "bwa-mem2 extracted BAM single $unit"
      bwa-mem2 mem -t "$THREADS" -R "$rg" "$REF" "$single" |
        samtools view -@ "$THREADS" -u -F 4 - |
        sort_bam "$out"
      bams+=("$out")
    fi
    if [[ "$singleton_count" -gt 0 ]]; then
      out="$tmp/${sample}.bam_extract_${unit}.singleton.bam"
      log "bwa-mem2 extracted BAM singleton $unit"
      bwa-mem2 mem -t "$THREADS" -R "$rg" "$REF" "$singleton" |
        samtools view -@ "$THREADS" -u -F 4 - |
        sort_bam "$out"
      bams+=("$out")
    fi
  done

  if [[ "${#bams[@]}" -eq 0 ]]; then
    log "no BAM-derived reads mapped"
    exit 1
  fi

  final="$manifest_bam"
  markdup_final "$final" "${bams[@]}"
  qc_bam "$final"
}

map_long_reads() {
  IFS=';' read -r -a reads <<< "$long_reads"
  final="$manifest_bam"
  log "minimap2 map-ont"
  minimap2 -t "$THREADS" -ax map-ont "$REF" "${reads[@]}" |
    samtools view -@ "$THREADS" -u -F 4 - |
    sort_bam "$final"
  qc_bam "$final"
}

align_assembly() {
  asm_sample_out="$ASM_OUTROOT/samples/$sample"
  mkdir -p "$asm_sample_out"
  final="$asm_sample_out/$sample.assembly_to_target.sorted.bam"
  paf="$asm_sample_out/$sample.assembly_to_target.paf"
  log "minimap2 asm5 assembly alignment"
  minimap2 -t "$THREADS" -x asm5 "$REF" "$assembly_fasta" > "$paf"
  minimap2 -t "$THREADS" -ax asm5 "$REF" "$assembly_fasta" | sort_bam "$final"
  samtools index -@ "$THREADS" "$final"
  samtools flagstat -@ "$THREADS" "$final" > "$asm_sample_out/$sample.flagstat.txt"
  samtools stats -@ "$THREADS" "$final" > "$asm_sample_out/$sample.samtools_stats.txt"
  samtools coverage "$final" > "$asm_sample_out/$sample.coverage.tsv"
  samtools idxstats "$final" > "$asm_sample_out/$sample.idxstats.txt"
}

{
  printf "id\t%s\n" "$row_id"
  printf "sample\t%s\n" "$sample"
  printf "input_kind\t%s\n" "$input_kind"
  printf "mapper\t%s\n" "$mapper"
  printf "deepvariant_model\t%s\n" "$deepvariant_model"
  printf "notes\t%s\n" "${notes:-}"
  printf "reference\t%s\n" "$REF"
  printf "started\t%s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
} > "$sample_out/$sample.mapping_info.tsv"

case "$input_kind" in
  short_reads)
    map_short_fastqs
    ;;
  bam_extract_short_reads)
    map_bam_extracted_reads
    ;;
  long_reads)
    map_long_reads
    ;;
  assembly)
    align_assembly
    ;;
  *)
    log "unsupported or missing input kind: $input_kind"
    exit 1
    ;;
esac

printf "finished\t%s\n" "$(date '+%Y-%m-%d %H:%M:%S')" >> "$sample_out/$sample.mapping_info.tsv"
touch "$done_marker"
rm -rf "$tmp"
log "complete"
