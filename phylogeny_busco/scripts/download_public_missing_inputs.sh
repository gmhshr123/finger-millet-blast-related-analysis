#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/scratch/finger_millet_blast_phylogeny}"
DEST="${DEST:-$BASE/reanalysis_inputs}"
MANIFEST_DIR="${MANIFEST_DIR:-$DEST/manifests}"
READS_TSV="${READS_TSV:-$MANIFEST_DIR/missing_public_read_downloads.tsv}"
ASSEMBLY_TSV="${ASSEMBLY_TSV:-$MANIFEST_DIR/missing_public_assembly_downloads.tsv}"
LOGDIR="$DEST/logs"

mkdir -p "$DEST/samples" "$LOGDIR"

DOWNLOAD_LOG="$LOGDIR/public_missing_downloaded.tsv"
SKIPPED_LOG="$LOGDIR/public_missing_skipped.tsv"
ERROR_LOG="$LOGDIR/public_missing_errors.tsv"

[[ -f "$DOWNLOAD_LOG" ]] || printf "timestamp\tsample\tkind\turl\tdestination\n" > "$DOWNLOAD_LOG"
[[ -f "$SKIPPED_LOG" ]] || printf "timestamp\tsample\tkind\turl\tdestination\treason\n" > "$SKIPPED_LOG"
[[ -f "$ERROR_LOG" ]] || printf "timestamp\tsample\tkind\turl\tdestination\treason\n" > "$ERROR_LOG"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

safe_name() {
  printf "%s" "$1" | tr -c "A-Za-z0-9._+-" "_"
}

download_one() {
  local sample="$1"
  local kind="$2"
  local url="$3"
  local dest_dir="$4"
  local dest done_marker

  url="${url//$'\r'/}"
  if [[ -z "$url" ]]; then
    return 0
  fi
  if [[ "$url" != http://* && "$url" != https://* && "$url" != ftp://* ]]; then
    printf "%s\t%s\t%s\t%s\t%s\tinvalid_url\n" "$(timestamp)" "$sample" "$kind" "$url" "$dest_dir" >> "$SKIPPED_LOG"
    return 0
  fi

  mkdir -p "$dest_dir"
  dest="$dest_dir/$(basename "$url")"
  done_marker="$dest.done"
  if [[ -s "$dest" && -f "$done_marker" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\talready_completed\n" "$(timestamp)" "$sample" "$kind" "$url" "$dest" >> "$SKIPPED_LOG"
    return 0
  fi

  if wget -c --tries=5 --timeout=60 --waitretry=10 -P "$dest_dir" "$url"; then
    touch "$done_marker"
    printf "%s\t%s\t%s\t%s\t%s\n" "$(timestamp)" "$sample" "$kind" "$url" "$dest" >> "$DOWNLOAD_LOG"
  else
    printf "%s\t%s\t%s\t%s\t%s\twget_failed\n" "$(timestamp)" "$sample" "$kind" "$url" "$dest" >> "$ERROR_LOG"
  fi
}

if [[ -f "$READS_TSV" ]]; then
  awk 'BEGIN { FS = "\t"; OFS = "\034" } NR > 1 { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' "$READS_TSV" |
  while IFS=$'\034' read -r sample isolate_id source_accession run_accession sample_accession study_accession platform instrument layout url_1 url_2 url_3 read_count base_count; do
    sample="${sample//$'\r'/}"
    isolate_id="${isolate_id//$'\r'/}"
    run_accession="${run_accession//$'\r'/}"
    platform="${platform//$'\r'/}"
    instrument="${instrument//$'\r'/}"
    layout="${layout//$'\r'/}"
    sample_dir="$DEST/samples/$(safe_name "$sample")"
    run_dir="$sample_dir/public_reads/$(safe_name "$run_accession")"
    mkdir -p "$run_dir"
    {
      printf "sample\t%s\n" "$sample"
      printf "isolate_id\t%s\n" "$isolate_id"
      printf "source_accession\t%s\n" "${source_accession//$'\r'/}"
      printf "run_accession\t%s\n" "$run_accession"
      printf "sample_accession\t%s\n" "${sample_accession//$'\r'/}"
      printf "study_accession\t%s\n" "${study_accession//$'\r'/}"
      printf "platform\t%s\n" "$platform"
      printf "instrument\t%s\n" "$instrument"
      printf "layout\t%s\n" "$layout"
      printf "read_count\t%s\n" "${read_count//$'\r'/}"
      printf "base_count\t%s\n" "${base_count//$'\r'/}"
    } > "$run_dir/run_info.tsv"
    download_one "$sample" "read_fastq" "$url_1" "$run_dir"
    download_one "$sample" "read_fastq" "$url_2" "$run_dir"
    download_one "$sample" "read_fastq" "$url_3" "$run_dir"
  done
fi

if [[ -f "$ASSEMBLY_TSV" ]]; then
  awk 'BEGIN { FS = "\t"; OFS = "\034" } NR > 1 { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' "$ASSEMBLY_TSV" |
  while IFS=$'\034' read -r sample isolate_id source_accession assembly_uid assembly_accession assembly_name organism ftp_path genomic_fna_ftp genomic_fna_url; do
    sample="${sample//$'\r'/}"
    case "$sample" in
      B51|G22|PH42)
        sample_dir="$DEST/samples/$(safe_name "$sample")"
        asm_dir="$sample_dir/assembly"
        mkdir -p "$asm_dir"
        {
          printf "sample\t%s\n" "$sample"
          printf "isolate_id\t%s\n" "${isolate_id//$'\r'/}"
          printf "source_accession\t%s\n" "${source_accession//$'\r'/}"
          printf "assembly_accession\t%s\n" "${assembly_accession//$'\r'/}"
          printf "assembly_name\t%s\n" "${assembly_name//$'\r'/}"
          printf "organism\t%s\n" "${organism//$'\r'/}"
        } > "$asm_dir/assembly_info.tsv"
        download_one "$sample" "assembly_fasta" "$genomic_fna_url" "$asm_dir"
        ;;
    esac
  done
fi

find "$DEST/samples" -maxdepth 5 -type f -printf "%s\t%p\n" | sort -k2,2 > "$LOGDIR/public_missing_files_with_sizes.tsv"
du -sh "$DEST/samples" > "$LOGDIR/public_missing_samples_total_size.txt"
