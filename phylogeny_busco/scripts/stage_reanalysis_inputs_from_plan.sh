#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/scratch/finger_millet_blast_phylogeny}"
PLAN="${PLAN:-$BASE/staging_plan_88_plus_CKF3589.tsv}"
DEST="${DEST:-$BASE/reanalysis_inputs}"
LOGDIR="$DEST/logs"

mkdir -p "$DEST/reference" "$DEST/samples" "$LOGDIR"

COPIED_LOG="$LOGDIR/staging_copied.tsv"
SKIPPED_LOG="$LOGDIR/staging_skipped.tsv"
MISSING_LOG="$LOGDIR/staging_missing.tsv"
ERROR_LOG="$LOGDIR/staging_errors.tsv"

printf "timestamp\tsample\tsource\tdestination\n" > "$COPIED_LOG"
printf "timestamp\tsample\tsource\tdestination\treason\n" > "$SKIPPED_LOG"
printf "timestamp\tsample\tisolate_id\tsection\tinput_method\n" > "$MISSING_LOG"
printf "timestamp\tsample\tsource\tdestination\treason\n" > "$ERROR_LOG"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

safe_name() {
  printf "%s" "$1" | tr -c "A-Za-z0-9._+-" "_"
}

copy_one() {
  local sample="$1"
  local src="$2"
  local dest_dir="$3"
  local base dest src_size dest_size

  if [[ -z "$src" ]]; then
    return 0
  fi
  if [[ ! -e "$src" ]]; then
    printf "%s\t%s\t%s\t%s\tmissing_source\n" "$(timestamp)" "$sample" "$src" "$dest_dir" >> "$ERROR_LOG"
    return 0
  fi

  mkdir -p "$dest_dir"
  base="$(basename "$src")"
  dest="$dest_dir/$base"
  src_size="$(stat -c%s "$src")"
  if [[ -e "$dest" ]]; then
    dest_size="$(stat -c%s "$dest")"
    if [[ "$src_size" == "$dest_size" ]]; then
      printf "%s\t%s\t%s\t%s\talready_present_same_size\n" "$(timestamp)" "$sample" "$src" "$dest" >> "$SKIPPED_LOG"
      return 0
    fi
  fi

  if rsync -a --partial "$src" "$dest_dir/"; then
    printf "%s\t%s\t%s\t%s\n" "$(timestamp)" "$sample" "$src" "$dest" >> "$COPIED_LOG"
  else
    printf "%s\t%s\t%s\t%s\trsync_failed\n" "$(timestamp)" "$sample" "$src" "$dest" >> "$ERROR_LOG"
  fi
}

copy_one "reference" "/path/to/lab_storage/E2Phylogenetic/trimmed/E2.single_copy.fasta" "$DEST/reference"
copy_one "reference" "/path/to/lab_storage/E2Phylogenetic/trimmed/E2.single_copy.fasta.fai" "$DEST/reference"
copy_one "reference" "/path/to/lab_storage/E2Phylogenetic/trimmed/E2.single_copy.dict" "$DEST/reference"

header="$(head -n 1 "$PLAN" | tr -d '\r')"
IFS=$'\t' read -r -a columns <<< "$header"
sample_idx=-1
isolate_idx=-1
section_idx=-1
input_idx=-1
paths_idx=-1
for i in "${!columns[@]}"; do
  case "${columns[$i]}" in
    sample) sample_idx="$i" ;;
    isolate_id) isolate_idx="$i" ;;
    section) section_idx="$i" ;;
    input_method) input_idx="$i" ;;
    source_paths) paths_idx="$i" ;;
  esac
done

if [[ "$sample_idx" -lt 0 || "$isolate_idx" -lt 0 || "$section_idx" -lt 0 || "$input_idx" -lt 0 || "$paths_idx" -lt 0 ]]; then
  printf "Required columns not found in %s\n" "$PLAN" >&2
  exit 1
fi

awk -v sample_col="$((sample_idx + 1))" \
    -v isolate_col="$((isolate_idx + 1))" \
    -v section_col="$((section_idx + 1))" \
    -v input_col="$((input_idx + 1))" \
    -v paths_col="$((paths_idx + 1))" \
    'BEGIN { FS = OFS = "\t" } NR > 1 { print $sample_col, $isolate_col, $section_col, $input_col, $paths_col }' "$PLAN" |
while IFS=$'\t' read -r sample isolate_id section input_method source_paths; do
  sample="${sample//$'\r'/}"
  isolate_id="${isolate_id//$'\r'/}"
  section="${section//$'\r'/}"
  input_method="${input_method//$'\r'/}"
  source_paths="${source_paths//$'\r'/}"
  sample_dir="$DEST/samples/$(safe_name "$sample")"
  mkdir -p "$sample_dir"
  printf "sample\t%s\nisolate_id\t%s\nsection\t%s\ninput_method\t%s\n" \
    "$sample" "$isolate_id" "$section" "$input_method" > "$sample_dir/sample_info.tsv"

  if [[ "$input_method" == "missing" || -z "${source_paths:-}" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\n" "$(timestamp)" "$sample" "$isolate_id" "$section" "$input_method" >> "$MISSING_LOG"
    continue
  fi

  IFS=';' read -r -a paths <<< "$source_paths"
  for src in "${paths[@]}"; do
    copy_one "$sample" "$src" "$sample_dir"
  done
done

find "$DEST" -maxdepth 3 -type f -printf "%s\t%p\n" | sort -k2,2 > "$LOGDIR/staged_files_with_sizes.tsv"
du -sh "$DEST" > "$LOGDIR/staged_total_size.txt"
