#!/bin/bash
#SBATCH --job-name=busco_226
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

set -euo pipefail

# Normal BUSCO genome-mode array run for many assembly FASTA files.
#
# Required input:
#   INPUT_DIR contains one FASTA per sample, for example:
#     assemble_E1.fasta
#     assemble_E12.fasta
#     assemble_UK6.fasta
#
# Recommended submission:
#   mkdir -p log
#   N=$(find /path/to/226strains -maxdepth 1 -name 'assemble_*.fasta' | wc -l)
#   sbatch --array=1-${N}%150 \
#     --export=ALL,INPUT_DIR=/path/to/226strains,OUTBASE=/path/to/results/busco_226_normal \
#     16_busco_226_normal_shareable.sbatch

PROJECT="${PROJECT:-$PWD}"
INPUT_DIR="${INPUT_DIR:-$PROJECT/226strains}"
OUTBASE="${OUTBASE:-$PROJECT/results/busco_226_normal}"
LINEAGE="${LINEAGE:-ascomycota_odb12}"
DOWNLOAD_PATH="${DOWNLOAD_PATH:-$PROJECT/busco_downloads}"
LIST_FILE="${LIST_FILE:-$OUTBASE/strains_list.txt}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
FASTA_GLOB="${FASTA_GLOB:-assemble_*.fasta}"
OFFLINE="${OFFLINE:-1}"

mkdir -p "$OUTBASE" "$DOWNLOAD_PATH" "$PROJECT/log"

if [[ ! -s "$LIST_FILE" ]]; then
  find "$INPUT_DIR" -maxdepth 1 -name "$FASTA_GLOB" | sort > "$LIST_FILE"
fi

TOTAL_FILES=$(wc -l < "$LIST_FILE")
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TOTAL_FILES < 1 )); then
  echo "ERROR: no FASTA files found in INPUT_DIR=$INPUT_DIR with FASTA_GLOB=$FASTA_GLOB" >&2
  exit 1
fi

if (( TASK_ID < 1 || TASK_ID > TOTAL_FILES )); then
  echo "ERROR: array task ID $TASK_ID is outside 1..$TOTAL_FILES" >&2
  exit 1
fi

FASTA=$(sed -n "${TASK_ID}p" "$LIST_FILE")
SAMPLE=$(basename "$FASTA")
SAMPLE="${SAMPLE%.fasta}"
SAMPLE="${SAMPLE%.fa}"
SAMPLE="${SAMPLE%.fna}"
RUN_NAME="${SAMPLE}_${LINEAGE}"
SUMMARY="$OUTBASE/$RUN_NAME/short_summary.specific.${LINEAGE}.${RUN_NAME}.txt"

echo "Started: $(date)"
echo "Task: $TASK_ID / $TOTAL_FILES"
echo "Input FASTA: $FASTA"
echo "Sample: $SAMPLE"
echo "Lineage: $LINEAGE"
echo "Output directory: $OUTBASE/$RUN_NAME"
echo "Threads: $THREADS"

if [[ ! -s "$FASTA" ]]; then
  echo "ERROR: input FASTA is missing or empty: $FASTA" >&2
  exit 1
fi

if [[ -s "$SUMMARY" ]]; then
  echo "Existing non-empty summary found; skipping: $SUMMARY"
  exit 0
fi

# Site-specific environment setup.
# Keep the first working option and delete/modify the others as needed.
module purge || true
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a || true

if ! command -v busco >/dev/null 2>&1; then
  echo "ERROR: busco command not found. Load a BUSCO module or activate a conda environment first." >&2
  exit 1
fi

busco --version

BUSCO_ARGS=(
  --in "$FASTA"
  --out "$RUN_NAME"
  --out_path "$OUTBASE"
  --lineage_dataset "$LINEAGE"
  --mode genome
  --cpu "$THREADS"
  --download_path "$DOWNLOAD_PATH"
  --force
)

if [[ "$OFFLINE" == "1" ]]; then
  BUSCO_ARGS+=(--offline)
fi

busco "${BUSCO_ARGS[@]}"

echo "Finished: $(date)"
