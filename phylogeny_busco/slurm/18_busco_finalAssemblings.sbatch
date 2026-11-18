#!/bin/bash
#SBATCH --job-name=fm_busco_final
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.err

set -euo pipefail

PROJECT="/path/to/scratch/finger_millet_blast_phylogeny"
INPUT_DIR="$PROJECT/finalAssemblings"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
LINEAGE="ascomycota_odb12"
OUTBASE="$PROJECT/results/busco_finalAssemblings_ascomycota_odb12"
DOWNLOAD_PATH="$PROJECT/busco_downloads"
LIST_FILE="$OUTBASE/strains_list.txt"

mkdir -p "$PROJECT/log" "$OUTBASE" "$DOWNLOAD_PATH"

if [[ ! -s "$LIST_FILE" ]]; then
  find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fasta" | sort > "$LIST_FILE"
fi

TOTAL_FILES=$(wc -l < "$LIST_FILE")
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

echo "Started: $(date)"
echo "Task ID: $TASK_ID / $TOTAL_FILES"
echo "Input directory: $INPUT_DIR"
echo "Output base: $OUTBASE"
echo "Lineage: $LINEAGE"
echo "Threads: $THREADS"

if (( TOTAL_FILES < 1 )); then
  echo "ERROR: no *.fasta files found in $INPUT_DIR" >&2
  exit 1
fi

if (( TASK_ID < 1 || TASK_ID > TOTAL_FILES )); then
  echo "ERROR: Slurm array task ID $TASK_ID is out of bounds (1 to $TOTAL_FILES)." >&2
  exit 1
fi

STRAIN_FILE=$(sed -n "${TASK_ID}p" "$LIST_FILE")
STRAIN_NAME=$(basename "$STRAIN_FILE" .fasta)
RUN_NAME="${STRAIN_NAME}_${LINEAGE}"
SUMMARY_FILE="$OUTBASE/$RUN_NAME/short_summary.specific.${LINEAGE}.${RUN_NAME}.txt"

echo "Assembly name: $STRAIN_NAME"
echo "FASTA: $STRAIN_FILE"
echo "Run name: $RUN_NAME"

if [[ ! -s "$STRAIN_FILE" ]]; then
  echo "ERROR: FASTA does not exist or is empty: $STRAIN_FILE" >&2
  exit 1
fi

if [[ -s "$SUMMARY_FILE" ]]; then
  echo "Existing BUSCO summary found; skipping: $SUMMARY_FILE"
  exit 0
fi

module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a
busco --version

busco \
  --in "$STRAIN_FILE" \
  --out "$RUN_NAME" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --offline \
  --force

echo "Finished BUSCO for $STRAIN_NAME: $(date)"
