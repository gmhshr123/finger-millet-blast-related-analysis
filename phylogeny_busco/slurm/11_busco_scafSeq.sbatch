#!/bin/bash
#SBATCH --job-name=fm_busco_scaf
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
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.err

set -euo pipefail

PROJECT="/path/to/scratch/finger_millet_blast_phylogeny"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
LINEAGE="ascomycota_odb12"
OUTBASE="$PROJECT/results/busco_scafSeq"
DOWNLOAD_PATH="$PROJECT/busco_downloads"

# Define the list of scaffolding files
SCAF_LIST="$PROJECT/scaf/scaf_list.txt"

# If the list file doesn't exist, create it by searching the scaf directory
if [[ ! -f "$SCAF_LIST" ]]; then
  echo "Generating scaffold list file at $SCAF_LIST..."
  mkdir -p "$PROJECT/scaf"
  find "$PROJECT/scaf" -name "*.scafSeq" | sort > "$SCAF_LIST"
fi

TOTAL_FILES=$(wc -l < "$SCAF_LIST")
echo "Total assemblies found: $TOTAL_FILES"

# Get the task ID (default to 1 for manual testing)
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TASK_ID < 1 || TASK_ID > TOTAL_FILES )); then
  echo "ERROR: Slurm array task ID $TASK_ID is out of bounds (1 to $TOTAL_FILES)." >&2
  exit 1
fi

# Extract the Nth file path
SCAF_FILE=$(sed -n "${TASK_ID}p" "$SCAF_LIST")
SCAF_NAME=$(basename "$SCAF_FILE" .scafSeq)

echo "Task ID: $TASK_ID"
echo "Assembly name: $SCAF_NAME"
echo "File path: $SCAF_FILE"

# Load BUSCO module
module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a

# Create output directories
mkdir -p "$OUTBASE" "$DOWNLOAD_PATH" "$PROJECT/log"

# Check if output already exists to allow safe rerunning/resuming
SUMMARY_FILE="$OUTBASE/${SCAF_NAME}_ascomycota_odb12/short_summary.specific.ascomycota_odb12.${SCAF_NAME}_ascomycota_odb12.txt"
if [[ -s "$SUMMARY_FILE" ]]; then
  echo "BUSCO summary file already exists and is non-empty at: $SUMMARY_FILE"
  echo "Skipping run for $SCAF_NAME."
  exit 0
fi

# Run BUSCO in genome mode
busco \
  --in "$SCAF_FILE" \
  --out "${SCAF_NAME}_ascomycota_odb12" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --offline \
  --force

echo "Finished BUSCO for $SCAF_NAME successfully."
