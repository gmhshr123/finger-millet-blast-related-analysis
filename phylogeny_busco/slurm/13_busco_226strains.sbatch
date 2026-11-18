#!/bin/bash
#SBATCH --job-name=fm_busco_226
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
OUTBASE="$PROJECT/results/busco_226strains"
DOWNLOAD_PATH="$PROJECT/busco_downloads"

# Define the list of assembly files
STRAINS_LIST="$PROJECT/226strains/strains_list.txt"

# If the list file doesn't exist, create it by searching the 226strains directory
if [[ ! -f "$STRAINS_LIST" ]]; then
  echo "Generating strains list file at $STRAINS_LIST..."
  mkdir -p "$PROJECT/226strains"
  find "$PROJECT/226strains" -name "assemble_*.fasta" | sort > "$STRAINS_LIST"
fi

TOTAL_FILES=$(wc -l < "$STRAINS_LIST")
echo "Total assemblies found: $TOTAL_FILES"

# Get the task ID (default to 1 for manual testing)
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TASK_ID < 1 || TASK_ID > TOTAL_FILES )); then
  echo "ERROR: Slurm array task ID $TASK_ID is out of bounds (1 to $TOTAL_FILES)." >&2
  exit 1
fi

# Extract the Nth file path
STRAIN_FILE=$(sed -n "${TASK_ID}p" "$STRAINS_LIST")
STRAIN_NAME=$(basename "$STRAIN_FILE" .fasta)

echo "Task ID: $TASK_ID"
echo "Assembly name: $STRAIN_NAME"
echo "File path: $STRAIN_FILE"

# Check if output already exists to allow safe rerunning/resuming
SUMMARY_FILE="$OUTBASE/${STRAIN_NAME}_ascomycota_odb12/short_summary.specific.ascomycota_odb12.${STRAIN_NAME}_ascomycota_odb12.txt"
if [[ -s "$SUMMARY_FILE" ]]; then
  echo "BUSCO summary file already exists and is non-empty at: $SUMMARY_FILE"
  echo "Skipping run for $STRAIN_NAME."
  exit 0
fi

# Load BUSCO module
module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a

# Create output directories
mkdir -p "$OUTBASE" "$DOWNLOAD_PATH" "$PROJECT/log"

# Run BUSCO in genome mode in offline mode
busco \
  --in "$STRAIN_FILE" \
  --out "${STRAIN_NAME}_ascomycota_odb12" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --offline \
  --force

echo "Finished BUSCO for $STRAIN_NAME successfully."
