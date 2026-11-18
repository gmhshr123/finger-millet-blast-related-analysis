#!/bin/bash
#SBATCH --job-name=fm_busco226_nsplit
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
RUN_BASE="$PROJECT/results/busco_226_validation_20260601"
OUTBASE="$RUN_BASE/nsplit_fasta_min200"
NSPLIT_BASE="$RUN_BASE/nsplit_inputs_min200"
DOWNLOAD_PATH="$PROJECT/busco_downloads"
STRAINS_LIST="$RUN_BASE/strains_list.txt"
SPLITTER="$PROJECT/scripts/split_fasta_at_n_runs.py"

mkdir -p "$RUN_BASE" "$OUTBASE" "$NSPLIT_BASE" "$DOWNLOAD_PATH" "$PROJECT/log"

if [[ ! -s "$STRAINS_LIST" ]]; then
  find "$PROJECT/226strains" -maxdepth 1 -name "assemble_*.fasta" | sort > "$STRAINS_LIST"
fi

TOTAL_FILES=$(wc -l < "$STRAINS_LIST")
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

if (( TASK_ID < 1 || TASK_ID > TOTAL_FILES )); then
  echo "ERROR: Slurm array task ID $TASK_ID is out of bounds (1 to $TOTAL_FILES)." >&2
  exit 1
fi

STRAIN_FILE=$(sed -n "${TASK_ID}p" "$STRAINS_LIST")
STRAIN_NAME=$(basename "$STRAIN_FILE" .fasta)
NSPLIT_FILE="$NSPLIT_BASE/${STRAIN_NAME}.nsplit.min200.fasta"
RUN_NAME="${STRAIN_NAME}_nsplit_min200_ascomycota_odb12"
SUMMARY_FILE="$OUTBASE/$RUN_NAME/short_summary.specific.ascomycota_odb12.$RUN_NAME.txt"

echo "Validation method: nsplit_fasta_min200"
echo "Task ID: $TASK_ID / $TOTAL_FILES"
echo "Assembly name: $STRAIN_NAME"
echo "Input: $STRAIN_FILE"
echo "N-split input: $NSPLIT_FILE"
echo "Output: $OUTBASE/$RUN_NAME"

if [[ -s "$SUMMARY_FILE" ]]; then
  echo "Existing summary found, skipping: $SUMMARY_FILE"
  exit 0
fi

python3 "$SPLITTER" --input "$STRAIN_FILE" --output "$NSPLIT_FILE" --min-len 200

module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a
busco --version

busco \
  --in "$NSPLIT_FILE" \
  --out "$RUN_NAME" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --offline \
  --force

echo "Finished nsplit_fasta_min200 BUSCO for $STRAIN_NAME"
