#!/bin/bash
#SBATCH --job-name=fm_busco_E2
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
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.err

set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
E2_GENOME="${E2_GENOME:-$PROJECT/reanalysis_inputs/reference/E2.genome.fasta}"
LINEAGE="${LINEAGE:-ascomycota_odb12}"
OUTBASE="${OUTBASE:-$PROJECT/results/busco_E2_ascomycota_odb12}"
RUN_NAME="${RUN_NAME:-E2_ascomycota_odb12}"
DOWNLOAD_PATH="${DOWNLOAD_PATH:-$PROJECT/busco_downloads}"

mkdir -p "$PROJECT/log" "$OUTBASE" "$DOWNLOAD_PATH"

LOG="$OUTBASE/busco_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date)"
echo "Project: $PROJECT"
echo "E2 genome: $E2_GENOME"
echo "Lineage: $LINEAGE"
echo "Output base: $OUTBASE"
echo "Threads: $THREADS"

if [[ ! -s "$E2_GENOME" ]]; then
  echo "ERROR: E2_GENOME does not exist or is empty: $E2_GENOME" >&2
  echo "Set E2_GENOME=/path/to/E2_fungal_genome.fasta when submitting this job." >&2
  exit 2
fi

module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a

busco --version
busco \
  --in "$E2_GENOME" \
  --out "$RUN_NAME" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --force

echo "Finished: $(date)"
