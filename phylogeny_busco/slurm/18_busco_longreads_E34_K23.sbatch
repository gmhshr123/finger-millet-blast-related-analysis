#!/bin/bash
#SBATCH --job-name=fm_busco_long
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
LINEAGE="ascomycota_odb12"
OUTBASE="$PROJECT/results/busco_longReads_E34_K23_ascomycota_odb12"
DOWNLOAD_PATH="$PROJECT/busco_downloads"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

FASTAS=(
  "$PROJECT/longReads/E34/E34.fasta"
  "$PROJECT/longReads/K23/K23.fasta"
)

FASTA="${FASTAS[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE=$(basename "$FASTA" .fasta)
RUN_NAME="${SAMPLE}_${LINEAGE}"

mkdir -p "$OUTBASE" "$DOWNLOAD_PATH" "$PROJECT/log"

module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a

busco \
  --in "$FASTA" \
  --out "$RUN_NAME" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --offline \
  --force
