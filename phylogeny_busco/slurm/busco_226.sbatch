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

INPUT_DIR="/path/to/226strains"
OUT_DIR="/path/to/busco_226_results"
DOWNLOAD_DIR="/path/to/busco_downloads"
LINEAGE="ascomycota_odb12"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUT_DIR" "$DOWNLOAD_DIR" log

LIST="$OUT_DIR/assemblies.list"
find "$INPUT_DIR" -maxdepth 1 -name "assemble_*.fasta" | sort > "$LIST"

FASTA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
SAMPLE=$(basename "$FASTA" .fasta)

module purge
module load BUSCO/6.0.0-conda

busco \
  --in "$FASTA" \
  --out "${SAMPLE}_${LINEAGE}" \
  --out_path "$OUT_DIR" \
  --lineage_dataset "$LINEAGE" \
  --mode genome \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_DIR" \
  --offline \
  --force
