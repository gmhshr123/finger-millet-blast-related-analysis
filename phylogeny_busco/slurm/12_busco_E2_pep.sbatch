#!/bin/bash
#SBATCH --job-name=fm_busco_E2_pep
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
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.err

set -euo pipefail

PROJECT="/path/to/scratch/finger_millet_blast_phylogeny"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
E2_PEP="$PROJECT/E2-annotation-Final/E2.pep.fasta"
LINEAGE="ascomycota_odb12"
OUTBASE="$PROJECT/results/busco_E2_pep_ascomycota_odb12"
RUN_NAME="E2_pep_ascomycota_odb12"
DOWNLOAD_PATH="$PROJECT/busco_downloads"

mkdir -p "$PROJECT/log" "$OUTBASE" "$DOWNLOAD_PATH"

LOG="$OUTBASE/busco_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date)"
echo "Project: $PROJECT"
echo "E2 proteins input: $E2_PEP"
echo "Lineage: $LINEAGE"
echo "Output base: $OUTBASE"
echo "Threads: $THREADS"

if [[ ! -s "$E2_PEP" ]]; then
  echo "ERROR: E2_PEP does not exist or is empty: $E2_PEP" >&2
  exit 2
fi

module purge
module load BUSCO/6.0.0-conda || module load BUSCO/5.8.3-foss-2023a

busco --version
busco \
  --in "$E2_PEP" \
  --out "$RUN_NAME" \
  --out_path "$OUTBASE" \
  --lineage_dataset "$LINEAGE" \
  --mode proteins \
  --cpu "$THREADS" \
  --download_path "$DOWNLOAD_PATH" \
  --offline \
  --force

echo "Finished E2 proteins BUSCO: $(date)"
