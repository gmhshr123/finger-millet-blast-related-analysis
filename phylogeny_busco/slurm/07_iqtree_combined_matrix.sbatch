#!/bin/bash
#SBATCH --job-name=fm_iqtree_99
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.err

set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
ALIGN="${ALIGN:-$PROJECT/results/reanalysis_phylogeny_inputs/read90_assembly9_MAFge0.05_missing30/finger_millet.read90_assembly9.MAFge0.05_missing30.phy}"
OUTDIR="${OUTDIR:-$PROJECT/results/iqtree_reanalysis_read90_assembly9_MAFge0.05_missing30}"
PREFIX="${PREFIX:-$OUTDIR/finger_millet.read90_assembly9.MAFge0.05_missing30}"

mkdir -p "$OUTDIR" "$PROJECT/log"

LOG="$OUTDIR/iqtree_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date)"
echo "Project: $PROJECT"
echo "Alignment: $ALIGN"
echo "Threads: $THREADS"
echo "Prefix: $PREFIX"

if command -v module >/dev/null 2>&1; then
  module load IQ-TREE/3.0.1-gompi-2024a || module load IQ-TREE || true
fi

if [[ ! -s "$ALIGN" ]]; then
  echo "ERROR: missing alignment: $ALIGN" >&2
  exit 1
fi

if command -v iqtree3 >/dev/null 2>&1; then
  IQTREE="$(command -v iqtree3)"
elif command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="$(command -v iqtree2)"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="$(command -v iqtree)"
else
  echo "ERROR: iqtree3/iqtree2/iqtree not found on PATH" >&2
  exit 1
fi

echo "IQ-TREE executable: $IQTREE"
"$IQTREE" --version || true

"$IQTREE" \
  -s "$ALIGN" \
  -m GTR+ASC \
  -B 1000 \
  -alrt 1000 \
  -T "$THREADS" \
  -redo \
  -pre "$PREFIX"

echo "Finished: $(date)"
