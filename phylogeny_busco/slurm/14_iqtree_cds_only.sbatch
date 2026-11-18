#!/bin/bash
#SBATCH --job-name=fm_iq_cds
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
INDIR="$PROJECT/results/reanalysis_phylogeny_inputs/CDSonly"
PREFIX="fm.CDSonly"
IQDIR="$PROJECT/results/iqtree_CDSonly"
IQPREFIX="$IQDIR/$PREFIX"

mkdir -p "$PROJECT/log" "$IQDIR"

LOG="$IQDIR/iqtree_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date)"
echo "Project: $PROJECT"
echo "Input: $INDIR/$PREFIX.phy"
echo "Tree: CDS-only SNP sensitivity analysis"

module load IQ-TREE/3.0.1-gompi-2024a

IQTREE="$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)"
echo "IQ-TREE executable: $IQTREE"
"$IQTREE" --version || true

"$IQTREE" \
  -s "$INDIR/$PREFIX.phy" \
  -m GTR+ASC \
  -B 1000 \
  -alrt 1000 \
  -T "$THREADS" \
  -redo \
  -pre "$IQPREFIX"

echo "Finished: $(date)"
