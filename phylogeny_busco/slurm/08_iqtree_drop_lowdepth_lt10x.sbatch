#!/bin/bash
#SBATCH --job-name=fm_iq_lt10
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
BASE="$PROJECT/results/reanalysis_phylogeny_inputs/read90_assembly9_MAFge0.03_missing30"
BASE_PREFIX="finger_millet.read90_assembly9.MAFge0.03_missing30"
OUTDIR="$PROJECT/results/reanalysis_phylogeny_inputs/drop_lowdepth_lt10x_MAFge0.03_missing30"
PREFIX="finger_millet.drop_lowdepth_lt10x.MAFge0.03_missing30"
IQDIR="$PROJECT/results/iqtree_drop_lowdepth_lt10x_MAFge0.03_missing30"
IQPREFIX="$IQDIR/$PREFIX"
READ_MANIFEST="$PROJECT/reanalysis_inputs/manifests/final_read_backed_samples.tsv"
REMOVE_LIST="$OUTDIR/$PREFIX.remove_samples.txt"

mkdir -p "$PROJECT/log" "$OUTDIR" "$IQDIR"

LOG="$IQDIR/iqtree_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date)"
echo "Project: $PROJECT"
echo "Subset: remove read-backed samples with mean depth <10x"

awk -F '\t' 'NR > 1 && ($9 + 0) < 10 {print $1}' "$READ_MANIFEST" | sort -u > "$REMOVE_LIST"
echo "Remove samples:"
cat "$REMOVE_LIST"

python3 "$PROJECT/scripts/filter_alignment_samples.py" \
  --fasta "$BASE/$BASE_PREFIX.fasta" \
  --sites "$BASE/$BASE_PREFIX.sites.tsv" \
  --sample-qc "$BASE/$BASE_PREFIX.sample_qc.tsv" \
  --manifest "$BASE/sample_manifest.tsv" \
  --remove-samples "$REMOVE_LIST" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"

module load IQ-TREE/3.0.1-gompi-2024a

IQTREE="$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)"
echo "IQ-TREE executable: $IQTREE"
"$IQTREE" --version || true

"$IQTREE" \
  -s "$OUTDIR/$PREFIX.phy" \
  -m GTR+ASC \
  -B 1000 \
  -alrt 1000 \
  -T "$THREADS" \
  -redo \
  -pre "$IQPREFIX"

echo "Finished: $(date)"
