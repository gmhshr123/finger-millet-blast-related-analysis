#!/bin/bash
#SBATCH --job-name=del_clip_edta
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%A_%a.err

# Clipped-read TE check (script 12) with the combined EDTA TE library, for the deletion calls
# and for the matched control intervals.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a SAMtools/1.21-GCC-13.3.0 BLAST+/2.16.0-gompi-2024a
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" meta/copy38.txt)
mkdir -p softclip_edta softclip_ctrl_edta
python scripts/12_softclip_te.py --strain "$S" --te-db te_lib/edta/EDTA --name-map te_lib/edta/name_map.tsv --out-dir softclip_edta
python scripts/12_softclip_te.py --strain "$S" --te-db te_lib/edta/EDTA --name-map te_lib/edta/name_map.tsv \
   --calls inputs/control_intervals.tsv --out-dir softclip_ctrl_edta
echo "[$(date)] finished $S"
