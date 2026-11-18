#!/bin/bash
#SBATCH --job-name=del_pairs
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%A_%a.err

# Read-pair evidence (is the DNA missing, or just unmappable?) for every depth-based call,
# and the same checks on matched control intervals where nothing is called deleted.
# Also runs the split-read TE check (script 12) on the control intervals, which gives the
# chance expectation for "clipped reads continue into a TE".
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a SAMtools/1.21-GCC-13.3.0 BLAST+/2.16.0-gompi-2024a
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" meta/copy38.txt)
mkdir -p pairs pairs_ctrl pairs_region softclip_ctrl
python scripts/14_read_pair_evidence.py --strain "$S" --out pairs/$S.pairs.tsv
python scripts/14_read_pair_evidence.py --strain "$S" --calls inputs/control_intervals.tsv --out pairs_ctrl/$S.pairs.tsv
# the same test on whole regions (fragments of one strain less than 1 kb apart, merged by 10_analyze.py)
python scripts/14_read_pair_evidence.py --strain "$S" --calls inputs/regions_strain_level.tsv --out pairs_region/$S.pairs.tsv
python scripts/12_softclip_te.py --strain "$S" --calls inputs/control_intervals.tsv --out-dir softclip_ctrl
echo "[$(date)] finished $S"
