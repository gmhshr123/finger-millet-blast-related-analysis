#!/bin/bash
#SBATCH --job-name=fm_map_mixed
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
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.err

set -eo pipefail

set +u
source ~/.bashrc || true
set -u
module load SAMtools/1.21-GCC-13.3.0 || true
module load minimap2/2.29-GCCcore-13.3.0 || true

PROJECT=/path/to/scratch/finger_millet_blast_phylogeny

for x in samtools bwa-mem2 minimap2; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: $x not found" >&2; exit 1; }
done

bash "$PROJECT/scripts/run_reanalysis_mapping_one.sh"
