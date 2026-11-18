#!/bin/bash
#SBATCH --job-name=rip_e2
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/rip_analysis
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/rip_analysis/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/rip_analysis/log/%x_%j.err

# RIP-type mutation analysis of full-length TE copies. GENOMES is a space-separated list of NAME=fasta.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a MAFFT/7.526-GCC-13.3.0-with-extensions
GENOMES=${GENOMES:-"E2=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/ref/E2.fasta"}
OUT=${OUT:-results_E2}
args=(); for g in $GENOMES; do args+=(--genome "$g"); done
mafft --version 2>&1 | head -1
python scripts/02_rip_family.py --copies inputs/full_length_copies.tsv "${args[@]}" --out "$OUT" --threads 8