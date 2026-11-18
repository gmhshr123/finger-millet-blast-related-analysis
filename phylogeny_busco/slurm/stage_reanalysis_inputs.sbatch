#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=fm_stage_inputs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu

set -euo pipefail

cd /path/to/scratch/finger_millet_blast_phylogeny
bash scripts/stage_reanalysis_inputs_from_plan.sh
