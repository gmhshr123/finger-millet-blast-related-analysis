#!/bin/bash
#SBATCH --job-name=del_bpseq
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
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.err

# Sequences at the breakpoints of every strain-level deleted region (script 28), for the 33 cohort strains.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a SAMtools/1.21-GCC-13.3.0 BLAST+/2.16.0-gompi-2024a
mkdir -p bpseq tmp/E2_blastdb
[ -s tmp/E2_blastdb/E2.nsq ] || makeblastdb -in ref/E2.fasta -dbtype nucl -out tmp/E2_blastdb/E2 > /dev/null
while read -r S; do
    [ -s "bpseq/$S.breakpoint_sequences.tsv" ] && continue
    python scripts/28_breakpoint_sequences.py --strain "$S"
done < meta/cohort33.txt
echo "[$(date)] finished"
