#!/bin/bash
#SBATCH --job-name=del_edta_prep
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.err

# Combined EDTA TE library (B71, E2, MZ5-1-6, P131): BLAST database for the clipped-read check,
# and a table that gives each library entry a family label (published name where one matches).
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a BLAST+/2.16.0-gompi-2024a

cd te_lib/edta
tr -d '\r' < 4G.EDTA.TElib.fa > EDTA.fa
md5sum 4G.EDTA.TElib.fa EDTA.fa > MD5SUMS.txt
makeblastdb -in EDTA.fa -dbtype nucl -out EDTA > /dev/null
blastn -task blastn -query EDTA.fa -db ../named_TEs -evalue 1e-10 -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" > EDTA_vs_named.tsv
grep '>' ../named_TEs.fa | sed 's/>//; s/\(_[0-9]*\)\?#.*//' | sort -u > named_families.txt
python ../../scripts/24_te_name_map.py EDTA.fa EDTA_vs_named.tsv name_map.tsv named_families.txt
