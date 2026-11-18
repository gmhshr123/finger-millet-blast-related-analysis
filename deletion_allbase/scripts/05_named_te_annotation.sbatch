#!/bin/bash
#SBATCH --job-name=del_te_annot
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.err

# Locate copies of the named M. oryzae transposable elements in E2, so that deletion
# breakpoints can be labelled by TE family (the existing E2.repeat.out only has
# anonymous RepeatModeler family names).
set -euo pipefail
module load RepeatMasker/4.1.9-foss-2023a

cd te_lib
mkdir -p rm_named
# -nolow / -no_is: only the supplied TE library; -s: sensitive; -a keeps alignments
RepeatMasker -pa 4 -s -nolow -no_is -lib named_TEs.fa -dir rm_named -gff ../ref/E2.fasta
awk 'NR>3 && $5!=""{print $5"\t"$6-1"\t"$7"\t"$10"\t"$1"\t"($9=="C"?"-":"+")"\t"$12"\t"$13"\t"$14"\t"$2}' rm_named/E2.fasta.out \
  | sort -k1,1 -k2,2n > ../ref/E2.named_TEs.bed
cut -f4 ../ref/E2.named_TEs.bed | sort | uniq -c
