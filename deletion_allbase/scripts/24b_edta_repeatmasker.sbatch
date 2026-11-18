#!/bin/bash
#SBATCH --job-name=del_edta_rm
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.err

# Locate copies of the combined EDTA TE library in E2 (same settings as 05 for the earlier library),
# then write a BED file in which each copy carries its family label from name_map.tsv.
set -euo pipefail
module load RepeatMasker/4.1.9-foss-2023a

cd te_lib/edta
mkdir -p rm
RepeatMasker -pa 4 -s -nolow -no_is -lib EDTA.fa -dir rm -gff ../../ref/E2.fasta
awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(FNR>1) lab[$1]=$2; next}
     {print $1,$2,$3,(($4 in lab)?lab[$4]:$4),$5,$6,$7,$8,$9,$10,$4}' name_map.tsv \
  <(awk 'NR>3 && $5!=""{print $5"\t"$6-1"\t"$7"\t"$10"\t"$1"\t"($9=="C"?"-":"+")"\t"$12"\t"$13"\t"$14"\t"$2}' rm/E2.fasta.out) \
  | sort -k1,1 -k2,2n > ../../ref/E2.edta_TEs.bed
cut -f4 ../../ref/E2.edta_TEs.bed | sort | uniq -c | sort -k1,1nr | head -40
tail -25 rm/E2.fasta.tbl
