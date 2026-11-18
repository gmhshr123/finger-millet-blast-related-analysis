#!/bin/bash
#SBATCH --job-name=del_prep_ref
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

# Freeze the E2 reference used by the original RefE2 alignments and build the BWA index.
set -euo pipefail
module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.21-GCC-13.3.0 BEDTools/2.31.1-GCC-13.3.0

SRC=/path/to/scratch/finger_millet_blast_phylogeny/E2-annotation-Final
mkdir -p ref
cp -n "$SRC"/E2.fasta "$SRC"/E2.gff3 "$SRC"/E2.repeat.out "$SRC"/E2-gene-pos.txt "$SRC"/E2_gap_pos.txt ref/

cd ref
samtools faidx E2.fasta
bwa index E2.fasta
cut -f1,2 E2.fasta.fai > E2.genome
# analysis is restricted to the seven pseudomolecules
grep '^Chr' E2.genome > E2.chr.genome
bedtools makewindows -g E2.chr.genome -w 100000 > E2.chr.100kb.bed
# assembly gaps: file is 1-based inclusive -> BED
awk 'NR>1{print $1"\t"$2-1"\t"$3}' E2_gap_pos.txt > E2.gaps.bed
# RepeatMasker .out -> BED (name = family, class in column 5)
awk 'NR>3 && $5!=""{print $5"\t"$6-1"\t"$7"\t"$10"\t"$11"\t"($9=="C"?"-":"+")}' E2.repeat.out \
  | sort -k1,1 -k2,2n > E2.repeats.bed
tr -d '\r' < E2-gene-pos.txt | awk 'BEGIN{OFS="\t"}{print $2,$3-1,$4,$1}' | sort -k1,1 -k2,2n > E2.genes.bed
md5sum E2.fasta E2.gff3 E2.repeat.out > MD5SUMS.txt
bwa 2>&1 | grep -i version > VERSIONS.txt || true
samtools --version | head -2 >> VERSIONS.txt
bedtools --version >> VERSIONS.txt
echo done
