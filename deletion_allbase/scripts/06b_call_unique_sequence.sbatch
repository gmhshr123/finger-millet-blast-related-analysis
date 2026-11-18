#!/bin/bash
#SBATCH --job-name=del_call2
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.err

# Primary all-base deletion calls: only non-repeat, callable sequence is used as evidence;
# repeats and gaps are bridged. Run for the BWA realignment (primary) and for the original
# Bowtie2 alignment (aligner check). Then score every interval set in every strain.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a

cut -f1-3 ref/E2.repeats.bed ref/E2.named_TEs.bed > ref/E2.all_repeats.bed
for al in bwa bowtie2; do
  python scripts/03_call_depth_deletions.py --aligner $al --tag $al.uniq --out calls_uniq \
     --exclude-bed ref/E2.all_repeats.bed --min-present-run 25 --threads 10 > log/call_uniq_$al.txt
done

awk 'BEGIN{OFS="\t"; print "id","chrom","start","end"} NR>1{k=$3":"$4"-"$5; if(!(k in s)){s[k]=1; print k,$3,$4,$5}}' \
  calls_uniq/depth_calls.bwa.uniq.primary.tsv > calls_uniq/depth_calls.intervals.tsv
for f in calls_uniq/depth_calls.intervals.tsv inputs/*.intervals.tsv; do
  b=$(basename "$f" .intervals.tsv)
  python scripts/07_interval_depth_stats.py --intervals "$f" --scorable-bed calls_uniq/scorable.bwa.uniq.bed \
     --out calls_uniq/$b.per_strain_depth.tsv.gz --threads 10
done
echo done