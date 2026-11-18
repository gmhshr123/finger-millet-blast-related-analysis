#!/bin/bash
#SBATCH --job-name=del_call
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

# All-base deletion calls for the 38 strains (BWA realignment = primary; original Bowtie2
# alignment = check that the result does not depend on the aligner), plus export of the
# Delly junction calls and scoring of the SNP-based candidates against the depth data.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a BCFtools/1.21-GCC-13.3.0

for al in bwa bowtie2; do
  python scripts/03_call_depth_deletions.py --aligner $al --threads 10 > log/call_$al.txt
done

# Delly calls -> one table (deletions only). DV/RV = read pairs / split reads supporting the variant.
mkdir -p sv_tsv
for s in $(cat meta/copy38.txt); do
  bcftools query -i 'SVTYPE="DEL"' \
    -f "$s\t%CHROM\t%POS\t%INFO/END\t%FILTER\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/PRECISE\t%INFO/HOMLEN\t[%GT\t%DR\t%DV\t%RR\t%RV]\n" \
    sv/$s.delly.bcf
done > sv_tsv/delly_deletions.tsv

# score supplied interval sets (SNP-gap runs, Data S10 rows) and the depth calls themselves in every strain
awk 'BEGIN{OFS="\t"; print "id","chrom","start","end"} NR>1{print "dc"NR-1,$3,$4,$5}' calls/depth_calls.bwa.primary.tsv > calls/depth_calls.bwa.primary.intervals.tsv
for f in calls/depth_calls.bwa.primary.intervals.tsv inputs/*.intervals.tsv; do
  [ -s "$f" ] || continue
  b=$(basename "$f" .intervals.tsv)
  python scripts/07_interval_depth_stats.py --intervals "$f" --out calls/$b.per_strain_depth.tsv.gz --threads 10
done
cat qc/*.depth_summary.tsv > calls/depth_summary_all.tsv
echo done
