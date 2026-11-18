#!/bin/bash
#SBATCH --job-name=del_call226
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%j.err

# All-base deletion calls for all 226 Data S1 strains, with exactly the rules of the 33-strain
# analysis (06b): BWA realignment, unique reads, repeats bridged, depth <= 10% of the strain's median.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a

# stop early, with a clear message, if any strain has no finished depth track
missing=0
while read -r S; do
  for q in q0 q20; do
    [ -s depth/$S.bwa.$q.bga.gz ] || { echo "missing depth track: $S $q"; missing=1; }
  done
done < meta/all226.txt
[ "$missing" -eq 0 ] || exit 1

mkdir -p calls_226 export_226
python scripts/03_call_depth_deletions.py --aligner bwa --tag bwa.uniq.all226 --out calls_226 --samples meta/all226.txt \
   --exclude-bed ref/E2.all_repeats.bed --min-present-run 25 --threads 10 > log/call_226.txt
cp calls_226/depth_calls.bwa.uniq.all226.primary.tsv calls_226/windows100kb.bwa.uniq.all226.tsv.gz \
   calls_226/sensitivity.bwa.uniq.all226.tsv calls_226/callable.bwa.uniq.all226.summary.txt export_226/
cat qc/*.depth_summary.tsv > export_226/depth_summary_all.tsv
echo "[$(date)] finished"
