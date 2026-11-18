#!/bin/bash
#SBATCH --job-name=fm_prep_ref
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.err

set -eo pipefail

set +u
source ~/.bashrc || true
set -u
module load SAMtools/1.21-GCC-13.3.0 || true

PROJECT=/path/to/scratch/finger_millet_blast_phylogeny
REF="$PROJECT/reanalysis_inputs/reference/E2.single_copy.fasta"

mkdir -p "$PROJECT"/{log,scripts,reanalysis_inputs/manifests,reanalysis_deepvariant/{bam,gvcf,vcf,qc,tmp},reanalysis_assembly_alignments}

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }
command -v bwa-mem2 >/dev/null 2>&1 || { echo "ERROR: bwa-mem2 not found after loading ~/.bashrc" >&2; exit 1; }

[[ -s "$REF" ]] || { echo "ERROR: missing reference: $REF" >&2; exit 1; }

samtools faidx "$REF"

DICT="${REF%.fasta}.dict"
if [[ ! -s "$DICT" ]]; then
  samtools dict -o "$DICT" "$REF"
fi

if [[ ! -s "${REF}.0123" ]]; then
  bwa-mem2 index "$REF"
fi

python3 "$PROJECT/scripts/build_reanalysis_mapping_manifest.py" \
  --base "$PROJECT" \
  --out "$PROJECT/reanalysis_inputs/manifests/mapping_manifest.tsv"

echo "[$(date)] Reference prep and manifest complete"
wc -l "$PROJECT/reanalysis_inputs/manifests/mapping_manifest.tsv"
