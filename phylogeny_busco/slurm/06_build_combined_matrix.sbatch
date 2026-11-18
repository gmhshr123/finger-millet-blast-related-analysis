#!/bin/bash
#SBATCH --job-name=fm_combined_matrix
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.err

set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
VCF="${VCF:-$PROJECT/reanalysis_deepvariant/filtered_vcf/finger_millet.final_read_backed.deepvariant.glnexus.biallelic_snps.MAFge0.05.missing_le0.30.vcf.gz}"
ASSEMBLY_MANIFEST="${ASSEMBLY_MANIFEST:-$PROJECT/reanalysis_inputs/manifests/final_assembly_backed_samples.tsv}"
REF="${REF:-$PROJECT/reanalysis_inputs/reference/E2.single_copy.fasta}"
OUTDIR="${OUTDIR:-$PROJECT/results/reanalysis_phylogeny_inputs/read90_assembly9_MAFge0.05_missing30}"
PREFIX="${PREFIX:-finger_millet.read90_assembly9.MAFge0.05_missing30}"

mkdir -p "$OUTDIR" "$PROJECT/log"

if command -v module >/dev/null 2>&1; then
  module load SAMtools/1.21-GCC-13.3.0 || module load SAMtools || true
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools is not available after module setup" >&2
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is not available" >&2
  exit 1
fi

[[ -s "$VCF" && -s "$VCF.tbi" ]] || { echo "ERROR: missing VCF or index: $VCF" >&2; exit 1; }
[[ -s "$ASSEMBLY_MANIFEST" ]] || { echo "ERROR: missing assembly manifest: $ASSEMBLY_MANIFEST" >&2; exit 1; }
[[ -s "$REF" && -s "$REF.fai" ]] || { echo "ERROR: missing reference or fai: $REF" >&2; exit 1; }
[[ -s "$PROJECT/scripts/build_combined_read_assembly_snp_matrix.py" ]] || { echo "ERROR: missing matrix builder script" >&2; exit 1; }

echo "Started: $(date)"
echo "VCF: $VCF"
echo "Assembly manifest: $ASSEMBLY_MANIFEST"
echo "OUTDIR: $OUTDIR"

python3 "$PROJECT/scripts/build_combined_read_assembly_snp_matrix.py" \
  --vcf "$VCF" \
  --assembly-manifest "$ASSEMBLY_MANIFEST" \
  --reference "$REF" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"

echo "Finished: $(date)"
cat "$OUTDIR/$PREFIX.summary.tsv"
