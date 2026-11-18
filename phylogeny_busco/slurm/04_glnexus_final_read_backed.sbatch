#!/bin/bash
#SBATCH --job-name=fm_glnexus_final
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%j.err

set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
LIST="${LIST:-$PROJECT/reanalysis_deepvariant/gvcf/final_read_backed_gvcf_list.txt}"
GLNEXUS_SIF="${GLNEXUS_SIF:-/path/to/home/glnexus/glnexus_v1.4.1.sif}"
OUT_DIR="$PROJECT/reanalysis_deepvariant/joint_vcf"
QC_DIR="$PROJECT/reanalysis_deepvariant/qc/glnexus"
TMP_BASE="$PROJECT/reanalysis_deepvariant/tmp/glnexus/${SLURM_JOB_ID:-manual}"
RUN_DIR="$TMP_BASE/work"

BCF_OUT="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.bcf"
VCF_OUT="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.vcf.gz"
BCF_TMP="$TMP_BASE/finger_millet.final_read_backed.deepvariant.glnexus.bcf"
VCF_TMP="$TMP_BASE/finger_millet.final_read_backed.deepvariant.glnexus.vcf.gz"
SAMPLES_OUT="$QC_DIR/final_read_backed.glnexus.samples.txt"
COUNTS_OUT="$QC_DIR/final_read_backed.glnexus.counts.txt"

mkdir -p "$OUT_DIR" "$QC_DIR" "$TMP_BASE" "$RUN_DIR" "$PROJECT/log"

if command -v module >/dev/null 2>&1; then
  module load BCFtools/1.21-GCC-13.3.0 || module load BCFtools || true
fi

if command -v apptainer >/dev/null 2>&1; then
  RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  RUNTIME="singularity"
else
  echo "ERROR: neither apptainer nor singularity is available" >&2
  exit 1
fi

if [[ ! -s "$LIST" ]]; then
  echo "ERROR: gVCF list is missing or empty: $LIST" >&2
  exit 1
fi

if [[ ! -s "$GLNEXUS_SIF" ]]; then
  echo "ERROR: GLnexus container is missing or empty: $GLNEXUS_SIF" >&2
  exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
  echo "ERROR: bcftools is not available after module setup" >&2
  exit 1
fi

mapfile -t GVCFS < "$LIST"
if [[ "${#GVCFS[@]}" -eq 0 ]]; then
  echo "ERROR: no gVCFs found in $LIST" >&2
  exit 1
fi

missing=0
for gvcf in "${GVCFS[@]}"; do
  if [[ ! -s "$gvcf" || ! -s "$gvcf.tbi" ]]; then
    echo "MISSING: $gvcf or $gvcf.tbi" >&2
    missing=$((missing + 1))
  fi
done

if [[ "$missing" -ne 0 ]]; then
  echo "ERROR: $missing gVCF inputs are missing data or indexes" >&2
  exit 1
fi

export TMPDIR="$TMP_BASE"
export APPTAINER_CACHEDIR="$TMP_BASE/cache"
export APPTAINER_TMPDIR="$TMP_BASE/tmp"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Runtime: $RUNTIME"
echo "gVCFs: ${#GVCFS[@]}"
echo "BCF: $BCF_OUT"
echo "VCF: $VCF_OUT"

pushd "$RUN_DIR" >/dev/null
"$RUNTIME" exec \
  --bind "$PROJECT:$PROJECT" \
  --bind "$TMP_BASE:$TMP_BASE" \
  "$GLNEXUS_SIF" \
  glnexus_cli \
  --config DeepVariantWGS \
  -t "${SLURM_CPUS_PER_TASK:-32}" \
  "${GVCFS[@]}" \
  > "$BCF_TMP"
popd >/dev/null

bcftools index --threads "${SLURM_CPUS_PER_TASK:-32}" "$BCF_TMP"
bcftools view --threads "${SLURM_CPUS_PER_TASK:-32}" -Oz -o "$VCF_TMP" "$BCF_TMP"
bcftools index --threads "${SLURM_CPUS_PER_TASK:-32}" -t "$VCF_TMP"

mv -f "$BCF_TMP" "$BCF_OUT"
mv -f "$BCF_TMP.csi" "$BCF_OUT.csi"
mv -f "$VCF_TMP" "$VCF_OUT"
mv -f "$VCF_TMP.tbi" "$VCF_OUT.tbi"

bcftools query -l "$VCF_OUT" > "$SAMPLES_OUT"
{
  echo -e "metric\tvalue"
  echo -e "samples\t$(wc -l < "$SAMPLES_OUT")"
  echo -e "records\t$(bcftools index -n "$VCF_OUT")"
} > "$COUNTS_OUT"

echo "Finished: $(date)"
cat "$COUNTS_OUT"
