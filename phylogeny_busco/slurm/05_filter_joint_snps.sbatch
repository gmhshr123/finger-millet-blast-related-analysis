#!/bin/bash
#SBATCH --job-name=fm_filter_snps
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

set -euo pipefail

PROJECT="${PROJECT:-/path/to/scratch/finger_millet_blast_phylogeny}"
IN="$PROJECT/reanalysis_deepvariant/joint_vcf/finger_millet.final_read_backed.deepvariant.glnexus.vcf.gz"
OUT_DIR="$PROJECT/reanalysis_deepvariant/filtered_vcf"
QC_DIR="$PROJECT/reanalysis_deepvariant/qc/filter_snps"

TAGGED_ALL="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.all_sites.with_AF_MAF_FMISS.vcf.gz"
SNPS_ALL="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.biallelic_snps.vcf.gz"
SNPS_MISS30="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.biallelic_snps.missing_le0.30.vcf.gz"
SNPS_MAF005_MISS20="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.biallelic_snps.MAFgt0.005.missing_le0.20.vcf.gz"
SNPS_MAF03_MISS30="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.biallelic_snps.MAFge0.03.missing_le0.30.vcf.gz"
SNPS_MAF05_MISS30="$OUT_DIR/finger_millet.final_read_backed.deepvariant.glnexus.biallelic_snps.MAFge0.05.missing_le0.30.vcf.gz"

mkdir -p "$OUT_DIR" "$QC_DIR" "$PROJECT/log"

if command -v module >/dev/null 2>&1; then
  module load BCFtools/1.21-GCC-13.3.0 || module load BCFtools || true
fi

if ! command -v bcftools >/dev/null 2>&1; then
  echo "ERROR: bcftools is not available after module setup" >&2
  exit 1
fi

[[ -s "$IN" && -s "$IN.tbi" ]] || { echo "ERROR: missing input VCF or index: $IN" >&2; exit 1; }

sample_qc() {
  local vcf="$1"
  local out="$2"
  bcftools query -f '[%SAMPLE\t%GT\t%DP\n]' "$vcf" |
    awk 'BEGIN{FS=OFS="\t"; print "sample","sites","called","missing","missing_pct","nonref","mean_dp_at_called_sites"}
      {
        sample=$1; gt=$2; dp=$3
        sites[sample]++
        if (gt == "." || gt == "./." || gt == ".|." || gt ~ /\./) {
          missing[sample]++
          next
        }
        called[sample]++
        if (gt ~ /[1-9]/) nonref[sample]++
        if (dp != "." && dp != "") {
          dp_sum[sample] += dp
          dp_n[sample]++
        }
      }
      END {
        for (sample in sites) {
          miss_pct = sites[sample] ? 100 * missing[sample] / sites[sample] : 0
          mean_dp = dp_n[sample] ? dp_sum[sample] / dp_n[sample] : 0
          printf "%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\n", sample, sites[sample], called[sample], missing[sample], miss_pct, nonref[sample], mean_dp
        }
      }' |
    sort -k1,1 > "$out"
}

echo "Started: $(date)"
echo "Input: $IN"

bcftools +fill-tags "$IN" --threads "${SLURM_CPUS_PER_TASK:-8}" -Oz -o "$TAGGED_ALL" -- -t AC,AN,AF,MAF,F_MISSING
bcftools index -t -f "$TAGGED_ALL"

# GLnexus emits unfiltered records as FILTER=".", so accept both "." and PASS.
bcftools view --threads "${SLURM_CPUS_PER_TASK:-8}" \
  -m2 -M2 -v snps -f '.,PASS' \
  -i 'REF~"^[ACGT]$" && ALT~"^[ACGT]$"' \
  "$TAGGED_ALL" -Oz -o "$SNPS_ALL"
bcftools index -t -f "$SNPS_ALL"

bcftools view --threads "${SLURM_CPUS_PER_TASK:-8}" \
  -i 'F_MISSING<=0.30' \
  "$SNPS_ALL" -Oz -o "$SNPS_MISS30"
bcftools index -t -f "$SNPS_MISS30"

bcftools view --threads "${SLURM_CPUS_PER_TASK:-8}" \
  -i 'MAF>0.005 && F_MISSING<=0.20' \
  "$SNPS_ALL" -Oz -o "$SNPS_MAF005_MISS20"
bcftools index -t -f "$SNPS_MAF005_MISS20"

bcftools view --threads "${SLURM_CPUS_PER_TASK:-8}" \
  -i 'MAF>=0.03 && F_MISSING<=0.30' \
  "$SNPS_ALL" -Oz -o "$SNPS_MAF03_MISS30"
bcftools index -t -f "$SNPS_MAF03_MISS30"

bcftools view --threads "${SLURM_CPUS_PER_TASK:-8}" \
  -i 'MAF>=0.05 && F_MISSING<=0.30' \
  "$SNPS_ALL" -Oz -o "$SNPS_MAF05_MISS30"
bcftools index -t -f "$SNPS_MAF05_MISS30"

{
  echo -e "file\trecords"
  for f in "$IN" "$TAGGED_ALL" "$SNPS_ALL" "$SNPS_MISS30" "$SNPS_MAF005_MISS20" "$SNPS_MAF03_MISS30" "$SNPS_MAF05_MISS30"; do
    echo -e "$(basename "$f")\t$(bcftools index -n "$f")"
  done
} > "$QC_DIR/filter_counts.tsv"

sample_qc "$SNPS_ALL" "$QC_DIR/sample_qc.biallelic_snps.tsv"
sample_qc "$SNPS_MISS30" "$QC_DIR/sample_qc.biallelic_snps.missing_le0.30.tsv"
sample_qc "$SNPS_MAF005_MISS20" "$QC_DIR/sample_qc.biallelic_snps.MAFgt0.005.missing_le0.20.tsv"
sample_qc "$SNPS_MAF03_MISS30" "$QC_DIR/sample_qc.biallelic_snps.MAFge0.03.missing_le0.30.tsv"
sample_qc "$SNPS_MAF05_MISS30" "$QC_DIR/sample_qc.biallelic_snps.MAFge0.05.missing_le0.30.tsv"

echo "Finished: $(date)"
cat "$QC_DIR/filter_counts.tsv"
