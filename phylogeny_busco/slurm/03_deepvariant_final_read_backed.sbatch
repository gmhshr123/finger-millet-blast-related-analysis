#!/bin/bash
#SBATCH --job-name=fm_dv_final
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.err

set -eo pipefail

set +u
source ~/.bashrc || true
set -u
module load SAMtools/1.21-GCC-13.3.0 || true

PROJECT=/path/to/scratch/finger_millet_blast_phylogeny
MANIFEST="$PROJECT/reanalysis_inputs/manifests/final_read_backed_samples.tsv"
REF="$PROJECT/reanalysis_inputs/reference/E2.single_copy.fasta"
DV_SIF="/path/to/home/deepvariant/deepvariant_1.9.0.sif"

mkdir -p "$PROJECT"/reanalysis_deepvariant/{gvcf,vcf,qc/deepvariant,tmp/deepvariant_cache,tmp/deepvariant_work,log}

line=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"; OFS="\034"} NR == id + 1 {for (i=1; i<=12; i++) printf "%s%s", $i, (i < 12 ? OFS : ORS); exit}' "$MANIFEST")
line=${line//$'\r'/}
[[ -n "$line" ]] || { echo "ERROR: no final read-backed manifest row for array id $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

IFS=$'\034' read -r sample input_kind mapper model bam out_gvcf out_vcf coverage_fraction mean_depth numreads depth_flag decision <<< "$line"
logdir="$PROJECT/reanalysis_deepvariant/qc/deepvariant/${sample}"

[[ "$decision" == "include_read_backed" ]] || { echo "ERROR: sample $sample is not marked include_read_backed" >&2; exit 1; }
[[ -n "$model" && -n "$bam" ]] || { echo "ERROR: missing model or BAM for $sample" >&2; exit 1; }
[[ -s "$REF" && -s "$REF.fai" ]] || { echo "ERROR: missing reference or fai: $REF" >&2; exit 1; }
[[ -s "$bam" ]] || { echo "ERROR: BAM missing: $bam" >&2; exit 1; }
[[ -s "${bam}.bai" ]] || samtools index -@ "$SLURM_CPUS_PER_TASK" "$bam"
[[ -s "$DV_SIF" ]] || { echo "ERROR: missing DeepVariant SIF: $DV_SIF" >&2; exit 1; }

if command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
elif command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
else
  echo "ERROR: neither apptainer nor singularity found on this node" >&2
  exit 1
fi

mkdir -p "$logdir"
export APPTAINER_CACHEDIR="$PROJECT/reanalysis_deepvariant/tmp/deepvariant_cache"
export APPTAINER_TMPDIR="$PROJECT/reanalysis_deepvariant/tmp/deepvariant_cache"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

work_tmp="$PROJECT/reanalysis_deepvariant/tmp/deepvariant_work/${sample}_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
intermediate_dir="$work_tmp/intermediate"
rm -rf "$work_tmp"
mkdir -p "$intermediate_dir"
cleanup_tmp() {
  rm -rf "$work_tmp"
}
trap cleanup_tmp EXIT

export TMPDIR="$work_tmp"
export TMP="$work_tmp"
export TEMP="$work_tmp"
export PARALLEL_TMPDIR="$work_tmp"
export APPTAINERENV_TMPDIR="$work_tmp"
export APPTAINERENV_TMP="$work_tmp"
export APPTAINERENV_TEMP="$work_tmp"
export APPTAINERENV_PARALLEL_TMPDIR="$work_tmp"
export SINGULARITYENV_TMPDIR="$work_tmp"
export SINGULARITYENV_TMP="$work_tmp"
export SINGULARITYENV_TEMP="$work_tmp"
export SINGULARITYENV_PARALLEL_TMPDIR="$work_tmp"

DV_NUM_SHARDS="${DV_NUM_SHARDS:-1}"
DEFAULT_EXTRA_ARGS='min_mapping_quality=1'
EXTRA_ARGS="${DV_EXTRA_ARGS:-$DEFAULT_EXTRA_ARGS}"
DV_DISABLE_SMALL_MODEL="${DV_DISABLE_SMALL_MODEL:-true}"

if [[ "${FORCE_RERUN:-0}" != "1" && -s "$out_vcf" && -s "${out_vcf}.tbi" && -s "$out_gvcf" && -s "${out_gvcf}.tbi" ]]; then
  echo "[$(date)] SKIP existing DeepVariant outputs for sample: $sample"
  exit 0
fi

echo "[$(date)] DeepVariant final read-backed sample=$sample input_kind=$input_kind model=$model mean_depth=${mean_depth} depth_flag=${depth_flag} runtime=$RUNTIME shards=$DV_NUM_SHARDS disable_small_model=$DV_DISABLE_SMALL_MODEL extra_args=$EXTRA_ARGS"
rm -f "$out_vcf" "${out_vcf}.tbi" "$out_gvcf" "${out_gvcf}.tbi"
"$RUNTIME" exec --bind "$PROJECT:$PROJECT" "$DV_SIF" run_deepvariant \
  --model_type="$model" \
  --ref="$REF" \
  --reads="$bam" \
  --output_vcf="$out_vcf" \
  --output_gvcf="$out_gvcf" \
  --intermediate_results_dir="$intermediate_dir" \
  --make_examples_extra_args="$EXTRA_ARGS" \
  --disable_small_model="$DV_DISABLE_SMALL_MODEL" \
  --num_shards="$DV_NUM_SHARDS" \
  --logging_dir="$logdir"

echo "[$(date)] DONE DeepVariant final read-backed sample: $sample"
