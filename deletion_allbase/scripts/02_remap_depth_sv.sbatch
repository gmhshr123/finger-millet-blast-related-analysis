#!/bin/bash
#SBATCH --job-name=del_remap
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/log/%x_%A_%a.err

# Per strain:
#  1. turn the original Bowtie2 BAM back into read pairs (all reads, mapped and unmapped)
#  2. realign with BWA-MEM, which reports split reads (soft clips / supplementary alignments)
#  3. write all-base depth tracks (all reads, and uniquely mapped reads with MAPQ>=20)
#     for both the new BWA alignment and the original Bowtie2 alignment
#  4. call structural variants with Delly (junction evidence from split reads and read pairs)
set -euo pipefail
module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.21-GCC-13.3.0 BEDTools/2.31.1-GCC-13.3.0 Delly/1.3.3-GCC-13.3.0 HTSlib/1.21-GCC-13.3.0

LIST=${LIST:-meta/copy38.txt}
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
T=${SLURM_CPUS_PER_TASK:-16}
REF=ref/E2.fasta
IN=bam_orig/$S.Gr.sorted.bam
mkdir -p bam_bwa depth qc sv tmp/$S
echo "[$(date)] sample $S"

OUT=bam_bwa/$S.bwa.markdup.bam
if [ ! -s "$OUT.bai" ]; then
  samtools collate -@4 -u -O "$IN" tmp/$S/collate \
    | samtools fastq -n - 2> qc/$S.fastq_export.log \
    | bwa mem -p -t "$T" -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA\tLB:$S" "$REF" - 2> qc/$S.bwa.log \
    | samtools fixmate -m -u - - \
    | samtools sort -@4 -m 2G -T tmp/$S/sort -u - \
    | samtools markdup -@4 -s -f qc/$S.markdup.txt - "$OUT"
  samtools index -@4 "$OUT"
fi
samtools flagstat -@4 "$IN"  > qc/$S.bowtie2.flagstat.txt
samtools flagstat -@4 "$OUT" > qc/$S.bwa.flagstat.txt

# all-base depth as run-length tracks (bedGraph with zero-depth runs kept).
# -F 0x704: drop unmapped, secondary, QC-fail and duplicate reads; -split: do not count
# bases inside an alignment gap (CIGAR D/N) as covered.
depth_track () {  # bam  minMAPQ  out
  samtools view -u -q "$2" -F 0x704 "$1" \
    | bedtools genomecov -ibam stdin -bga -split \
    | bgzip -@4 > "$3"
}
depth_track "$OUT" 0  depth/$S.bwa.q0.bga.gz
depth_track "$OUT" 20 depth/$S.bwa.q20.bga.gz
depth_track "$IN"  0  depth/$S.bowtie2.q0.bga.gz
depth_track "$IN"  20 depth/$S.bowtie2.q20.bga.gz

# genome-wide depth summary on the seven chromosomes (length-weighted median and mean)
for trk in bwa.q0 bwa.q20 bowtie2.q0 bowtie2.q20; do
  zcat depth/$S.$trk.bga.gz | awk -v s="$S" -v t="$trk" '$1 ~ /^Chr/ {l=$3-$2; h[$4]+=l; n+=l; sum+=l*$4; if($4==0) z+=l}
    END{half=n/2; c=0; m=asorti(h,k,"@ind_num_asc"); for(i=1;i<=m;i++){c+=h[k[i]]; if(c>=half){med=k[i]; break}}
        printf "%s\t%s\t%d\t%.2f\t%d\t%d\n", s,t,n,sum/n,med,z}'
done > qc/$S.depth_summary.tsv

delly call -g "$REF" -o sv/$S.delly.bcf "$OUT" > qc/$S.delly.log 2>&1

rm -rf tmp/$S
echo "[$(date)] finished $S"
