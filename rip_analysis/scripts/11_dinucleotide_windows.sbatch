#!/bin/bash
#SBATCH --job-name=rip_dinuc
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24GB
#SBATCH --time=6-22:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/rip_analysis
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/rip_analysis/log/%x_%j.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/rip_analysis/log/%x_%j.err
# 1. All 16 dinucleotides in sliding windows of the four genomes (The RIPper's window, then RIPCaL's scan window).
# 2. RIPCaL's dinucleotide-frequency mode (-type index) on the copies of every TE group.
set -euo pipefail
module load SciPy-bundle/2024.05-gfbf-2024a
E2=/path/to/scratch/finger_millet_blast_phylogeny/deletion_allbase/ref/E2.fasta
G="--genome E2=$E2 --genome MZ5=genomes/MZ5.renamed.fa --genome B71=genomes/B71.renamed.fa --genome PoP131=genomes/PoP131.renamed.fa"
python scripts/11_dinucleotide_windows.py $G --out results_dinuc --window 1000 --step 500
python scripts/11_dinucleotide_windows.py $G --out results_dinuc/scan300 --window 300 --step 150 --ripcal-scan

B=$PWD/ripcal
PERL="$B/env/bin/perl -I$B/stub"
out=results_dinuc/ripcal_dinuc_groups.tsv
: > "$out"
for f in "$B"/groups/*.fa; do
  g=$(basename "$f" .fa)
  d="$B/runs/$g/index"; mkdir -p "$d"
  sed '/^>/! s/-//g' "$f" > "$d/input.fa"            # RIPCaL's dinucleotide mode takes unaligned sequences
  ( cd "$d" && rm -f ./*_dinuc.txt && $PERL "$B/src/ripcal2_install/perl/ripcal" -c -seq input.fa -type index > ripcal.log 2>&1 ) || echo "FAILED $g"
  t=$(ls "$d"/*_dinuc.txt | head -1)
  [ -s "$out" ] || { printf 'group\t'; head -1 "$t" | cut -f2-; } > "$out"
  printf '%s\t' "$g" >> "$out"; sed -n 2p "$t" | cut -f2- >> "$out"
done
wc -l "$out"
