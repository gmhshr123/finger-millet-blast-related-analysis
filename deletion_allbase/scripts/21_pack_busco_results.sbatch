#!/bin/bash
#SBATCH --job-name=pack_busco
#SBATCH --partition=batch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=6-22:00:00
#SBATCH --array=1-6
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.edu
#SBATCH --chdir=/path/to/scratch/finger_millet_blast_phylogeny/results
#SBATCH --output=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.out
#SBATCH --error=/path/to/scratch/finger_millet_blast_phylogeny/log/%x_%A_%a.err

# Pack one BUSCO result folder (millions of small files) into a single .tar.gz, to get the
# scratch file count back under quota. Steps, in this order:
#   1. copy the small summary files (short_summary*, full_table.tsv, missing_busco_list.tsv and
#      any loose files at the top of the folder) to <folder>_summaries/, so they stay readable;
#   2. write busco_archives/<folder>.tar.gz;
#   3. check the archive: gzip integrity, and the list of files (with sizes) and of symbolic
#      links inside the archive must be identical to what is on disk;
#   4. only if every check passes, remove the original folder.
# Safe to re-run: a folder that is already packed and removed is skipped, and a half-written
# archive (job killed) is thrown away and rebuilt.
# To unpack again:  cd results && tar -I pigz -xf busco_archives/<folder>.tar.gz
set -euo pipefail
FOLDERS=(busco_226_validation_20260601 busco_final_ragtag_ascomycota_odb12 busco_scafSeq
         busco_chro_final_Mum_v2_ascomycota_odb12 busco_226strains busco_finalAssemblings_ascomycota_odb12)
# FOLDER / ARCHIVE_DIR / KEEP_ORIGINAL=1 can be set by hand to test the script on another folder
F=${FOLDER:-${FOLDERS[$((${SLURM_ARRAY_TASK_ID:-1} - 1))]}}
A=${ARCHIVE_DIR:-busco_archives}
mkdir -p $A
echo "[$(date)] folder $F"

if [ ! -d "$F" ]; then
  [ -f "$A/$F.verified" ] && { echo "already packed and removed"; exit 0; }
  echo "folder not found and no verified archive: nothing done"; exit 1
fi

# 1. keep the summaries readable
mkdir -p "${F}_summaries"
{ find "$F" -type f \( -name 'short_summary*' -o -name 'full_table.tsv' -o -name 'missing_busco_list.tsv' \)
  find "$F" -maxdepth 1 -type f; } | sort -u > "$A/$F.summary_files.txt"
rsync -a --files-from="$A/$F.summary_files.txt" . "${F}_summaries/"
echo "[$(date)] summaries copied: $(wc -l < "$A/$F.summary_files.txt") files"

# 2. archive (a leftover archive without the .verified marker is incomplete: rebuild it)
if [ ! -f "$A/$F.verified" ]; then
  rm -f "$A/$F.tar.gz"
  tar -cf - "$F" | pigz -p "$SLURM_CPUS_PER_TASK" > "$A/$F.tar.gz.part"
  mv "$A/$F.tar.gz.part" "$A/$F.tar.gz"
  echo "[$(date)] archive written: $(du -h "$A/$F.tar.gz" | cut -f1)"

  # 3. verify
  pigz -t "$A/$F.tar.gz"
  find "$F" -type f -printf '%p\t%s\n' | LC_ALL=C sort > "$A/$F.disk_files.txt"
  find "$F" -type l -printf '%p\n'     | LC_ALL=C sort > "$A/$F.disk_links.txt"
  tar -I pigz -tvf "$A/$F.tar.gz" > "$A/$F.tar_listing.txt"
  awk '$1 ~ /^-/ {print $6 "\t" $3}' "$A/$F.tar_listing.txt" | LC_ALL=C sort > "$A/$F.tar_files.txt"
  awk '$1 ~ /^l/ {print $6}'         "$A/$F.tar_listing.txt" | LC_ALL=C sort > "$A/$F.tar_links.txt"
  nd=$(wc -l < "$A/$F.disk_files.txt"); nt=$(wc -l < "$A/$F.tar_files.txt")
  echo "files on disk: $nd   files in archive: $nt   links on disk: $(wc -l < "$A/$F.disk_links.txt")"
  cmp "$A/$F.disk_files.txt" "$A/$F.tar_files.txt" || { echo "FILE LIST MISMATCH - original folder kept"; exit 1; }
  cmp "$A/$F.disk_links.txt" "$A/$F.tar_links.txt" || { echo "LINK LIST MISMATCH - original folder kept"; exit 1; }
  [ "$nd" -gt 0 ] || { echo "no files found - original folder kept"; exit 1; }
  echo "$nd files, checked $(date)" > "$A/$F.verified"
  rm -f "$A/$F.tar_listing.txt" "$A/$F.tar_files.txt" "$A/$F.tar_links.txt" "$A/$F.disk_links.txt"
  pigz -f "$A/$F.disk_files.txt"          # kept as the table of contents of the archive
  echo "[$(date)] archive verified"
fi

# 4. remove the original
[ "${KEEP_ORIGINAL:-0}" = 1 ] && { echo "KEEP_ORIGINAL=1: original folder kept"; exit 0; }
rm -rf "$F"
echo "[$(date)] original folder removed; finished $F"
