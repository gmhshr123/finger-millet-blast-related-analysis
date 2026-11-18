#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/scratch/finger_millet_blast_phylogeny}"

download_candidate() {
  local sample="$1"
  local accession_prefix="$2"
  local url_base="$3"
  local biosample="$4"
  local host="$5"
  local out="$BASE/reanalysis_inputs/samples/$sample/assembly"

  mkdir -p "$out"
  cd "$out"
  wget -c -nv "$url_base/${accession_prefix}_genomic.fna.gz"
  wget -c -nv "$url_base/${accession_prefix}_assembly_report.txt"
  wget -c -nv "$url_base/${accession_prefix}_assembly_stats.txt"
  wget -c -nv "$url_base/md5checksums.txt"

  : > md5checksums.selected.txt
  for f in "./${accession_prefix}_genomic.fna.gz" "./${accession_prefix}_assembly_report.txt" "./${accession_prefix}_assembly_stats.txt"; do
    grep -F "$f" md5checksums.txt >> md5checksums.selected.txt
  done
  md5sum -c md5checksums.selected.txt
  gzip -t "${accession_prefix}_genomic.fna.gz"
  touch "${accession_prefix}_genomic.fna.gz.done"

  {
    printf "field\tvalue\n"
    printf "sample\t%s\n" "$sample"
    printf "role\tcandidate_crabgrass_outgroup\n"
    printf "organism\tPyricularia grisea / Magnaporthe grisea\n"
    printf "host\t%s\n" "$host"
    printf "host_common_name\tcrabgrass\n"
    printf "biosample\t%s\n" "$biosample"
    printf "source_url\t%s\n" "$url_base"
    printf "input_type\tassembly_fasta\n"
  } > "$BASE/reanalysis_inputs/samples/$sample/sample_info.tsv"
}

download_candidate "DsLIZ" \
  "GCA_002925245.1_ASM292524v1" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/925/245/GCA_002925245.1_ASM292524v1" \
  "SAMN08009550" \
  "Digitaria sanguinalis"

download_candidate "VO107" \
  "GCA_002924675.1_ASM292467v1" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/924/675/GCA_002924675.1_ASM292467v1" \
  "SAMN08009577" \
  "Digitaria sanguinalis"

download_candidate "DS0505" \
  "GCA_001548815.1_ASM154881v1" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/548/815/GCA_001548815.1_ASM154881v1" \
  "SAMN04299554" \
  "Digitaria sanguinalis"

download_candidate "DS9461" \
  "GCA_001548795.1_DS9461" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/548/795/GCA_001548795.1_DS9461" \
  "SAMN04318446" \
  "Digitaria sanguinalis"

download_candidate "W98-20" \
  "GCA_002866905.1_Mg_W98-20_1.0" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/866/905/GCA_002866905.1_Mg_W98-20_1.0" \
  "SAMN07818931" \
  "Digitaria sanguinalis"
