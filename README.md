# Finger millet blast: related analyses

Scripts for analyses of finger millet-infecting *Magnaporthe oryzae* (MoE) strains from eastern Africa. The sequencing data are in NCBI BioProject
[PRJNA1279063](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1279063). The reference is the E2 genome
assembly (GenBank GCA_059469275.1).

This repository contains the analysis and plotting code, and the interactive deletion figures. No sequence
data or alignments are included.

## Interactive figures

The interactive figures and the deletion database are in `docs/` and can be viewed at
https://gmhshr123.github.io/finger-millet-blast-related-analysis/

| Page | Content |
|---|---|
| `heatmap_226_strains_all_base.html` | All-base deletions in 226 strains, per window |
| `heatmap_226_strains_SNP_based.html` | SNP-based deletions in 226 strains, per window |
| `deletion_tracks_33_strains.html` | Deletion tracks and genome features along E2 in 33 Ethiopian strains |
| `deletion_database.html` | All all-base deletions in 33 Ethiopian strains, with breakpoints and read evidence |

## Contents

| Folder | Analysis |
|---|---|
| `phylogeny_busco/` | Phylogeny of MoE and other *M. oryzae* lineages, and BUSCO completeness of the assemblies |
| `deletion_allbase/` | Deletions detected from read depth at all bases, compared with the SNP-based deletion calls; transposable elements (TEs) at deletion breakpoints |
| `rip_analysis/` | Repeat-induced point mutation (RIP) in full-length TE copies and genome-wide RIP indices |
| `fosbury_ltr/` | Identity of the two long terminal repeats (LTRs), target-site duplications and LTR indels of full-length Fosbury elements |

### `phylogeny_busco/`

Reads were mapped to E2 and variants were called per sample with DeepVariant and jointly with GLnexus.
SNPs were filtered and combined with SNPs projected from public assemblies of strains without
reads. Maximum-likelihood trees were inferred with IQ-TREE from the combined SNP matrix, with
sensitivity trees after removing low-depth samples. BUSCO was run on the assemblies of the 226 strains
and on other assemblies, and the pattern of missing BUSCO genes was analysed. `scripts/` holds the
Python and shell scripts, `slurm/` the job files in run order.

### `deletion_allbase/`

Scripts are numbered in run order. The main steps are:

1. Realignment of each strain with BWA-MEM, per-base depth of all and of uniquely mapped reads
   (mapping quality 20 or higher), and Delly deletion calls (`01`, `02`).
2. Deletion calls from depth: a callable, non-repeat base is called deleted when its depth of
   uniquely mapped reads is at most 5% of the strain median; calls are merged within a strain
   (1 kb) and across strains (80% reciprocal overlap) (`03`, `06`, `10`).
3. Reproduction of the SNP-based deletion calls from the SNP matrix, and comparison with the
   all-base calls (`04`, `10`, `13`).
4. Read evidence at each call: discordant and concordant read pairs, split reads, and soft-clipped
   reads aligned to a TE library; matched control intervals give the background rate of each kind of
   evidence (`12`, `14`, `14a`, `16`, `17`).
5. TE annotation of E2 with RepeatMasker (public TE sequences, and the combined EDTA library of four
   genomes), TEs at breakpoints, and comparison with the manually checked breakpoints in Data S9 (`05`,
   `18`, `24`-`25`, `28`-`29`).
6. Deletion calls in all 226 strains, gene and effector statistics, and the genome-wide figures (`19`,
   `20`, `27`, `34`-`38`, `40`-`43`, `45`-`48`).

Note: the figure scripts `45`-`48` read the per-window values from the interactive pages in `docs/`
(`heatmap_226_strains_all_base.html` and `deletion_tracks_33_strains.html`). The scripts that built those
pages are not included.

### `rip_analysis/`

Full-length copies of LTR retrotransposon families and of Pot2 were extracted from four assemblies (E2,
MZ5-1-6, B71, P131), aligned with MAFFT, split into subfamilies, and compared with a majority-rule
consensus. Substitutions are counted by type and dinucleotide context (`01`-`03`). RIPCaL 2.0 was run on
the same alignments (`06`, `07`, `10`). Composite, product and substrate RIP indices were computed in
sliding windows along each genome with the settings of The RIPper (`05`, `09`, `11`, `12`, `14`, `16`). Scripts `13` and `16` make the figures.

### `fosbury_ltr/`

Checks of the full-length Fosbury elements in Data S8: identity between the two LTRs of each element,
target-site duplications (TSDs) next to the elements, and indels between the two LTRs. Sequences are
retrieved from NCBI with E-utilities.

## Running the scripts

The cluster jobs were written for a Slurm cluster with environment modules. Paths in the scripts must be
changed for another system: `/path/to/scratch` stands for the project folder on the cluster, and `YOUR_ACCOUNT` and the e-mail
address in the job files are placeholders. Most local Python scripts are run from the
parent folder of `deletion_allbase/` or `rip_analysis/`.

Main software (versions as loaded on the cluster):

| Software | Version |
|---|---|
| BWA | 0.7.18 |
| minimap2 | 2.29 |
| SAMtools, HTSlib, BCFtools | 1.21 |
| BEDTools | 2.31.1 |
| Delly | 1.3.3 |
| DeepVariant, GLnexus | container images, see `phylogeny_busco/slurm/` |
| IQ-TREE | 3.0.1 |
| BUSCO | 5.8.3 and 6.0.0 |
| RepeatMasker | 4.1.9 |
| BLAST+ | 2.16.0 |
| MAFFT | 7.526 |
| SeqKit | 2.9.0 |
| RIPCaL | 2.0 |
| Python | 3.12 with pandas, NumPy, SciPy, matplotlib, openpyxl |
