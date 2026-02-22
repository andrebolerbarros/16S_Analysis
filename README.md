# 16S Amplicon Sequencing Analysis Pipeline

This repository presents an R-based pipeline for processing and analysing 16S rRNA gene amplicon sequencing data. The pipeline covers the complete workflow from raw reads to ecological statistics and differential abundance testing.

## Repository Structure

    16s_Analysis/
    ├── README.md                          # This file
    ├── data/                              # Data directory (raw files not tracked by git)
    │   └── README.md                      # Expected data layout and metadata format
    ├── docs/                              # Documentation and setup notebooks
    └── workflow/
        ├── README.md                      # Workflow overview
        ├── 97_packages/
        │   └── rFigaro/                   # Local rFigaro package source
        ├── 98_notebooks/                  # Analysis notebooks (run in order)
        │   ├── 00_setup.qmd
        │   ├── 01_primer_trimming.qmd
        │   ├── 02_figaro_trim_optimisation.qmd
        │   ├── 03_dada2_processing.qmd
        │   ├── 04_controls_qc.qmd
        │   ├── 05_phyloseq_assembly.qmd
        │   ├── 06_diversity_analysis.qmd
        │   └── 07_differential_abundance.qmd
        └── 99_scripts/                    # Relevant scripts & functions


## Dependencies

| Package | Purpose |
| --- | --- |
| dada2 | ASV denoising |
| phyloseq | Microbiome data container & ecology |
| decontam | Contaminant detection |
| ANCOMBC | Differential abundance (ANCOM-BC2) |
| breakaway | Bias-corrected alpha diversity |
| DivNet | Network-based diversity estimation |
| vegan | Multivariate ecology (PERMANOVA etc.) |
| DirichletReg | Dirichlet regression |
| microbiome | Phyloseq utilities |
| ggplot2 | Plotting |
| dplyr / tidyr | Data wrangling |
| rFigaro | Trim parameter optimisation (local) |
| DECIPHER | Taxonomic classification |
| phangorn | Phylogenetic tree |

## Quick Start

1. Place raw FASTQ files in `data/raw/`
2. Place sample metadata in `data/metadata/metadata.csv`
3. Run notebooks in `workflow/98_notebooks/` from `01` → `07`
4. Render reports with Quarto; HTML outputs appear alongside each `.qmd`

## Key Parameters to Set

Each notebook has a clearly marked **User Parameters** chunk at the top. The most important are:

| Parameter | Notebook | Description |
| --- | --- | --- |
| `primer_fwd` / `primer_rev` | 01  | Primer sequences |
| `amplicon_length` | 02  | Expected amplicon length (without primers) |
| `binned_quality` | 03  | `TRUE` for NovaSeq/NextSeq (2-channel), `FALSE` for MiSeq |
| `taxonomy_db` | 05  | Path to SILVA/GTDB reference database |
| `group_var` | 06, 07 | Metadata column for group comparisons |

## References

* Callahan BJ *et al.* (2016) DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581–583. doi:10.1038/nmeth.3869
* McMurdie PJ & Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLoS ONE* 8(4):e61217.
* Microbiome Orchestrating Multiomics Analysis (OMA) Book – https://microbiome.github.io/OMA/
* Willis A & Bunge J (2015) Estimating diversity via frequency ratios. *Biometrics* 71:1042–1049.
* Willis AD & Martin BD (2020) DivNet: Estimating diversity in networked communities. *Peer Community J* 1:e27.
* Lin H & Peddada SD (2020) Analysis of compositions of microbiomes with bias correction. *Nature Communications* 11:3514.
* Cao Y *et al.* (2023) ANCOM-BC2: Analysis of Compositions of Microbiomes with Bias Correction 2. *F1000Research* 13:369.
* Davis NM *et al.* (2018) Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome* 6:226.