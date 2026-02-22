## 16s Analysis - Workflow

### Step 0 – Setup

Install required packages (see `workflow/98_notebooks/00_setup.qmd`).

### Step 1 – Primer Trimming (`workflow/98_notebooks/01_primer_trimming.qmd`)

Removes primer sequences from raw reads using **cutadapt** (called via R). Generates quality reports and a trimmed FASTQ directory used by all downstream steps.

### Step 2 – Trim Parameter Optimisation (`workflow/98_notebooks/02_figaro_trim_optimisation.qmd`)

Runs **rFigaro** on the trimmed reads to determine optimal `truncLen` and `maxEE` parameters for DADA2. Accounts for whether the sequencing platform uses **continuous** or **binned** quality scores (e.g., NovaSeq 6000 / NextSeq).

### Step 3 – DADA2 Processing (`workflow/98_notebooks/03_dada2_processing.qmd`)

Full DADA2 pipeline:

* Quality filtering and truncation
* Error learning (with platform-aware settings for binned vs continuous Q-scores)
* Denoising (sample inference)
* Paired-end read merging
* Chimera removal
* Read-tracking summary

### Step 4 – Controls QC (`workflow/98_notebooks/04_controls_qc.qmd`)

* **Negative controls**: contamination detection using `decontam` (prevalence and frequency methods)
* **Positive controls** (mock community): accuracy assessment, expected vs observed composition
* Contaminant ASV removal

### Step 5 – phyloseq Assembly (`workflow/98_notebooks/05_phyloseq_assembly.qmd`)

* Merge DADA2 outputs with sample metadata
* Taxonomic annotation (SILVA/GTDB)
* Optional: phylogenetic tree construction (DECIPHER + phangorn)
* Filtering (prevalence, minimum counts)
* Object saved as `data/phyloseq/ps.rds`

### Step 6 – Diversity Analysis (`workflow/98_notebooks/06_diversity_analysis.qmd`)

**Alpha diversity**

* Observed richness, Shannon, Simpson, Faith's PD
* Bias-corrected richness using **breakaway** (Willis & Bunge 2015) and **DivNet** (Willis & Martin 2020) for Shannon/Simpson
* Statistical testing with `betta()` (breakaway) for group comparisons

**Beta diversity**

* Bray-Curtis, weighted & unweighted UniFrac, Aitchison (CLR-transformed)
* PCoA / NMDS ordination
* PERMANOVA (`adonis2`), PERMDISP

### Step 7 – Differential Abundance (`workflow/98_notebooks/07_differential_abundance.qmd`)

* **ANCOM-BC2** (`ANCOMBC` package): composition-aware differential testing with bias correction, multiple comparison control, and sensitivity/specificity analysis
* **Dirichlet regression** (`DirichletReg`) as an alternative compositional approach
* Volcano plots, heatmaps, and bar charts of significant taxa
