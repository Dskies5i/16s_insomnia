# 🧬 Preprocessing Workflow (FASTQ → Unified OTU Table & Master Metadata)

This directory contains the modular bioinformatics pipeline to download, import, trim primers, perform *denoising* (DADA2), unify via clustering (VSEARCH 99%), visualize, and structure metadata for 16S rRNA insomnia datasets. The final product is ready for integration and analysis in R (`phyloseq`).

---

## 📁 1. Project Structure and Scripts

```text
16s_insomnia/
├── 0_preprocessing/
│   ├── config/
│   │   └── datasets_config.tsv             <-- Dynamic parameters (primers, truncation, layout)
│   ├── script_drafts_ideas/                <-- Drafts and experimental tests
│   ├── utils/                              <-- Support scripts and auxiliary utilities
│   └── scripts/
│       ├── 00_install_sra_tools.sh         <-- SRA tools installation
│       ├── 01_download_silva_classifier.sh <-- SILVA classifier download (Full-Length)
│       ├── 01_download_study_sras_w_metadata.sh <-- Raw FASTQ and SRA metadata download
│       ├── 01_install_qiime2.sh            <-- QIIME 2 environment setup
│       ├── 02_metadata_qiime2_format.sh    <-- Base metadata formatting for QIIME 2
│       ├── 03_merge_study_w_downloaded_metadata.py <-- Enrich and harmonize clinical metadata (ISI, BMI, etc.)
│       ├── 04_qiime_import_and_trim.sh     <-- FASTQ import and primer removal (Cutadapt)
│       ├── 04_viz_quality_plot.sh          <-- .qzv generation for quality inspection
│       ├── 05_dada2_denoise_datasets.sh    <-- Denoising and error correction with DADA2 (dynamic Paired/Single)
│       ├── 05b_visualize_denoise_stats.sh  <-- Statistical summary of DADA2 process
│       ├── 06_merge_cluster_export.sh      <-- Merge, VSEARCH 99% clustering, SILVA taxonomy, and R export
│       ├── 06b_visualize_unified_table.sh  <-- Generates summary visualization of unified table (.qzv)
│       └── 06D_build_master_metadata.py    <-- Builds final master metadata table (master_metadata_unified.tsv)
├── data/
│   ├── pre_processed/
│   │   ├── metadata/                       <-- Formatted (*_q2.tsv) and enriched (*_format_enriched.tsv) metadata
│   │   └── unified/                        <-- Final matrices, FASTA, and master_metadata_unified.tsv for R
│   └── raw/                                <-- Raw FASTQs per study and supplementary tables
└── resources/                              <-- SILVA classifiers and external resources
```

ENVIRONMENT: conda activate qiime2-amplicon-2024.5
---

## 🚀 2. Execution Order, Descriptions, and Output Files

### Phase A: Data and Metadata Acquisition

#### 1. Base Tools and Environment Installation

- **Script:** `bash 0_preprocessing/scripts/00_install_sra_tools.sh` and `bash 0_preprocessing/scripts/01_install_qiime2.sh`
- **Description:** Configures the Conda virtual environment with the necessary QIIME 2 dependencies (amplicon version) and SRA Toolkit tools for genomic data download.
- **Output Files:** Active Conda environment (`qiime2-amplicon-2024.5`) and configured binaries.

#### 2. FASTQ Sequence and SRA Metadata Download

- **Script:** `bash 0_preprocessing/scripts/01_download_study_sras_w_metadata.sh`
- **Description:** Connects to the SRA repository using the defined BioProjects to download raw `.fastq.gz` files for each study and their associated metadata.
- **Output Files:**
  - `data/raw/<study_id>/fastq/*.fastq.gz`
  - Raw SRA metadata tables.

#### 3. Metadata Formatting and Enrichment

- **Scripts:**
  - `bash 0_preprocessing/scripts/02_metadata_qiime2_format.sh`
  - `python3 0_preprocessing/scripts/03_merge_study_w_downloaded_metadata.py`
- **Description:** Adapts initial metadata to QIIME 2 formal standards (injecting `#q2:types` header) and harmonizes key clinical variables (such as ISI scales, BMI, and normalization of study groups into `Insomnia` vs `Control`).
- **Output Files:**
  - `data/pre_processed/metadata/*_q2.tsv`
  - `data/pre_processed/metadata/*_format_enriched.tsv`

---

### Phase B: Bioinformatics Processing in QIIME 2

#### 4. Taxonomic Classifier SILVA Download

- **Script:** `bash 0_preprocessing/scripts/01_download_silva_classifier.sh`
- **Description:** Downloads the Naive Bayes classifier trained with the SILVA database (99% identity) for subsequent taxonomic assignment.
- **Output Files:**
  - `resources/classifiers/silva-138-99-nb-classifier.qza`

#### 5. Read Import and Primer Trimming (Cutadapt)

- **Script:** `bash 0_preprocessing/scripts/04_qiime_import_and_trim.sh`
- **Description:** Generates manifest files, imports FASTQ sequences into native QIIME 2 objects (`.qza`), and removes primers specific to each amplified region (V3-V4 or V4) using Cutadapt.
- **Output Files:**
  - `data/pre_processed/qiime2_workspace/<study_id>/demux.qza`
  - `data/pre_processed/qiime2_workspace/<study_id>/demux_trimmed.qza`

#### 6. Quality Plot Generation

- **Script:** `bash 0_preprocessing/scripts/04_viz_quality_plot.sh`
- **Description:** Generates interactive visual summaries of per-base quality (Phred scores) for each study before filtering.
- **Output Files:**
  - `data/pre_processed/qiime2_workspace/<study_id>/demux_summary.qzv` *(for inspection at view.qiime2.org)*

#### 7. Denoising and Error Correction with DADA2

- **Script:** `bash 0_preprocessing/scripts/05_dada2_denoise_datasets.sh`
- **Description:** Runs the DADA2 algorithm dynamically (automatically detecting whether the design is *Paired-End* or *Single-End* based on configuration) using the truncation parameters (`trunc_f`, `trunc_r`) defined in `datasets_config.tsv`.
- **Output Files:**
  - `data/pre_processed/qiime2_workspace/<study_id>/table.qza` (ASV count table)
  - `data/pre_processed/qiime2_workspace/<study_id>/rep_seqs.qza` (Representative sequences)
  - `data/pre_processed/qiime2_workspace/<study_id>/denoising_stats.qza` (Retention statistics)

#### 8. Denoising Statistics Visualization

- **Script:** `bash 0_preprocessing/scripts/05b_visualize_denoise_stats.sh`
- **Description:** Visualizes performance and percentage of non-chimeric reads retained per sample after DADA2.
- **Output Files:**
  - `.qzv` files for filtering and denoising statistics.

---

### Phase C: Global Integration, Clustering, and R Export

#### 9. Merge, VSEARCH 99% Clustering, Taxonomy, and Export

- **Script:** `bash 0_preprocessing/scripts/06_merge_cluster_export.sh`
- **Description:** Merges all individual tables, groups representative sequences via *de novo* clustering at 99% identity with VSEARCH (creating a unified feature vocabulary), assigns global taxonomy with SILVA, and exports clean files readable by R.
- **Output Files:**
  - `data/pre_processed/unified/unified_otu_table.tsv` and `.csv`
  - `data/pre_processed/unified/unified_rep_seqs.fasta`
  - `data/pre_processed/unified/taxonomy.tsv`

#### 10. Unified Table Visualization (Optional)

- **Script:** `bash 0_preprocessing/scripts/06b_visualize_unified_table.sh`
- **Description:** Generates a visual summary of the final unified artifact for sequencing depth validation.
- **Output Files:**
  - `.qzv` summary artifact of the global unified table.

#### 11. Master Metadata Table Construction

- **Script:** `python3 0_preprocessing/scripts/06D_build_master_metadata.py`
- **Description:** Intelligently consolidates all individual enriched metadata into a single clean master matrix, unifying sample identifiers with count tables.
- **Output Files:**
  - `data/pre_processed/unified/master_metadata_unified.tsv` *(definitive file for merging with `phyloseq` in R)*

---

## 🆔 3. SRA Identifiers by Dataset

| Dataset | BioProject ID | Region | Strategy | Status |
| --- | --- | --- | --- | --- |
| **liub_2019** | `PRJNA527914` | V3-V4 | Paired-End | Active |
| **chen_2022** | `PRJNA882369` | V3-V4 | Single-End | Active |
| **haimov_2022** | `PRJNA843882` | V4 | Paired-End | Active |

---

## 🛠 Troubleshooting

### Rcpp / DADA2 Error on macOS (Apple Silicon M1/M2/M3)

If `05_dada2_denoise_datasets.sh` fails due to Rcpp shared objects:

```bash
conda activate qiime2-amplicon-2024.5
Rscript -e "install.packages('Rcpp', repos='https://cloud.r-project.org')"
Rscript -e "library(dada2); cat('✅ DADA2 OK\n')"
```
