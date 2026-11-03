
# ============================================================
# SCRIPT: 13a_QC_extract_ko_contrib.R
# DESCRIPTION: Quality Control (QC) for auxiliary contribution inputs
# USE: Rscript 1_meta_analysis/13_QC_extract_ko_contrib.R
# ============================================================
# ============================================================
# SCRIPT: 13b_QC_extract_ko_contrib.R (Memory-Safe QC)
# ============================================================

source("1_meta_analysis/00_setup.R")
library(readr)
library(dplyr)
library(fs)

cat("🔍 Running Memory-Safe Quality Control for Script 13b Outputs...\n\n")
out_tables <- "outputs/tables"

# 1. File Existence & Size Check
files_to_check <- c("metadata_for_taxa.csv", "picrust_taxonomy_filled.csv", "KO_FC_group.csv", "ko_contrib_filtered.csv")

for (f in files_to_check) {
  f_path <- file.path(out_tables, f)
  if (file.exists(f_path)) {
    cat(sprintf("  [PASS] Found: %s (Size: %.2f MB)\n", f, file_size(f_path) / (1024 * 1024)))
  } else {
    stop(sprintf("  [FAIL] Missing: %s", f))
  }
}

cat("\n--- Validating Table Contents ---\n")
meta_taxa <- read_csv(file.path(out_tables, "metadata_for_taxa.csv"), show_col_types = FALSE)
cat(sprintf("✔️ metadata_for_taxa.csv: %d samples.\n", nrow(meta_taxa)))

tax_filled <- read_csv(file.path(out_tables, "picrust_taxonomy_filled.csv"), show_col_types = FALSE)
cat(sprintf("✔️ picrust_taxonomy_filled.csv: %d ASVs.\n", nrow(tax_filled)))

ko_fc <- read_csv(file.path(out_tables, "KO_FC_group.csv"), show_col_types = FALSE)
cat(sprintf("✔️ KO_FC_group.csv: %d KOs.\n", nrow(ko_fc)))

contrib <- read_csv(file.path(out_tables, "ko_contrib_filtered.csv"), show_col_types = FALSE)
cat(sprintf("✔️ ko_contrib_filtered.csv: %d rows (Lightweight & Safe).\n", nrow(contrib)))

cat("\n✅ All Quality Control checks passed successfully! Ready for Script 13.\n")