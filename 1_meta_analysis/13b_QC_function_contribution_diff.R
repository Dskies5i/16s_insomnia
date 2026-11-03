# ============================================================
# SCRIPT: 13b_QC_function_contribution_diff.R
# DESCRIPTION: Quality Control for Script 13 Functional Contribution Outputs
# USAGE: Rscript 1_meta_analysis/13b_QC_function_contribution_diff.R
# ============================================================

source("1_meta_analysis/00_setup.R")
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(fs)

cat("🔍 Running Quality Control for Script 13 (Functional Contribution)...\n\n")

out_contrib_dir <- "outputs/function_contrib"

# 1. Check if output directory exists
if (!dir_exists(out_contrib_dir)) {
  stop(sprintf("❌ QC Failed: Directory '%s' does not exist. Did Script 13 execute completely?", out_contrib_dir))
}

# 2. List all generated files
generated_files <- dir_ls(out_contrib_dir, recurse = TRUE, type = "file")

if (length(generated_files) == 0) {
  stop("❌ QC Failed: The output directory 'outputs/function_contrib' is empty!")
}

cat(sprintf("✔️ Found %d generated files in '%s':\n", length(generated_files), out_contrib_dir))

for (f in generated_files) {
  rel_path <- path_rel(f, start = out_contrib_dir)
  cat(sprintf("  [PASS] %s (Size: %.2f KB)\n", rel_path, file_size(f) / 1024))
}

cat("\n✅ All Script 13 QC checks passed successfully!\n")
