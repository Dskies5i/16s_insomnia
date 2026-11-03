# ============================================================
# NOMBRE: 06_QC_boruta_pca_heatmap_outputs.R
# DESCRIPCIÓN: Verificación de integridad de los outputs de 06_boruta_pca_heatmap.R
# USO: Rscript 1_meta_analysis/06_QC_boruta_pca_heatmap_outputs.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(fs); library(readr)
})

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- file.path(cfg$outputs$dir, "boruta")

cat("🔍 Verificando integridad de outputs de 06_boruta_pca_heatmap.R\n\n")

expected_files <- c(
  "boruta_pca.png",
  "boruta_heatmap_samples.png"
)

all_ok <- TRUE
for (f in expected_files) {
  f_path <- file.path(out_dir, f)
  if (!file_exists(f_path)) {
    cat(sprintf("❌ [FALTA] %s\n", f)); all_ok <- FALSE
  } else {
    fs_size <- file_size(f_path)
    if (fs_size == 0) { cat(sprintf("⚠️ [VACÍO] %s\n", f)); all_ok <- FALSE }
    else { cat(sprintf("✅ [OK] %s (%s)\n", f, as.character(fs_size))) }
  }
}

cat("\n------------------------------------------------------------\n")
cat("🧪 Verificando que los PNG son imágenes válidas (>5KB)...\n")

for (f in expected_files) {
  f_path <- file.path(out_dir, f)
  if (file_exists(f_path)) {
    sz <- file_size(f_path)
    if (sz < fs_bytes("5KB")) {
      cat(sprintf("⚠️ %s tiene solo %s — puede estar corrupto.\n",
                  f, as.character(sz)))
      all_ok <- FALSE
    }
  }
}

# Optionally, cross-check that 05 output exists (06 depends on it)
sel_path <- file.path(out_dir, "boruta_selected_otu.csv")
if (!file_exists(sel_path)) {
  cat("❌ boruta_selected_otu.csv ausente — 06 no pudo haberse ejecutado.\n")
  all_ok <- FALSE
} else {
  cat("✅ boruta_selected_otu.csv (input de 06) presente\n")
}

cat("\n============================================================\n")
if (all_ok) { 
  cat("🎉 06 — TODAS LAS VERIFICACIONES PASARON.\n")
} else { 
  cat("🚨 06 — Problemas de integridad.\n"); stop("Revisa arriba.") 
  }
cat("============================================================\n")