# ============================================================
# NOMBRE: 04_QC_fb_ratio_outputs.R
# DESCRIPCIÓN: Verificación de integridad de los outputs de 04_fb_ratio.R
# USO: Rscript 1_meta_analysis/04_QC_fb_ratio_outputs.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(fs); library(readr)
})

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- cfg$outputs$dir

cat("🔍 Verificando integridad de outputs de 04_fb_ratio.R\n\n")

expected_files <- c(
  "tables/fb_ratio.csv",
  "figures/fb_ratio.png",
  "logs/fb_ratio_stats.txt"
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
cat("🧪 Validando fb_ratio.csv...\n")

tryCatch({
  fb <- readr::read_csv(file.path(out_dir, "tables/fb_ratio.csv"),
                        show_col_types = FALSE)
  required <- c("Sample", "FB_ratio", "Group")
  missing  <- setdiff(required, colnames(fb))
  if (length(missing) > 0) {
    cat("❌ Columnas faltantes:", paste(missing, collapse = ", "), "\n")
    all_ok <- FALSE
  } else {
    cat(sprintf("✅ fb_ratio.csv: %d muestras\n", nrow(fb)))
    
    na_ratio <- sum(is.na(fb$FB_ratio))
    cat(sprintf("   NA en FB_ratio: %d\n", na_ratio))
    if (na_ratio == nrow(fb)) {
      cat("❌ Todos los FB_ratio son NA — ¡fallo total!\n"); all_ok <- FALSE
    }
    
    tab <- table(fb$Group)
    cat("   Group counts: ",
        paste(names(tab), "=", as.integer(tab), collapse = ", "), "\n")
    
    # Range sanity
    cat(sprintf("   FB_ratio rango: [%.3f, %.3f], mediana=%.3f\n",
                min(fb$FB_ratio, na.rm = TRUE),
                max(fb$FB_ratio, na.rm = TRUE),
                median(fb$FB_ratio, na.rm = TRUE)))
  }
}, error = function(e) {
  cat("❌ Error leyendo fb_ratio.csv:", e$message, "\n"); all_ok <<- FALSE
})

# Stats file must contain a t-test
tryCatch({
  st <- readLines(file.path(out_dir, "logs/fb_ratio_stats.txt"))
  if (!any(grepl("t = ", st)) || !any(grepl("p-value", st))) {
    cat("⚠️ fb_ratio_stats.txt no contiene 't =' o 'p-value'.\n")
  } else {
    cat("✅ fb_ratio_stats.txt contiene t-test\n")
  }
}, error = function(e) {
  cat("❌ Error leyendo fb_ratio_stats.txt:", e$message, "\n"); all_ok <<- FALSE
})

cat("\n============================================================\n")
if (all_ok) { 
  cat("🎉 04 — TODAS LAS VERIFICACIONES PASARON.\n")
} else { 
  cat("🚨 04 — Problemas de integridad.\n"); stop("Revisa arriba.") 
  }
cat("============================================================\n")