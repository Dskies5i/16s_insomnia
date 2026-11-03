# ============================================================
# NOMBRE: 02_QC_alpha_diversity_outputs.R
# DESCRIPCIÓN: Verificación de integridad de los outputs de 02_alpha_diversity.R
# USO: Rscript 1_meta_analysis/02_QC_alpha_diversity_outputs.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(fs); library(readr)
})

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- cfg$outputs$dir

cat("🔍 Verificando integridad de outputs de 02_alpha_diversity.R\n\n")

expected_files <- c(
  "tables/alpha_diversity.csv",
  "figures/alpha_shannon_simpson.png",
  "logs/alpha_diversity_stats.txt"    # or whatever your 02 writes; adjust
)

all_ok <- TRUE
for (f in expected_files) {
  f_path <- file.path(out_dir, f)
  if (!file_exists(f_path)) {
    cat(sprintf("❌ [FALTA] %s no existe.\n", f)); all_ok <- FALSE
  } else {
    fs_size <- file_size(f_path)
    if (fs_size == 0) { cat(sprintf("⚠️ [VACÍO] %s\n", f)); all_ok <- FALSE }
    else { cat(sprintf("✅ [OK] %s (%s)\n", f, as.character(fs_size))) }
  }
}

cat("\n------------------------------------------------------------\n")
cat("🧪 Validando contenido de alpha_diversity.csv...\n")

tryCatch({
  alpha <- readr::read_csv(file.path(out_dir, "tables/alpha_diversity.csv"),
                           show_col_types = FALSE)
  required_cols <- c("Sample", "Shannon", "Simpson", "Group")
  missing <- setdiff(required_cols, colnames(alpha))
  if (length(missing) > 0) {
    cat("❌ Columnas faltantes:", paste(missing, collapse = ", "), "\n")
    all_ok <- FALSE
  } else {
    cat(sprintf("✅ alpha_diversity.csv: %d muestras, %d columnas\n",
                nrow(alpha), ncol(alpha)))
    
    # No NAs in Shannon / Simpson
    na_sh  <- sum(is.na(alpha$Shannon))
    na_si  <- sum(is.na(alpha$Simpson))
    if (na_sh > 0 || na_si > 0) {
      cat(sprintf("⚠️ NA detectados — Shannon: %d, Simpson: %d\n", na_sh, na_si))
      all_ok <- FALSE
    } else {
      cat("   ✅ Sin NAs en Shannon / Simpson\n")
    }
    
    # Group counts
    tab <- table(alpha$Group)
    cat("   Group counts: ",
        paste(names(tab), "=", as.integer(tab), collapse = ", "), "\n")
  }
}, error = function(e) {
  cat("❌ Error leyendo alpha_diversity.csv:", e$message, "\n")
  all_ok <<- FALSE
})

cat("\n============================================================\n")
if (all_ok) {
  cat("🎉 02 — TODAS LAS VERIFICACIONES PASARON.\n")
} else { 
  cat("🚨 02 — Problemas de integridad.\n")
  stop("Revisa arriba.") 
}
cat("============================================================\n")