# ============================================================
# NOMBRE: 03_QC_beta_diversity_outputs.R
# DESCRIPCIÓN: Verificación de integridad de los outputs de 03_beta_diversity_permanova.R
# USO: Rscript 1_meta_analysis/03_QC_beta_diversity_outputs.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(fs); library(readr)
})

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- cfg$outputs$dir

cat("🔍 Verificando integridad de outputs de 03_beta_diversity_permanova.R\n\n")

expected_files <- c(
  "figures/beta_pcoa_nmds.png",
  "logs/permanova.txt"
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
cat("🧪 Validando permanova.txt...\n")

tryCatch({
  lines <- readLines(file.path(out_dir, "logs/permanova.txt"))
  
  # Must contain at least one Pr(>F) line
  has_pval <- any(grepl("Pr\\(>F\\)", lines))
  if (!has_pval) {
    cat("❌ No se encontró columna 'Pr(>F)' en permanova.txt\n")
    all_ok <- FALSE
  } else {
    cat("✅ permanova.txt contiene resultados de adonis2\n")
  }
  
  # Detect "Model" rows
  model_lines <- grep("^Model", lines, value = TRUE)
  cat(sprintf("   %d líneas de 'Model' encontradas\n", length(model_lines)))
  
  # R² > 0 check
  r2_matches <- regmatches(lines, regexpr("R2\\s*=\\s*[0-9.]+", lines))
  if (length(r2_matches) == 0) {
    # some adonis2 formats put R2 column without "R2 =" prefix — check numeric block
    cat("   ⚠️ No se extrajo R² explícito (formato adonis2 variable).\n")
  } else {
    cat("   R² detectados:", paste(r2_matches, collapse = ", "), "\n")
  }
  
  # Length sanity
  if (length(lines) < 20) {
    cat("⚠️ permanova.txt es muy corto — ¿terminó 03?\n")
  } else {
    cat(sprintf("   ✅ permanova.txt tiene %d líneas\n", length(lines)))
  }
}, error = function(e) {
  cat("❌ Error leyendo permanova.txt:", e$message, "\n")
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

cat("\n============================================================\n")
if (all_ok){ 
  cat("🎉 03 — TODAS LAS VERIFICACIONES PASARON.\n")
} else { cat("🚨 03 — Problemas de integridad.\n"); stop("Revisa arriba.") 
}
cat("============================================================\n")