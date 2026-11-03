# ============================================================
# NOMBRE: 05_QC_boruta_outputs.R
# DESCRIPCIÓN: Verificación de integridad de los outputs de 05_boruta.R
# USO: Rscript 1_meta_analysis/05_QC_boruta_outputs.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(fs); library(readr)
})

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- file.path(cfg$outputs$dir, "boruta")

cat("🔍 Verificando integridad de outputs de 05_boruta.R\n\n")

expected_files <- c(
  "boruta_model.rds",
  "boruta_selected_otu.csv",
  "boruta_selected_otus_with_taxonomy.csv",
  "boruta_importance_clean.png"
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
cat("🧪 Cargando modelo Boruta...\n")

tryCatch({
  library(Boruta)
  res <- readRDS(file.path(out_dir, "boruta_model.rds"))
  stopifnot(inherits(res, "Boruta"))
  cat("✅ boruta_model.rds cargado\n")
  
  # Decision summary
  dec <- res$finalDecision
  tab <- table(dec)
  cat("   Decisiones: ",
      paste(names(tab), "=", as.integer(tab), collapse = ", "), "\n")
  
  n_confirmed <- sum(dec == "Confirmed")
  if (n_confirmed == 0) {
    cat("⚠️ Cero features 'Confirmed' — puede ser esperado o fallo.\n")
  } else {
    cat(sprintf("   ✅ %d ASVs confirmados\n", n_confirmed))
  }
}, error = function(e) {
  cat("❌ Error cargando boruta_model.rds:", e$message, "\n"); all_ok <<- FALSE
})

cat("\n------------------------------------------------------------\n")
cat("🧪 Validando tablas de ASVs seleccionados...\n")

tryCatch({
  sel <- readr::read_csv(file.path(out_dir, "boruta_selected_otu.csv"),
                         show_col_types = FALSE)
  cat(sprintf("✅ boruta_selected_otu.csv: %d muestras × %d columnas\n",
              nrow(sel), ncol(sel)))
  
  if (ncol(sel) <= 1) {
    cat("❌ No hay columnas de OTUs (solo 'Sample').\n"); all_ok <- FALSE
  } else {
    # All feature columns must be prefixed
    feature_cols <- setdiff(colnames(sel), "Sample")
    bad <- !grepl("^OTU_", feature_cols)
    if (any(bad)) {
      cat(sprintf("⚠️ %d columnas sin prefijo OTU_: %s\n",
                  sum(bad), paste(head(feature_cols[bad], 3), collapse = ", ")))
    } else {
      cat("   ✅ Todas las columnas tienen prefijo 'OTU_'\n")
    }
  }
}, error = function(e) {
  cat("❌ Error leyendo boruta_selected_otu.csv:", e$message, "\n"); all_ok <<- FALSE
})

tryCatch({
  annot <- readr::read_csv(file.path(out_dir, "boruta_selected_otus_with_taxonomy.csv"),
                           show_col_types = FALSE)
  cat(sprintf("✅ boruta_selected_otus_with_taxonomy.csv: %d filas × %d cols\n",
              nrow(annot), ncol(annot)))
  
  required <- c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  missing  <- setdiff(required, colnames(annot))
  if (length(missing) > 0) {
    cat("❌ Columnas de taxonomía faltantes:", paste(missing, collapse = ", "), "\n")
    all_ok <- FALSE
  }
  
  # OTU_ID should have NO prefix
  if (any(grepl("^OTU_", annot$OTU_ID))) {
    cat("❌ OTU_ID contiene prefijo 'OTU_' — join falló.\n"); all_ok <- FALSE
  } else {
    cat("   ✅ OTU_ID sin prefijo (join consistente)\n")
  }
  
  # Genus not all NA
  n_genus <- sum(!is.na(annot$Genus) & annot$Genus != "")
  if (n_genus == 0) {
    cat("⚠️ Columna Genus está vacía — taxonomía no se anotó.\n")
    all_ok <- FALSE
  } else {
    cat(sprintf("   ✅ %d/%d ASVs con Genus asignado\n", n_genus, nrow(annot)))
  }
}, error = function(e) {
  cat("❌ Error leyendo boruta_selected_otus_with_taxonomy.csv:", e$message, "\n")
  all_ok <<- FALSE
})


cat("\n============================================================\n")
if (all_ok) {
  cat("🎉 05 — TODAS LAS VERIFICACIONES PASARON.\n")
} else { 
  cat("🚨 05 — Problemas de integridad.\n"); 
  stop("Revisa arriba.")
}
cat("============================================================\n")
