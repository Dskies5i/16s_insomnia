# ============================================================
# NOMBRE: 08_QC_lefse.R
# DESCRIPCIÓN: Verificación de integridad + diagnóstico de join de LEfSe
# USO: Rscript 1_meta_analysis/08_QC_lefse.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(fs); library(readr); library(dplyr)
})

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- file.path(cfg$outputs$dir, "lefse")
qc_dir  <- file.path(out_dir, "qc")
dir_create(qc_dir)

cat("==========================================================\n")
cat("QC LEfSe — rank_filter =", cfg$differential$lefse$rank_filter, "\n")
cat("==========================================================\n\n")

# ---- 1. Existencia + tamaño ----
expected_files <- c(
  "lefse_result.csv",
  "lefse_sig_with_otu.csv",
  "lefse_barplot.png"
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

cat("\n----------------------------------------------------------\n")
cat("🧪 Validando contenido de lefse_result.csv y lefse_sig_with_otu.csv\n")

# ---- 2. Contenido del resultado LEfSe ----
res <- NULL
tryCatch({
  res <- readr::read_csv(file.path(out_dir, "lefse_result.csv"),
                         show_col_types = FALSE)
  cat(sprintf("✅ lefse_result.csv: %d filas × %d columnas\n",
              nrow(res), ncol(res)))
  cat("   Columnas:", paste(colnames(res), collapse = ", "), "\n")
  
  # Must have pvalue + ef_lda
  req <- c("pvalue", "ef_lda")
  miss <- setdiff(req, colnames(res))
  if (length(miss) > 0) {
    cat("❌ Columnas faltantes en resultado LEfSe:", paste(miss, collapse = ", "), "\n")
    all_ok <- FALSE
  } else {
    cat(sprintf("   ✅ %d/%d filas con pvalue <= %.3f\n",
                sum(res$pvalue <= cfg$differential$lefse$p_max),
                nrow(res), cfg$differential$lefse$p_max))
    cat(sprintf("   ✅ %d/%d filas con ef_lda >= %.1f\n",
                sum(res$ef_lda >= cfg$differential$lefse$lda_min),
                nrow(res), cfg$differential$lefse$lda_min))
  }
}, error = function(e) {
  cat("❌ Error leyendo lefse_result.csv:", e$message, "\n"); all_ok <<- FALSE
})

# ---- 3. Contenido del resultado enriquecido ----
sig <- NULL
tryCatch({
  sig <- readr::read_csv(file.path(out_dir, "lefse_sig_with_otu.csv"),
                         show_col_types = FALSE)
  cat(sprintf("\n✅ lefse_sig_with_otu.csv: %d filas × %d columnas\n",
              nrow(sig), ncol(sig)))
  cat("   Columnas:", paste(colnames(sig), collapse = ", "), "\n")
  
  # No duplicated features
  n_dup <- sum(duplicated(sig$feature))
  if (n_dup > 0) {
    cat(sprintf("⚠️ %d features DUPLICADAS — fan-out de join.\n", n_dup))
  } else {
    cat("   ✅ Sin features duplicadas\n")
  }
}, error = function(e) {
  cat("❌ Error leyendo lefse_sig_with_otu.csv:", e$message, "\n"); all_ok <<- FALSE
})

# ---- 4. Cobertura taxonómica ----
cat("\n----------------------------------------------------------\n")
cat("🧪 Cobertura taxonómica en lefse_sig_with_otu.csv\n")

if (!is.null(sig)) {
  rank_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rank_cols <- rank_cols[rank_cols %in% colnames(sig)]
  
  coverage <- sapply(rank_cols, function(rc) {
    sum(!is.na(sig[[rc]]) & sig[[rc]] != "" & sig[[rc]] != "Unassigned")
  })
  
  for (rc in rank_cols) {
    cat(sprintf("   %-8s : %d / %d con dato\n", rc, coverage[rc], nrow(sig)))
  }
  
  rank_filter <- cfg$differential$lefse$rank_filter
  if (rank_filter %in% names(coverage)) {
    pct <- round(100 * coverage[rank_filter] / nrow(sig), 1)
    cat(sprintf("\n   %s = %.1f%% cobertura\n", rank_filter, pct))
    if (pct < 90) {
      cat(sprintf("   ⚠️ Cobertura < 90%% para '%s' — revisa el join.\n", rank_filter))
    } else {
      cat(sprintf("   ✅ Cobertura aceptable (>90%%)\n"))
    }
  }
}# close: if (!is.null(sig))

# ---- 5. Reporte ----
report <- list(
  rank_filter = cfg$differential$lefse$rank_filter,
  n_res       = if (!is.null(res)) nrow(res) else 0L,
  n_sig       = if (!is.null(sig)) nrow(sig) else 0L,
  n_dup       = if (!is.null(sig)) sum(duplicated(sig$feature)) else NA_integer_,
  coverage    = if (!is.null(sig) && exists("coverage")) coverage else NULL,
  all_ok      = all_ok
)

saveRDS(report, file.path(qc_dir, "lefse_qc_report.rds"))
cat("\n📝 Reporte QC guardado en:", file.path(qc_dir, "lefse_qc_report.rds"), "\n")

# ---- 6. Veredicto final ----
cat("\n============================================================\n")
if (all_ok) {
  cat("🎉 08b — TODAS LAS VERIFICACIONES PASARON.\n")
} else {
  cat("🚨 08b — Se encontraron problemas de integridad.\n");
  stop("Revisa los errores anteriores.")
}
cat("============================================================\n")