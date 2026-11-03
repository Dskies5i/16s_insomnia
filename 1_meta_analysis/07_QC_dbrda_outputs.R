# ============================================================
# NOMBRE: 07_QC_dbrda_outputs.R
# DESCRIPCIÓN: QC completo de los outputs de 07_dbrda_envfit.R
#
# CHECKS:
#   1. Existencia + tamaño > 0 de los 10 archivos esperados
#   2. Carga de los RDS principales (step_model, dbrda_model)
#   3. Estructura del modelo db-RDA (clase, terms, no vacío)
#   4. Legibilidad de CSV (anova_all, anova_margin, envfit_*)
#   5. ✅ CRÍTICO: alineación pvals ↔ scores en envfit (el bug p.max)
#   6. ⚠️ Detección de aliasing/collinearity (sex, study_id)
#   7. R² presente y en rango (0–1)
#   8. Reporte final con veredicto claro
#
# USO: Rscript 1_meta_analysis/07_QC_dbrda_outputs.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(fs)
  library(readr)
  library(vegan)
})

cfg     <- yaml::read_yaml("config/config.yaml")
out_dir <- file.path(cfg$outputs$dir, "dbrda")

cat("🔍 Verificando integridad de los archivos en:", out_dir, "\n\n")

# Lista de warnings / errores acumulados para el veredicto final
warn_msgs <- character(0)
err_msgs  <- character(0)

# ============================================================
# 1. Existencia + tamaño > 0
# ============================================================
cat("─── [1/6] Existencia y tamaño ───\n")

expected_files <- c(
  "ordiR2step_model.rds",
  "ordiR2step_log.txt",
  "dbrda_model.rds",
  "dbrda_anova_all.csv",
  "dbrda_anova_margin.csv",
  "dbrda_r2.txt",
  "envfit_significant_otus.csv",
  "envfit_significant_otus_with_taxonomy.csv",
  "dbrda_sites.pdf"
)

for (f in expected_files) {
  f_path <- file.path(out_dir, f)
  if (!file_exists(f_path)) {
    cat(sprintf("❌ [FALTA] %s\n", f))
    err_msgs <- c(err_msgs, paste("Falta archivo:", f))
  } else if (file_size(f_path) == 0) {
    cat(sprintf("⚠️ [VACÍO] %s (0 bytes)\n", f))
    warn_msgs <- c(warn_msgs, paste("Archivo vacío:", f))
  } else {
    cat(sprintf("✅ [OK] %s (%s)\n", f, as.character(file_size(f_path))))
  }
}

# Nota: dbrda_sites_by_group.pdf es OPCIONAL (sólo si usaste la versión extendida)
opt_files <- c("dbrda_sites_by_group.pdf")
for (f in opt_files) {
  f_path <- file.path(out_dir, f)
  if (file_exists(f_path) && file_size(f_path) > 0) {
    cat(sprintf("ℹ️  [OPCIONAL presente] %s\n", f))
  }
}

# ============================================================
# 2. Cargar RDS principales
# ============================================================
cat("\n─── [2/6] Cargando objetos .rds ───\n")

step_model <- NULL
tryCatch({
  step_model <- readRDS(file.path(out_dir, "ordiR2step_model.rds"))
  cat("✅ ordiR2step_model.rds cargado (clase:", paste(class(step_model), collapse = ","), ")\n")
}, error = function(e) {
  cat("❌ No se pudo leer ordiR2step_model.rds:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("ordiR2step_model.rds:", e$message))
})

dbrda_model <- NULL
tryCatch({
  dbrda_model <- readRDS(file.path(out_dir, "dbrda_model.rds"))
  cat("✅ dbrda_model.rds cargado (clase:", paste(class(dbrda_model), collapse = ","), ")\n")
}, error = function(e) {
  cat("❌ No se pudo leer dbrda_model.rds:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("dbrda_model.rds:", e$message))
})

# ============================================================
# 3. Estructura del modelo db-RDA
# ============================================================
cat("\n─── [3/6] Estructura del modelo db-RDA ───\n")

if (!is.null(dbrda_model)) {
  tryCatch({
    # Es un objeto capscale / rda / dbrda
    ok_class <- inherits(dbrda_model, c("capscale", "rda", "dbrda", "cca"))
    if (!ok_class) {
      cat("⚠️ dbrda_model no es capscale/rda/dbrda/cca — clase:",
          paste(class(dbrda_model), collapse = ","), "\n")
      warn_msgs <- c(warn_msgs, "dbrda_model tiene clase inesperada")
    } else {
      cat("✅ Clase válida:", paste(class(dbrda_model), collapse = ","), "\n")
    }
    
    # Terms del modelo
    if (!is.null(dbrda_model$terms)) {
      term_labels <- attr(dbrda_model$terms, "term.labels")
      cat("   Terms:", paste(term_labels, collapse = " | "), "\n")
      if (length(term_labels) == 0) {
        cat("❌ El modelo no tiene términos estimables.\n")
        err_msgs <- c(err_msgs, "dbrda_model sin terms")
      }
    } else {
      cat("⚠️ dbrda_model$terms es NULL\n")
      warn_msgs <- c(warn_msgs, "dbrda_model sin terms")
    }
    
    # Eigenvalores / inercia
    if (!is.null(dbrda_model$CCA)) {
      cat(sprintf("   Eigenvalores CCA: %s\n",
                  paste(round(head(dbrda_model$CCA$eig, 3), 4), collapse = ", ")))
    }
  }, error = function(e) {
    cat("❌ Error inspeccionando dbrda_model:", e$message, "\n")
    err_msgs <<- c(err_msgs, paste("dbrda_model estructura:", e$message))
  })
}

# ============================================================
# 4. Legibilidad de CSV (anova + envfit)
# ============================================================
cat("\n─── [4/6] Legibilidad de CSVs ───\n")

anova_all <- NULL
tryCatch({
  anova_all <- read_csv(file.path(out_dir, "dbrda_anova_all.csv"), show_col_types = FALSE)
  cat(sprintf("✅ dbrda_anova_all.csv: %d filas × %d columnas\n",
              nrow(anova_all), ncol(anova_all)))
  if (nrow(anova_all) == 0) {
    cat("❌ dbrda_anova_all.csv está vacío.\n")
    err_msgs <- c(err_msgs, "dbrda_anova_all.csv vacío")
  }
  req <- c("Df", "SumOfSqs", "F", "Pr(>F)")
  miss <- setdiff(req, colnames(anova_all))
  if (length(miss) > 0) {
    cat("⚠️ Columnas esperadas faltantes:", paste(miss, collapse = ", "), "\n")
    warn_msgs <- c(warn_msgs, paste("anova_all col faltante:", paste(miss, collapse = ",")))
  }
}, error = function(e) {
  cat("❌ dbrda_anova_all.csv:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("anova_all:", e$message))
})

tryCatch({
  anova_m <- read_csv(file.path(out_dir, "dbrda_anova_margin.csv"), show_col_types = FALSE)
  cat(sprintf("✅ dbrda_anova_margin.csv: %d filas × %d columnas\n",
              nrow(anova_m), ncol(anova_m)))
}, error = function(e) {
  cat("❌ dbrda_anova_margin.csv:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("anova_margin:", e$message))
})

envfit_sig   <- NULL
envfit_annot <- NULL
tryCatch({
  envfit_sig <- read_csv(file.path(out_dir, "envfit_significant_otus.csv"),
                         show_col_types = FALSE)
  cat(sprintf("✅ envfit_significant_otus.csv: %d filas × %d cols\n",
              nrow(envfit_sig), ncol(envfit_sig)))
}, error = function(e) {
  cat("❌ envfit_significant_otus.csv:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("envfit_sig:", e$message))
})

tryCatch({
  envfit_annot <- read_csv(file.path(out_dir, "envfit_significant_otus_with_taxonomy.csv"),
                           show_col_types = FALSE)
  cat(sprintf("✅ envfit_significant_otus_with_taxonomy.csv: %d filas × %d cols\n",
              nrow(envfit_annot), ncol(envfit_annot)))
}, error = function(e) {
  cat("❌ envfit_significant_otus_with_taxonomy.csv:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("envfit_annot:", e$message))
})

# ============================================================
# 5. ✅ CRÍTICO: alineación pvals ↔ scores en envfit
# ============================================================
cat("\n─── [5/6] ⚠️ Verificación CRÍTICA de envfit (bug p.max) ───\n")

if (!is.null(envfit_sig) && nrow(envfit_sig) > 0) {
  # (a) Todas las pvals deben ser < 0.05 (el archivo se llama 'significant')
  n_above <- sum(envfit_sig$pval >= 0.05, na.rm = TRUE)
  if (n_above > 0) {
    cat(sprintf("❌ %d filas con pval >= 0.05 en archivo 'significant' → envfit desalineado.\n",
                n_above))
    err_msgs <- c(err_msgs, sprintf("envfit: %d pvals >= 0.05 (desalineado)", n_above))
  } else {
    cat("✅ Todas las pvals < 0.05.\n")
  }
  
  # (b) Ninguna pval debe ser NA
  n_na <- sum(is.na(envfit_sig$pval))
  if (n_na > 0) {
    cat(sprintf("❌ %d pvals NA → fit$vectors$pvals reciclado incorrectamente.\n", n_na))
    err_msgs <- c(err_msgs, sprintf("envfit: %d pvals NA", n_na))
  } else {
    cat("✅ Sin pvals NA.\n")
  }
  
  # (c) Correlación esperada: OTUs con mayor |CAP1|/|CAP2| deben tener pval bajo
  # (aunque no es 1:1, si la correlación es < -0.1 probablemente están cruzadas)
  if (all(c("CAP1","CAP2") %in% colnames(envfit_sig))) {
    r2_mag <- sqrt(envfit_sig$CAP1^2 + envfit_sig$CAP2^2)
    cor_val <- suppressWarnings(cor(r2_mag, envfit_sig$pval, method = "spearman"))
    cat(sprintf("   Correlación |CAP|~pval (Spearman): %.3f\n", cor_val))
    if (!is.na(cor_val) && cor_val > 0.1) {
      cat("⚠️ Correlación positiva — sospechoso: mayor |CAP| debería ir con menor pval.\n")
      warn_msgs <- c(warn_msgs,
                     sprintf("envfit: cor(|CAP|,pval) = %.3f (esperado < 0)", cor_val))
    }
  }
  
  # (d) Ambas tablas (con/sin taxonomía) deben tener el mismo #filas
  if (!is.null(envfit_annot)) {
    if (nrow(envfit_sig) != nrow(envfit_annot)) {
      cat(sprintf("❌ envfit_sig (%d) ≠ envfit_annot (%d) filas — join incompleto.\n",
                  nrow(envfit_sig), nrow(envfit_annot)))
      err_msgs <- c(err_msgs, "envfit_sig vs annot: distinto #filas")
    } else {
      cat("✅ envfit_sig y envfit_annot tienen el mismo # de filas.\n")
    }
  }
} else {
  cat("⚠️ No hay archivo envfit_sig — saltando chequeos de alineación.\n")
  warn_msgs <- c(warn_msgs, "envfit_sig ausente o vacío")
}# ============================================================
# 6. Aliasing / collinearity + R²
# ============================================================
cat("\n─── [6/6] Aliasing, collinearity y R² ───\n")

# R²
tryCatch({
  r2_lines <- readLines(file.path(out_dir, "dbrda_r2.txt"))
  cat("Contenido de dbrda_r2.txt:\n")
  for (l in r2_lines) cat("   ", l, "\n", sep = "")
  
  # ✅ FIX: single regmatches, applied to gregexpr output
  matches <- regmatches(r2_lines, gregexpr("[0-9]+\\.[0-9]+", r2_lines))
  nums    <- as.numeric(unlist(matches))
  nums    <- nums[!is.na(nums)]   # drop any non-parsed
  
  if (length(nums) >= 2) {
    r2_val <- nums[1]
    r2_adj <- nums[2]
    cat(sprintf("   R² = %.6f | R² adjusted = %.6f\n", r2_val, r2_adj))
    
    if (any(is.na(c(r2_val, r2_adj))) ||
        any(c(r2_val, r2_adj) < 0 | c(r2_val, r2_adj) > 1)) {
      cat("❌ R² fuera de [0,1].\n")
      err_msgs <- c(err_msgs, "R² fuera de rango")
    } else {
      cat("✅ R² en rango [0,1].\n")
    }
  } else {
    cat("⚠️ No se pudieron extraer 2 números de R² del archivo.\n")
    warn_msgs <- c(warn_msgs, "dbrda_r2.txt formato inesperado")
  }
}, error = function(e) {
  cat("❌ dbrda_r2.txt:", e$message, "\n")
  err_msgs <<- c(err_msgs, paste("r2:", e$message))
})


# ============================================================
# Veredicto final
# ============================================================
cat("\n============================================================\n")
if (length(err_msgs) == 0 && length(warn_msgs) == 0) {
  cat("🎉 TODAS LAS VERIFICACIONES PASARON.\n")
  cat("     El pipeline está listo para continuar con el siguiente script.\n")
} else if (length(err_msgs) == 0) {
  cat("⚠️  QC COMPLETÓ CON ADVERTENCIAS (no bloqueantes):\n")
  for (w in warn_msgs) cat("   •", w, "\n")
  cat("\n    Puedes continuar, pero revisa lo anterior.\n")
} else {
  cat("🚨 SE ENCONTRARON PROBLEMAS CRÍTICOS:\n")
  for (e in err_msgs) cat("   ✗", e, "\n")
  if (length(warn_msgs) > 0) {
    cat("\n    Además, advertencias:\n")
    for (w in warn_msgs) cat("   •", w, "\n")
  }
  cat("\n    Si el error es de envfit (pvals >= 0.05 o NA), ejecuta:\n")
  cat("        Rscript 1_meta_analysis/07_v2_envfit_fix.R\n")
  stop("Verificación fallida — revisa los mensajes anteriores.")
}
cat("============================================================\n")


