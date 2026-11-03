# ============================================================
# utils_helper_functions.R
# Funciones auxiliares para el análisis funcional (PICRUSt2 y Correlaciones)
# ============================================================

library(readr)
library(dplyr)
library(phyloseq)

# 1) Leer la lista de KOs significativos desde un archivo de bordes o tabla
make_ko_list <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  # Asegurar que devuelva una columna ko_function o similar
  if (!"ko_function" %in% colnames(df)) {
    if ("Target" %in% colnames(df)) {
      df$ko_function <- df$Target
    } else if ("KO" %in% colnames(df)) {
      df$ko_function <- df$KO
    } else {
      df$ko_function <- df[[1]]
    }
  }
  return(df)
}

# 2) Subseleccionar la matriz de KOs de PICRUSt2
subset_ko_matrix <- function(ko_path, sig_ko) {
  ko_df <- read.delim(ko_path, header = TRUE, row.names = 1, check.names = FALSE)
  # Filtrar filas cuyos nombres coincidan con sig_ko (ignorando prefijos si es necesario)
  common_kos <- intersect(rownames(ko_df), sig_ko)
  if (length(common_kos) == 0) {
    # Intentar limpiar prefijos "ko:" si los hay
    rownames(ko_df) <- sub("^ko:", "", rownames(ko_df))
    sig_ko_clean <- sub("^ko:", "", sig_ko)
    common_kos <- intersect(rownames(ko_df), sig_ko_clean)
  }
  if (length(common_kos) == 0) {
    stop("No se encontraron KOs coincidentes entre la matriz PICRUSt2 y la lista significativa.")
  }
  as.matrix(ko_df[common_kos, , drop = FALSE])
}

# 3) Leer taxonomía, OTUs y metadatos desde el objeto phyloseq filtrado o archivos
read_tax_otu_meta <- function() {
  ps_path <- file.path(OUT_DIR, "models/physeq_filtered.rds")
  if (file.exists(ps_path)) {
    ps <- readRDS(ps_path)
    list(
      otu = as.matrix(otu_table(ps)),
      tax = as.matrix(tax_table(ps)),
      meta = as.data.frame(sample_data(ps))
    )
  } else {
    stop("No se encontró el objeto phyloseq filtrado en: ", ps_path)
  }
}

# 4) Crear matriz de taxones objetivo (por Género o Especie)
make_target_taxa_matrix <- function(tax_mat, otu_mat) {
  # Agrupar por Genus si está disponible, o usar OTUs directamente
  if ("Genus" %in% colnames(tax_mat)) {
    df <- as.data.frame(otu_mat)
    df$Genus <- tax_mat[rownames(df), "Genus"]
    df <- df %>% filter(!is.na(Genus) & Genus != "" & Genus != "Unassigned")
    # Agrupar por media o suma por género
    df_agg <- df %>% group_by(Genus) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE))
    mat <- as.matrix(df_agg %>% tibble::column_to_rownames("Genus"))
    return(mat)
  } else {
    return(otu_mat)
  }
}

# 5) Transformación logarítmica con pseudocuento
log_transform <- function(mat, base = 2, pseudocount = 1) {
  if (base == 2) {
    log2(mat + pseudocount)
  } else if (base == 10) {
    log10(mat + pseudocount)
  } else {
    log(mat + pseudocount, base = base)
  }
}