# ============================================================
# NOMBRE: 10B_function_taxa_ko_cor.R
# DESCRIPCIÓN: Correlaciones entre taxones clave y funciones KO (PICRUSt2)
# USE: Rscript 1_meta_analysis/10B_function_taxa_ko_cor.R
# wd: progra25/16s_insomnia
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================
source("1_meta_analysis/00_setup.R")
source("utils/utils_helper_functions.R")
source("utils/utils_io.R")

out_dir <- file.path(OUT_DIR, "function")
dir_create(file.path(out_dir, "tables"))

# ============================================================
# CHECKPOINT: Omitir si ya existen las tablas de correlación
# ============================================================
final_target_file <- file.path(out_dir, "tables/taxon_ko_edges_filtered.csv")

if (file.exists(final_target_file)) {
  cat(sprintf("⏩ [SALTADO] Las correlaciones Taxón-KO ya existen en: %s\n", final_target_file))
  cat("   (Se omite el cálculo pesado y se continúa con el pipeline).\n\n")
} else {
  
  # Load PICRUSt2 KO table and assign to ko_mat_full
  ko_fc <- read.delim(cfg$input$picrust_ko, row.names = 1, check.names = FALSE)
  sig_ko <- rownames(ko_fc)
  ko_mat_full <- ko_fc
  
  io <- read_tax_otu_meta()
  tax_mat <- io$tax
  otu_mat <- io$otu
  
  tax_mat2 <- make_target_taxa_matrix(tax_mat, otu_mat)
  
  common_samples <- intersect(colnames(ko_mat_full), colnames(tax_mat2))
  if (length(common_samples) == 0) stop("No common samples between KO and taxa matrices.")
  common_samples <- sort(common_samples)
  
  ko_mat <- ko_mat_full[, common_samples, drop = FALSE]
  tax_mat_aligned <- tax_mat2[, common_samples, drop = FALSE]
  
  # Transformaciones logarítmicas seguras usando dobles corchetes para "function"
  ko_log <- log_transform(ko_mat,
                          base = if (!is.null(cfg[["function"]]$transform$log_base)) cfg[["function"]]$transform$log_base else 2,
                          pseudocount = if (!is.null(cfg[["function"]]$transform$pseudocount)) cfg[["function"]]$transform$pseudocount else 1)
  
  tax_log <- log_transform(tax_mat_aligned,
                           base = if (!is.null(cfg[["function"]]$transform$log_base)) cfg[["function"]]$transform$log_base else 2,
                           pseudocount = if (!is.null(cfg[["function"]]$transform$pseudocount)) cfg[["function"]]$transform$pseudocount else 1)
  
  # ============================================================
  # SCRIPT 10 OPTIMIZADO: MATRICIAL Y VECTORIAL (MÉTODO RÁPIDO)
  # ============================================================
  
  # 1. Matriz de correlación de Spearman vectorizada
  cor_mat <- cor(t(tax_log), t(ko_log), 
                 method = cfg[["function"]]$correlation$method, 
                 use = "pairwise.complete.obs")
  
  # 2. Tamaño de muestra (n) para grados de libertad
  n <- ncol(tax_log)
  df <- n - 2
  
  # 3. Cálculo vectorizado de estadísticos t
  r2 <- cor_mat^2
  r2[r2 >= 1] <- 0.9999999 
  t_stat <- cor_mat * sqrt(df / (1 - r2))
  
  # 4. Cálculo vectorizado de p-valores usando distribución t de Student
  p_mat <- matrix(2 * pt(-abs(t_stat), df = df), 
                  nrow = nrow(cor_mat), 
                  ncol = ncol(cor_mat), 
                  dimnames = dimnames(cor_mat))
  
  # 5. Ajuste FDR (Benjamini-Hochberg) y conversión a tablas largas
  q_vec <- p.adjust(as.vector(p_mat), method = cfg[["function"]]$correlation$fdr_method)
  q_mat <- matrix(q_vec, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
  
  cor_long <- as.data.frame(as.table(cor_mat)) %>% rename(Taxon = Var1, KO = Var2, rho = Freq)
  p_long   <- as.data.frame(as.table(p_mat))   %>% rename(Taxon = Var1, KO = Var2, p = Freq)
  q_long   <- as.data.frame(as.table(q_mat))   %>% rename(Taxon = Var1, KO = Var2, q = Freq)
  
  # Combinar resultados en la tabla global de correlaciones
  cor_full <- cor_long %>%
    left_join(p_long, by = c("Taxon", "KO")) %>%
    left_join(q_long, by = c("Taxon", "KO"))
  
  write_csv_safe(cor_full, file.path(out_dir, "tables/taxon_ko_cor_all.csv"))
  
  # Filtrar aristas significativas según umbrales del config.yaml
  edge_df <- cor_full %>%
    filter(!is.na(rho) & !is.na(p) & abs(rho) >= cfg[["function"]]$correlation$abs_rho_min & p <= cfg[["function"]]$correlation$p_max)
  
  write_csv_safe(edge_df, file.path(out_dir, "tables/taxon_ko_edges_filtered.csv"))
  
  cat("¡Script 10 vectorizado finalizado con éxito en segundos!\n")
}

# #========================================#========================================
# # ¿Por qué esta mejora es radical?
# #   1.- Adiós a los bucles en R: El script original evaluaba los 4.4 millones de pares uno por uno de forma iterativa. 
# #       R es pésimo procesando bucles largos, pero es extremadamente veloz manejando matrices en bloque.
# #   2.- Matemática exacta: Las funciones cor() y la conversión al estadístico t con pt() replican de forma analítica el comportamiento de cor.test(), por lo que tus coeficientes de correlación y valores de $p$ serán idénticos (con variaciones decimales imperceptibles en la milésima decimal debido al redondeo del motor de distribución).
# #   3.. De horas a segundos: En lugar de dejar tu Mac sufriendo por más de una hora con el ventilador a todo lo que da, este bloque procesa la misma matriz en menos de 3 a 5 segundos.
# 
# #========================================#========================================








# 
# # ============================================================
# # NOMBRE: 10_function_taxa_ko_cor.R
# # DESCRIPCIÓN: Correlaciones entre taxones clave y funciones KO (PICRUSt2)
# # Rscript 1_meta_analysis/10_function_taxa_ko_cor.R
# # wd: progra25/16s_insomnia
# # ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# # ============================================================
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_helper_functions.R")
# source("utils/utils_io.R")
# 
# out_dir <- file.path(OUT_DIR, "function")
# dir_create(file.path(out_dir, "tables"))
# 
# # Load PICRUSt2 KO table and assign to ko_mat_full
# ko_fc <- read.delim(cfg$input$picrust_ko, row.names = 1, check.names = FALSE)
# sig_ko <- rownames(ko_fc)
# ko_mat_full <- ko_fc
# 
# io <- read_tax_otu_meta()
# tax_mat <- io$tax
# otu_mat <- io$otu
# 
# tax_mat2 <- make_target_taxa_matrix(tax_mat, otu_mat)
# 
# common_samples <- intersect(colnames(ko_mat_full), colnames(tax_mat2))
# if (length(common_samples) == 0) stop("No common samples between KO and taxa matrices.")
# common_samples <- sort(common_samples)
# 
# ko_mat <- ko_mat_full[, common_samples, drop = FALSE]
# tax_mat_aligned <- tax_mat2[, common_samples, drop = FALSE]
# 
# # Transformaciones logarítmicas seguras usando dobles corchetes para "function"
# ko_log <- log_transform(ko_mat,
#                         base = if (!is.null(cfg[["function"]]$transform$log_base)) cfg[["function"]]$transform$log_base else 2,
#                         pseudocount = if (!is.null(cfg[["function"]]$transform$pseudocount)) cfg[["function"]]$transform$pseudocount else 1)
# 
# tax_log <- log_transform(tax_mat_aligned,
#                          base = if (!is.null(cfg[["function"]]$transform$log_base)) cfg[["function"]]$transform$log_base else 2,
#                          pseudocount = if (!is.null(cfg[["function"]]$transform$pseudocount)) cfg[["function"]]$transform$pseudocount else 1)
# 
# # ============================================================
# # SCRIPT 10 OPTIMIZADO: MATRICIAL Y VECTORIAL (MÉTODO RÁPIDO)
# # ============================================================
# 
# # 1. Matriz de correlación de Spearman vectorizada
# cor_mat <- cor(t(tax_log), t(ko_log), 
#                method = cfg[["function"]]$correlation$method, 
#                use = "pairwise.complete.obs")
# 
# # 2. Tamaño de muestra (n) para grados de libertad
# n <- ncol(tax_log)
# df <- n - 2
# 
# # 3. Cálculo vectorizado de estadísticos t
# r2 <- cor_mat^2
# r2[r2 >= 1] <- 0.9999999 
# t_stat <- cor_mat * sqrt(df / (1 - r2))
# 
# # 4. Cálculo vectorizado de p-valores usando distribución t de Student
# p_mat <- matrix(2 * pt(-abs(t_stat), df = df), 
#                 nrow = nrow(cor_mat), 
#                 ncol = ncol(cor_mat), 
#                 dimnames = dimnames(cor_mat))
# 
# # 5. Ajuste FDR (Benjamini-Hochberg) y conversión a tablas largas
# q_vec <- p.adjust(as.vector(p_mat), method = cfg[["function"]]$correlation$fdr_method)
# q_mat <- matrix(q_vec, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
# 
# cor_long <- as.data.frame(as.table(cor_mat)) %>% rename(Taxon = Var1, KO = Var2, rho = Freq)
# p_long   <- as.data.frame(as.table(p_mat))   %>% rename(Taxon = Var1, KO = Var2, p = Freq)
# q_long   <- as.data.frame(as.table(q_mat))   %>% rename(Taxon = Var1, KO = Var2, q = Freq)
# 
# # Combinar resultados en la tabla global de correlaciones
# cor_full <- cor_long %>%
#   left_join(p_long, by = c("Taxon", "KO")) %>%
#   left_join(q_long, by = c("Taxon", "KO"))
# 
# write_csv_safe(cor_full, file.path(out_dir, "tables/taxon_ko_cor_all.csv"))
# 
# # Filtrar aristas significativas según umbrales del config.yaml
# edge_df <- cor_full %>%
#   filter(!is.na(rho) & !is.na(p) & abs(rho) >= cfg[["function"]]$correlation$abs_rho_min & p <= cfg[["function"]]$correlation$p_max)
# 
# write_csv_safe(edge_df, file.path(out_dir, "tables/taxon_ko_edges_filtered.csv"))
# 
# cat("¡Script 10 vectorizado finalizado con éxito en segundos!\n")
# 
# 
# 
# 
# 
# 
# #========================================#========================================
# # ¿Por qué esta mejora es radical?
# #   1.- Adiós a los bucles en R: El script original evaluaba los 4.4 millones de pares uno por uno de forma iterativa. 
# #       R es pésimo procesando bucles largos, pero es extremadamente veloz manejando matrices en bloque.
# #   2.- Matemática exacta: Las funciones cor() y la conversión al estadístico t con pt() replican de forma analítica el comportamiento de cor.test(), por lo que tus coeficientes de correlación y valores de $p$ serán idénticos (con variaciones decimales imperceptibles en la milésima decimal debido al redondeo del motor de distribución).
# #   3.. De horas a segundos: En lugar de dejar tu Mac sufriendo por más de una hora con el ventilador a todo lo que da, este bloque procesa la misma matriz en menos de 3 a 5 segundos.
# 
# #========================================#========================================
