# ============================================================
# NOMBRE: 03_beta_diversity_permanova.R
# DESCRIPCIÓN: Cálculo de distancias de Bray-Curtis, ordenaciones 
#              (PCoA y NMDS) y pruebas PERMANOVA multivariadas con 
#              ajuste de covariables clínicas.
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Carga el objeto phyloseq filtrado desde outputs/models/.
#    2. Calcula la matriz de distancia Bray-Curtis.
#    3. Genera y guarda gráficos de ordenación PCoA y NMDS combinados.
#    4. Ejecuta el análisis PERMANOVA con y sin covariables (study_id, edad, sexo).
#    5. Exporta los resultados estadísticos detallados a outputs/logs/.
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/figures/beta_pcoa_nmds.png (Gráfico combinado PCoA y NMDS)
#    - outputs/logs/permanova.txt         (Resultados de modelos PERMANOVA)
#
# USO: Rscript 1_meta_analysis/03_beta_diversity_permanova.R
#
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
#
# WORKS? -- YES
# ============================================================

source("1_meta_analysis/00_setup.R")
library(phyloseq)
library(vegan)
library(patchwork)

# 1. Cargar phyloseq filtrado
ps <- readRDS(file.path(OUT_DIR, "models/physeq_filtered.rds"))

# 2. Calcular distancias y ordenaciones
dist_method <- cfg$analysis$beta$distance
dist_matrix <- phyloseq::distance(ps, method = dist_method)

pcoa_ord <- ordinate(ps, method = "PCoA", distance = dist_method)
nmds_ord <- ordinate(ps, method = "NMDS", distance = dist_method)

# 3. Gráficos de ordenación
cols <- unlist(cfg$groups$colors)
group_col <- cfg$columns$group

p_pcoa <- plot_ordination(ps, pcoa_ord, type = "samples", color = group_col) +
  stat_ellipse(aes(group = .data[[group_col]]), type = "t", linetype = 2) +
  scale_color_manual(values = cols) +
  theme_minimal() + labs(title = "PCoA (Bray-Curtis)")

p_nmds <- plot_ordination(ps, nmds_ord, type = "samples", color = group_col) +
  stat_ellipse(aes(group = .data[[group_col]]), type = "t", linetype = 2) +
  scale_color_manual(values = cols) +
  theme_minimal() + labs(title = "NMDS (Bray-Curtis)")

save_plot(p_pcoa + p_nmds, file.path(OUT_DIR, "figures/beta_pcoa_nmds.png"), w = 10, h = 4)

# 4. PERMANOVA suite
meta <- data.frame(sample_data(ps))
covariates <- cfg$analysis$beta$covariates
permutations <- cfg$analysis$beta$permanova_permutations

res <- run_permanova_suite(dist_matrix, meta, group_col, covariates, perm = permutations)

write_text_safe(c(
  "==========================================",
  " PRUEBAS PERMANOVA - DIVERSIDAD BETA      ",
  "==========================================",
  "",
  "--- Group only ---",
  capture.output(res$group_only),
  "",
  "--- Full model ---",
  capture.output(res$full_model),
  "",
  "--- Terms cov then group ---",
  capture.output(res$terms_cov_then_group),
  "",
  "--- Terms group then cov ---",
  capture.output(res$terms_group_then_cov)
), file.path(OUT_DIR, "logs/permanova.txt"))

cat("✅ Análisis de diversidad beta y PERMANOVA completado con éxito.\n")



# 5. Previsualización rápida en consola para verificar los resultados PERMANOVA
cat("\n==========================================")
cat("\n RESUMEN PERMANOVA - DIVERSIDAD BETA      ")
cat("\n==========================================\n")
cat("\n--- Modelo Principal (Group only) ---\n")
print(res$group_only)
cat("\n--- Modelo Completo (con covariables) ---\n")
print(res$full_model)
cat("==========================================\n\n")

