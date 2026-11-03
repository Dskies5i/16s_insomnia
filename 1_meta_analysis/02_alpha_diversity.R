# ============================================================
# NOMBRE: 02_alpha_diversity.R
# DESCRIPCIÓN: Calcula la diversidad alfa (Shannon, Simpson, etc.) 
#              a partir del objeto phyloseq filtrado, genera gráficos 
#              de violín/caja comparativos entre grupos y ejecuta pruebas 
#              estadísticas (t-test).
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Carga el objeto phyloseq filtrado desde outputs/models/.
#    2. Estima los índices de diversidad alfa configurados en config.yaml.
#    3. Combina los valores de diversidad con los metadatos clínicos.
#    4. Genera y guarda gráficos combinados (Shannon y Simpson) en outputs/figures/.
#    5. Realiza pruebas estadísticas t-test y exporta los logs a outputs/logs/.
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/tables/alpha_diversity.csv        (Tabla completa con métricas de diversidad por muestra)
#    - outputs/figures/alpha_shannon_simpson.png (Gráfico combinado en formato violín/boxplots)
#    - outputs/logs/alpha_stats.txt              (Resultados detallados de las pruebas estadísticas t-test)
#
# USO: Rscript 1_meta_analysis/02_alpha_diversity.R
#
# PARA ESTE SCRIPT EN TERMINAL: Rscript -e "install.packages(c('patchwork', 'ggpubr', 'Boruta', 'pheatmap', 'RColorBrewer'), repos='https://cloud.r-project.org')"
# 
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
#
# WORKS? -- YES 
#============================================================
source("1_meta_analysis/00_setup.R")
library(phyloseq)
library(patchwork)

# 1. Cargar el objeto phyloseq filtrado
ps_path <- file.path(OUT_DIR, "models/physeq_filtered.rds")
if (!file.exists(ps_path)) {
  stop("❌ ERROR: No se encontró el archivo physeq_filtered.rds. Ejecuta primero 01_build_phyloseq.R")
}
ps <- readRDS(ps_path)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# 2. Estimar diversidad alfa con base en las medidas del archivo de configuración
alpha_measures <- cfg$analysis$alpha$measures
alpha_tab <- estimate_richness(ps, measures = alpha_measures)
alpha_tab$Sample <- rownames(alpha_tab)

# 3. Incorporar metadatos (Grupo de estudio)
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("Sample")
group_col <- cfg$columns$group

alpha_tab <- alpha_tab %>%  
  left_join(meta %>% select(Sample, Group = all_of(group_col)), by = "Sample")

# Forzar el orden de los niveles de los grupos según la configuración
groups_levels <- cfg$groups$levels
alpha_tab$Group <- factor(alpha_tab$Group, levels = groups_levels)

# Guardar tabla completa de diversidad alfa
write_csv_safe(alpha_tab, file.path(OUT_DIR, "tables/alpha_diversity.csv"))

# 4. Generación de gráficos de violín y cajas
cols <- unlist(cfg$groups$colors)

p_shannon <- ggplot(alpha_tab, aes(x = Group, y = Shannon, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "white") +
  scale_fill_manual(values = cols) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold")) +
  labs(title = "Shannon Diversity", x = "", y = "Shannon Index")

p_simpson <- ggplot(alpha_tab, aes(x = Group, y = Simpson, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "white") +
  scale_fill_manual(values = cols) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold")) +
  labs(title = "Simpson Diversity", x = "", y = "Simpson Index")

# Guardar gráfico combinado mediante patchwork
save_plot(p_shannon + p_simpson, file.path(OUT_DIR, "figures/alpha_shannon_simpson.png"), w = 9, h = 4)

# 5. Ejecución de pruebas estadísticas (t-test) y guardado de logs
sh_t <- t.test(Shannon ~ Group, data = alpha_tab)
si_t <- t.test(Simpson ~ Group, data = alpha_tab)

stats_output <- c(
  "==========================================",
  " RESUMEN DE PRUEBAS ESTADÍSTICAS - DIVERSIDAD ALFA",
  "==========================================",
  "",
  "--- Shannon Index (t-test) ---",
  capture.output(print(sh_t)),
  "",
  "--- Simpson Index (t-test) ---",
  capture.output(print(si_t)),
  "=========================================="
)

write_text_safe(stats_output, file.path(OUT_DIR, "logs/alpha_diversity_stats.txt"))

cat("✅ Análisis de diversidad alfa completado con éxito. Resultados guardados en outputs/\n")

# 6. Previsualización rápida en consola
cat("\n==========================================")
cat("\n RESUMEN DE DIVERSIDAD ALFA (Primeras filas)")
cat("\n==========================================\n")
print(head(alpha_tab, 5))
cat("---------------------\nResultados rápidos de las pruebas t-test:\n")
print(sh_t)
print(si_t)
cat("==========================================\n\n")


