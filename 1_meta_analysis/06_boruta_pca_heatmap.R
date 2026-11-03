# ============================================================
# NOMBRE: 06_boruta_pca_heatmap.R
# DESCRIPCIÓN: Análisis de componentes principales (PCA) y mapas 
#              de calor (Heatmap Z-score) basados exclusivamente en 
#              los ASVs seleccionados por Boruta.
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Carga los ASVs seleccionados por Boruta y el objeto phyloseq.
#    2. Ejecuta un PCA y genera el gráfico de dispersión por grupo.
#    3. Genera un Heatmap escalado (Z-score) ordenado por el grupo Control.
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/boruta/boruta_pca.png             (Gráfico PCA con ASVs clave)
#    - outputs/boruta/boruta_heatmap_samples.png (Heatmap de expresión/abundancia)
#
# USO: Rscript 1_meta_analysis/06_boruta_pca_heatmap.R
#
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
#
# WORKS? -- YES!!! (CONFIRMED)
# ============================================================

source("1_meta_analysis/00_setup.R")
library(tibble)
library(ggplot2)

out_dir <- file.path(OUT_DIR, "boruta")
sel_file <- file.path(out_dir, "boruta_selected_otu.csv")

if (!file.exists(sel_file)) {
  stop("❌ No se encontró boruta_selected_otu.csv. Ejecuta primero 05_boruta.R")
}

sel <- readr::read_csv(sel_file, show_col_types = FALSE) %>%
  tibble::column_to_rownames("Sample")

#should be numeric — but if Boruta's output CSV includes the Group or covariate columns (it shouldn't, but check), prcomp will error. Add a guard:
# sel <- sel[, sapply(sel, is.numeric), drop = FALSE]

ps <- readRDS(file.path(OUT_DIR, "models/physeq_filtered.rds"))
meta <- as.data.frame(sample_data(ps))
group_var <- cfg$columns$group

# 1. PCA
pca_res <- prcomp(sel, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Sample <- rownames(pca_df)

raw_group_pca <- meta[rownames(pca_df), group_var]
if (is.data.frame(raw_group_pca) || inherits(raw_group_pca, "sample_data")) {
  raw_group_pca <- raw_group_pca[[1]]
}
pca_df$Group <- as.factor(as.character(raw_group_pca))

color_values <- unlist(cfg$groups$colors)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA using Boruta-selected ASVs",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
  ) +
  scale_color_manual(values = color_values)

save_plot(p_pca, file.path(out_dir, "boruta_pca.png"), w = 6, h = 5)

# 2. Heatmap
required <- c("pheatmap", "RColorBrewer")
miss <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Faltan paquetes para el heatmap: ", paste(miss, collapse = ", "))
}

library(pheatmap)
library(RColorBrewer)

otu_sd <- apply(sel, 2, sd, na.rm = TRUE)
sel2 <- sel[, otu_sd != 0, drop = FALSE]

otu_z <- t(scale(sel2))
otu_z[!is.finite(otu_z)] <- NA
otu_z <- otu_z[complete.cases(otu_z), , drop = FALSE]

sample_ids <- colnames(otu_z)

# Safe extraction for annotations
raw_group_hm <- meta[sample_ids, group_var]
if (is.data.frame(raw_group_hm) || inherits(raw_group_hm, "sample_data")) {
  raw_group_hm <- raw_group_hm[[1]]
}

annotation_col <- data.frame(Group = factor(as.character(raw_group_hm)))
rownames(annotation_col) <- sample_ids

sample_order <- sample_ids[order(annotation_col$Group)]
otu_z2 <- otu_z[, sample_order, drop = FALSE]
ann2 <- annotation_col[sample_order, , drop = FALSE]

ctrl_name <- cfg$groups$levels[1]
ctrl_samples <- rownames(ann2)[ann2$Group == ctrl_name]
if (length(ctrl_samples) > 0) {
  ctrl_means <- rowMeans(otu_z2[, ctrl_samples, drop = FALSE], na.rm = TRUE)
  otu_order <- order(ctrl_means, decreasing = TRUE)
  otu_z3 <- otu_z2[otu_order, , drop = FALSE]
} else {
  otu_z3 <- otu_z2
}

ann_colors <- list(Group = color_values)
breaks_list <- seq(-0.7, 0.7, length.out = 101)

png(file.path(out_dir, "boruta_heatmap_samples.png"), width = 1800, height = 1200, res = 200)
pheatmap(
  otu_z3,
  annotation_col = ann2,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  scale = "none",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = breaks_list,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  main = "Heatmap of Boruta-selected ASVs (Z-score)"
)
dev.off()

cat("✅ PCA y Heatmap de Boruta generados con éxito en outputs/boruta/.\n\n")
