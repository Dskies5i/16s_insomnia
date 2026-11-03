
# ============================================================
# NOMBRE: 05_boruta.R
# DESCRIPCIÓN: Selección de características (features) importantes 
#              mediante el algoritmo de Machine Learning Boruta 
#              para discriminar entre grupos de estudio.
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Prepara las matrices de metadatos y abundancia de ASVs.
#    2. Ejecuta el algoritmo Boruta incorporando covariables.
#    3. Extrae los ASVs relevantes identificados como importantes.
#    4. Anota taxonómicamente los ASVs seleccionados y genera gráfico.
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/boruta/boruta_model.rds                  (Modelo Boruta guardado)
#    - outputs/boruta/boruta_selected_otu.csv           (Matriz de ASVs seleccionados)
#    - outputs/boruta/boruta_selected_otus_with_taxonomy.csv (Taxonomía de ASVs clave)
#    - outputs/boruta/boruta_importance.png             (Gráfico de importancia de variables)
#
# USO: Rscript 1_meta_analysis/05_boruta.R
#
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
#
# WORKS? -- YES !!!! 
#
# REQUIRED INSTALL: conda install -c conda-forge r-boruta -y 
#  (IN TERMINAL - INSIDE ENVIRONMENT)
# ============================================================

source("1_meta_analysis/00_setup.R")
library(Boruta)
library(tibble)

out_dir <- file.path(OUT_DIR, "boruta")
dir_create(out_dir)

# Cargar phyloseq filtrado
ps <- readRDS(file.path(OUT_DIR, "models/physeq_filtered.rds"))

# Extraer OTUs y metadatos limpios
otu_mat <- as.data.frame(otu_table(ps), check.names = FALSE)
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
X <- as.data.frame(t(otu_mat))

# Aplicar prefijo seguro para evitar errores en fórmulas de R si los IDs son numéricos
otu_prefix <- cfg$differential$boruta$otu_prefix
colnames(X) <- paste0(otu_prefix, colnames(X))

meta <- as.data.frame(sample_data(ps))

target <- cfg$columns$group
covs <- cfg$differential$boruta$covariates

# Ajustar tipos de datos
meta[[target]] <- as.factor(meta[[target]])
valid_covs <- intersect(covs, colnames(meta))
for (cc in valid_covs) {
  meta[[cc]] <- as.factor(meta[[cc]])
}

Y <- meta[[target]]
X_full <- cbind(meta[, valid_covs, drop = FALSE], X)

df_boruta <- data.frame(Group = Y, X_full)
colnames(df_boruta)[1] <- "Group"

# Ejecutar Boruta
boruta_res <- Boruta(
  Group ~ .,
  data = df_boruta,
  doTrace = cfg$differential$boruta$doTrace,
  maxRuns = cfg$differential$boruta$maxRuns
)

saveRDS(boruta_res, file.path(out_dir, "boruta_model.rds"))

# Extracción de variables seleccionadas
final_vars <- getSelectedAttributes(boruta_res, withTentative = FALSE)
final_otus <- final_vars[grepl(paste0("^", otu_prefix), final_vars)]

if(length(final_otus) > 0) {
  selected_otu_table <- X[, final_otus, drop = FALSE]
  write_csv_safe(
    selected_otu_table %>% tibble::rownames_to_column("Sample"),
    file.path(out_dir, "boruta_selected_otu.csv")
  )
  
  # Anotación taxonómica (quitando el prefijo temporal para hacer el join)
  taxonomy <- as.data.frame(tax_table(ps)) %>% tibble::rownames_to_column("OTU_ID")
  #boruta_otus_df <- data.frame(OTU_ID = gsub(otu_prefix, "", final_otus))
  #AFTER
  boruta_otus_df <- data.frame(OTU_ID = sub(paste0("^", otu_prefix), "", final_otus))
  boruta_annot <- dplyr::left_join(boruta_otus_df, taxonomy, by = "OTU_ID") %>%
    dplyr::select(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)
  
  write_csv_safe(boruta_annot, file.path(out_dir, "boruta_selected_otus_with_taxonomy.csv"))
} else {
  warning("⚠️ Boruta no seleccionó ningún ASV con los parámetros actuales.")
}

cat("✅ Selección de características con Boruta finalizada con éxito.\n")
cat(sprintf("\nASVs seleccionados como importantes: %d\n\n", length(final_otus)))




# ============================================================
# VISUALIZACIÓN MEJORADA DE BORUTA (ggplot2 horizontal)
# ============================================================

library(ggplot2)

# 1. Extraer importancia de todas las decisiones
imp_mat <- attStats(boruta_res)
imp_mat$Feature <- rownames(imp_mat)

# 2. Filtrar solo confirmados, tentativas y sombras de referencia
imp_filtered <- imp_mat[imp_mat$decision %in% c("Confirmed", "Tentative") | 
                          grepl("^shadow", imp_mat$Feature), ]

# 3. Mapear con la taxonomía para nombres legibles
taxonomy <- as.data.frame(tax_table(ps)) %>% tibble::rownames_to_column("OTU_ID")
imp_filtered$OTU_ID <- gsub(otu_prefix, "", imp_filtered$Feature)

imp_filtered <- dplyr::left_join(imp_filtered, taxonomy, by = "OTU_ID")

# Crear etiqueta legible (Género/Especie o la sombra)
imp_filtered$Label <- ifelse(grepl("^shadow", imp_filtered$Feature), 
                             imp_filtered$Feature,
                             ifelse(!is.na(imp_filtered$Genus) & imp_filtered$Genus != "g__other",
                                    paste0(imp_filtered$Genus, " (", substr(imp_filtered$OTU_ID, 1, 6), ")"),
                                    paste0("ASV_", substr(imp_filtered$OTU_ID, 1, 8))))

# Ordenar por importancia media
imp_filtered$Label <- factor(imp_filtered$Label, levels = imp_filtered$Label[order(imp_filtered$meanImp)])

# 4. Graficar con ggplot2
p_boruta <- ggplot(imp_filtered, aes(x = meanImp, y = Label, fill = decision)) +
  geom_col(alpha = 0.85, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Confirmed" = "#2b5c8f", 
                               "Tentative" = "#e6ab02", 
                               "Shadow"    = "#7570b3",
                               "Rejected"  = "#d95f02")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Boruta Feature Importance (Confirmed ASVs)",
    x = "Mean Importance Score",
    y = "Taxon / Feature",
    fill = "Decision"
  ) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_blank()
  )

# Guardar figura en alta resolución
save_plot(p_boruta, file.path(out_dir, "boruta_importance_clean.png"), w = 8, h = 7)

# Abrir automáticamente en macOS
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", file.path(out_dir, "boruta_importance_clean.png"))
}

cat("✅ Gráfico limpio de importancia guardado en outputs/boruta/boruta_importance_clean.png\n")



# 
# ##============================================================
# source("1_meta_analysis/00_setup.R")
# library(Boruta)
# library(tibble)
# 
# out_dir <- file.path(OUT_DIR, "boruta")
# dir_create(out_dir)
# 
# # Cargar phyloseq filtrado
# ps <- readRDS(file.path(OUT_DIR, "models/physeq_filtered.rds"))
# 
# # Extraer OTUs y metadatos limpios
# otu_mat <- as.data.frame(otu_table(ps))
# if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
# X <- as.data.frame(t(otu_mat))
# 
# meta <- as.data.frame(sample_data(ps))
# 
# target <- cfg$columns$group
# covs <- cfg$analysis$beta$covariates
# 
# # Ajustar tipos de datos
# meta[[target]] <- as.factor(meta[[target]])
# for (cc in covs) {
#   if (cc %in% colnames(meta)) meta[[cc]] <- as.factor(meta[[cc]])
# }
# 
# Y <- meta[[target]]
# valid_covs <- intersect(covs, colnames(meta))
# X_full <- cbind(meta[, valid_covs, drop = FALSE], X)
# 
# df_boruta <- data.frame(Group = Y, X_full)
# colnames(df_boruta)[1] <- "Group"
# 
# # Ejecutar Boruta
# boruta_res <- Boruta(
#   Group ~ .,
#   data = df_boruta,
#   doTrace = 2,
#   maxRuns = 100
# )
# 
# saveRDS(boruta_res, file.path(out_dir, "boruta_model.rds"))
# 
# # Extracción de variables seleccionadas
# final_vars <- getSelectedAttributes(boruta_res, withTentative = FALSE)
# final_otus <- final_vars[final_vars %in% colnames(X)]
# 
# if(length(final_otus) > 0) {
#   selected_otu_table <- X[, final_otus, drop = FALSE]
#   write_csv_safe(
#     selected_otu_table %>% tibble::rownames_to_column("Sample"),
#     file.path(out_dir, "boruta_selected_otu.csv")
#   )
#   
#   # Anotación taxonómica
#   taxonomy <- as.data.frame(tax_table(ps)) %>% tibble::rownames_to_column("OTU_ID")
#   boruta_otus_df <- data.frame(OTU_ID = final_otus)
#   boruta_annot <- dplyr::left_join(boruta_otus_df, taxonomy, by = "OTU_ID") %>%
#     dplyr::select(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)
#   
#   write_csv_safe(boruta_annot, file.path(out_dir, "boruta_selected_otus_with_taxonomy.csv"))
# } else {
#   warning("⚠️ Boruta no seleccionó ningún ASV con los parámetros actuales.")
# }
# 
# # Gráfico de importancia
# png(file.path(out_dir, "boruta_importance.png"), width = 1800, height = 900, res = 200)
# plot(boruta_res, las = 2, cex.axis = 0.6, main = "Boruta Feature Importance")
# dev.off()
# 
# cat("✅ Selección de características con Boruta finalizada con éxito.\n")
# 
# 
# cat("\n==========================================")
# cat("\n RESUMEN DE SELECCIÓN BORUTA            ")
# cat("\n==========================================\n")
# print(boruta_res)
# cat(sprintf("\nASVs seleccionados como importantes: %d\n", length(final_otus)))
# if(length(final_otus) > 0) {
#   print(head(boruta_annot, 5))
# }
# cat("==========================================\n\n")