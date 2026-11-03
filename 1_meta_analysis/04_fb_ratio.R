# ============================================================
# NOMBRE: 04_fb_ratio.R
# DESCRIPCIÓN: Cálculo de la razón Firmicutes (Bacillota) a 
#              Bacteroidota por muestra, pruebas estadísticas y 
#              generación de boxplots comparativos.
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Carga el phyloseq y extrae perfiles taxonómicos a nivel Phylum.
#    2. Calcula la razón F:B para cada muestra de manera robusta.
#    3. Genera un gráfico de cajas con valores de dispersión (jitter).
#    4. Realiza t-test y exporta tanto la tabla de valores como los logs.
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/tables/fb_ratio.csv          (Valores de la razón F:B por muestra)
#    - outputs/figures/fb_ratio.png         (Gráfico de cajas comparativo)
#    - outputs/logs/fb_ratio_stats.txt      (Resultados estadísticos de la prueba t)
#
# USO: Rscript 1_meta_analysis/04_fb_ratio.R
#
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
#
# WORKS? -- ???
# ============================================================
# EN TERMINAL: Rscript -e "install.packages('ggpubr', repos='https://cloud.r-project.org')"
# OR: conda install -c conda-forge r-ggpubr -y


source("1_meta_analysis/00_setup.R")
library(phyloseq)
library(ggpubr)

ps <- readRDS(file.path(OUT_DIR, "models/physeq_filtered.rds"))

otu <- as.data.frame(otu_table(ps), check.names = FALSE)
if (!taxa_are_rows(ps)) otu <- t(otu)
tax <- as.data.frame(tax_table(ps), check.names = FALSE)
meta <- as.data.frame(sample_data(ps))

ph_a <- cfg$analysis$fb_ratio$phylum_a  # "Firmicutes"
ph_b <- cfg$analysis$fb_ratio$phylum_b  # "Bacteroidota"

# Cálculo de la razón F:B por muestra
bf_ratio <- apply(otu, 2, function(x) {
  b <- sum(x[which(tax[, "Phylum"] == ph_b)], na.rm = TRUE)
  f <- sum(x[which(tax[, "Phylum"] == ph_a)], na.rm = TRUE)
  ifelse(b == 0, NA, f / b)
})

# Alineación limpia de muestras entre OTU table y metadatos
common_samples <- intersect(names(bf_ratio), rownames(meta))
group_col <- cfg$columns$group

# Extracción robusta garantizando vector plano de caracteres
group_vals <- meta[common_samples, group_col]
if (is.data.frame(group_vals)) group_vals <- group_vals[[1]]
group_vals <- trimws(as.character(group_vals))

df <- data.frame(
  Sample = common_samples,
  FB_ratio = bf_ratio[common_samples],
  Group_raw = group_vals,
  stringsAsFactors = FALSE
)

# Mapeo explícito tanto para códigos numéricos ("1"/"2") como texto ("Control"/"Insomnia")
df$Group <- ifelse(df$Group_raw == "1" | tolower(df$Group_raw) == "control", "Control",
                   ifelse(df$Group_raw == "2" | tolower(df$Group_raw) == "insomnia", "Insomnia", NA))

# Filtrado de NAs y asignación de niveles de factor
df <- df[!is.na(df$FB_ratio) & !is.na(df$Group), ]
df$Group <- factor(df$Group, levels = cfg$groups$levels)

cat("\nConteo de muestras por grupo (mapeo correcto):\n")
print(table(df$Group))

write_csv_safe(df, file.path(OUT_DIR, "tables/fb_ratio.csv"))

# Gráfico
cols <- unlist(cfg$groups$colors)

p <- ggplot(df, aes(x = Group, y = FB_ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.format",
                     #label.y = max(df$FB_ratio, na.rm = TRUE) * 1.05) +
                     label.y = max(df$FB_ratio, na.rm = TRUE) * 1.10) +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  labs(title = "Firmicutes:Bacteroidota ratio", y = "F:B ratio")

save_plot(p, file.path(OUT_DIR, "figures/fb_ratio.png"), w = 6, h = 4)

# Estadísticas
st <- t.test(FB_ratio ~ Group, data = df)
write_text_safe(capture.output(st), file.path(OUT_DIR, "logs/fb_ratio_stats.txt"))

cat("\n✅ Análisis de ratio Firmicutes:Bacteroidota completado con éxito.\n")

cat("\n==========================================")
cat("\n RESUMEN DE LA RAZÓN FIRMICUTES:BACTEROIDOTA")
cat("\n==========================================\n")
print(summary(df$FB_ratio))
cat("\nResultados de la prueba t:\n")
print(st)
cat("==========================================\n\n")

