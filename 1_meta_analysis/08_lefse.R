# ============================================================
# NOMBRE: 08_lefse.R
# DESCRIPCIÓN: Identificación de biomarcadores diferenciales con 
#              LEfSe (microbiomeMarker) y gráficos de barras firmados.
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Carga el objeto phyloseq filtrado y normaliza las abundancias 
#       a proporciones relativas.
#    2. Ejecuta el algoritmo LEfSe (run_lefse) para identificar 
#       biomarcadores taxonómicos diferenciales entre los grupos de estudio.
#    3. Enriquece los resultados con las rutas taxonómicas completas 
#       del objeto original mediante unión por hoja taxonómica.
#    4. Filtra los biomarcadores basándose en los umbrales de rango taxonómico,
#       valor de p, puntuación LDA mínima y patrones a excluir definidos en config.yaml.
#    5. Genera un gráfico de barras con puntuaciones LDA firmadas y colores 
#       personalizados por grupo.
#
# ARCHIVOS DE ENTRADA (INPUTS):
#    - outputs/models/physeq_filtered.rds      (Objeto phyloseq filtrado)
#    - config/config.yaml                      (Archivo de configuración central)
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/lefse/lefse_result.csv          (Resultados completos del análisis LEfSe)
#    - outputs/lefse/lefse_sig_with_otu.csv    (Biomarcadores significativos filtrados con taxonomía)
#    - outputs/lefse/lefse_barplot.png         (Gráfico de barras con puntuaciones LDA firmadas)
#
# USO: Rscript 1_meta_analysis/08_lefse.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================


source("1_meta_analysis/00_setup.R")
source("utils/utils_phyloseq.R")
source("utils/utils_io.R")
source("utils/utils_plot.R")      # ← ADDED: save_plot lives here

library(phyloseq)
library(tibble)
library(dplyr)
library(ggplot2)

cfg <- yaml::read_yaml("config/config.yaml")
out_dir <- file.path(cfg$outputs$dir, "lefse")
dir_create(out_dir)

required <- c("microbiomeMarker")
miss <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) stop("Falta el paquete obligatorio: microbiomeMarker")
library(microbiomeMarker)

ps_path <- file.path(cfg$outputs$dir, "models/physeq_filtered.rds")
if (!file.exists(ps_path)) {
  stop("❌ No se encontró physeq_filtered.rds en: ", ps_path)
}
ps <- readRDS(ps_path)

p_comp <- transform_sample_counts(ps, function(x) x / sum(x))

# ---- LEFSe params ----
grp          <- cfg$differential$lefse$group
rank_filter  <- cfg$differential$lefse$rank_filter
p_max        <- cfg$differential$lefse$p_max
lda_min      <- cfg$differential$lefse$lda_min
drop_pattern <- cfg$differential$lefse$drop_tax_pattern

cat("🔄 Ejecutando LEfSe (puede tardar varios minutos)...\n")
res_marker <- run_lefse(
  p_comp,
  group      = grp,
  taxa_rank  = rank_filter,
  kw_cutoff  = p_max,
  lda_cutoff = lda_min
)

# Extract plain data.frame from S4 marker_table
res <- as.data.frame(res_marker@marker_table)
class(res) <- "data.frame"
res$rank <- rank_filter

write_csv_safe(res, file.path(out_dir, "lefse_result.csv"))

# ---- Taxonomy annotations ----
tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("OTU_ID")

rank_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_cols <- rank_cols[rank_cols %in% colnames(tax_df)]

tax_df$tax_path <- apply(tax_df[, rank_cols, drop = FALSE], 1, function(x) {
  x[is.na(x)] <- "Unassigned"
  paste(trimws(x), collapse = "|")
})

# ---- Dynamic column names (microbiomeMarker version-safe) ----
col_p   <- if ("pvalue"       %in% colnames(res)) "pvalue"       else "p.value"
col_lda <- if ("ef_lda"       %in% colnames(res)) "ef_lda"       else "LDAscore"
col_tax <- if ("feature"      %in% colnames(res)) "feature"      else "tax"
col_dir <- if ("enrich_group" %in% colnames(res)) "enrich_group" else "direction"

# ---- Filter ----
lefse_filtered <- res %>%
  dplyr::filter(.data[[col_p]] <= p_max, .data[[col_lda]] >= lda_min)

if (!is.null(drop_pattern) && nzchar(drop_pattern)) {
  lefse_filtered <- lefse_filtered %>%
    dplyr::filter(is.na(.data[[col_tax]]) | !grepl(drop_pattern, .data[[col_tax]]))
}

# ---- Strip rank prefix on BOTH sides (double underscore — your SILVA format) ----
strip_rank_prefix <- function(x) sub("^[a-z]__", "", as.character(x))

lefse_filtered$tax_leaf <- strip_rank_prefix(lefse_filtered[[col_tax]])
tax_df$tax_leaf         <- strip_rank_prefix(tax_df[[rank_filter]])

# Deduplicate tax_df on leaf to avoid join fan-out (multiple ASVs → same species)
tax_df_leaf <- tax_df %>%
  dplyr::select(OTU_ID, dplyr::all_of(rank_cols), tax_leaf, tax_path) %>%
  dplyr::distinct(tax_leaf, .keep_all = TRUE)

lefse_otu <- dplyr::left_join(lefse_filtered, tax_df_leaf, by = "tax_leaf")

cat(sprintf("NAs en %s tras join: %d / %d\n",
            rank_filter, sum(is.na(lefse_otu[[rank_filter]])), nrow(lefse_otu)))

write_csv_safe(lefse_otu, file.path(out_dir, "lefse_sig_with_otu.csv"))

# ---- Plot ----
if (nrow(lefse_otu) > 0) {
  pos_direction <- cfg$groups$levels[2]
  lefse_otu[[col_dir]] <- as.character(lefse_otu[[col_dir]])
  
  if (!tolower(pos_direction) %in% tolower(unique(lefse_otu[[col_dir]]))) {
    warning("No markers for positive group '", pos_direction, "'. Skipping signed barplot.")
  } else {
    lefse_otu <- lefse_otu %>%
      mutate(
        LDAscore   = .data[[col_lda]],
        tax_name   = .data[[col_tax]],
        direction  = .data[[col_dir]],
        LDA_signed = ifelse(tolower(direction) == tolower(pos_direction),
                            LDAscore, -LDAscore)
      ) %>%
      arrange(LDA_signed) %>%
      mutate(tax_name = factor(tax_name, levels = tax_name))
    
    p <- ggplot(lefse_otu, aes(x = LDA_signed, y = tax_name, fill = direction)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      theme_minimal(base_size = 12) +
      labs(
        title    = "LEfSe Biomarkers (Filtered)",
        subtitle = sprintf("Rank: %s | LDA > %s | p < %s", rank_filter, lda_min, p_max),
        x = "Signed LDA Score", y = "Taxon"
      ) +
      scale_fill_manual(values = unlist(cfg$groups$colors)) +
      theme(axis.text.y = element_text(size = 9))
    
    save_plot(p, file.path(out_dir, "lefse_barplot.png"),
              w = 9, h = max(6, nrow(lefse_otu) * 0.25))
    
    cat(sprintf("✅ %d biomarcadores significativos. Gráfico guardado.\n", nrow(lefse_otu)))
  }
} else {
  warning("⚠️ No hay biomarcadores significativos con los filtros actuales.")
}

cat("✅ Análisis LEfSe finalizado en outputs/lefse/.\n")









# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================



# ============================================================

# source("1_meta_analysis/00_setup.R")
# source("utils/utils_phyloseq.R")
# source("utils/utils_io.R")
# 
# library(phyloseq)
# library(tibble)
# library(dplyr)
# library(ggplot2)
# 
# # 1. Cargar configuración y definir rutas de forma robusta
# cfg <- yaml::read_yaml("config/config.yaml")
# out_dir <- file.path(cfg$outputs$dir, "lefse")
# dir_create(out_dir)
# 
# # 2. Verificar disponibilidad de microbiomeMarker
# required <- c("microbiomeMarker")
# miss <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
# if (length(miss) > 0) stop("Falta el paquete obligatorio: microbiomeMarker")
# 
# library(microbiomeMarker)
# 
# # 3. Cargar objeto phyloseq filtrado
# ps_path <- file.path(cfg$outputs$dir, "models/physeq_filtered.rds")
# if (!file.exists(ps_path)) {
#   stop("❌ No se encontró el archivo physeq_filtered.rds en: ", ps_path, 
#        "\nAsegúrate de ejecutar el script de filtrado previo.")
# }
# 
# ps <- readRDS(ps_path)
# 
# # Normalizar a proporciones/abundancia relativa para LEfSe
# p_comp <- transform_sample_counts(ps, function(x) x / sum(x))
# 
# # 4. Ejecutar análisis LEfSe
# grp <- cfg$differential$lefse$group
# rank_filter <- cfg$differential$lefse$rank_filter
# p_max <- cfg$differential$lefse$p_max
# lda_min <- cfg$differential$lefse$lda_min
# drop_pattern <- cfg$differential$lefse$drop_tax_pattern
# 
# cat("🔄 Ejecutando análisis LEfSe (esto puede tomar unos minutos)...\n")
# res_marker <- run_lefse(
#   p_comp,
#   group = grp,
#   taxa_rank = rank_filter,
#   kw_cutoff = p_max,
#   lda_cutoff = lda_min
# )
# 
# # ✅ Correct: extract the plain data.frame from the S4 slot
# res <- as.data.frame(res_marker@marker_table)
# class(res) <- "data.frame"   # ✅ Key: completely strip the marker_table class
# res$rank <- rank_filter
# 
# write_csv_safe(res, file.path(out_dir, "lefse_result.csv"))
# 
# # 5. Preparar anotaciones taxonómicas
# tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("OTU_ID")
# 
# rank_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# rank_cols <- rank_cols[rank_cols %in% colnames(tax_df)]
# 
# tax_df$tax_path <- apply(tax_df[, rank_cols, drop = FALSE], 1, function(x) {
#   x[is.na(x)] <- "Unassigned"
#   paste(trimws(x), collapse = "|")
# })
# 
# # 6. Mapear columnas dinámicas, filtrar y hacer join por hoja taxonómica
# col_p   <- if ("pvalue"  %in% colnames(res)) "pvalue"  else "p.value"
# col_lda <- if ("ef_lda"  %in% colnames(res)) "ef_lda"  else "LDAscore"
# col_tax <- if ("feature" %in% colnames(res)) "feature" else "tax"
# col_dir <- if ("enrich_group" %in% colnames(res)) "enrich_group" else "direction"
# 
# lefse_filtered <- res %>%
#   dplyr::filter(.data[[col_p]] <= p_max, .data[[col_lda]] >= lda_min)
# 
# if (!is.null(drop_pattern) && nzchar(drop_pattern)) {
#   lefse_filtered <- lefse_filtered %>% 
#     #OLDER:dplyr::filter(!grepl(drop_pattern, .data[[col_tax]]))
#     dplyr::filter(is.na(.data[[col_tax]]) | !grepl(drop_pattern, .data[[col_tax]]))
#   }
# 
# # Limpiar prefijos en la salida de LEfSe para hacer join limpio con tax_df
# lefse_filtered$tax_leaf <- sub("^[a-z]__", "", lefse_filtered[[col_tax]])
# #OLDER tax_df$tax_leaf <- tax_df[[rank_filter]]
# tax_df$tax_leaf <- sub("^[a-z]__", "", tax_df[[rank_filter]])
# 
# lefse_otu <- dplyr::left_join(lefse_filtered, tax_df, by = "tax_leaf")
# cat(sprintf("NAs en %s tras join: %d / %s\n", rank_filter, sum(is.na(lefse_otu[[rank_filter]])), nrow(lefse_otu)))
# 
# write_csv_safe(lefse_otu, file.path(out_dir, "lefse_sig_with_otu.csv"))
# 
# # 7. Generar gráfico si hay resultados significativos
# if (nrow(lefse_otu) > 0) {
#   pos_direction <- cfg$groups$levels[2] 
#   stopifnot(pos_direction %in% unique(lefse_otu[[col_dir]]))
#   
#   lefse_otu <- lefse_otu %>%
#     mutate(
#       LDAscore   = .data[[col_lda]],
#       tax_name   = .data[[col_tax]],
#       direction  = .data[[col_dir]],
#       LDA_signed = ifelse(direction == pos_direction, LDAscore, -LDAscore)
#     ) %>%
#     arrange(LDA_signed) %>%
#     mutate(tax_name = factor(tax_name, levels = tax_name))
#   
#   p <- ggplot(lefse_otu, aes(x = LDA_signed, y = tax_name, fill = direction)) +
#     geom_bar(stat = "identity", alpha = 0.8) +
#     theme_minimal(base_size = 12) +
#     labs(
#       title = "LEfSe Biomarkers (Filtered)", 
#       subtitle = sprintf("Rank: %s | LDA > %s | p < %s", rank_filter, lda_min, p_max),
#       x = "Signed LDA Score", 
#       y = "Taxon"
#     ) +
#     scale_fill_manual(values = unlist(cfg$groups$colors)) +
#     theme(axis.text.y = element_text(size = 9))
#   
#   if (exists("save_plot")) {
#     save_plot(p, file.path(out_dir, "lefse_barplot.png"), w = 9, h = max(6, nrow(lefse_otu) * 0.25))
#   } else {
#     ggsave(file.path(out_dir, "lefse_barplot.png"), p, width = 9, height = max(6, nrow(lefse_otu) * 0.25))
#   }
#   
#   cat(sprintf("✅ Gráfico LEfSe guardado. Se encontraron %d biomarcadores significativos.\n", nrow(lefse_otu)))
# } else {
#   warning("⚠️ No se encontraron biomarcadores significativos con los filtros de LEfSe actuales.")
# }
# 
# cat("✅ Análisis LEfSe finalizado con éxito en outputs/lefse/.\n")
# 
# 
