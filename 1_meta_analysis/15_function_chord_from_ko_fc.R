# ============================================================
# NOMBRE: 15_function_chord_from_ko_fc.R
# DESCRIPCIÓN: Gráfico de acorde basado en cambios de pliegue (FC) de KOs
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# USAGE: Rscript 1_meta_analysis/15_function_chord_from_ko_fc.R
# ============================================================

source("1_meta_analysis/00_setup.R")
source("utils/utils_io.R")
source("utils/utils_helper_functions.R")

library(circlize)
library(readr)
library(dplyr)
library(fs)

cat("Running Script 15: Functional Chord Diagram from KO Fold Change...\n")

# 1. Definir directorios de salida
out_dir <- file.path(OUT_DIR, "function")
dir_create(file.path(out_dir, "figures"))

# 2. Cargar datos de contribución diferencial o fold change de KOs
fc_path <- "outputs/tables/KO_FC_group.csv"
if (!file.exists(fc_path)) {
  stop("❌ Error: KO_FC_group.csv no encontrado. Ejecuta primero el script 13a.")
}

fc_df <- read_csv(fc_path, show_col_types = FALSE)

# 3. Preparar datos para el diagrama de acorde
# Filtrar o seleccionar las características más relevantes si es necesario
if (nrow(fc_df) > 40) {
  message(sprintf("⚠️ Limitando el diagrama de acorde a los primeros 40 KOs (de %d) para claridad visual.", nrow(fc_df)))
  fc_df <- fc_df %>% arrange(desc(abs(fc))) %>% head(40)
}

# Construir tabla de enlaces/sectores para circlize
chord_df <- fc_df %>%
  mutate(
    value = abs(fc),
    direction = ifelse(fc > 0, "Insomnia", "Control"),
    edge_col = ifelse(fc > 0, "#d95f02", "#2b5c8f")
  ) %>%
  select(from = KO, to = direction, value, edge_col)

# 4. Configurar paletas y colores del circlize
sectors <- unique(c(chord_df$from, chord_df$to))
grid_cols <- setNames(rep("#bdbdbd", length(sectors)), sectors)
grid_cols["Insomnia"] <- "#d95f02"
grid_cols["Control"] <- "#2b5c8f"

# 5. Generar PDF y PNG del diagrama de acorde
pdf_path <- file.path(out_dir, "figures/chord_ko_fc.pdf")
png_path <- file.path(out_dir, "figures/chord_ko_fc.png")

# PNG de alta resolución
png(png_path, width = 3600, height = 3600, res = 300)
circos.clear()
circos.par(start.degree = 90, gap.after = 5)
chordDiagram(
  x = chord_df[, c("from", "to", "value")],
  grid.col = grid_cols,
  col = chord_df$edge_col,
  transparency = 0.4,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1)
)
circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    nm <- get.cell.meta.data("sector.index")
    xc <- get.cell.meta.data("xcenter")
    circos.text(xc, y = 0.5, labels = nm, facing = "clockwise", niceFacing = TRUE, cex = 0.6, font = 2)
  }, 
  bg.border = NA
)
circos.clear()
dev.off()

cat(sprintf("✅ Gráfico de acorde generado con éxito en: %s\n", png_path))






# # ============================================================
# # NOMBRE: 15_function_chord_from_ko_fc.R
# # DESCRIPCIÓN: Diagrama de acorde basado en cambio de pliegue funcional (FC)
# # ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# # USAGE: Rscript 1_meta_analysis/15_function_chord_from_ko_fc.R
# # ============================================================
# 
# source("1_meta_analysis/00_setup.R")
# source("analysis/_common_function.R")
# source("utils/utils_io.R")
# library(circlize)
# library(readr)
# library(dplyr)
# 
# out_dir <- file.path(OUT_DIR, "function_fc_chord")
# dir_create(out_dir)
# 
# ko_list <- readr::read_csv(cfg$function_contrib$inputs$ko_fc_group, show_col_types = FALSE)
# 
# sig_ko <- ko_list %>%
#   dplyr::filter(Strength %in% c("Upregulated", "Downregulated")) %>%
#   dplyr::pull(ko_function) %>%
#   unique()
# 
# if (length(sig_ko) == 0) stop("No KO with Strength Upregulated/Downregulated.")
# 
# ko_mat_full <- subset_ko_matrix(cfg$function$inputs$picrust_ko, sig_ko)
# 
# io <- read_tax_otu_meta()
# tax_mat <- io$tax
# otu_mat <- io$otu
# tax_mat2 <- make_target_taxa_matrix(tax_mat, otu_mat)
# 
# common_samples <- intersect(colnames(ko_mat_full), colnames(tax_mat2))
# if (length(common_samples) == 0) stop("No common samples.")
# common_samples <- sort(common_samples)
# 
# ko_log <- log_transform(ko_mat_full[, common_samples, drop = FALSE],
#                         base = cfg$function$transform$log_base,
#                         pseudocount = cfg$function$transform$pseudocount)
# tax_log <- log_transform(tax_mat2[, common_samples, drop = FALSE],
#                          base = cfg$function$transform$log_base,
#                          pseudocount = cfg$function$transform$pseudocount)
# 
# cor_mat <- cor(t(tax_log), t(ko_log), method = "spearman", use = "pairwise.complete.obs")
# 
# tax_list <- rownames(tax_log)
# ko_list_names <- rownames(ko_log)
# p_mat <- matrix(NA_real_, nrow = length(tax_list), ncol = length(ko_list_names),
#                 dimnames = list(tax_list, ko_list_names))
# 
# for (i in tax_list) {
#   x <- as.numeric(tax_log[i, ])
#   for (j in ko_list_names) {
#     y <- as.numeric(ko_log[j, ])
#     ok <- complete.cases(x, y)
#     if (sum(ok) >= 3) p_mat[i, j] <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")$p.value)
#   }
# }
# 
# q_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
#                 nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
# 
# cor_long <- as.data.frame(as.table(cor_mat)) %>% rename(Taxon = Var1, KO = Var2, rho = Freq)
# p_long   <- as.data.frame(as.table(p_mat))   %>% rename(Taxon = Var1, KO = Var2, p = Freq)
# q_long   <- as.data.frame(as.table(q_mat))   %>% rename(Taxon = Var1, KO = Var2, q = Freq)
# 
# edge_df <- cor_long %>%
#   left_join(p_long, by = c("Taxon","KO")) %>%
#   left_join(q_long, by = c("Taxon","KO")) %>%
#   filter(!is.na(rho), p < 0.05, abs(rho) >= 0.5)
# 
# write_csv_safe(edge_df, file.path(out_dir, "edges_filtered.csv"))
# 
# if (nrow(edge_df) == 0) stop("No edges after filtering.")
# 
# chord_df <- edge_df %>%
#   mutate(value = abs(rho),
#          edge_col = ifelse(rho > 0, cfg$function$chord$edge_pos, cfg$function$chord$edge_neg)) %>%
#   select(group_from = Taxon, group_to = KO, value, edge_col)
# 
# sectors <- unique(c(chord_df$group_from, chord_df$group_to))
# taxon_levels <- intersect(rownames(tax_mat2), sectors)
# 
# grid_cols <- setNames(rep(cfg$function$chord$ko_color, length(sectors)), sectors)
# tax_colors <- cfg$function$chord$taxa_palette
# if (length(taxon_levels) > length(tax_colors)) stop("Increase taxa_palette.")
# grid_cols[taxon_levels] <- tax_colors[seq_along(taxon_levels)]
# 
# sector_order <- sectors
# is_taxa <- sector_order %in% taxon_levels
# change_pos <- which(diff(is_taxa) != 0)
# if (length(sector_order) > 1 && is_taxa[length(is_taxa)] != is_taxa[1]) change_pos <- c(change_pos, length(is_taxa))
# 
# gap.after <- rep(cfg$function$chord$gap_small, length(sector_order))
# gap.after[change_pos] <- cfg$function$chord$gap_big
# 
# pdf(file.path(out_dir, "chord_fc_based.pdf"), width = 10, height = 8)
# circos.clear()
# circos.par(start.degree = cfg$function$chord$start_degree, gap.after = gap.after)
# 
# chordDiagram(
#   x = chord_df[, c("group_from","group_to","value")],
#   grid.col = grid_cols[sector_order],
#   col = chord_df$edge_col,
#   transparency = cfg$function$chord$transparency,
#   annotationTrack = "grid",
#   order = sector_order,
#   preAllocateTracks = list(track.height = 0.08)
# )
# 
# circos.trackPlotRegion(
#   track.index = 1,
#   panel.fun = function(x, y) {
#     nm <- get.cell.meta.data("sector.index")
#     xc <- get.cell.meta.data("xcenter")
#     circos.text(x = xc, y = 0.5, labels = nm,
#                 facing = "clockwise", niceFacing = TRUE,
#                 adj = c(0.5, 0.5),
#                 cex = cfg$function$chord$label_cex, font = 2)
#   },
#   bg.border = NA
# )
# 
# circos.clear()
# dev.off()
# cat("✅ Diagrama de acorde por cambio de pliegue generado con éxito.\n")



