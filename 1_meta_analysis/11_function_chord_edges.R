# 
# ============================================================
# NAME: 11_function_chord_edges.R
# DESCRIPTION: Gráfico de acorde (Chord Diagram) optimizado para interacciones Taxón-KO
# USAGE: Rscript 1_meta_analysis/11_function_chord_edges.R
# ============================================================

source("1_meta_analysis/00_setup.R")
source("utils/utils_io.R")
library(circlize)
library(readr)
library(dplyr)

out_dir <- file.path(OUT_DIR, "function")
dir_create(file.path(out_dir, "figures"))

edge_path <- file.path(out_dir, "tables/taxon_ko_edges_filtered.csv")
if (!file.exists(edge_path)) {
  stop("❌ Error: filtered edges table not found. Run script 10 first.")
}

edge_df <- read_csv(edge_path, show_col_types = FALSE)
if (nrow(edge_df) == 0) {
  stop("❌ Error: No edges found after filtering. Adjust thresholds in config.yaml.")
}

# Memory & Readability Protection: Keep top 40 strongest/most significant edges for the chord diagram
if (nrow(edge_df) > 40) {
  message(sprintf("ℹ️ Limiting chord diagram to top 40 strongest edges out of %d for clarity and memory safety.", nrow(edge_df)))
  edge_df <- edge_df %>%
    arrange(p, desc(abs(rho))) %>%
    head(40)
}

# Safe config extractors with fallbacks
cfg_chord <- if (!is.null(cfg[["function"]]$chord)) cfg[["function"]]$chord else list()
edge_pos     <- if (!is.null(cfg_chord$edge_pos)) cfg_chord$edge_pos else "#2b8cbe"
edge_neg     <- if (!is.null(cfg_chord$edge_neg)) cfg_chord$edge_neg else "#de2d26"
ko_color     <- if (!is.null(cfg_chord$ko_color)) cfg_chord$ko_color else "#bdbdbd"
taxa_palette <- if (!is.null(cfg_chord$taxa_palette)) cfg_chord$taxa_palette else c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#e6ab02")
gap_small    <- if (!is.null(cfg_chord$gap_small)) cfg_chord$gap_small else 2
gap_big      <- if (!is.null(cfg_chord$gap_big)) cfg_chord$gap_big else 8
start_degree <- if (!is.null(cfg_chord$start_degree)) cfg_chord$start_degree else 90
transparency <- if (!is.null(cfg_chord$transparency)) cfg_chord$transparency else 0.4
label_cex    <- if (!is.null(cfg_chord$label_cex)) cfg_chord$label_cex else 0.7

chord_df <- edge_df %>%
  mutate(
    value = abs(rho),
    edge_col = ifelse(rho > 0, edge_pos, edge_neg)
  ) %>%
  select(group_from = Taxon, group_to = KO, value, edge_col)

sectors <- unique(c(chord_df$group_from, chord_df$group_to))

io <- read_tax_otu_meta()
tax_mat <- as.data.frame(io$tax)

# Robust taxonomy sector matching across all available columns in tax table
taxon_levels <- unique(unlist(lapply(tax_mat, function(col) {
  intersect(unique(as.character(col)), sectors)
})))

grid_cols <- setNames(rep(ko_color, length(sectors)), sectors)
if (length(taxon_levels) > length(taxa_palette)) {
  taxa_palette <- colorRampPalette(taxa_palette)(length(taxon_levels))
}
if (length(taxon_levels) > 0) {
  grid_cols[taxon_levels] <- taxa_palette[seq_along(taxon_levels)]
}

sector_order <- sectors
is_taxa <- sector_order %in% taxon_levels
change_pos <- which(diff(as.integer(is_taxa)) != 0)
if (length(sector_order) > 1 && is_taxa[length(is_taxa)] != is_taxa[1]) {
  change_pos <- c(change_pos, length(sector_order))
}

gap.after <- rep(gap_small, length(sector_order))
if (length(change_pos) > 0) {
  gap.after[change_pos] <- gap_big
}

# ==============================================================================
# 1. GENERAR PNG DE ALTA RESOLUCIÓN (3600x3600, Lienzo Cuadrado)
# ==============================================================================
png(file.path(out_dir, "figures/chord_taxon_ko.png"), width = 3600, height = 3600, res = 300)

circos.clear()
circos.par(start.degree = start_degree, gap.after = gap.after)

chordDiagram(
  x = chord_df[, c("group_from", "group_to", "value")],
  grid.col = grid_cols[sector_order],
  col = chord_df$edge_col,
  transparency = transparency,
  annotationTrack = "grid",
  order = sector_order,
  preAllocateTracks = list(track.height = 0.08)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    nm <- get.cell.meta.data("sector.index")
    xc <- get.cell.meta.data("xcenter")
    # Truncar nombres extremadamente largos para evitar solapamiento
    nm_display <- if(nchar(nm) > 25) paste0(substr(nm, 1, 22), "...") else nm
    
    circos.text(x = xc, y = 0.5, labels = nm_display,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0.5, 0.5),
                cex = min(label_cex, 0.55), font = 2)
  },
  bg.border = NA
)

circos.clear()
dev.off()

# ==============================================================================
# 2. GENERAR PDF VECTORIAL CUADRADO (12x12 pulgadas)
# ==============================================================================
pdf(file.path(out_dir, "figures/chord_taxon_ko.pdf"), width = 12, height = 12)

circos.clear()
circos.par(start.degree = start_degree, gap.after = gap.after)

chordDiagram(
  x = chord_df[, c("group_from", "group_to", "value")],
  grid.col = grid_cols[sector_order],
  col = chord_df$edge_col,
  transparency = transparency,
  annotationTrack = "grid",
  order = sector_order,
  preAllocateTracks = list(track.height = 0.08)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    nm <- get.cell.meta.data("sector.index")
    xc <- get.cell.meta.data("xcenter")
    nm_display <- if(nchar(nm) > 25) paste0(substr(nm, 1, 22), "...") else nm
    
    circos.text(x = xc, y = 0.5, labels = nm_display,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0.5, 0.5),
                cex = min(label_cex, 0.55), font = 2)
  },
  bg.border = NA
)

circos.clear()
dev.off()

cat("✅ Gráficos de acorde mejorados generados con éxito en PNG y PDF (outputs/function/figures/)\n")






# 
# # # ============================================================
# # # NAME: 11_function_chord_edges.R
# # # DESCRIPTION: Gráfico de acorde (Chord Diagram) para interacciones Taxón-KO
# # # USAGE: Rscript 1_meta_analysis/11_function_chord_edges.R
# # # ============================================================
# 
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_io.R")
# library(circlize)
# library(readr)
# library(dplyr)
# 
# out_dir <- file.path(OUT_DIR, "function")
# dir_create(file.path(out_dir, "figures"))
# 
# edge_path <- file.path(out_dir, "tables/taxon_ko_edges_filtered.csv")
# if (!file.exists(edge_path)) {
#   stop("❌ Error: filtered edges table not found. Run script 10 first.")
# }
# 
# edge_df <- read_csv(edge_path, show_col_types = FALSE)
# if (nrow(edge_df) == 0) {
#   stop("❌ Error: No edges found after filtering. Adjust thresholds in config.yaml.")
# }
# 
# # Memory & Readability Protection: Keep top 40 strongest/most significant edges for the chord diagram
# if (nrow(edge_df) > 40) {
#   message(sprintf("ℹ️ Limiting chord diagram to top 40 strongest edges out of %d for clarity and memory safety.", nrow(edge_df)))
#   edge_df <- edge_df %>%
#     arrange(p, desc(abs(rho))) %>%
#     head(40)
# }
# 
# # Safe config extractors with fallbacks
# cfg_chord <- if (!is.null(cfg[["function"]]$chord)) cfg[["function"]]$chord else list()
# edge_pos     <- if (!is.null(cfg_chord$edge_pos)) cfg_chord$edge_pos else "#2b8cbe"
# edge_neg     <- if (!is.null(cfg_chord$edge_neg)) cfg_chord$edge_neg else "#de2d26"
# ko_color     <- if (!is.null(cfg_chord$ko_color)) cfg_chord$ko_color else "#bdbdbd"
# taxa_palette <- if (!is.null(cfg_chord$taxa_palette)) cfg_chord$taxa_palette else c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#e6ab02")
# gap_small    <- if (!is.null(cfg_chord$gap_small)) cfg_chord$gap_small else 2
# gap_big      <- if (!is.null(cfg_chord$gap_big)) cfg_chord$gap_big else 8
# start_degree <- if (!is.null(cfg_chord$start_degree)) cfg_chord$start_degree else 90
# transparency <- if (!is.null(cfg_chord$transparency)) cfg_chord$transparency else 0.4
# label_cex    <- if (!is.null(cfg_chord$label_cex)) cfg_chord$label_cex else 0.7
# 
# chord_df <- edge_df %>%
#   mutate(
#     value = abs(rho),
#     edge_col = ifelse(rho > 0, edge_pos, edge_neg)
#   ) %>%
#   select(group_from = Taxon, group_to = KO, value, edge_col)
# 
# sectors <- unique(c(chord_df$group_from, chord_df$group_to))
# 
# io <- read_tax_otu_meta()
# tax_mat <- as.data.frame(io$tax)
# 
# # Robust taxonomy sector matching across all available columns in tax table
# taxon_levels <- unique(unlist(lapply(tax_mat, function(col) {
#   intersect(unique(as.character(col)), sectors)
# })))
# 
# grid_cols <- setNames(rep(ko_color, length(sectors)), sectors)
# if (length(taxon_levels) > length(taxa_palette)) {
#   taxa_palette <- colorRampPalette(taxa_palette)(length(taxon_levels))
# }
# if (length(taxon_levels) > 0) {
#   grid_cols[taxon_levels] <- taxa_palette[seq_along(taxon_levels)]
# }
# 
# sector_order <- sectors
# is_taxa <- sector_order %in% taxon_levels
# change_pos <- which(diff(as.integer(is_taxa)) != 0)
# if (length(sector_order) > 1 && is_taxa[length(is_taxa)] != is_taxa[1]) {
#   change_pos <- c(change_pos, length(is_taxa))
# }
# 
# gap.after <- rep(gap_small, length(sector_order))
# if (length(change_pos) > 0) {
#   gap.after[change_pos] <- gap_big
# }
# 
# pdf(file.path(out_dir, "figures/chord_taxon_ko.pdf"), width = 10, height = 8)
# circos.clear()
# circos.par(start.degree = start_degree, gap.after = gap.after)
# 
# chordDiagram(
#   x = chord_df[, c("group_from", "group_to", "value")],
#   grid.col = grid_cols[sector_order],
#   col = chord_df$edge_col,
#   transparency = transparency,
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
#                 cex = label_cex, font = 2)
#   },
#   bg.border = NA
# )
# 
# circos.clear()
# dev.off()
# cat("✅ Gráfico de acorde generado con éxito en outputs/function/figures/chord_taxon_ko.pdf\n")
