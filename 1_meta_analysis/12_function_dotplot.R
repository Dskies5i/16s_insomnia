# ============================================================
# NOMBRE: 12_function_dotplot.R
# DESCRIPCIÓN: Dotplot de correlaciones significativas Taxón-KO (Figura 5B)
# USAGE: Rscript 1_meta_analysis/12_function_dotplot.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

source("1_meta_analysis/00_setup.R")
library(readr)
library(dplyr)
library(ggplot2)
library(fs)

out_dir <- file.path(OUT_DIR, "function")
dir_create(file.path(out_dir, "figures"))

edge_path <- file.path(out_dir, "tables/taxon_ko_edges_filtered.csv")
if (!file.exists(edge_path)) {
  stop("❌ Error: filtered edges table not found. Run script 10 first.")
}

edge_df <- read_csv(edge_path, show_col_types = FALSE)
if (nrow(edge_df) == 0) {
  stop("❌ Error: No edges after filtering.")
}

# Check if node mapping file exists before trying to read it
node_file <- if (!is.null(cfg[["function"]]$inputs$ko_nodes)) cfg[["function"]]$inputs$ko_nodes else "outputs/function/tables/cytoscape_nodes.csv"

if (file.exists(node_file)) {
  nodes <- read_csv(node_file, show_col_types = FALSE)
  map <- if (!is.null(cfg[["function"]]$dotplot$ko_direction_map)) cfg[["function"]]$dotplot$ko_direction_map else c("Up" = "insomnia", "Down" = "control")
  
  nodes2 <- nodes %>%
    mutate(direction = dplyr::case_when(
      Strength %in% names(map) ~ unname(map[Strength]),
      TRUE ~ "control"
    )) %>%
    select(ID, direction)
  
  edge_df2 <- edge_df %>%
    mutate(ko_id = sub("^ko:", "", KO)) %>%
    left_join(nodes2, by = c("ko_id" = "ID")) %>%
    mutate(
      direction = ifelse(is.na(direction), "control", direction),
      abs_rho = abs(rho),
      direction = factor(direction, levels = c("control", "insomnia"))
    )
} else {
  message("ℹ️ Note: ko_nodes file not found. Inferring direction from rho sign for dotplot.")
  edge_df2 <- edge_df %>%
    mutate(
      abs_rho = abs(rho),
      direction = factor(ifelse(rho > 0, "insomnia", "control"), levels = c("control", "insomnia"))
    )
}

if (nrow(edge_df2) > 200) {
  message("ℹ️ Limiting dotplot to top 200 strongest associations for visual clarity.")
  edge_df2 <- edge_df2 %>% arrange(desc(abs_rho)) %>% head(200)
}

ko_levels <- edge_df2 %>%
  mutate(KO_clean = sub("^ko:", "", KO),
         KO_num = suppressWarnings(as.numeric(sub("^K", "", KO_clean)))) %>%
  arrange(KO_num) %>% distinct(KO) %>% pull(KO)

tax_levels <- sort(unique(edge_df2$Taxon))

plot_df <- edge_df2 %>%
  mutate(
    KO = factor(KO, levels = ko_levels),
    Taxon = factor(Taxon, levels = tax_levels)
  )

colors_map <- if (!is.null(cfg[["function"]]$dotplot$colors)) unlist(cfg[["function"]]$dotplot$colors) else c("control" = "#2b8cbe", "insomnia" = "#de2d26")

# CORRECCIÓN AQUÍ: Uso seguro de cfg[["function"]] en lugar de cfg$function
size_range <- if (!is.null(cfg[["function"]]$dotplot$size_range)) unlist(cfg[["function"]]$dotplot$size_range) else c(1, 4)

p <- ggplot(plot_df, aes(x = KO, y = Taxon)) +
  geom_point(aes(size = abs_rho, color = direction), alpha = 0.9) +
  scale_color_manual(values = colors_map) +
  scale_size_continuous(range = size_range) +
  labs(
    x = "KO",
    y = "Taxon",
    color = "KO direction",
    size = "|Spearman rho|",
    title = "Significant taxon–KO correlations (Figure 5B)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

plot_path <- file.path(out_dir, "figures/figure5b_taxon_ko_dotplot.png")
ggsave(plot_path, p, width = 11, height = 6)
cat(sprintf("✅ Dotplot (Figure 5B) successfully generated at: %s\n", plot_path))


