# ============================================================
# NOMBRE: 13b_function_contribution_diff.R
# DESCRIPCIÓN: Análisis de contribución diferencial de funciones metabólicas
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# USAGE: Rscript 1_meta_analysis/13b_function_contribution_diff.R
# ============================================================
# ============================================================
# SCRIPT: 13b_function_contribution_diff.R
# DESCRIPTION: Differential functional contribution module (Fixed join)
# ============================================================

source("1_meta_analysis/00_setup.R")
source("utils/utils_io.R")
library(dplyr)
library(tidyr)
library(readr)
library(fs)

cat("🚀 Running Script 13: Differential Functional Contribution Analysis...\n")

out_contrib_dir <- "outputs/function_contrib"
fs::dir_create(out_contrib_dir)

# 1. Load required tables
contrib_path <- "outputs/tables/ko_contrib_filtered.csv"
meta_path <- "outputs/tables/metadata_for_taxa.csv"
fc_path <- "outputs/tables/KO_FC_group.csv"

if (!file.exists(contrib_path) || !file.exists(meta_path) || !file.exists(fc_path)) {
  stop("❌ Error: Missing required input tables in outputs/tables/. Run script 13b first.")
}

contrib_df <- read_csv(contrib_path, show_col_types = FALSE)
meta_taxa <- read_csv(meta_path, show_col_types = FALSE)
ko_fc <- read_csv(fc_path, show_col_types = FALSE)

# 2. Perform differential contribution calculations
cat("📊 Processing contribution shifts between groups...\n")

# Ensure metadata has lowercase 'sample' for proper joining
if ("Sample" %in% colnames(meta_taxa) && !"sample" %in% colnames(meta_taxa)) {
  meta_taxa <- meta_taxa %>% rename(sample = Sample)
}

# Merge contribution data with metadata groups
contrib_merged <- contrib_df %>%
  left_join(meta_taxa, by = "sample")

# Calculate mean contribution per group for each KO and taxon
contrib_summary <- contrib_merged %>%
  group_by(KO, taxon, Group) %>%
  summarise(mean_contrib = mean(norm_taxon_function_contrib, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_contrib, values_fill = 0)

# Compute contribution change (Insomnia vs Control)
if (all(c("Insomnia", "Control") %in% colnames(contrib_summary))) {
  contrib_summary <- contrib_summary %>%
    mutate(contrib_diff = Insomnia - Control)
} else {
  cols <- setdiff(colnames(contrib_summary), c("KO", "taxon"))
  if (length(cols) >= 2) {
    contrib_summary <- contrib_summary %>%
      mutate(contrib_diff = .data[[cols[2]]] - .data[[cols[1]]])
  } else {
    contrib_summary <- contrib_summary %>% mutate(contrib_diff = 0)
  }
}

# 3. Save final results
output_file <- file.path(out_contrib_dir, "differential_functional_contributions.csv")
write_csv_safe(contrib_summary, output_file)

cat(sprintf("✅ Success! Differential functional contribution table saved to: %s\n", output_file))







# # ============================================================
# # NOMBRE: 13_function_contribution_diff.R
# # DESCRIPCIÓN: Análisis de contribución diferencial de funciones metabólicas
# # ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# # ============================================================
# 
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_io.R")
# source("utils/utils_plot.R")
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(readr)
# 
# out_dir <- cfg$function_contrib$output$dir
# dir_create(file.path(out_dir, "tables"))
# dir_create(file.path(out_dir, "figures"))
# 
# metadata <- readr::read_csv(cfg$function_contrib$inputs$metadata_for_taxa, show_col_types = FALSE)
# taxonomy <- readr::read_csv(cfg$function_contrib$inputs$taxonomy_filled, show_col_types = FALSE)
# ko_contrib <- readr::read_csv(cfg$function_contrib$inputs$ko_contrib_filtered, show_col_types = FALSE)
# 
# metadata_simple <- metadata %>% dplyr::select(Sample, Set_for_taxa)
# 
# ko_contrib2 <- ko_contrib %>%
#   left_join(metadata_simple, by = c("sample" = "Sample")) %>%
#   mutate(Set_taxon = paste(Set_for_taxa, taxon, sep = "_"))
# 
# taxonomy_long <- taxonomy %>%
#   select(Set_taxon, Genus, Species) %>%
#   pivot_longer(cols = c(Genus, Species), names_to = "Level", values_to = "Taxa_name") %>%
#   filter(!is.na(Taxa_name) & Taxa_name != "")
# 
# ko_contrib_annotated <- ko_contrib2 %>%
#   left_join(taxonomy_long, by = "Set_taxon")
# 
# target_taxa <- cfg$function_contrib$targets$taxa
# ko_contrib_filtered <- ko_contrib_annotated %>%
#   filter(Taxa_name %in% target_taxa)
# 
# metadata2 <- metadata %>% rename(sample = Sample)
# ko_merged <- ko_contrib_filtered %>%
#   left_join(metadata2, by = "sample")
# 
# ko_summary <- ko_merged %>%
#   group_by(KO, Taxa_name, Group) %>%
#   summarise(mean_contrib = mean(norm_taxon_function_contrib, na.rm = TRUE), .groups = "drop")
# 
# ko_diff <- ko_summary %>%
#   pivot_wider(names_from = Group, values_from = mean_contrib) %>%
#   mutate(diff = Insomnia - Ctrl) %>%
#   filter(!is.na(diff))
# 
# write_csv_safe(ko_diff, file.path(out_dir, "tables/ko_diff_by_taxa.csv"))
# 
# base_size <- cfg$function_contrib$plot$base_size
# 
# p1 <- ggplot(ko_diff, aes(x = reorder(KO, diff), y = diff, fill = Taxa_name)) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.7) +
#   coord_flip() +
#   labs(title = "KO contribution difference (Insomnia - Ctrl)",
#        x = "KO", y = "Mean contribution difference", fill = "Taxa") +
#   theme_minimal(base_size = base_size)
# 
# save_plot(p1, file.path(out_dir, "figures/ko_diff_overall.png"), w = 10, h = 7)
# 
# p2 <- ggplot(ko_diff, aes(x = reorder(KO, diff), y = diff, fill = KO)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   coord_flip() +
#   facet_wrap(~ Taxa_name, scales = "free_y") +
#   labs(title = "KO contribution difference by taxa",
#        x = "KO", y = "Insomnia - Ctrl") +
#   theme_minimal(base_size = base_size)
# 
# save_plot(p2, file.path(out_dir, "figures/ko_diff_by_taxa_facet.png"), w = 12, h = 8)
# 
# plot_df_long <- ko_diff %>%
#   pivot_longer(cols = c(Ctrl, Insomnia), names_to = "Group", values_to = "mean_contrib")
# 
# p3 <- ggplot(plot_df_long, aes(x = reorder(KO, diff), y = mean_contrib, fill = Group)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
#   facet_wrap(~ Taxa_name, scales = "free_y") +
#   coord_flip() +
#   labs(title = "Group-wise KO contribution",
#        x = "KO", y = "Mean contribution") +
#   theme_minimal(base_size = base_size)
# 
# save_plot(p3, file.path(out_dir, "figures/ko_contrib_groupwise.png"), w = 12, h = 8)
# cat("✅ Contribución funcional calculada y graficada con éxito.\n")