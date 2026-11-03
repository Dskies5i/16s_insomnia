# MEMORY EFFICIENT METHOD
# ============================================================
# SCRIPT: 13a_extract_ko_contrib.R  ---  (Pair-Restricted & Lightweight)
# USE: Rscript 1_meta_analysis/13a_extract_ko_contrib.R
# ============================================================

source("1_meta_analysis/00_setup.R")
source("utils/utils_io.R")
library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)

cat("⏳ Creating lightweight functional contribution tables based on filtered edges...\n")
out_tables <- "outputs/tables"
fs::dir_create(out_tables)

# 1. Load Phyloseq and Edge Data
ps <- readRDS("outputs/models/physeq_filtered.rds")
edges_path <- "outputs/function/tables/taxon_ko_edges_filtered.csv"
if (!file.exists(edges_path)) stop("❌ Error: taxon_ko_edges_filtered.csv not found. Run script 10 first!")

edges <- read_csv(edges_path, show_col_types = FALSE)

# 2. Metadata
meta <- as.data.frame(as(sample_data(ps), "data.frame")) %>% tibble::rownames_to_column("Sample")
group_col <- intersect(c("Group", "group"), colnames(meta))[1]

meta_taxa <- meta %>% 
  dplyr::rename(Group = all_of(group_col)) %>%
  dplyr::select(Sample, Group) %>% 
  dplyr::mutate(Set_for_taxa = "Set1")

write_csv_safe(meta_taxa, file.path(out_tables, "metadata_for_taxa.csv"))

# 3. Taxonomy
tax_df <- as.data.frame(as(tax_table(ps), "matrix")) %>% tibble::rownames_to_column("ASV")
tax_filled <- tax_df %>%
  dplyr::mutate(
    Set_taxon = paste("Set1", Genus, sep = "_"),
    Taxa_name = ifelse(is.na(Genus) | Genus == "", "g__Unassigned", Genus)
  )

write_csv_safe(tax_filled, file.path(out_tables, "picrust_taxonomy_filled.csv"))

# 4. KO Fold Change
ko_full <- read_tsv("outputs/tables/merged_KO_metagenome.tsv", show_col_types = FALSE)
colnames(ko_full)[1] <- "KO"

ko_long <- ko_full %>%
  pivot_longer(cols = -KO, names_to = "Sample", values_to = "Abundance") %>%
  left_join(meta_taxa, by = "Sample") %>%
  group_by(KO, Group) %>%
  summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_ab) %>%
  mutate(
    fc = Insomnia - Control,
    Strength = case_when(fc > 0 ~ "Upregulated", fc < 0 ~ "Downregulated", TRUE ~ "Unchanged"),
    ko_function = KO
  )

write_csv_safe(ko_long, file.path(out_tables, "KO_FC_group.csv"))

# 5. Filtered Contribution (Restricted strictly to significant edges in `edges`)
otu_mat <- as.data.frame(otu_table(ps))
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

otu_long <- as.data.frame(otu_mat) %>%
  tibble::rownames_to_column("ASV") %>%
  pivot_longer(cols = -ASV, names_to = "sample", values_to = "count") %>%
  left_join(tax_filled %>% select(ASV, Genus), by = "ASV") %>%
  filter(!is.na(Genus) & Genus != "")

# Aggregate counts per sample and genus
taxon_sample_counts <- otu_long %>%
  group_by(sample, Genus) %>%
  summarise(norm_taxon_function_contrib = sum(count), .groups = "drop") %>%
  rename(taxon = Genus)

# Join ONLY with the valid Taxon-KO pairs present in edges
# edges typically has columns like Taxon and KO
edge_pairs <- edges %>% select(Taxon, KO) %>% distinct()

contrib_filtered <- edge_pairs %>%
  tidyr::crossing(sample = unique(taxon_sample_counts$sample)) %>%
  left_join(taxon_sample_counts, by = c("sample" = "sample", "Taxon" = "taxon")) %>%
  mutate(
    taxon = Taxon,
    norm_taxon_function_contrib = ifelse(is.na(norm_taxon_function_contrib), 0, norm_taxon_function_contrib),
    Set_for_taxa = "Set1"
  ) %>%
  select(sample, KO, taxon, norm_taxon_function_contrib, Set_for_taxa)

write_csv_safe(contrib_filtered, file.path(out_tables, "ko_contrib_filtered.csv"))
cat("✅ Successfully generated lightweight ko_contrib_filtered.csv!\n")





# 
# # ============================================================
# # SCRIPT: 13b_extract_ko_contrib.R (Memory-Optimized Version)
# # ============================================================
# ### BIG FILE!!!
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_io.R")
# library(dplyr)
# library(tidyr)
# library(readr)
# library(phyloseq)
# 
# cat("⏳ Creating functional contribution input tables (optimized)...\n")
# out_tables <- "outputs/tables"
# fs::dir_create(out_tables)
# 
# # 1. Load Phyloseq and Edge Data
# ps <- readRDS("outputs/models/physeq_filtered.rds")
# edges_path <- "outputs/function/tables/taxon_ko_edges_filtered.csv"
# edges <- read_csv(edges_path, show_col_types = FALSE)
# target_kos <- unique(edges$KO)
# 
# # 2. Metadata
# meta <- as.data.frame(as(sample_data(ps), "data.frame")) %>% tibble::rownames_to_column("Sample")
# group_col <- intersect(c("Group", "group"), colnames(meta))[1]
# 
# meta_taxa <- meta %>% 
#   dplyr::rename(Group = all_of(group_col)) %>%
#   dplyr::select(Sample, Group) %>% 
#   dplyr::mutate(Set_for_taxa = "Set1")
# 
# write_csv_safe(meta_taxa, file.path(out_tables, "metadata_for_taxa.csv"))
# 
# # 3. Taxonomy
# tax_df <- as.data.frame(as(tax_table(ps), "matrix")) %>% tibble::rownames_to_column("ASV")
# tax_filled <- tax_df %>%
#   dplyr::mutate(
#     Set_taxon = paste("Set1", Genus, sep = "_"),
#     Taxa_name = ifelse(is.na(Genus) | Genus == "", "g__Unassigned", Genus)
#   )
# 
# write_csv_safe(tax_filled, file.path(out_tables, "picrust_taxonomy_filled.csv"))
# 
# # 4. KO Fold Change
# ko_full <- read_tsv("outputs/tables/merged_KO_metagenome.tsv", show_col_types = FALSE)
# colnames(ko_full)[1] <- "KO"
# 
# ko_long <- ko_full %>%
#   pivot_longer(cols = -KO, names_to = "Sample", values_to = "Abundance") %>%
#   left_join(meta_taxa, by = "Sample") %>%
#   group_by(KO, Group) %>%
#   summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from = Group, values_from = mean_ab) %>%
#   mutate(
#     fc = Insomnia - Control,
#     Strength = case_when(fc > 0 ~ "Upregulated", fc < 0 ~ "Downregulated", TRUE ~ "Unchanged"),
#     ko_function = KO
#   )
# 
# write_csv_safe(ko_long, file.path(out_tables, "KO_FC_group.csv"))
# 
# # 5. Filtered Contribution (Memory Efficient)
# otu_mat <- as.data.frame(otu_table(ps))
# if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
# 
# otu_long <- as.data.frame(otu_mat) %>%
#   tibble::rownames_to_column("ASV") %>%
#   pivot_longer(cols = -ASV, names_to = "sample", values_to = "count") %>%
#   left_join(tax_filled %>% select(ASV, Genus), by = "ASV") %>%
#   filter(!is.na(Genus) & Genus != "")
# 
# # Aggregate directly by sample and Genus without heavy expand_grid explosion
# contrib_df <- otu_long %>%
#   group_by(sample, Genus) %>%
#   summarise(norm_taxon_function_contrib = sum(count), .groups = "drop") %>%
#   rename(taxon = Genus) %>%
#   mutate(Set_for_taxa = "Set1")
# 
# # Cross join only for target KOs present in edges
# contrib_df <- tidyr::crossing(
#   sample = unique(contrib_df$sample),
#   KO = target_kos,
#   taxon = unique(contrib_df$taxon)
# ) %>%
#   left_join(contrib_df, by = c("sample", "taxon")) %>%
#   mutate(
#     norm_taxon_function_contrib = ifelse(is.na(norm_taxon_function_contrib), 0, norm_taxon_function_contrib),
#     Set_for_taxa = "Set1"
#   )
# 
# write_csv_safe(contrib_df, file.path(out_tables, "ko_contrib_filtered.csv"))
# cat("✅ ¡Tablas generadas exitosamente sin saturar la memoria!\n")
# 
# 
# 
# # # ============================================================
# # # RAM HEAVY!!!!
# # # SCRIPT: 13b_extract_ko_contrib.R
# # # DESCRIPTION: Generates downstream contribution tables safely
# # # ============================================================
# # source("1_meta_analysis/00_setup.R")
# # source("utils/utils_io.R")
# # library(dplyr)
# # library(tidyr)
# # library(readr)
# # library(phyloseq)
# # 
# # cat("⏳ Creating functional contribution input tables...\n")
# # out_tables <- "outputs/tables"
# # fs::dir_create(out_tables)
# # 
# # ps <- readRDS("outputs/models/physeq_filtered.rds")
# # edges_path <- "outputs/function/tables/taxon_ko_edges_filtered.csv"
# # edges <- read_csv(edges_path, show_col_types = FALSE)
# # target_kos <- unique(edges$KO)
# # 
# # meta <- as.data.frame(as(sample_data(ps), "data.frame")) %>% tibble::rownames_to_column("Sample")
# # group_col <- intersect(c("Group", "group"), colnames(meta))[1]
# # 
# # meta_taxa <- meta %>% 
# #   dplyr::rename(Group = all_of(group_col)) %>%
# #   dplyr::select(Sample, Group) %>% 
# #   dplyr::mutate(Set_for_taxa = "Set1")
# # 
# # write_csv_safe(meta_taxa, file.path(out_tables, "metadata_for_taxa.csv"))
# # 
# # tax_df <- as.data.frame(as(tax_table(ps), "matrix")) %>% tibble::rownames_to_column("ASV")
# # tax_filled <- tax_df %>%
# #   dplyr::mutate(
# #     Set_taxon = paste("Set1", Genus, sep = "_"),
# #     Taxa_name = ifelse(is.na(Genus) | Genus == "", "g__Unassigned", Genus)
# #   )
# # 
# # write_csv_safe(tax_filled, file.path(out_tables, "picrust_taxonomy_filled.csv"))
# # 
# # ko_full <- read_tsv("outputs/tables/merged_KO_metagenome.tsv", show_col_types = FALSE)
# # colnames(ko_full)[1] <- "KO"
# # 
# # ko_long <- ko_full %>%
# #   pivot_longer(cols = -KO, names_to = "Sample", values_to = "Abundance") %>%
# #   left_join(meta_taxa, by = "Sample") %>%
# #   group_by(KO, Group) %>%
# #   summarise(mean_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
# #   pivot_wider(names_from = Group, values_from = mean_ab) %>%
# #   mutate(
# #     fc = Insomnia - Control,
# #     Strength = case_when(fc > 0 ~ "Upregulated", fc < 0 ~ "Downregulated", TRUE ~ "Unchanged"),
# #     ko_function = KO
# #   )
# # 
# # write_csv_safe(ko_long, file.path(out_tables, "KO_FC_group.csv"))
# # 
# # otu_mat <- as.data.frame(otu_table(ps))
# # if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
# # 
# # otu_long <- as.data.frame(otu_mat) %>%
# #   tibble::rownames_to_column("ASV") %>%
# #   pivot_longer(cols = -ASV, names_to = "sample", values_to = "count") %>%
# #   left_join(tax_filled %>% select(ASV, Genus), by = "ASV") %>%
# #   filter(!is.na(Genus) & Genus != "")
# # 
# # contrib_df <- expand_grid(
# #   sample = unique(otu_long$sample),
# #   KO = target_kos,
# #   taxon = unique(otu_long$Genus)
# # ) %>%
# #   left_join(
# #     otu_long %>% group_by(sample, Genus) %>% summarise(norm_taxon_function_contrib = sum(count), .groups = "drop"),
# #     by = c("sample" = "sample", "taxon" = "Genus")
# #   ) %>%
# #   mutate(
# #     norm_taxon_function_contrib = ifelse(is.na(norm_taxon_function_contrib), 0, norm_taxon_function_contrib),
# #     Set_for_taxa = "Set1"
# #   )
# # 
# # write_csv_safe(contrib_df, file.path(out_tables, "ko_contrib_filtered.csv"))
# # cat("✅ 完了しました！\n")