# ============================================================
# NOMBRE: 14_export_ko_list_for_cluego.R
# DESCRIPCIÓN: Exportación de lista de genes KO para análisis con CluEgo
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# USAGE: Rscript 1_meta_analysis/14_export_ko_list_for_cluego.R
# ============================================================

source("1_meta_analysis/00_setup.R")
library(readr)

out_dir <- cfg[["function_contrib"]]$output$dir
dir_create(file.path(out_dir, "network"))

edges <- readr::read_csv("outputs/function/tables/taxon_ko_edges_filtered.csv", show_col_types = FALSE)
genes <- unique(edges$KO) 

write.table(
  genes,
  file = file.path(out_dir, "network/cluego_gene_list_KO.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
cat("✅ Lista de genes KO exportada para CluEgo.\n")




# # ============================================================
# # NOMBRE: 14_export_ko_list_for_cluego.R
# # DESCRIPCIÓN: Exportación de lista de genes KO para análisis con CluEgo
# # ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# # ============================================================
# 
# source("1_meta_analysis/00_setup.R")
# library(readr)
# 
# out_dir <- cfg$function_contrib$output$dir
# dir_create(file.path(out_dir, "network"))
# 
# 
# # Change this:
# # edges <- readr::read_csv(cfg$function_contrib$inputs$cytoscape_edges, show_col_types = FALSE)
# # genes <- unique(edges$Source)
# 
# # To this:
# edges <- readr::read_csv("outputs/function/tables/taxon_ko_edges_filtered.csv", show_col_types = FALSE)
# genes <- unique(edges$KO) 
# 
# write.table(
#   genes,
#   file = file.path(out_dir, "network/cluego_gene_list_KO.txt"),
#   row.names = FALSE, col.names = FALSE, quote = FALSE
# )
# cat("✅ Lista de genes KO exportada para CluEgo.\n")