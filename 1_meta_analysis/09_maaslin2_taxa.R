# ============================================================
# NOMBRE: 09_maaslin2_taxa.R
# DESCRIPCIÓN: Análisis de asociaciones diferenciales con MaAsLin2
# USO Rscript 1_meta_analysis/09_maaslin2_taxa.R
#
# INSTALL: conda install -c bioconda bioconductor-maaslin2
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================
source("1_meta_analysis/00_setup.R")
source("utils/utils_phyloseq.R")
source("utils/utils_io.R")
library(Maaslin2)
library(phyloseq)
library(dplyr)

out_dir <- file.path(OUT_DIR, "maaslin2")
dir_create(out_dir)

ps <- readRDS(file.path(OUT_DIR, "models/physeq_filtered.rds"))

otu_mat <- as.data.frame(otu_table(ps), check.names = FALSE)
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
otu <- as.data.frame(t(otu_mat))

otu_prefix <- cfg$differential$boruta$otu_prefix
colnames(otu) <- paste0(otu_prefix, colnames(otu))

# Build a pristine, flat data frame containing ONLY required fixed effects
raw_meta <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
fixed_vars <- cfg$differential$maaslin2$fixed_effects

meta_clean <- data.frame(row.names = rownames(raw_meta))
for (v in fixed_vars) {
  if (v %in% colnames(raw_meta)) {
    col_data <- raw_meta[[v]]
    if (is.data.frame(col_data)) col_data <- col_data[[1]]
    meta_clean[[v]] <- unlist(col_data)
  }
}

# Explicitly format fixed effects variables
group_col <- cfg$columns$group
if (group_col %in% colnames(meta_clean)) meta_clean[[group_col]] <- as.factor(meta_clean[[group_col]])
if ("study_id" %in% colnames(meta_clean)) meta_clean$study_id <- as.factor(meta_clean$study_id)
if ("age" %in% colnames(meta_clean)) meta_clean$age <- as.numeric(as.character(meta_clean$age))

meta <- meta_clean
otu <- otu[rownames(meta), , drop = FALSE]



fit <- Maaslin2(
  input_data = otu,
  input_metadata = meta,
  output = file.path(out_dir, cfg$differential$maaslin2$output_dir),
  fixed_effects = cfg$differential$maaslin2$fixed_effects,
  plot_heatmap = FALSE,  # Disable internal plotting to prevent ggplot2 version crash
  plot_scatter = FALSE   # Disable internal scatter plots
)


# Unir con taxonomía
sig_path <- file.path(out_dir, cfg$differential$maaslin2$output_dir, "significant_results.tsv")
if (file.exists(sig_path)) {
  sig_res <- read.delim(sig_path, sep = "\t", check.names = FALSE)
  sig_res$OTU_ID <- gsub(otu_prefix, "", sig_res$feature)
  
  taxonomy <- as.data.frame(tax_table(ps)) %>% tibble::rownames_to_column("OTU_ID")
  sig_merged <- dplyr::left_join(sig_res, taxonomy, by = "OTU_ID") %>%
    dplyr::select(feature, metadata, value, coef, qval, Phylum, Family, Genus, Species)
  
  write_csv_safe(sig_merged, file.path(out_dir, "maaslin2_sig_with_taxonomy.csv"))
}
cat("✅ Análisis MaAsLin2 completado con éxito.\n")


#IMPROVEMENTS FROM PARK S.
#Pristine Metadata Extraction: Added safe extraction wrappers to prevent R's internal xtfrm crash caused by phyloseq's internal metadata column wrappers.

#Plotting Bypass (plot_heatmap = FALSE, plot_scatter = FALSE): Disabled Maaslin2's internal plotting engine to prevent crashes caused by deprecated label syntax in modern versions of ggplot2 (v3.4.0+).
