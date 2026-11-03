library(phyloseq)
library(dplyr)
library(tibble)

build_phyloseq <- function() {
  # Read paths and columns directly from the global config (cfg) using the correct YAML keys
  otu_path <- cfg$input$otu
  tax_path <- cfg$input$taxonomy
  meta_path <- cfg$input$metadata
  cols <- cfg$columns
  
  otu_df <- read.delim(otu_path, header = TRUE, check.names = FALSE, row.names = NULL)
  tax_df <- read.delim(tax_path, header = TRUE, check.names = FALSE, row.names = NULL)
  meta_df <- read.delim(meta_path, header = TRUE, check.names = FALSE, row.names = NULL)
  
  # Asignar rownames usando la primera columna de la tabla OTU
  otu_id_col <- colnames(otu_df)[1]
  otu_mat <- otu_df %>% tibble::column_to_rownames(otu_id_col) %>% as.matrix()
  
  # Asignar rownames usando la primera columna de la taxonomía
  tax_id_col <- colnames(tax_df)[1]
  tax_df <- tax_df %>% tibble::column_to_rownames(tax_id_col)
  
  # Procesar el string taxonómico de QIIME 2
  tax_col_name <- grep("Taxon|tax|taxonomy", colnames(tax_df), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(tax_col_name)) {
    tax_split <- strsplit(as.character(tax_df[[tax_col_name]]), ";\\s*")
    max_len <- max(sapply(tax_split, length))
    
    tax_mat <- do.call(rbind, lapply(tax_split, function(x) {
      length(x) <- max_len
      x[is.na(x)] <- "Unassigned"
      return(x)
    }))
    
    colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(tax_mat)]
    rownames(tax_mat) <- rownames(tax_df)
    
    tax_mat <- gsub("^[a-z]__", "", tax_mat)
    tax_mat[tax_mat == "" | tax_mat == "unidentified" | tax_mat == "Ambiguous"] <- NA
    tax_mat <- as.matrix(tax_mat)
  } else {
    tax_mat <- as.matrix(tax_df)
  }
  
  sample_id_col <- cols$sample_id
  meta_df <- meta_df %>% tibble::column_to_rownames(sample_id_col)
  
  if (rownames(meta_df)[1] == "#q2:types" || grepl("categorical|numeric", meta_df[1, 1])) {
    meta_df <- meta_df[-1, ]
  }
  
  common_samples <- intersect(colnames(otu_mat), rownames(meta_df))
  common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
  
  otu_mat <- otu_mat[common_taxa, common_samples]
  meta_df <- meta_df[common_samples, ]
  tax_mat <- tax_mat[common_taxa, ]
  
  ps <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat),
    sample_data(meta_df)
  )
  
  
  
  return(ps)
}

# ------------------------------------------------------------
# Fija los niveles (orden) del factor de grupo en sample_data.
# Necesario para que las comparaciones Control-vs-Insomnia y las
# paletas de color se mantengan consistentes en todo el pipeline.
# ------------------------------------------------------------
enforce_group_levels <- function(ps, group_col, levels) {
  sd <- data.frame(sample_data(ps))
  if (!group_col %in% colnames(sd)) {
    stop("Columna de grupo '", group_col, "' no existe en sample_data.")
  }
  sd[[group_col]] <- factor(sd[[group_col]], levels = levels)
  sample_data(ps) <- sample_data(sd)
  ps
}
################################


# 
# 
# build_phyloseq <- function() {
#   # Read paths and columns directly from the global config (cfg) using the correct YAML keys
#   otu_path <- cfg$input$otu
#   tax_path <- cfg$input$taxonomy
#   meta_path <- cfg$input$metadata
#   cols <- cfg$columns
#   
#   otu_df <- read.delim(otu_path, header = TRUE, check.names = FALSE, row.names = NULL)
#   tax_df <- read.delim(tax_path, header = TRUE, check.names = FALSE, row.names = NULL)
#   meta_df <- read.delim(meta_path, header = TRUE, check.names = FALSE, row.names = NULL)
#   
#   # Asignar rownames usando la primera columna de la tabla OTU
#   otu_id_col <- colnames(otu_df)[1]
#   otu_mat <- otu_df %>% tibble::column_to_rownames(otu_id_col) %>% as.matrix()
#   
#   # Asignar rownames usando la primera columna de la taxonomía
#   tax_id_col <- colnames(tax_df)[1]
#   tax_df <- tax_df %>% tibble::column_to_rownames(tax_id_col)
#   
#   # Procesar el string taxonómico de QIIME 2
#   tax_col_name <- grep("Taxon|tax|taxonomy", colnames(tax_df), ignore.case = TRUE, value = TRUE)[1]
#   if (!is.na(tax_col_name)) {
#     tax_split <- strsplit(as.character(tax_df[[tax_col_name]]), ";\\s*")
#     max_len <- max(sapply(tax_split, length))
#     
#     tax_mat <- do.call(rbind, lapply(tax_split, function(x) {
#       length(x) <- max_len
#       x[is.na(x)] <- "Unassigned"
#       return(x)
#     }))
#     
#     colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(tax_mat)]
#     rownames(tax_mat) <- rownames(tax_df)
#     
#     tax_mat <- gsub("^[a-z]__", "", tax_mat)
#     tax_mat[tax_mat == "" | tax_mat == "unidentified" | tax_mat == "Ambiguous"] <- NA
#     tax_mat <- as.matrix(tax_mat)
#   } else {
#     tax_mat <- as.matrix(tax_df)
#   }
#   
#   sample_id_col <- cols$sample_id
#   meta_df <- meta_df %>% tibble::column_to_rownames(sample_id_col)
#   
#   if (rownames(meta_df)[1] == "#q2:types" || grepl("categorical|numeric", meta_df[1, 1])) {
#     meta_df <- meta_df[-1, ]
#   }
#   
#   common_samples <- intersect(colnames(otu_mat), rownames(meta_df))
#   common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
#   
#   otu_mat <- otu_mat[common_taxa, common_samples]
#   meta_df <- meta_df[common_samples, ]
#   tax_mat <- tax_mat[common_taxa, ]
#   
#   ps <- phyloseq(
#     otu_table(otu_mat, taxa_are_rows = TRUE),
#     tax_table(tax_mat),
#     sample_data(meta_df)
#   )
#   
#   return(ps)
# }




