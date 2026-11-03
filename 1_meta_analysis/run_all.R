# ============================================================
# NOMBRE: run_all.R
# DESCRIPCIÓN: Ejecución secuencial de todo el pipeline de metaanálisis 16S
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# USAGE: Rscript 1_meta_analysis/run_all.R
# ============================================================

scripts <- c(
  "1_meta_analysis/01_build_phyloseq.R",
  "1_meta_analysis/02_alpha_diversity.R",
  "1_meta_analysis/02_QC_alpha_diversity_outputs.R",
  "1_meta_analysis/03_beta_diversity_permanova.R",
  "1_meta_analysis/03_QC_beta_diversity_outputs.R",
  "1_meta_analysis/04_fb_ratio.R",
  "1_meta_analysis/04_QC_fb_ratio_outputs.R",
  "1_meta_analysis/05_boruta.R",
  "1_meta_analysis/05_QC_boruta_outputs.R",
  "1_meta_analysis/06_boruta_pca_heatmap.R",
  "1_meta_analysis/06_QC_boruta_pca_heatmap_outputs.R",
  "1_meta_analysis/07_dbrda_envfit_fix.R",
  "1_meta_analysis/07_QC_dbrda_outputs.R",
  "1_meta_analysis/08_lefse.R",
  "1_meta_analysis/08_QC_lefse.R",
  "1_meta_analysis/09_maaslin2_taxa.R",
  "1_meta_analysis/10A_run_picrust2.sh",
  "1_meta_analysis/10B_function_taxa_ko_cor.R",
  "1_meta_analysis/10B_QC_function_taxa_ko.R",
  "1_meta_analysis/11_function_chord_edges.R",
  "1_meta_analysis/12_function_dotplot.R",
  "1_meta_analysis/13a_extract_ko_contrib.R",
  "1_meta_analysis/13a_QC_extract_ko_contrib.R",
  "1_meta_analysis/13b_function_contribution_diff.R",
  "1_meta_analysis/13b_QC_function_contribution_diff.R",
  "1_meta_analysis/14_export_ko_list_for_cluego.R",
  "1_meta_analysis/15_function_chord_from_ko_fc.R"
)

for (s in scripts) {
  message("\n==================================================")
  message("🚀 Running: ", s)
  message("==================================================\n")
  
  if (grepl("\\.sh$", s)) {
    # Execute shell scripts using system()
    status <- system(paste("bash", s))
    if (status != 0) {
      stop("❌ Error executing shell script: ", s)
    }
  } else {
    # Execute R scripts using source()
    source(s)
  }
}

message("\n✅ ¡Pipeline completo ejecutado con éxito de principio a fin!\n")