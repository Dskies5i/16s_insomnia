### Github base
#- https://doi.org/10.1016/j.crmicr.2026.100577
#- https://github.com/Sun980610/26_Insomnia_meta-analysis

### Crear carpetas principales
dir.create("setup")
dir.create("analysis")
dir.create("R")
dir.create("outputs")
dir.create("data")
dir.create("config")
dir.create("preprocessing")

### Crear subcarpetas dentro de outputs
dir.create("outputs/tables")
dir.create("outputs/figures")
dir.create("outputs/models")
dir.create("outputs/logs")

### Crear subcarpetas dentro de data
dir.create("data/raw")
dir.create("data/processed")

### -----------------------------

# Crear archivos vacíos en analysis/
file.create("analysis/00_setup.R")
file.create("analysis/01_build_phyloseq.R")
file.create("analysis/02_alpha_diversity.R")
file.create("analysis/03_beta_diversity_permanova.R")
file.create("analysis/04_fb_ratio.R")
file.create("analysis/05_boruta.R")
file.create("analysis/06_boruta_pca_heatmap.R")
file.create("analysis/07_dbrda_envfit.R")
file.create("analysis/08_lefse.R")
file.create("analysis/09_maaslin2_taxa.R")
file.create("analysis/10_function_taxa_ko_cor.R")
file.create("analysis/11_function_chord_edges.R")
file.create("analysis/12_function_dotplot_fig5b.R")
file.create("analysis/13_function_contribution_diff.R")
file.create("analysis/14_export_ko_list_for_cluego.R")
file.create("analysis/15_function_chord_from_ko_fc.R")
file.create("analysis/run_all.R")
file.create("config/config.yaml")

# Crear archivos en R/
file.create("R/utils_io.R")
file.create("R/utils_phyloseq.R")
file.create("R/utils_plot.R")
file.create("R/utils_stats.R")


# Crear archivos en Setup/
