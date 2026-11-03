# ============================================================
# NAME: 10B_QC_function_taxa_ko.R
# DESCRIPTION: Quality Control and Summary Diagnostics for 
#              Taxon-KO Functional Correlations (Script 10)
# USAGE: Rscript 1_meta_analysis/10B_QC_function_taxa_ko.R
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(fs)
})

cat("\n════════════════════════════════════════════════════════════\n")
cat("  Starting QC & Diagnostics: Taxon-KO Correlations\n")
cat("════════════════════════════════════════════════════════════\n\n")

# 1. Define paths
out_dir <- "outputs/function"
fig_dir <- file.path(out_dir, "figures")
dir_create(fig_dir)

all_cor_path <- file.path(out_dir, "tables/taxon_ko_cor_all.csv")
filtered_edges_path <- file.path(out_dir, "tables/taxon_ko_edges_filtered.csv")

# 2. Check if output files exist
if (!file.exists(all_cor_path) || !file.exists(filtered_edges_path)) {
  stop("❌ Error: Output tables not found! Make sure script 10 finished successfully first.")
}

# 3. Load data
cat("📂 Loading correlation tables...\n")
cor_full <- read_csv(all_cor_path, show_col_types = FALSE)
edge_df  <- read_csv(filtered_edges_path, show_col_types = FALSE)

# 4. Print Summary Diagnostics
cat("\n------------------------------------------------------------\n")
cat("📊 SUMMARY STATISTICS\n")
cat("------------------------------------------------------------\n")
cat(sprintf("• Total Taxon-KO pairs evaluated: %d\n", nrow(cor_full)))
cat(sprintf("• Unique Taxa analyzed          : %d\n", n_distinct(cor_full$Taxon)))
cat(sprintf("• Unique KO functions analyzed  : %d\n", n_distinct(cor_full$KO)))
cat(sprintf("• Filtered significant edges    : %d (%.2f%% of total)\n", 
            nrow(edge_df), (nrow(edge_df) / nrow(cor_full)) * 100))

if (nrow(edge_df) > 0) {
  cat(sprintf("• Correlation range (rho) in edges: [%.3f, %.3f]\n", 
              min(edge_df$rho, na.rm = TRUE), max(edge_df$rho, na.rm = TRUE)))
  cat(sprintf("• Mean absolute rho in edges    : %.3f\n", 
              mean(abs(edge_df$rho), na.rm = TRUE)))
  
  cat("\n🔍 Top 5 Strongest Significant Associations:\n")
  top_edges <- edge_df %>% arrange(p) %>% head(5)
  print(top_edges %>% select(Taxon, KO, rho, p, q))
} else {
  cat("⚠️ Warning: 0 edges passed the current filtering thresholds in config.yaml.\n")
}

# 5. Generate Quick QC Visualization
cat("\n🎨 Generating QC distribution plot...\n")
qc_plot <- ggplot(cor_full, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "QC: Distribution of Taxon-KO Correlation Coefficients",
    subtitle = "Dashed red lines show typical filter cutoffs (abs(rho) >= 0.5)",
    x = "Spearman Correlation (rho)",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 12)

plot_path <- file.path(fig_dir, "qc_cor_distribution.pdf")
ggsave(plot_path, qc_plot, width = 7, height = 5)
cat(sprintf("✅ QC plot successfully saved to: %s\n", plot_path))

cat("\n════════════════════════════════════════════════════════════\n")
cat("  QC Check Complete! All systems nominal.\n")
cat("════════════════════════════════════════════════════════════\n\n")


