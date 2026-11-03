#utils_stats.R
library(vegan)

run_permanova_suite <- function(dist, metadata, group_col, covariates, perm = 999) {
  f1 <- as.formula(paste0("dist ~ ", group_col))
  res1 <- vegan::adonis2(f1, data = metadata, permutations = perm)
  
  f2 <- as.formula(paste0("dist ~ ", paste(c(group_col, covariates), collapse = " + ")))
  res2 <- vegan::adonis2(f2, data = metadata, permutations = perm)
  
  f_terms_a <- as.formula(paste0("dist ~ ", paste(c(covariates, group_col), collapse = " + ")))
  res_terms_a <- vegan::adonis2(f_terms_a, data = metadata, by = "terms", permutations = perm)
  
  f_terms_b <- as.formula(paste0("dist ~ ", paste(c(group_col, covariates), collapse = " + ")))
  res_terms_b <- vegan::adonis2(f_terms_b, data = metadata, by = "terms", permutations = perm)
  
  list(
    group_only = res1,
    full_model = res2,
    terms_cov_then_group = res_terms_a,
    terms_group_then_cov = res_terms_b
  )
}