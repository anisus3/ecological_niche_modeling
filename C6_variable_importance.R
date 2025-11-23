# =============================================================================
# Script purpose:
# Extract variable importance from fitted models, aggregate mean and SD by
# algorithm and variable, save the results to 'environmental_significance.csv',
# print the table, and list ensemble members via get_built_models(ens).
# =============================================================================

#---------- Libraries ----------
library(dplyr)

#---------- 1. Compute variable importance (mean & SD by algo × variable) ----------
vi <- get_variables_importance(models)
head(vi)

algos <- unique(vi$algo)
variables <- unique(vi$expl.var)

variable_importance <- data.frame()

for (m in algos) {
  for (v in variables) {
    sub <- vi[vi$algo == m & vi$expl.var == v, "var.imp", drop = TRUE]
    mean_imp <- mean(sub, na.rm = TRUE)
    sd_imp   <- sd(sub, na.rm = TRUE)
    variable_importance <- rbind(
      variable_importance,
      data.frame(algo = m, expl.var = v,
                 mean_import = mean_imp, sd_import = sd_imp)
    )
  }
}

#---------- 2. Save results to CSV (project root) ----------
out_csv <- file.path(getwd(), "environmental_significance.csv")
write.csv(variable_importance, out_csv, row.names = FALSE)

#---------- 3. Print table and inspect ensemble composition ----------
variable_importance

# Ensemble overview and individual model weights
get_built_models(ens)
