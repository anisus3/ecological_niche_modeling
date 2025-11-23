# =============================================================================
# Script purpose:
# Build an ensemble model (EMwmean) from previously trained single models
# using TSS and ROC for weighting. Stores the ensemble in 'ens'.
# =============================================================================

#---------- Libraries ----------
library(biomod2)

#---------- 1. Preconditions ----------
stopifnot(exists("models"))

#---------- 2. Ensemble modeling ----------
ens <- BIOMOD_EnsembleModeling(
  bm.mod               = models,
  models.chosen        = "all",
  em.by                = "all",
  em.algo              = "EMwmean",       # weighted mean
  metric.select        = c("TSS","ROC"),
  metric.select.thresh = c(0, 0),
  metric.eval          = c("TSS","ROC"),
  EMwmean.decay        = "proportional"
)

#---------- Note on warnings ----------
# You may see warnings like:
#   In coords.roc(...): 'transpose=TRUE' is deprecated and will be removed.
# These come from internal ROC code and can be safely ignored.
