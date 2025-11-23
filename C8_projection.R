# =============================================================================
# Script purpose:
# Project trained single models (biomod2) onto the current climate stack 'env_r',
# then compute the ensemble forecasting from these projections. Both new and
# legacy argument names are provided for compatibility. Parallel execution uses
# nb.cpu = 8; outputs are written to the default biomod2 directories.
# =============================================================================

#---------- Libraries ----------
library(biomod2)
library(raster)

#---------- 1. Preconditions ----------
stopifnot(exists("models"), exists("ens"), exists("env_r"))

#---------- 2. Projection to current climate ----------
proj_cur <- BIOMOD_Projection(
  bm.mod           = models,          # new argument name
  modeling.output  = models,          # legacy name for compatibility
  new.env          = env_r,
  proj.id          = "current",
  proj.name        = "current",       # duplicate legacy name to pass checks
  models.chosen    = "all",
  selected.models  = "all",           # legacy alias to avoid 'unused' warnings
  CV.keep.models   = TRUE,
  build.clamping.mask = TRUE,
  metric.binary    = c("TSS","ROC"),
  binary.meth      = c("TSS","ROC"),  # keep both variants
  do.stack         = TRUE,
  compress         = "xz",
  nb.cpu           = 8
)

#---------- 3. Ensemble forecasting ----------
ens_cur <- BIOMOD_EnsembleForecasting(
  bm.em              = ens,           # new argument name
  EM.output          = ens,           # legacy name for compatibility
  bm.proj            = proj_cur,
  projection.output  = proj_cur,      # legacy alias
  metric.binary      = c("TSS","ROC"),
  binary.meth        = c("TSS","ROC"),
  compress           = "xz",
  nb.cpu             = 8
)

message("✅ Projection complete. Rasters saved in ./outputs/biomod/Arion_pooled_run/")
