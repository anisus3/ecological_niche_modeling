# =============================================================================
# Script purpose:
# Define user modeling options for GBM, RF, and MAXENT (GLM uses defaults),
# then run BIOMOD_Modeling with random CV (10 reps, 70/30), evaluate via TSS/ROC,
# compute variable importance, and use prevalence = 0.5 on previously formatted data.
# =============================================================================

#---------- Libraries ----------
library(biomod2)

#---------- 1. Preconditions ----------
stopifnot(exists("formatted"))

#---------- 2. User-defined hyperparameters ----------
user.val <- list(
  GBM.binary.gbm.gbm = list(
    for_all_datasets = list(
      n.trees = 2500,
      interaction.depth = 3,
      shrinkage = 0.01,
      bag.fraction = 0.5,
      train.fraction = 1
    )
  ),
  RF.binary.randomForest.randomForest = list(
    for_all_datasets = list(
      ntree = 1500,
      mtry  = NULL
    )
  ),
  MAXENT.binary.MAXENT.MAXENT = list(
    for_all_datasets = list(
    )
  )
  # IMPORTANT: no GLM block here — let it use defaults and auto-build the formula
)

#---------- 3. Build modeling options ----------
opt <- bm_ModelingOptions(
  data.type = "binary",
  models    = c("GLM","GBM","RF","MAXENT"),
  strategy  = "user.defined",
  user.val  = user.val,
  bm.format = formatted  # provide data access so GLM can auto-build its formula
)

#---------- 4. Run BIOMOD modeling ----------
models <- BIOMOD_Modeling(
  bm.format         = formatted,
  modeling.id       = "Arion_pooled_run",
  models            = c("GLM","GBM","RF","MAXENT"),
  CV.strategy       = "random",
  CV.nb.rep         = 10,
  CV.perc           = 0.7,
  CV.do.full.models = FALSE,
  OPT.user          = opt,             # OPT.strategy is embedded in 'opt'
  metric.eval       = c("TSS","ROC"),  # in 4.x "ROC" is valid; "AUCroc" also ok
  var.import        = 3,
  prevalence        = 0.5,
  nb.cpu            = 8,
  do.progress       = TRUE
)
