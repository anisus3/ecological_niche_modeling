# ==============================================================
# Regional permutation importance for BIOMOD ensemble – GRID version
# ==============================================================

library(raster)
library(dplyr)
library(purrr)
library(biomod2)

# --------------------------------------------------------------
# 0. Region settings
# --------------------------------------------------------------
REGION_NAME <- "Asia"

# Main region bounds (must match PD script)
GRID_LAT_MIN <- 25
GRID_LAT_MAX <- 63
GRID_LON_MIN <- 40
GRID_LON_MAX <- 143

# Grid splitting: 3 rows (latitude), 4 cols (longitude)
N_LAT_SPLIT <- 3
N_LON_SPLIT <- 4

# Explicit mapping: long WC2.1 names → short bio* names
name_map <- c(
  "wc2.1_30s_bio_1"  = "bio1",
  "wc2.1_30s_bio_2"  = "bio2",
  "wc2.1_30s_bio_7"  = "bio7",
  "wc2.1_30s_bio_12" = "bio12",
  "wc2.1_30s_bio_14" = "bio14"
)

# env_r_regional – predictor stack for the full region (same as PD script)
nm <- names(env_r_regional)
mask <- nm %in% names(name_map)
nm[mask] <- name_map[nm[mask]]
names(env_r_regional) <- nm

cat("Predictor names after renaming:\n  ",
    paste(names(env_r_regional), collapse = ", "), "\n", sep = "")

# --------------------------------------------------------------
# 1. Create output directory for PI
# --------------------------------------------------------------
dir.create("outputs/PI", recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------
# 2. Make subregion boundaries (same as PD)
# --------------------------------------------------------------

lat_breaks <- seq(GRID_LAT_MIN, GRID_LAT_MAX,
                  length.out = N_LAT_SPLIT + 1)
lon_breaks <- seq(GRID_LON_MIN, GRID_LON_MAX,
                  length.out = N_LON_SPLIT + 1)

# Function to generate extent for (i=lat band, j=lon band)
subregion_extent <- function(i, j) {
  extent(
    lon_breaks[j], lon_breaks[j + 1],
    lat_breaks[i], lat_breaks[i + 1]
  )
}

# --------------------------------------------------------------
# 3. Prepare BIOMOD list of models (same as PD)
# --------------------------------------------------------------

if (inherits(models, "BIOMOD.models.out")) {
  model_names <- BIOMOD_LoadModels(models)
  models_list <- mget(model_names, envir = .GlobalEnv)
} else if (is.list(models)) {
  models_list <- models
} else {
  stop("Unknown model format in 'models'")
}

models_list <- Filter(Negate(is.null), models_list)

# --------------------------------------------------------------
# 4. Ensemble prediction with safety wrapper (same as PD)
# --------------------------------------------------------------

predict_ensemble <- function(newdata, models_list) {
  pred_list <- lapply(models_list, function(m) {
    out <- try(predict(m, newdata = newdata), silent = TRUE)
    if (inherits(out, "try-error")) {
      return(rep(NA_real_, nrow(newdata)))
    } else {
      return(as.numeric(out))
    }
  })
  
  pred_mat <- do.call(cbind, pred_list)
  rowMeans(pred_mat, na.rm = TRUE)
}

# --------------------------------------------------------------
# 5. Permutation importance for one subregion
#    Metric: 1 - cor(base_pred, perm_pred)
# --------------------------------------------------------------

compute_pi_for_region <- function(env_sub, models_list,
                                  pred_cols,
                                  n_rep = 10) {
  # base ensemble prediction
  base_pred <- predict_ensemble(env_sub[, pred_cols, drop = FALSE],
                                models_list)
  base_pred <- as.numeric(base_pred)
  
  # if all NA or zero variance, return NA immediately
  if (all(is.na(base_pred)) || sd(base_pred, na.rm = TRUE) == 0) {
    return(
      data.frame(
        variable   = pred_cols,
        importance = NA_real_,
        sd         = NA_real_
      )
    )
  }
  
  pi_list <- lapply(pred_cols, function(var) {
    drops <- numeric(n_rep)
    
    for (k in seq_len(n_rep)) {
      env_perm <- env_sub
      env_perm[[var]] <- sample(env_perm[[var]])
      
      perm_pred <- predict_ensemble(env_perm[, pred_cols, drop = FALSE],
                                    models_list)
      perm_pred <- as.numeric(perm_pred)
      
      drops[k] <- 1 - cor(base_pred, perm_pred, use = "complete.obs")
    }
    
    data.frame(
      variable   = var,
      importance = mean(drops, na.rm = TRUE),
      sd         = sd(drops,   na.rm = TRUE)
    )
  })
  
  bind_rows(pi_list) |>
    arrange(desc(importance))
}

# ==============================================================
# 6. LOOP over all subregions (numbered 1..12, same as PD)
# ==============================================================

region_id <- 1   # counter 1..12

for (i in N_LAT_SPLIT:1) {       # latitude bands (north -> south)
  for (j in 1:N_LON_SPLIT) {     # longitude bands
    
    cat("Processing subregion", region_id, "...\n")
    
    # ---- define the extent ----
    sub_ext <- subregion_extent(i, j)
    
    # ---- crop raster to subregion ----
    env_sub_r <- crop(env_r_regional, sub_ext, snap = "out")
    
    # skip empty tiles (e.g. ocean/void areas)
    if (is.null(env_sub_r) || ncell(env_sub_r) == 0) {
      cat("  Skipping region", region_id, "(empty)\n")
      region_id <- region_id + 1
      next
    }
    
    # convert to data frame
    env_df <- as.data.frame(env_sub_r, xy = FALSE, na.rm = TRUE)
    
    if (nrow(env_df) < 100) {
      cat("  Skipping region", region_id, "(too few pixels)\n")
      region_id <- region_id + 1
      next
    }
    
    # ---- variable names ----
    pred_cols <- names(env_sub_r)
    
    # ---- sample pixels (same as PD) ----
    set.seed(123)
    N_SAMPLES <- min(10000, nrow(env_df))
    idx <- sample(seq_len(nrow(env_df)), N_SAMPLES)
    env_sub <- env_df[idx, pred_cols, drop = FALSE]
    
    # ---- compute PI for this region ----
    pi_region <- compute_pi_for_region(
      env_sub    = env_sub,
      models_list = models_list,
      pred_cols  = pred_cols,
      n_rep      = 10   # adjust if needed
    )
    
    # ---- save ----
    out_file <- sprintf("outputs/PI/PI_%s_region_%02d.csv",
                        REGION_NAME, region_id)
    write.csv(pi_region, out_file, row.names = FALSE)
    
    cat("  ✔ Saved:", out_file, "\n")
    
    region_id <- region_id + 1
  }
}

cat("\n✔ All subregions completed.\n")
