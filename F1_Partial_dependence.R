# ==============================================================
# Regional partial dependence for BIOMOD ensemble – GRID version
# ==============================================================

library(raster)
library(dplyr)
library(purrr)
library(biomod2)
library(ggplot2)

# --------------------------------------------------------------
# 0. Region settings
# --------------------------------------------------------------
REGION_NAME <- "Asia"

# Main region bounds
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

nm <- names(env_r_regional)
mask <- nm %in% names(name_map)
nm[mask] <- name_map[nm[mask]]
names(env_r_regional) <- nm

cat("Predictor names after renaming:\n  ",
    paste(names(env_r_regional), collapse = ", "), "\n", sep = "")

# --------------------------------------------------------------
# 1. Create output directory
# --------------------------------------------------------------
dir.create("outputs/PD", recursive = TRUE, showWarnings = FALSE)


# --------------------------------------------------------------
# 2. Make subregion boundaries
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
# 3. Prepare BIOMOD list of models
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
# 4. Ensemble prediction with safety wrapper
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
# 5. Partial dependence for one variable
# --------------------------------------------------------------

compute_pd_for_var <- function(var_name, env_sub, models_list,
                               n_grid = 25) {
  
  x_vals <- env_sub[[var_name]]
  
  x_min <- quantile(x_vals, 0.01, na.rm = TRUE)
  x_max <- quantile(x_vals, 0.99, na.rm = TRUE)
  
  grid_vals <- seq(x_min, x_max, length.out = n_grid)
  
  pd_results <- lapply(grid_vals, function(gval) {
    newdata <- env_sub
    newdata[[var_name]] <- gval
    
    preds <- predict_ensemble(newdata, models_list)
    
    data.frame(
      variable = var_name,
      x_value  = gval,
      pred_mean = mean(preds, na.rm = TRUE),
      pred_sd   = sd(preds, na.rm = TRUE)
    )
  })
  
  bind_rows(pd_results)
}




# ==============================================================
# 6. LOOP over all subregions
# ==============================================================

region_id <- 1   # counter 1..12

# IMPORTANT: iterate latitude bands north → south
for (i in N_LAT_SPLIT:1) {       # latitude bands (north -> south)
  for (j in 1:N_LON_SPLIT) {     # longitude bands (west -> east)
    
    # ---- define the extent ----
    sub_ext <- subregion_extent(i, j)
    
    # upper-right corner coordinates (lon = xmax, lat = ymax)
    ur_lon <- xmax(sub_ext)
    ur_lat <- ymax(sub_ext)
    
    cat(
      sprintf(
        "Processing subregion %d ... | upper-right corner: lon = %.4f, lat = %.4f\n",
        region_id, ur_lon, ur_lat
      )
    )
    
    # ---- crop raster to subregion ----
    env_sub_r <- crop(env_r_regional, sub_ext, snap = "out")
    
    # skip empty tiles (e.g. ocean/void areas)
    if (ncell(env_sub_r) == 0) {
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
    
    # ---- sample pixels ----
    set.seed(123)
    N_TARGET   <- 100000
    N_SAMPLES  <- min(N_TARGET, nrow(env_df))
    cat("  Sampling", N_SAMPLES, "pixels out of", nrow(env_df), "\n")
    
    idx <- sample(seq_len(nrow(env_df)), N_SAMPLES)
    env_sub <- env_df[idx, pred_cols, drop = FALSE]
    
    # ---- compute PD for this region ----
    pd_region <- map_dfr(
      pred_cols,
      ~ compute_pd_for_var(.x,
                           env_sub     = env_sub,
                           models_list = models_list,
                           n_grid      = 25)
    )
    
    # ---- save ----
    out_file <- sprintf("outputs/PD/PD_%s_region_%02d.csv",
                        REGION_NAME, region_id)
    write.csv(pd_region, out_file, row.names = FALSE)
    
    cat("  ✔ Saved:", out_file, "\n")
    
    region_id <- region_id + 1
  }
}

cat("\n✔ All subregions completed.\n")
