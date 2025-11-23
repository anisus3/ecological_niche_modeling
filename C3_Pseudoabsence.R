# =============================================================================
# Script purpose:
# Generate pseudo-absence points inside the background mask (bg_mask).
# Returns dataframe 'pa' with columns decimal_longitude, decimal_latitude, resp=0.
# =============================================================================

library(terra)
library(dplyr)

stopifnot(exists("bg_mask"))

# Parameters
N_PA  <- 12000   # number of pseudo-absence points
SEED  <- 42      # RNG seed
set.seed(SEED)

# Non-NA cells of background mask
bg_cells <- which(!is.na(values(bg_mask)))
stopifnot(length(bg_cells) > 0)

# Random sample
n_draw <- min(N_PA, length(bg_cells))
pa_cells <- sample(bg_cells, size = n_draw)

# Convert to coordinates
pa <- terra::xyFromCell(bg_mask, pa_cells) %>%
  as.data.frame() %>%
  rename(decimal_longitude = x, decimal_latitude = y) %>%
  mutate(resp = 0L)

# Result is stored in variable `pa`
cat("✔ Pseudo-absences generated:", n_draw, "points\n")
