# =============================================================================
# Script purpose:
# Print the list of predictor layers in 'env_r', combine presence (occ) and
# pseudo-absence (pa) coordinates into a single table with response, and
# prepare data for 'biomod2' using BIOMOD_FormatingData (no PA generation here).
# One occurrence per raster cell is enforced via 'filter.raster = TRUE'.
# =============================================================================

#---------- Libraries ----------
library(biomod2); library(dplyr)

#---------- 1. Inspect predictor stack ----------
cat("Слои в env_r:\n", paste(names(env_r), collapse = ", "), "\n")

#---------- 2. Combine presence and pseudo-absence ----------
all_xy <- bind_rows(
  occ %>% transmute(x = decimal_longitude, y = decimal_latitude, resp = 1L),
  pa  %>% transmute(x = decimal_longitude, y = decimal_latitude, resp = 0L)
)

#---------- 3. Prepare BIOMOD data ----------
resp.name <- "Arion_pooled"
resp.var  <- all_xy$resp
resp.xy   <- as.matrix(all_xy[,c("x","y")])

formatted <- BIOMOD_FormatingData(
  resp.var   = resp.var,
  expl.var   = env_r,      # SpatRaster / RasterStack
  resp.xy    = resp.xy,
  resp.name  = resp.name,
  PA.nb.rep  = 0,
  filter.raster = TRUE     # one point per raster cell
)
