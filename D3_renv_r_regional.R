# =============================================================================
# Script purpose:
# Prepare predictor stack for SDM: load climate (WorldClim) and add anthropogenic
# layers (C3S cropland strict/mosaic, urban) as 0–1 fractions on the same grid.
# Align rasters to the climate template if needed, save a combined GeoTIFF,
# and create a RasterStack 'env_r_regional' compatible with biomod2 workflows.
# =============================================================================

#---------- Libraries ----------
library(terra)
library(raster)

#---------- Settings ----------
WC_RES_ARCMIN <- 0.5      # WorldClim resolution (arc-min, for labels)
THIN_KM       <- 25       # thinning radius (only for reporting)
REGION_NAME   <- "Asia" # <<< имя региона: "Asia", "Europe" и т.п.

dir.create("outputs", showWarnings = FALSE)

#---------- 1. Input paths ----------
occ_file <- sprintf("GBIF_thinned_%dkm.csv", THIN_KM)

# Климатический стек из предыдущего скрипта
clim_file <- file.path(
  "env_data",
  sprintf("WorldClim_v2.1_selected_%gmin_%s.tif", WC_RES_ARCMIN, REGION_NAME)
)

# C3S-слои из скрипта derivation C3S_LC2022_..._frac_..._%s.tif
c3s_dir <- file.path("env_data", "Landscape_use", "out_c3s_derived")
fn_cropland_strict <- file.path(
  c3s_dir,
  sprintf("C3S_LC2022_cropland_strict_frac_%sarcmin_%s.tif", WC_RES_ARCMIN, REGION_NAME)
)
fn_cropland_mosaic <- file.path(
  c3s_dir,
  sprintf("C3S_LC2022_cropland_mosaic_frac_%sarcmin_%s.tif", WC_RES_ARCMIN, REGION_NAME)
)
fn_urban <- file.path(
  c3s_dir,
  sprintf("C3S_LC2022_urban_frac_%sarcmin_%s.tif", WC_RES_ARCMIN, REGION_NAME)
)

# Basic existence checks
stopifnot(file.exists(clim_file))
stopifnot(file.exists(fn_cropland_strict),
          file.exists(fn_cropland_mosaic),
          file.exists(fn_urban))
if (!file.exists(occ_file)) {
  message("⚠ Occurrence file not found: ", occ_file, " (only reported in summary).")
}

#---------- 2. Load climate stack ----------
env_terra_regional <- terra::rast(clim_file)
cat("✔ Climate source: ", clim_file, "\n",
    "  Layers: ", paste(names(env_terra_regional), collapse = ", "), "\n",
    sep = "")

template <- env_terra_regional[[1]]

#---------- 3. Function to align rasters ----------
align_to_template <- function(r, template) {
  if (terra::compareGeom(r, template, stopOnError = FALSE)) {
    r
  } else {
    terra::resample(r, template, method = "bilinear")
  }
}

#---------- 4. Load anthropogenic layers (percent → fraction) ----------
message("\n🏙️🌾  Adding anthropogenic layers (C3S: cropland_strict, cropland_mosaic, urban)…")

r_cropland_strict <- terra::rast(fn_cropland_strict)
r_cropland_mosaic <- terra::rast(fn_cropland_mosaic)
r_urban           <- terra::rast(fn_urban)

# Align to climate template if needed (continuous → bilinear)
r_cropland_strict <- align_to_template(r_cropland_strict, template)
r_cropland_mosaic <- align_to_template(r_cropland_mosaic, template)
r_urban           <- align_to_template(r_urban,           template)

# Convert percent (0–100) → fractions (0–1) and clamp
r_cropland_strict <- terra::clamp(r_cropland_strict / 100, 0, 1, values = TRUE)
names(r_cropland_strict) <- "anthro_cropland_strict"

r_cropland_mosaic <- terra::clamp(r_cropland_mosaic / 100, 0, 1, values = TRUE)
names(r_cropland_mosaic) <- "anthro_cropland_mosaic"

r_urban <- terra::clamp(r_urban / 100, 0, 1, values = TRUE)
names(r_urban) <- "anthro_urban"

#---------- 5. Combine climate + anthropogenic ----------
env_terra_regional <- c(
  env_terra_regional,
  r_cropland_strict,
  r_cropland_mosaic,
  r_urban
)

stopifnot(terra::compareGeom(env_terra_regional, template, stopOnError = FALSE))
stopifnot(terra::same.crs(env_terra_regional, template))
message("✔ Geometry and CRS match the climate template.")

#---------- 6. Export for biomod2 ----------
out_tif_regional <- file.path(
  "outputs",
  sprintf("env_%gmin_with_anthro_%s.tif", WC_RES_ARCMIN, REGION_NAME)
)

terra::writeRaster(
  env_terra_regional,
  out_tif_regional,
  overwrite = TRUE,
  filetype  = "GTiff"
)

env_r_regional <- raster::stack(out_tif_regional)

cat(
  "\n— Data prepared for modeling (region: ", REGION_NAME, ") —\n",
  "  Occurrences: ", occ_file, "\n",
  "  Climate:     ", clim_file, "\n",
  "  Anthro:      ", paste(
    basename(c(fn_cropland_strict, fn_cropland_mosaic, fn_urban)),
    collapse = ", "
  ), "\n",
  "  Resolution:  ", WC_RES_ARCMIN, " arc-min (~",
  round(WC_RES_ARCMIN * 1.85, 1), " km)\n",
  "  Thinning:    ", THIN_KM, " km\n",
  "  Layers:      ", paste(names(env_r_regional), collapse = ", "), "\n",
  "  Output:      ", out_tif_regional, "\n",
  sep = ""
)
