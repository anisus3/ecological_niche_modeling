# =============================================================================
# Script purpose:
# Prepare predictor stack for SDM: load climate (WorldClim) and add anthropogenic
# layers (C3S cropland strict/mosaic, urban) as 0–1 fractions on the same grid.
# Align rasters to the climate template if needed, save a combined GeoTIFF,
# and create a RasterStack 'env_r' compatible with biomod2 workflows.
# =============================================================================

#---------- Libraries ----------
library(terra)
library(raster)

#---------- Functions ----------
# Align a raster to a template (same grid/extent/CRS). Uses bilinear for continuous data.
align_to_template <- function(r, template) {
  if (terra::compareGeom(r, template, stopOnError = FALSE)) {
    r
  } else {
    terra::resample(r, template, method = "bilinear")
  }
}

#---------- 1. Parameters ----------
WC_RES_ARCMIN <- 0.5     # WorldClim resolution (arc-min, used for file labels)
THIN_KM       <- 25      # thinning radius used in occurrence preprocessing (for reporting)
dir.create("outputs", showWarnings = FALSE)

#---------- 2. Input paths ----------
occ_file  <- sprintf("GBIF_thinned_%dkm.csv", THIN_KM)
clim_file <- file.path("env_data", "WorldClim_selected_Asia.tiff")

c3s_dir <- file.path("env_data", "Landscape_use", "out_c3s_derived")
fn_cropland_strict <- file.path(c3s_dir, sprintf("C3S_LC2022_cropland_strict_frac_%sarcmin.tif", WC_RES_ARCMIN))
fn_cropland_mosaic <- file.path(c3s_dir, sprintf("C3S_LC2022_cropland_mosaic_frac_%sarcmin.tif", WC_RES_ARCMIN))
fn_urban           <- file.path(c3s_dir, sprintf("C3S_LC2022_urban_frac_%sarcmin.tif",           WC_RES_ARCMIN))

# Basic existence checks
stopifnot(file.exists(clim_file))
stopifnot(file.exists(fn_cropland_strict), file.exists(fn_cropland_mosaic), file.exists(fn_urban))
if (!file.exists(occ_file)) message("⚠ Occurrence file not found: ", occ_file, " (only reported in summary).")

#---------- 3. Load climate stack ----------
env_terra <- terra::rast(clim_file)
cat("✔ Climate source: ", clim_file, "\n",
    "  Layers: ", paste(names(env_terra), collapse = ", "), "\n", sep = "")

template <- env_terra[[1]]

#---------- 4. Load anthropogenic layers (percent → fraction) ----------
message("\n🏙️🌾  Adding anthropogenic layers (C3S: cropland_strict, cropland_mosaic, urban)…")

r_cropland_strict <- terra::rast(fn_cropland_strict)
r_cropland_mosaic <- terra::rast(fn_cropland_mosaic)
r_urban           <- terra::rast(fn_urban)

# Align to climate template if needed (continuous data → bilinear)
r_cropland_strict <- align_to_template(r_cropland_strict, template)
r_cropland_mosaic <- align_to_template(r_cropland_mosaic, template)
r_urban           <- align_to_template(r_urban,           template)

# Convert percent to 0–1 fractions and clamp to [0,1]
r_cropland_strict <- terra::clamp(r_cropland_strict / 100, 0, 1, values = TRUE); names(r_cropland_strict) <- "anthro_cropland_strict"
r_cropland_mosaic <- terra::clamp(r_cropland_mosaic / 100, 0, 1, values = TRUE); names(r_cropland_mosaic) <- "anthro_cropland_mosaic"
r_urban           <- terra::clamp(r_urban           / 100, 0, 1, values = TRUE); names(r_urban)           <- "anthro_urban"

#---------- 5. Combine climate + anthropogenic ----------
env_terra <- c(env_terra, r_cropland_strict, r_cropland_mosaic, r_urban)

stopifnot(terra::compareGeom(env_terra, template, stopOnError = FALSE))
stopifnot(terra::same.crs(env_terra, template))
message("✔ Geometry and CRS match the climate template.")

#---------- 6. Export for biomod2 ----------
out_tif  <- file.path("outputs", sprintf("env_%gmin_with_anthro_Asia.tif", WC_RES_ARCMIN))
terra::writeRaster(env_terra, out_tif, overwrite = TRUE, filetype = "GTiff")

env_r <- raster::stack(out_tif)

cat(
  "\n— Data prepared for modeling —\n",
  "  Occurrences: ", occ_file, "\n",
  "  Climate:     ", clim_file, "\n",
  "  Anthro:      ", paste(basename(c(fn_cropland_strict, fn_cropland_mosaic, fn_urban)), collapse = ", "), "\n",
  "  Resolution:  ", WC_RES_ARCMIN, " arc-min (~", round(WC_RES_ARCMIN * 1.85, 1), " km)\n",
  "  Thinning:    ", THIN_KM, " km\n",
  "  Layers:      ", paste(names(env_r), collapse = ", "), "\n",
  "  Output:      ", out_tif, "\n",
  sep = ""
)
