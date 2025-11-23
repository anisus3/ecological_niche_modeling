# =============================================================================
# Script purpose:
# Load WorldClim bioclim layers (v2.1) at given resolution, crop them to the
# predefined region, keep only selected variables (bio1,2,7,12,14),
# and save as a single GeoTIFF stack with region name in the filename.
# Paths/data:
# - Downloads WorldClim to `env_data/` (created if missing).
# - Uses hardcoded ROI extent in this script (adjust as needed).
# - Outputs GeoTIFF: `env_data/WorldClim_v2.1_selected_<WC_RES_ARCMIN>min_<REGION_NAME>.tif`.
# =============================================================================

#---------- Libraries ----------
library(terra)
library(geodata)

#---------- User-defined parameters ----------
REGION_NAME <- "Asia"        # <--- ИМЯ РЕГИОНА: "Asia", "Europe"
WC_RES_ARCMIN <- 0.5         # WorldClim resolution in arc-minutes

# Bounding box (нужный регион)

roi <- terra::ext(40, 143, 25, 63)

#roi <- terra::ext(6, 32, 40, 52)

# Пример для Европы (если понадобится):
# roi <- terra::ext(6, 32, 40, 52)

#---------- Output directory ----------
outdir <- "env_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#---------- 1. Load WorldClim bioclim (v2.1) ----------
bio <- geodata::worldclim_global(var = "bio", res = WC_RES_ARCMIN, path = outdir)
stopifnot(terra::nlyr(bio) > 0)

message(sprintf("✔ WorldClim v2.1 loaded: %g arc-min (~%.1f km), layers: %d",
                WC_RES_ARCMIN, WC_RES_ARCMIN * 1.85, terra::nlyr(bio)))
message("  Layer names: ", paste(names(bio), collapse = ", "))

# wc2.1_<tag>_bio_<id>
res_tag <- if (WC_RES_ARCMIN == 0.5) {
  "30s"
} else {
  sprintf("%gm", WC_RES_ARCMIN)
}

# Needed variables
bio_ids    <- c(1, 2, 7, 12, 14)
keep_names <- sprintf("wc2.1_%s_bio_%d", res_tag, bio_ids)

# Validate presence
missing_layers <- setdiff(keep_names, names(bio))
if (length(missing_layers) > 0) {
  stop("Не найдены слои: ",
       paste(missing_layers, collapse = ", "),
       "\nДоступные имена: ",
       paste(names(bio), collapse = ", "))
}

#---------- 2. Crop to region ----------
bio_crop <- terra::crop(bio, roi)
message("✔ Cropped to region: ", REGION_NAME)

# Select only needed layers
bio_sel <- bio_crop[[keep_names]]
message("✔ Selected layers: ", paste(names(bio_sel), collapse = ", "))

#---------- 3. Save final stack ----------
stack_path <- file.path(
  outdir,
  sprintf("WorldClim_v2.1_selected_%gmin_%s.tif", WC_RES_ARCMIN, REGION_NAME)
)

terra::writeRaster(
  bio_sel,
  stack_path,
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6")
)

message("✔ Final stack saved: ", stack_path)
message("  Layers: ", paste(names(bio_sel), collapse = ", "))
