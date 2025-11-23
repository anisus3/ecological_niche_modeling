# =============================================================================
# Script purpose:
# Build a background mask by computing an MCP around occurrences, buffering it
# in an azimuthal equidistant CRS by BUF_KM, projecting back to the raster CRS,
# and cropping/masking the full predictor stack (climate + anthropogenic).
# Saves: raster mask (GeoTIFF), buffered polygon (GPKG), and masked stack (GeoTIFF).
# =============================================================================

#---------- Libraries ----------
library(terra)

#---------- Helpers ----------
# Build azimuthal equidistant proj string centered at (lat0, lon0)
aeqd_string <- function(lat0, lon0) {
  sprintf("+proj=aeqd +lat_0=%.6f +lon_0=%.6f +datum=WGS84 +units=m +no_defs", lat0, lon0)
}

#---------- 1. Preconditions & tags ----------
stopifnot(exists("occ"), exists("env_terra"), inherits(env_terra, "SpatRaster"))
stopifnot(exists("WC_RES_ARCMIN"), exists("THIN_KM"))
wc_tag <- sprintf("%g", WC_RES_ARCMIN)   # e.g., 2.5 → "2.5", 10 → "10"

#---------- 2. Buffer parameter (km → m) ----------
if (!exists("BUF_KM")) BUF_KM <- 200
BUF_M <- BUF_KM * 1000

#---------- 3. Occurrence points (WGS84) ----------
pts_ll <- terra::vect(occ, geom = c("decimal_longitude","decimal_latitude"), crs = "EPSG:4326")

#---------- 4. MCP (convex hull) in lon/lat ----------
mcp_ll <- terra::convHull(pts_ll)

#---------- 5. Buffer in local AEQD (meters), then project back to raster CRS ----------
r_template <- env_terra[[1]]
cent <- terra::centroids(mcp_ll)
lon0 <- terra::geom(cent)[1, "x"]; lat0 <- terra::geom(cent)[1, "y"]
aeqd <- aeqd_string(lat0, lon0)

mcp_m       <- terra::project(mcp_ll, aeqd)
mcp_buf_m   <- terra::buffer(mcp_m, width = BUF_M)
mcp_buf_prj <- terra::project(mcp_buf_m, terra::crs(r_template))

#---------- 6. Crop + mask (robust: rasterize polygon, then mask) ----------
r_crop <- terra::crop(r_template, terra::ext(mcp_buf_prj))
mask_r <- terra::rasterize(mcp_buf_prj, r_crop, field = 1)
bg_mask <- terra::mask(r_crop, mask_r)

#---------- 7. Apply mask to ALL predictors (climate + anthropogenic) ----------
env_masked <- terra::mask(terra::crop(env_terra, terra::ext(bg_mask)), mask_r)

#---------- 8. Save outputs ----------
dir.create("outputs", showWarnings = FALSE)

mask_out <- sprintf("outputs/background_mask_%dkm_%smin_%dkm.tif", BUF_KM, wc_tag, THIN_KM)
buf_out  <- sprintf("outputs/background_buffer_%dkm_%smin_%dkm.gpkg", BUF_KM, wc_tag, THIN_KM)
env_out  <- sprintf("outputs/env_all_masked_%dkm_%smin_%dkm.tif",   BUF_KM, wc_tag, THIN_KM)

try(terra::writeRaster(bg_mask, mask_out, overwrite = TRUE), silent = TRUE)
try(terra::writeVector(mcp_buf_prj, buf_out, filetype = "GPKG", overwrite = TRUE), silent = TRUE)
try(terra::writeRaster(env_masked, env_out, overwrite = TRUE), silent = TRUE)

cat(
  "\n— Background mask —\n",
  "Buffer around MCP: ", BUF_KM, " km\n",
  "WorldClim: ", wc_tag, " arc-min; thinning: ", THIN_KM, " km\n",
  "Files:\n  • ", mask_out, "\n  • ", buf_out, "\n  • ", env_out, "\n",
  sep = ""
)
