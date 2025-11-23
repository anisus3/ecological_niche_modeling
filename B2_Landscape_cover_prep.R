# =============================================================================
# Script purpose:
# Derive land-cover fractions (0–100%) from C3S Land Cover 2022 on the
# WorldClim grid (WGS84). For selected class codes and groups, build binary
# masks, resample by area-averaging to the target WorldClim raster, save each
# fraction as GeoTIFF, and also write a combined stack.
# Paths/data:
# - Input NetCDF: `env_data/Landscape_use/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc`.
# - Requires a WorldClim template `bio` SpatRaster already in memory (from climate script).
# - Outputs written to `env_data/Landscape_use/out_c3s_derived/` as per-variable and stacked GeoTIFFs.
# =============================================================================

#---------- Libraries ----------
library(terra)
library(parallel)

#---------- Constants & paths ----------
datadir <- "env_data/Landscape_use"
ncfile  <- file.path(datadir, "C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc")
outdir  <- file.path(datadir, "out_c3s_derived")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Target output resolution label (arc-min) for filenames only
RES_ARCMIN <- 0.5

# Selected codes and logical groups (C3S LCCS classes)
SELECTED_CODES <- c(10, 11, 12, 20, 30, 40, 190)
GROUPS <- list(
  cropland_strict = c(10, 11, 12, 20),
  cropland_mosaic = c(30, 40),
  cropland_total  = c(10, 11, 12, 20, 30, 40),
  urban           = c(190)
)

#---------- Terra threading ----------
terraOptions(
  threads  = max(1, detectCores(logical = TRUE) - 1),
  todisk   = TRUE,
  progress = 1
)
message("Threads for terra: ", terraOptions()$threads)

#---------- 1. Read C3S LCCS and align CRS/grid ----------
stopifnot(file.exists(ncfile))
r_lccs <- rast(ncfile, subds = "lccs_class")  # categorical class layer
message("✔ Read C3S file: ", ncfile)

# Expect a WorldClim template 'bio' (SpatRaster) to be present in the session
stopifnot(exists("bio"), inherits(bio, "SpatRaster"))
target_r   <- bio[[1]]          # template grid/CRS (WGS84)
target_crs <- crs(target_r)

# If input has missing CRS, assume WGS84 like the product specification
if (is.na(crs(r_lccs)) || nchar(crs(r_lccs)) == 0) {
  message("⚠️  Input C3S layer has no CRS. Assuming WGS84: ", target_crs)
  crs(r_lccs) <- target_crs
}

# If CRS differ, project with nearest-neighbor to preserve categories
if (!terra::same.crs(r_lccs, target_r)) {
  message("↻ CRS differs. Reprojecting C3S → target CRS (WGS84)…")
  r_lccs <- project(r_lccs, target_r, method = "near")
  message("   ✔ Now: ", crs(r_lccs))
}

message(
  "WorldClim grid: ",
  paste0(res(target_r)[1] * 60, "×", res(target_r)[2] * 60, " arc-min; CRS = ", crs(target_r))
)

#---------- 2. Helper functions ----------
# Binary mask for one or multiple class codes (categorical → 0/1 numeric)
make_binary <- function(r_codes, codes) {
  r_bin <- ifel(r_codes %in% codes, 1, 0)
  levels(r_bin) <- NULL
  as.numeric(r_bin)
}

# Fraction (0..1) on target grid via area-weighted average
fraction_to_target <- function(r_bin, target_r) {
  resample(r_bin, target_r, method = "average")
}

# Save fraction (0..1) as percent (0..100) GeoTIFF
write_fraction_percent <- function(r_frac01, basename_stub) {
  r_pct <- r_frac01 * 100
  nm    <- sprintf("%s_%sarcmin.tif", basename_stub, RES_ARCMIN)
  fpath <- file.path(outdir, nm)
  writeRaster(
    r_pct, fpath, overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"),
    datatype = "FLT4S"
  )
  message("   ✔ Saved: ", basename(fpath))
  list(path = fpath, rast = r_pct)
}

#---------- 3. Fractions for individual class codes ----------
out_list  <- list()
out_paths <- character(0)

for (code in SELECTED_CODES) {
  nm_stub <- sprintf("C3S_LC2022_code%03d_frac", code)
  message("→ Computing fraction for code: ", code)
  
  r_bin <- make_binary(r_lccs, code)
  stopifnot(terra::same.crs(r_bin, target_r))
  
  r_frac01 <- fraction_to_target(r_bin, target_r)  # 0..1
  saved    <- write_fraction_percent(r_frac01, nm_stub)
  
  out_list[[sprintf("code%03d", code)]] <- saved$rast  # percent
  out_paths <- c(out_paths, saved$path)
}

#---------- 4. Fractions for groups (e.g., cropland, urban) ----------
for (g in names(GROUPS)) {
  codes  <- GROUPS[[g]]
  nm_stub <- sprintf("C3S_LC2022_%s_frac", g)
  message("→ Group: ", g, " (", paste(codes, collapse = ","), ")")
  
  r_bin <- make_binary(r_lccs, codes)
  stopifnot(terra::same.crs(r_bin, target_r))
  
  r_frac01 <- fraction_to_target(r_bin, target_r)  # 0..1
  saved    <- write_fraction_percent(r_frac01, nm_stub)
  
  out_list[[g]] <- saved$rast                      # percent
  out_paths <- c(out_paths, saved$path)
}

#---------- 5. Write combined stack (percent) on WorldClim grid ----------
r_stack <- rast(out_list)
names(r_stack) <- names(out_list)

f_stack <- file.path(outdir, sprintf("C3S_LC2022_selected_frac_%sarcmin_stack.tif", RES_ARCMIN))
writeRaster(r_stack, f_stack, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"))

message("\n✅ Done. Saved to: ", outdir)
message("📄 Individual files:"); print(basename(out_paths))
message("📦 Combined stack: ", basename(f_stack))
