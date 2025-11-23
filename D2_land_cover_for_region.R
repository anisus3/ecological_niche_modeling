# =============================================================================
# Script purpose:
# Read C3S Land Cover 2022, align its CRS to a WorldClim template, compute
# land-cover fractions (0–100%) for selected LCCS code groups on the template
# grid, crop to a given region, and save each layer and a combined stack.
# Paths/data:
# - Input NetCDF: `env_data/Landscape_use/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc`.
# - Climate template from previous script: `env_data/WorldClim_v2.1_selected_<RES_ARCMIN>min_<REGION_NAME>.tif`.
# - Outputs (per group + stack) to `env_data/Landscape_use/out_c3s_derived/`.
# =============================================================================

#---------- Libraries ----------
library(terra)

#---------- Settings ----------
RES_ARCMIN  <- 0.5          # target resolution label (arc-min) for filenames
REGION_NAME <- "Asia"     # имя региона (должно совпадать с климатическим скриптом)

#---------- Paths ----------
datadir <- "env_data/Landscape_use"
ncfile  <- file.path(datadir, "C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc")
outdir  <- file.path(datadir, "out_c3s_derived")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# КЛИМАТИЧЕСКИЙ ШАБЛОН из предыдущего скрипта (WorldClim)
# Файл создаётся климатическим скриптом в папке "env_data":
# env_data/WorldClim_v2.1_selected_<RES_ARCMIN>min_<REGION_NAME>.tif
template_path <- sprintf(
  "env_data/WorldClim_v2.1_selected_%gmin_%s.tif",
  RES_ARCMIN, REGION_NAME
)
stopifnot(file.exists(template_path))
template <- terra::rast(template_path)[[1]]  # single layer as template

# C3S LCCS code groups (v2.1.1)
GROUPS <- list(
  cropland_strict = c(10,11,12,20),
  cropland_mosaic = c(30,40),
  urban           = c(190)
)

# bbox для REGION_NAME
bbox_region <- terra::ext(40, 143, 25, 63)
# для Asia можно будет тут задать другой прямоугольник

#---------- 1. Read C3S LCCS and ensure CRS ----------
stopifnot(file.exists(ncfile))
r_lccs <- try(terra::rast(ncfile, subds = "lccs_class"), silent = TRUE)
if (inherits(r_lccs, "try-error")) r_lccs <- terra::rast(ncfile)
message("✔ Read C3S LCCS: ", ncfile)

# Ensure CRS exists (WGS84); if missing, copy from template
if (is.na(terra::crs(r_lccs)) || nchar(terra::crs(r_lccs)) == 0) {
  terra::crs(r_lccs) <- terra::crs(template)
}
# If CRS differs, reproject to template CRS (preserve native resolution; nearest for categorical)
if (!terra::same.crs(r_lccs, template)) {
  message("🔁 CRS differs. Projecting C3S → template CRS (near)…")
  r_lccs <- terra::project(
    r_lccs,
    terra::crs(template),
    method = "near",
    res    = terra::res(r_lccs)
  )
}

#---------- 2. Helper functions ----------
# (1) Binary mask (0/1) on native grid (ensure numeric, not factor)
make_binary <- function(r_codes, codes) {
  r_bin <- terra::ifel(r_codes %in% codes, 1, 0)
  levels(r_bin) <- NULL
  as.numeric(r_bin)
}
# (2) Fraction on the template grid: resample average of 0/1 → [0..1]
fraction_to_template <- function(r_bin, template) {
  terra::resample(r_bin, template, method = "average")
}

#---------- 3. Derive fractions for groups and save ----------
out_list <- list()
for (g in names(GROUPS)) {
  codes <- GROUPS[[g]]
  message("→ Group: ", g, " (codes: ", paste(codes, collapse = ","), ")")
  
  r_bin    <- make_binary(r_lccs, codes)
  r_frac01 <- fraction_to_template(r_bin, template)                  # 0..1
  r_pct    <- terra::clamp(r_frac01 * 100, 0, 100, values = TRUE)    # 0..100 %
  
  # Crop to REGION_NAME after alignment to template
  r_pct <- terra::crop(r_pct, bbox_region, snap = "out")
  
  nm_tif <- file.path(
    outdir,
    sprintf("C3S_LC2022_%s_frac_%sarcmin_%s.tif", g, RES_ARCMIN, REGION_NAME)
  )
  terra::writeRaster(
    r_pct, nm_tif, overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE","ZLEVEL=6"),
    datatype = "FLT4S"
  )
  message("   ✔ Written: ", nm_tif)
  
  out_list[[g]] <- r_pct
}

#---------- 4. Combined stack (strict + mosaic + urban) ----------
r_stack <- terra::rast(out_list)
names(r_stack) <- names(out_list)

stack_tif <- file.path(
  outdir,
  sprintf("C3S_LC2022_selected_frac_%sarcmin_stack_%s.tif", RES_ARCMIN, REGION_NAME)
)
terra::writeRaster(
  r_stack, stack_tif, overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE","ZLEVEL=6"),
  datatype = "FLT4S"
)

# Consistency checks with climate template
stopifnot(terra::same.crs(r_stack, template))
stopifnot(terra::compareGeom(r_stack, template, stopOnError = FALSE))

message(
  "\n✅ Done for region: ", REGION_NAME,
  "\n   ", file.path(outdir, sprintf("C3S_LC2022_cropland_strict_frac_%sarcmin_%s.tif", RES_ARCMIN, REGION_NAME)),
  "\n   ", file.path(outdir, sprintf("C3S_LC2022_cropland_mosaic_frac_%sarcmin_%s.tif", RES_ARCMIN, REGION_NAME)),
  "\n   ", file.path(outdir, sprintf("C3S_LC2022_urban_frac_%sarcmin_%s.tif", RES_ARCMIN, REGION_NAME)),
  "\n   ", stack_tif
)
