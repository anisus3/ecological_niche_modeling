# =============================================================================
# Script purpose:
# Load WorldClim v1.4 bioclim layers at a chosen arc-minute resolution,
# (optionally) crop by the extent of thinned occurrence points, compute a
# Spearman correlation matrix from random samples, identify highly correlated
# pairs (|r| ≥ threshold), perform a greedy de-correlation variable selection,
# and save the final subset of layers as a GeoTIFF.
# Paths/data:
# - Downloads WorldClim to `env_data/` (created if missing).
# - Optional cropping extent from `GBIF_thinned_<THIN_KM>km.csv` in the project root.
# - Output GeoTIFF written to `env_data/WorldClim_v1.4_selected_*vars_*min_<THIN_KM>km.tif`.
# =============================================================================

#---------- Libraries ----------
library(geodata)
library(terra)
library(dplyr)

#---------- Constants ----------
WC_RES_ARCMIN <- 2.5          # WorldClim resolution (arc-minutes)
THIN_KM       <- 25            # thinning radius used previously (km), for file naming
occ_file      <- sprintf("GBIF_thinned_%dkm.csv", THIN_KM)
outdir        <- "env_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#---------- 1. Load WorldClim bioclim ----------
bio <- geodata::worldclim_global(var = "bio", res = WC_RES_ARCMIN, path = outdir)
stopifnot(terra::nlyr(bio) > 0)
message(sprintf("✔ WorldClim v1.4 loaded: %g arc-min (~%.1f km), layers: %d",
                WC_RES_ARCMIN, WC_RES_ARCMIN * 1.85, terra::nlyr(bio)))

#---------- 2. Optional crop by occurrence extent ----------
if (file.exists(occ_file)) {
  occ <- read.csv(occ_file) %>%
    filter(!is.na(decimal_latitude), !is.na(decimal_longitude))
  if (nrow(occ) > 0) {
    occ_sp <- terra::vect(occ, geom = c("decimal_longitude","decimal_latitude"), crs = "EPSG:4326")
    bio <- terra::crop(bio, terra::ext(occ_sp))
    message(sprintf("Cropped climate rasters to the extent of: %s", occ_file))
  } else {
    message(sprintf("No valid coordinates in %s — using global extent", occ_file))
  }
} else {
  message(sprintf("File %s not found — using global extent", occ_file))
}

#---------- 3. Sample values & compute Spearman correlation ----------
set.seed(42)
# Note: if the cropped area is small, reduce 'size' to avoid errors
sample_points <- terra::spatSample(bio, size = 100000, method = "random", as.points = TRUE)
bio_values    <- terra::extract(bio, sample_points, df = TRUE)[, -1]

cor_m <- stats::cor(bio_values, use = "pairwise.complete.obs", method = "spearman")

#---------- 4. Find highly correlated pairs ----------
cor_threshold <- 0.75   # absolute Spearman threshold (tune as needed)

cor_df <- as.data.frame(as.table(cor_m)) |>
  dplyr::rename(var1 = Var1, var2 = Var2, r = Freq) |>
  dplyr::filter(var1 != var2) |>
  dplyr::mutate(abs_r = abs(r)) |>
  dplyr::filter(abs_r >= cor_threshold) |>
  dplyr::arrange(desc(abs_r))

# remove duplicates of A–B vs B–A
cor_df <- cor_df[!duplicated(t(apply(cor_df[, c("var1","var2")], 1, sort))), ]

cat(sprintf("\nPairs with |r| ≥ %.2f:\n", cor_threshold))
if (nrow(cor_df) == 0) {
  cat("   — none —\n")
} else {
  print(cor_df, row.names = FALSE)
}

#---------- 5. Greedy selection of low-collinearity variables ----------
# Build unique list of correlated pairs (as matrix with two columns)
if (nrow(cor_df) > 0) {
  pairs <- unique(t(apply(cor_df[, c("var1","var2")], 1, sort)))
  colnames(pairs) <- c("a","b")
} else {
  pairs <- matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("a","b")))
}

drop <- character()
if (nrow(pairs) > 0) {
  for (i in seq_len(nrow(pairs))) {
    a <- pairs[i, "a"]; b <- pairs[i, "b"]
    if (a %in% drop || b %in% drop) next
    # greedy rule: drop the second variable (you may replace by an importance-based rule later)
    drop <- c(drop, b)
  }
}

target_vars <- setdiff(names(bio), drop)
target_ids  <- which(names(bio) %in% target_vars)
target_ids  <- sort(unique(target_ids))

cat(
  "Selected variable indices: ", paste(target_ids, collapse = ", "), "\n",
  "Selected variable names: ", paste(target_vars, collapse = ", "), "\n",
  sep = ""
)

#---------- 6. Save final stack ----------
bio_sel <- bio[[target_ids]]
# ensure simple names like "bio1", "bio2", ...
names(bio_sel) <- names(bio)[target_ids]

out_tif <- file.path(
  outdir,
  sprintf("WorldClim_v1.4_selected_%dvars_%gmin_%dkm.tif",
          length(target_ids), WC_RES_ARCMIN, THIN_KM)
)

terra::writeRaster(bio_sel, out_tif, overwrite = TRUE)

cat(
  "\n✔ Final stack saved to: ", out_tif, "\n",
  "Layers: ", paste(names(bio_sel), collapse = ", "), "\n",
  sep = ""
)
