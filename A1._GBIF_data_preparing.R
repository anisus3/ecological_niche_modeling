# =============================================================================
# Script purpose:
# Reads all GBIF/iNat occurrence CSVs from 'GBIF_manual' with auto-delimiter
# detection, cleans & merges records, then performs geodesic thinning so that
# points are at least THIN_KM kilometers apart (greedy, block-wise).
# Outputs a thinned CSV and a GeoPackage (GPKG) layer and prints a short report.
# Before running, make sure a directory with CSV files exists (default:
# `GBIF_manual`), and each file contains columns:
# species, scientificName, decimalLatitude, decimalLongitude.
# =============================================================================

#---------- Libraries ----------
library(dplyr)
library(readr)
library(sf)
library(purrr)
library(tibble)

#---------- Constants ----------
# Generalization parameter (in kilometers)
THIN_KM <- 25
THIN_M  <- THIN_KM * 1000

# Input folder with CSVs
indir <- "GBIF_manual"

#---------- Functions ----------
# Auto-detect GBIF/iNat CSV delimiter and coerce key columns
read_gbif_csv <- function(path){
  loc <- locale(encoding = "UTF-8")
  try_delims <- c(",", "\t", ";")
  cols_keep <- c(
    "species","scientificName","decimalLatitude","decimalLongitude",
    "coordinateUncertaintyInMeters","basisOfRecord","year",
    "datasetKey","occurrenceID","countryCode"
  )
  for (del in try_delims) {
    df <- try(
      read_delim(
        path, delim = del, locale = loc, show_col_types = FALSE,
        guess_max = 100000, col_types = cols(.default = col_character()),
        col_select = any_of(cols_keep)
      ),
      silent = TRUE
    )
    if (!inherits(df, "try-error") &&
        any(c("decimalLatitude","decimalLongitude") %in% names(df))) {
      
      # ensure all expected columns exist
      missing <- setdiff(cols_keep, names(df))
      if (length(missing) > 0) df[missing] <- NA_character_
      
      message(sprintf("✔ %s — read with delimiter '%s'", basename(path), del))
      
      # typing/clean-up + safe species formation
      df <- df %>%
        mutate(
          species = coalesce(na_if(species, ""), na_if(scientificName, "")),
          decimalLatitude  = as.numeric(gsub(",", ".", decimalLatitude)),
          decimalLongitude = as.numeric(gsub(",", ".", decimalLongitude)),
          year = suppressWarnings(as.integer(year)),
          coordinateUncertaintyInMeters =
            suppressWarnings(as.numeric(gsub(",", ".", coordinateUncertaintyInMeters)))
        )
      return(df)
    }
  }
  stop("Could not detect a valid delimiter or required columns in file: ", path)
}

# Greedy block-wise geodesic thinning (≥ dist_m meters)
thin_block_greedy <- function(block_df, kept_sf, dist_m){
  if (nrow(block_df) == 0) return(block_df[0, , drop = FALSE])
  
  block_sf <- st_as_sf(
    block_df, coords = c("decimal_longitude","decimal_latitude"),
    crs = 4326, remove = FALSE
  )
  
  # discard points too close to already kept points
  if (!is.null(kept_sf) && nrow(kept_sf) > 0) {
    too_close_list <- st_is_within_distance(block_sf, kept_sf, dist = dist_m)
    far_mask <- lengths(too_close_list) == 0L
    block_sf <- block_sf[far_mask, , drop = FALSE]
    if (nrow(block_sf) == 0) {
      return(
        block_sf[0, c("species","decimal_latitude","decimal_longitude"), drop = FALSE] %>%
          st_drop_geometry()
      )
    }
  }
  
  keep_idx <- integer(0)
  if (nrow(block_sf) > 0) {
    for (i in seq_len(nrow(block_sf))) {
      if (length(keep_idx) == 0) {
        keep_idx <- c(keep_idx, i)
      } else {
        cand_sf   <- block_sf[i, ]
        chosen_sf <- block_sf[keep_idx, ]
        close_to_chosen <- any(st_is_within_distance(cand_sf, chosen_sf, dist = dist_m)[[1]])
        if (!close_to_chosen) keep_idx <- c(keep_idx, i)
      }
    }
  }
  
  kept_block_sf <- block_sf[keep_idx, , drop = FALSE]
  kept_block_sf %>%
    st_drop_geometry() %>%
    dplyr::select(species, decimal_latitude, decimal_longitude)
}

#---------- 1. Discover & read inputs ----------
files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(files) >= 1)

d_list <- files %>% map(read_gbif_csv)

#---------- 2. Merge & clean ----------
occ_all <- bind_rows(d_list)

occ_raw <- occ_all %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  distinct(decimalLatitude, decimalLongitude, species, .keep_all = TRUE) %>%
  transmute(
    species,
    decimal_latitude  = decimalLatitude,
    decimal_longitude = decimalLongitude,
    year,
    coordinateUncertaintyInMeters,
    basisOfRecord,
    datasetKey, occurrenceID, countryCode
  )

#---------- 3. Shuffle for stable greedy behavior ----------
set.seed(42)
occ_small <- occ_raw %>%
  dplyr::select(species, decimal_latitude, decimal_longitude) %>%
  { .[sample.int(nrow(.)), , drop = FALSE] }

#---------- 4. Block configuration ----------
block_size <- 5000L
dist_m     <- THIN_M

n <- nrow(occ_small)
idx_starts <- seq(1L, n, by = block_size)
idx_ends   <- pmin(idx_starts + block_size - 1L, n)

kept_df <- tibble(species = character(),
                  decimal_latitude = numeric(),
                  decimal_longitude = numeric())
kept_sf <- NULL

#---------- 5. Block-wise thinning ----------
for (i in seq_along(idx_starts)) {
  s <- idx_starts[i]; e <- idx_ends[i]
  message(sprintf("Block %d/%d: rows %d–%d (threshold: %d km)",
                  i, length(idx_starts), s, e, THIN_KM))
  block_df <- occ_small[s:e, , drop = FALSE]
  kept_add <- thin_block_greedy(block_df, kept_sf, dist_m)
  if (nrow(kept_add) > 0) {
    kept_df <- bind_rows(kept_df, kept_add)
    new_sf <- st_as_sf(_
                       
