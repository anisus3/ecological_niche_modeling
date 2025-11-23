# =============================================================================
# Script purpose:
# Remove records whose latitudes match a list of suspicious values (within a
# numeric tolerance). Logs removed rows to a separate CSV, creates a backup of
# the original file, and overwrites the input CSV with the cleaned data.
# =============================================================================

#---------- Libraries ----------
library(readr)
library(dplyr)
library(stringr)

#---------- Functions ----------
# Vectorized latitude flag: TRUE if any |x - bad_val| <= tol
is_bad_lat <- function(x, bad_vals, tol){
  if (length(bad_vals) == 0) return(rep(FALSE, length(x)))
  rowSums(abs(outer(x, bad_vals, `-`)) <= tol) > 0
}

#---------- 1. Parameters ----------
in_csv  <- "GBIF_thinned_25km.csv"   # input CSV to clean (e.g., "GBIF_raw_combined.csv")
tol     <- 1e-5                      # numeric tolerance for latitude comparison
bad_lat <- c(
  78.065,
  33.7176,
  0       # add more here (e.g., 0, 90, -90, etc.)
)

#---------- 2. Read & basic checks ----------
stopifnot(file.exists(in_csv))
df <- readr::read_csv(in_csv, show_col_types = FALSE)

# Ensure required columns exist
stopifnot(all(c("decimal_latitude","decimal_longitude") %in% names(df)))

# Coerce to numeric (safeguard)
df <- df %>%
  mutate(
    decimal_latitude  = as.numeric(decimal_latitude),
    decimal_longitude = as.numeric(decimal_longitude)
  )

#---------- 3. Flag & split ----------
flag_bad <- is_bad_lat(df$decimal_latitude, bad_lat, tol)
removed  <- df[flag_bad, , drop = FALSE]
kept     <- df[!flag_bad, , drop = FALSE]

message(sprintf("Records to be removed: %d (of %d)", nrow(removed), nrow(df)))

#---------- 4. Write removal log (if any) ----------
if (nrow(removed) > 0) {
  log_path <- str_replace(in_csv, "\\.csv$", "_removed_by_lat_log.csv")
  readr::write_csv(removed, log_path)
  message("Removal log written: ", log_path)
}

#---------- 5. Backup & overwrite cleaned CSV ----------
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
backup_path <- sprintf("%s.bak_%s", in_csv, ts)
ok <- file.copy(in_csv, backup_path, overwrite = FALSE)
if (ok) {
  message("Backup created: ", backup_path)
} else {
  message("Warning: backup not created (a backup with the same name may already exist).")
}

readr::write_csv(kept, in_csv)
message("Done. Cleaned file overwritten: ", in_csv)
