library(ggplot2)
library(dplyr)
library(purrr)
library(tools)

# Folder with CSV files
pd_dir <- "outputs/PD"

# All files PD_*.csv
pd_files <- list.files(pd_dir, pattern = "^PD_.*\\.csv$", full.names = TRUE)

# Read files
pd_all <- map_dfr(pd_files, function(f) {
  df <- read.csv(f)
  df$region <- file_path_sans_ext(basename(f))
  df
})

# Order of variables
var_order <- c(
  "bio1", "bio2", "bio7", "bio12", "bio14",
  "anthro_cropland_mosaic",
  "anthro_cropland_strict",
  "anthro_urban"
)

pd_all <- map_dfr(pd_files, function(f) {
  df <- read.csv(f)
  df$region <- file_path_sans_ext(basename(f))
  df$source_file <- basename(f)   # added the source file name
  df
})

# --- SET FILE INDEX FOR COLOR ---
region_info <- data.frame(
  region     = unique(pd_all$region),
  file_index = seq_along(unique(pd_all$region)),
  stringsAsFactors = FALSE
)

pd_all <- pd_all %>%
  left_join(region_info, by = "region")

# --- SAVE THE MERGED PD TABLE FOR ALL REGIONS ---
out_csv <- file.path(pd_dir, "PD_all_regions.csv")
write.csv(pd_all, out_csv, row.names = FALSE)
cat("✔ The table with all PD values has been saved to:", out_csv, "\n")

# --- COLORS IN BLOCKS OF 4 FILES ---
color_scheme <- c("green3", "red3", "blue3", "black")  # cycle of 4 colors

region_colors <- sapply(region_info$file_index, function(i) {
  block_id <- (i - 1) %/% 4
  color_scheme[(block_id %% 4) + 1]
})

names(region_colors) <- region_info$region

# --- KEEP THE LINETYPES UNCHANGED ---
unique_regions <- region_info$region
region_linetypes <- rep(c("solid", "dashed", "twodash", "dotted"),
                        length.out = length(unique_regions))
names(region_linetypes) <- unique_regions

# --- Output PDF ---
out_pdf <- file.path(pd_dir, "PD_all_regions_comparison.pdf")

pdf(out_pdf, width = 46, height = 31)

ggplot(pd_all, aes(x = x_value,
                   y = pred_mean,
                   colour = region,
                   linetype = region)) +
  geom_line(size = 2) +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  scale_colour_manual(
    values = region_colors,
    name   = "Region"
  ) +
  scale_linetype_manual(
    values = region_linetypes,
    name   = "Region"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(3.5, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 1.1),
    strip.text       = element_text(size = 40),
    axis.text        = element_text(size = 40),
    axis.title       = element_text(size = 20)
  )

dev.off()

cat("✔ Plot saved to:", out_pdf, "\n")
