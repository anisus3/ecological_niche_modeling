library(ggplot2)
library(dplyr)
library(purrr)
library(tools)

# Folder with PD CSVs
pd_dir <- "outputs/PD"

# All files PD_*.csv
pd_files <- list.files(pd_dir, pattern = "^PD_.*\\.csv$", full.names = TRUE)

# Table with file info and linetype grouping
region_info <- data.frame(
  region     = file_path_sans_ext(basename(pd_files)),
  file_index = seq_along(pd_files),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # groups of 4 files: 1–4, 5–8, 9–12, ...
    block_id = (file_index - 1) %/% 4,
    # cycle block_id to 1,2,3,1,2,3,... for linetype pattern
    linetype_group = (block_id %% 3) + 1,
    linetype_group = factor(linetype_group, levels = 1:3)
  )

# Read all files and add region/subregion name from filename
pd_all <- map_dfr(pd_files, function(f) {
  df <- read.csv(f)
  region_name <- file_path_sans_ext(basename(f))
  df$region <- region_name
  df
})

# Variable order for facets
var_order <- c(
  "bio1",
  "bio2",
  "bio7",
  "bio12",
  "bio14",
  "anthro_cropland_mosaic",
  "anthro_cropland_strict",
  "anthro_urban"
)

# Add linetype info and set variable order
pd_all <- pd_all %>%
  left_join(region_info, by = "region") %>%
  mutate(variable = factor(variable, levels = var_order))

out_pdf <- file.path(pd_dir, "PD_all_regions_comparison.pdf")

pdf(out_pdf, width = 12, height = 8)
ggplot(pd_all, aes(x = x_value,
                   y = pred_mean,
                   colour = region,
                   linetype = linetype_group)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotdash"),
    guide = "none"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "lines")
  )
dev.off()

cat("✔ Plot saved to file:", out_pdf, "\n")
