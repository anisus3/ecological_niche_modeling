library(ggplot2)
library(dplyr)
library(purrr)
library(tools)

# Folder with permutation importance CSVs
pi_dir <- "outputs/PI"

# All PI_*.csv files
pi_files <- list.files(pi_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(pi_files) == 0) {
  stop("No PI_*.csv files found in 'outputs/PI'")
}

# --- Read all PI files ---
pi_all <- map_dfr(pi_files, function(f) {
  df <- read.csv(f)
  
  # filename without path/extension, e.g., "PI_Asia_region_01"
  region_label <- file_path_sans_ext(basename(f))
  
  df$region_label <- region_label
  df
})

# Cast to factors and fix level order
pi_all <- pi_all |>
  mutate(
    importance    = as.numeric(importance),
    variable      = as.factor(variable),
    region_label  = factor(region_label, levels = unique(region_label))
  )

# --- Plot, color by region_label ---
out_pdf <- file.path(pi_dir, "PI_all_variables_by_region_barplots.pdf")

pdf(out_pdf, width = 12, height = 8)

ggplot(pi_all, aes(x = region_label, y = importance, fill = region_label)) +
  geom_col() +
  facet_wrap(~ variable, ncol = 4) +   # 4×2 для 8 переменных
  theme_bw() +
  labs(
    x = "Region (file name)",
    y = "Permutation importance",
    fill = "Region (file)",
    title = "Regional permutation importance by predictor"
  ) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(face = "bold"),
    legend.position   = "bottom",
    legend.key.width  = unit(1.5, "lines")
  )

dev.off()

cat("✔ Plot saved to file:", out_pdf, "\n")

# ==========================================================
# Heatmap: variables × regions (by PI), with original order
# ==========================================================

# variable order — same as env_r_regional
var_order <- names(env_r_regional)

# aggregate
pi_mat <- pi_all |>
  group_by(variable, region_label) |>
  summarise(importance = mean(importance, na.rm = TRUE), .groups = "drop") |>
  mutate(
    variable     = factor(variable, levels = var_order),
    region_label = factor(region_label, levels = unique(as.character(region_label)))
  )

heatmap_pdf <- file.path(pi_dir, "PI_heatmap_variables_by_region.pdf")

pdf(heatmap_pdf, width = 10, height = 6)

ggplot(pi_mat, aes(x = region_label, y = variable, fill = importance)) +
  geom_tile(color = "grey90") +
  scale_y_discrete(limits = rev(var_order)) +
  
  # --- Custom red-white palette ---
  scale_fill_gradientn(
    colours = c("white", "#FFDFD0" ,"#FFB89C", "#FF5D3C", "#DE0003", "#710006"),
    na.value = "white",
    name = "PI"
  ) +
  
  theme_bw() +
  labs(
    x = "Region (file name)",
    y = "Predictor",
    title = "Regional permutation importance (heatmap)"
  ) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90")
  )


dev.off()

cat("✔ Heatmap saved to file:", heatmap_pdf, "\n")
