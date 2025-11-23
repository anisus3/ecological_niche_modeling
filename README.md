# Overview
Scripts for ecological niche modeling (ENM) and projections using BIOMOD2. The pipeline prepares occurrence data, assembles climatic and land-cover predictors (WorldClim, Sentinel-3/C3S Land Cover), fits ensemble species distribution models, projects them to chosen regions, and quantifies predictor influence via permutation importance (PI) and partial dependence (PD).

# Data
- Occurrences: GBIF/iNaturalist CSVs placed in `GBIF_manual/`, cleaned and thinned to `GBIF_thinned_<km>.csv` and GPKG.
- Climate: WorldClim v1.4/v2.1 bioclim rasters downloaded to `env_data/`.
- Land cover: C3S Land Cover 2022 (Sentinel-3) NetCDF in `env_data/Landscape_use/`, derived fraction rasters in `env_data/Landscape_use/out_c3s_derived/`.
- Outputs: model objects, rasters, and plots in `outputs/` and BIOMOD-generated subfolders.

# Methods
- Preprocess occurrences (delimiter auto-detect, cleaning, geodesic thinning).
- Prepare predictors: climate subsets and cropped stacks, land-cover fractions aligned to the climate grid.
- Modeling: BIOMOD2 single/ensemble models, background/pseudoabsence generation, projection, ensemble metrics (AUC/TSS).
- Interpretation (optional): regional PI (permutation drop in correlation) and PD curves on a 3×4 spatial grid; graphics for PI/PD across regions.

# Folder structure
- `A*.R` – occurrence ingestion, cleaning, thinning.
- `B*.R` – climate and land-cover predictor preparation (WorldClim, C3S/Sentinel-3).
- `C*.R` – modeling workflow (background, formatting, fitting, projection, ensemble).
- `D*.R` – optional region-specific predictor stacks.
- `F*.R` – optional PI/PD calculation and plotting.
- `env_data/` – downloaded/derived rasters; `outputs/` – model exports and diagnostics.

# How to reproduce
1) Set R libs (or `renv` via `D3_renv_r_regional.R` if used).  
2) Block A: run occurrence prep (`A1._GBIF_data_preparing.R`, `A2. Cleaning_from_bad_data.R`) to produce thinned CSV/GPKG.  
3) Block B: build predictor stacks (`B1_Climate_layer_preparing.R`, `B2_Landscape_cover_prep.R`).  
4) Block C: run modeling sequence (`C1`–`C9`) to fit, project, and visualize ensembles.  
5) Optional Block D: region-focused stacks (`D1`–`D3`).  
6) Optional Block F: compute PI/PD on subregions and plot (`F2_Par_Depen_Asia_blocks.R`, `F2_Permutation_importance...R`, `F3_PD_Graphics.R`, `F4_PI_graphs.R`).  
Adjust region names/paths inside scripts as needed before execution.

# Requirements
- R (≥4.x recommended). Core packages: `biomod2`, `terra`, `raster`, `geodata`, `sf`, `readr`, `dplyr`, `purrr`, `ggplot2`, `tibble`.  
- For land cover: C3S NetCDF file `C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc` placed under `env_data/Landscape_use/`.  
- Sufficient disk space for WorldClim downloads and derived rasters.

# Citation
- Cite BIOMOD2 (Thuiller et al.) for modeling, WorldClim for climate layers, and C3S Land Cover (Sentinel-3) for land cover products.  
- Include GBIF/iNaturalist data DOIs for occurrence records used. Add your species/region-specific citation as appropriate in publications.
