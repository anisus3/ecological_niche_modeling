# =============================================================================
# Script purpose:
# Visualize BIOMOD2 ensemble outputs: habitat suitability map, binary map
# (TSS threshold), and clamping mask. Save the main suitability map as
# high-resolution PNG. Optionally overlay occurrence points.
# =============================================================================

#---------- Libraries ----------
library(raster)

#---------- 1. Main projection (habitat suitability) ----------
r <- raster("Arion.pooled/proj_current/proj_current_Arion.pooled_ensemble.tif")
plot(r, main = "Arion pooled – ensemble habitat suitability")

#---------- 2. Save high-resolution PNG ----------
png("Arion_pooled_ensemble_urban.png", width = 10000, height = 6000, res = 800)
plot(
  r,
  main = "Arion pooled – ensemble habitat suitability",
  col  = rev(terrain.colors(100))  # optional: replace with viridis or heat.colors
)
# Optionally overlay occurrence points:
# points(occ$decimal_longitude, occ$decimal_latitude, pch = 21, bg = "red", cex = 0.6)
dev.off()

# Overlay points interactively if needed
points(occ$decimal_longitude, occ$decimal_latitude, pch = 21, bg = "red", cex = 0.6)

#---------- 3. Binary map by TSS ----------
r_bin <- raster("Arion.pooled/proj_current/proj_current_Arion.pooled_ensemble_TSSbin.tif")
plot(r_bin, main = "Binary map (TSS threshold)")

#---------- 4. Clamping mask ----------
clamp <- raster("Arion.pooled/proj_current/proj_current_ClampingMask.tif")
plot(clamp, main = "Clamping mask")
