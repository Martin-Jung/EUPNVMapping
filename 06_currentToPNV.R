# --- --- #
# This script calculates the change of land from current to PNV estimate
# using a reclassified 2018 MAES layer as baseline
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")

# Set run cores for parallization
options('ibis.nthread' = cores)
options('ibis.runparallel' = ifelse(cores>1, FALSE, FALSE) ) # Set to FALSE for slower, but memory proof runs

# -------------------- #
# Get reference raster for the given grain
background <- terra::rast( "data/background_ref.tif")

## Legend
# 1 - Woodland and forest
# 2 - Heathland and shrub
# 3 - Grassland
# 4 - Sparsely vegetated areas
# 5 - Cropland
# 6 - Urban
# 7 - Wetlands
# 8 - Rivers and lakes
# 9 - Marine inlets and transitional waters
# 11 - Pasture
ll <- list.files("corine/",full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
ll <- ll[grep("Corine_", ll)]
ras <- rast(ll)
names(ras) <- stringr::str_remove(tools::file_path_sans_ext( basename(ll) ), "Corine_2018_")
names(ras) <- c("Woodland.and.forest", "Pasture", "Heathland.and.shrub", "Grassland", "Sparsely.vegetated.areas",
                "Cropland", "Urban", "Wetlands", "Rivers.and.lakes", "Marine.inlets.and.transitional.waters")

# Align with background
if(!terra::compareGeom(ras, background, stopOnError = FALSE)){
  if(sf::st_crs(ras) != sf::st_crs(background)) ras <- terra::project(ras, background)
  # Crop and mask
  ras <- terra::crop(ras, background)
  ras <- terra::extend(ras, background)
  ras <- terra::mask(ras, background)
}

# Since we have a equal area, projection, the area is just the grid size.
# But multiply regardless to get to km2 
ras <- (ras / 10000) * terra::cellSize(ras, unit = "km")

plot(ras$Pasture)

# --- #
# Get the ensemble predictions
path_ensemble <- "./ensemble/Masked/"
ll <- list.files(path_ensemble, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
#ll <- ll[grep("Threshold", ll)]
ll <- ll[grep("Masked", ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

#### Summarize directly from the predictions ####
# Since those are probability layer, we multiply directly with area

results <- data.frame()
for(lyr in ll){
  preds <- terra::rast(lyr)
  
  new_med <- ras * preds[['q50']]
  new_low <- ras * preds[['q05']]
  new_high <- ras * preds[['q95']]
  df <- cbind(
    terra::global(new_med, "sum",na.rm = TRUE),
    terra::global(new_low, "sum",na.rm = TRUE),
    terra::global(new_high, "sum",na.rm = TRUE)
  ) |> as.data.frame() 
  names(df) <- c("q50", "q05", "q95")
  df <- df |> tibble::rownames_to_column("variable")
}
# TODO: Figure out what to actually summarize. Imo it makes only sense over
# managed areas, e.g. cropland, pasture, urban, plantation.

#### Process the thresholded layers ####
# Load each layer, then assess the amount of area that requires transitioning 
# from current (Corine) to potential

results <- data.frame()

for(lyr in ll){
  message(basename(lyr))
  
  tr <- rast(lyr)
  # tr[tr<1] <- NA
  
  # Calculate zonal statistics for each
  for(i in 1:terra::nlyr(tr)){
    z <- terra::zonal(ras,z = tr[[i]], fun = "sum", na.rm = TRUE, as.raster = FALSE)
    z <- tidyr::pivot_longer(z,cols = -1)
    names(z)[1] <- "lyr"
    z$lyr <- names(tr[[i]])
    z$pnv <- stringr::str_remove(tools::file_path_sans_ext(basename(lyr)), "Threshold__")
    # Also add the total amount of land per layer
    z$pnv_km2 <- terra::global(tr[[i]], "sum", na.rm = TRUE)[,1] / (1000)
    results <- dplyr::bind_rows(results, z)
    rm(z)
  }
}

# Save the output
saveRDS(results, "results/Summary_current2PNV_km2.rds")
