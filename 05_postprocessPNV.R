# --- --- #
# This script post-processes the PNV layers by removing potentially unsuitable land
# as well areas with an unlikely transition such as intensively cultivated urban land.
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
background <- terra::rast( paste0(path_background, "background_ref.tif") )

# Check results
ll <- list.files(path_output, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

#### Make a combined mask ####

# Background as basis
bb <- background

#Load latest CORINE layer 
#clc <- terra::rast(paste0(path_rawdata, "CORINE/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))
# NOTE: Subset outsourced, load preprocessed layer
clc_mask <- terra::rast("/media/martin/AAB4A0AFB4A08005/CorineMask_2018_clundefined.tif")

# Align with background
if(!terra::compareGeom(bb, clc_mask, stopOnError = FALSE)){
  if(sf::st_crs(bb) != sf::st_crs(clc_mask)) clc_mask <- terra::project(clc_mask, bb)
  # Crop and mask
  clc_mask <- terra::crop(clc_mask, bb)
  clc_mask <- terra::mask(clc_mask, bb)
}

# Any grid cells with more than 50% water bodies or urban surface -> mask out from bb
m <- clc_mask
clc_mask[clc_mask <=5000] <- NA
clc_mask[clc_mask>0] <- 1

bb <- terra::mask(bb, clc_mask, inverse = TRUE)

#### Now remove from all tif files the masked areas
for(f in ll){
  print(basename(f))
  
  ofd <- paste0(dirname(f), "/", "Masked/")
  dir.create(ofd, showWarnings = FALSE)
  ofname <- paste0(ofd, basename(f))
  if(file.exists(ofname)) next()
  
  # Mask all potential natural vegetation layers
  ras <- terra::rast(f)
  ras <- terra::crop(ras, bb)
  ras <- terra::extend(ras, bb)
  # Now mask
  ras <- terra::mask(ras, bb)
  
  # Save output
  writeRaster(ras, ofname)
}
