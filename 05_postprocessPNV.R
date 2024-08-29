# --- --- #
# This script post-processes the PNV layers by 
# 1) calculating an ensemble among the two Bayesian predictions
# 2) removing potentially unsuitable land as well areas with an unlikely transition 
# such as intensively cultivated urban land.
# --- --- #

# Load in functions and packages 
library(assertthat)
library(terra)
library(sf)
library(ibis.iSDM)
source("00_Functions.R")


path_output = "./PNVHabitats__ClimateRun/"

# -------------------- #
# Get reference raster for the given grain
background <- terra::rast( "data/background_ref.tif") 

# Check results
ll <- list.files(path_output, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("Prediction", ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

# And validation statistics
ll_val <- list.files(path_output, full.names = TRUE, recursive = TRUE)
ll_val <- ll_val[has_extension(ll_val,"rds")]
ll_val <- ll_val[grep("Fold", ll_val)]

# Function to read in validation results
vals <-
  ll_val %>% 
  purrr::map_df(~readRDS(.) |> dplyr::mutate(filename = basename(.)))
# Format filename
vals <- tidyr::separate(vals, filename, c("Valid", "variable", "model", "FoldRepeat"),sep = "__") |> 
  # Drop NA (old validation results)
  tidyr::drop_na(FoldRepeat)

# --------------- #
#### Take each prediction and make an ensemble ####
assertthat::assert_that(exists("ll"))

# Path ensemble
path_ensemble <- "./ensemble"
dir.create(path_ensemble, showWarnings = FALSE)

# Create ensembles
for(vars in unique(vals$variable)){
  ras <- terra::rast( grep(vars,ll,ignore.case = TRUE,value = TRUE) )
  w <- vals |> dplyr::filter(variable == vars) |> 
    # Aggregate per model 
    dplyr::filter(metric == "f1") |> dplyr::group_by(model) |> 
    dplyr::summarise(value = mean(value)) |> dplyr::pull(value)
  
  new <- emptyraster(ras)
  for(n in unique(names(ras))){
    message(vars, n)
    out <- terra::weighted.mean(ras[[which(names(ras)==n)]],w = w)
    names(out) <- n
    new <- c(new, out)
  }
  writeRaster(new, paste0(path_ensemble, "/Ensemble__",vars,".tif"), datatype = "FLT4S")
  rm(new)
}
gc(verbose = TRUE)

#### Make a combined mask ####

# Background as basis
bb <- background

#Load latest CORINE layer 
#clc <- terra::rast(paste0(path_rawdata, "CORINE/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))
# NOTE: Subset outsourced, load preprocessed layer
# clc_mask <- terra::rast("/media/martin/AAB4A0AFB4A08005/CorineMask_2018_clundefined.tif")
clc_mask <- terra::rast("corine/CorineMask_2018_clundefined.tif")

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

# New list using ensemble only
ll <- list.files(path_ensemble,full.names = TRUE)
ll <- ll[has_extension(ll, 'tif')]

# Now remove from all tif files the masked areas
for(f in ll){
  print(basename(f))
  
  ofd <- paste0(dirname(f), "/", "Masked/")
  dir.create(ofd, showWarnings = FALSE)
  ofname <- paste0(ofd, tools::file_path_sans_ext(basename(f)), "__masked.tif")
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
gc()
