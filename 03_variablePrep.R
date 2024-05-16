library(terra)
library(stars)

# We are using primarily the data already prepared through dataclima, 
# however add a range of variables related to future climate

# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")

# path_habitat <- "/mnt/pdrive/bec_data/200_processeddata/dataclima/habitat_occurrence/"
# Path variables to input data. This one is expected to point to a 

# Set run cores for parallization
options('ibis.nthread' = cores)
options('ibis.runparallel' = ifelse(cores>1, FALSE, FALSE) ) # Set to FALSE for slower, but memory proof runs

#### New background layer for alignment ####

if(!file.exists("data/background_ref.tif")){
  # Make new background layer from previous extent
  ll <- list.files(path_presentcovs, full.names = TRUE)
  background <- terra::rast(ll[10])
  # Set to 1
  background[!is.na(background)] <- 1
  
  # Near resample to new reference grid
  refs <- terra::rast( paste0(path_background, "ReferenceGrid_Europe_bin_1000m.tif") )
  refs <- terra::crop(refs, background)
  background <- terra::resample(background, refs)
  names(background) <- "background"
  
  # Remove irrelevant countries, in particular Turkey as not all variables extend there anyway.
  # Clip with global GADM in this case (excluding Turkey in the process)
  nuts <- sf::st_read(path_nuts, layer = "EU_28plus", quiet = TRUE) %>% dataclima:::laea_transform()
    # dplyr::filter(CNTR_CODE != "TR") #%>% dplyr::filter(CNTR_CODE != "NO", CNTR_CODE != "IS")
  background <- terra::crop(background, nuts)
  background <- terra::mask(background, nuts)
  
  # Save the output
  terra::writeRaster(background, "data/background_ref.tif")
  # Also project to WGS84 for later
  background_wgs84 <- terra::project(background, terra::crs("WGS84") )
  terra::writeRaster(background_wgs84, "data/background_refWGS84.tif")
}

# Load backgrounds here
background <- terra::rast("data/background_ref.tif")
# Also project to WGS84 for later
background_wgs84 <- terra::rast( "data/background_refWGS84.tif")

# Path output for covariates
path_covs <- paste0("/media/martin/AAB4A0AFB4A08005/PNV_covs")
dir.create(path_covs,showWarnings = FALSE)

# ------------------ #
#### Alignment of all current and future climate variables ####

# Using CHELSA climatologies here #
# Contemporary first
path_chelsa <- "/mnt/pdrive/watxene/CHELSA/Climatologies/1981-2010"
# Get all files here
ll <- list.files(path_chelsa,full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "tif")]
# Get only for bio|pr|tas|tasmax|tasmin as those are the only ones projected
ll <- ll[basename(dirname(ll)) %in% c("bio")]

# Contemporary first
ll_contemp <- ll
# Loop through each file, crop, mask and then project
for(f in ll_contemp){
  
  timeperiod <- basename(dirname(dirname(f)))
  filename <- basename(f)
  
  message(timeperiod, " - ",filename)
  
  ofname <- paste0(path_covs, "/", filename)
  if(file.exists(ofname)) next()
  
  ras <- terra::rast(f)
  # Crop
  ras <- terra::crop(ras, background_wgs84)
  writeRaster(ras, ofname)
  ras <- terra::rast(ofname) # Load again
  
  if(!compareGeom(ras, background_wgs84, stopOnError=F)){
    if(length(grep("kg0|kg1|kg2|kg3|kg4|kg5",basename(f)))>0){
      ras <- terra::resample(ras, background_wgs84, method = "near")
    } else {
      ras <- terra::resample(ras, background_wgs84, method = "bilinear")
    }
  }
  # Mask
  ras <- terra::mask(ras, background_wgs84)
  
  # Project
  ras <- terra::project(ras, terra::crs(background))
  # plot(ras)
  # Save output
  writeRaster(ras, ofname,overwrite = TRUE)
}

# Same for future layers
path_chelsa <- c("/mnt/pdrive/watxene/CHELSA/Climatologies/2011-2040/",
                 "/mnt/pdrive/watxene/CHELSA/Climatologies/2041-2070/",
                 "/mnt/pdrive/watxene/CHELSA/Climatologies/2071-2100/")
# Get all files here
ll <- list.files(path_chelsa,full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "tif")]
# Get only for bio|pr|tas|tasmax|tasmin as those are the only ones projected
ll <- ll[basename(dirname(ll)) %in% c("bio")]

# Loop through each file, crop, mask and then project
for(f in ll){
  
  ssp <- basename(dirname(dirname(f)))
  gcm <- basename(dirname(dirname(dirname(f))))
  filename <- basename(f)
  
  message(ssp, " - ",gcm, " - ", filename)
  
  od <- paste0(path_covs,"/",ssp, "/", gcm)
  dir.create(od,recursive = TRUE,showWarnings = FALSE)
  ofname <- paste0(od, "/", filename)
  if(file.exists(ofname)) next()
  
  ras <- terra::rast(f)
  # Crop
  ras <- terra::crop(ras, background_wgs84)
  writeRaster(ras, ofname)
  ras <- terra::rast(ofname) # Load again
  
  if(!compareGeom(ras, background_wgs84, stopOnError=F)){
    if(length(grep("kg0|kg1|kg2|kg3|kg4|kg5",basename(f)))>0){
      ras <- terra::resample(ras, background_wgs84, method = "near")
    } else {
      ras <- terra::resample(ras, background_wgs84, method = "bilinear")
    }
  }
  # Mask
  ras <- terra::mask(ras, background_wgs84)
  
  # Project
  ras <- terra::project(ras, terra::crs(background))
  # plot(ras)
  # Save output
  writeRaster(ras, ofname,overwrite = TRUE)
}
stop("Done")

#### PNV other variables ###
# Make a folder for all the other variables useful or related to PN
# Place it also on the data drive

# Everything DEM and Distance coast
paths <- c("/mnt/pdrive/bec_data/200_processeddata/dataclima/predictors/1000",
           "/mnt/pdrive/bec_data/200_processeddata/dataclima/pnv_baselines/1000"
                 )
# Get all files here
ll <- list.files(paths,full.names = TRUE,recursive = TRUE)
ll <- ll[grep("Distance|EUDEM|Lithology|PNV|EuroVegMap|EU_ecologicalRegions|Wetland", ll)]
ll <- ll[has_extension(ll, "tif")]

# Now align with the new background and save them all
for(f in ll){
  print(basename(f))
  od <- paste0(path_covs,"/","OtherVariables")
  dir.create(od,recursive = TRUE,showWarnings = FALSE)
  ofname <- paste0(od, "/", basename(f))
  if(file.exists(ofname)) next()
  
  ras <- terra::rast(f)
  
  # Crop
  ras <- terra::crop(ras, background)
  writeRaster(ras, ofname)
  ras <- terra::rast(ofname) # Load again
  
  if(!compareGeom(ras, background, stopOnError=F)){
    if(length(grep("kg0|kg1|kg2|kg3|kg4|kg5",basename(f)))>0){
      ras <- terra::resample(ras, background, method = "near")
    } else {
      ras <- terra::resample(ras, background, method = "bilinear")
    }
  }
  # Set ranges to 0 and then mask again
  ras[is.na(ras)] <- 0
  
  # Mask
  ras <- terra::mask(ras, background)
  
  # Project
  ras <- terra::project(ras, terra::crs(background))
  # plot(ras)
  # Save output
  writeRaster(ras, ofname,overwrite = TRUE)
  
}
