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

# Path ensemble
path_ensemble <- "./ensemble"
dir.create(path_ensemble, showWarnings = FALSE)

path_output = paste0("/mnt/hdrive/","PNVHabitats__ClimateRun/")

# Adjust directory
terra::terraOptions(tempdir = "/media/martin/AAB4A0AFB4A08005/tmp/")

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
clc_mask <- terra::rast("/media/martin/AAB4A0AFB4A08005/CorineMask_2018_clundefined.tif")
# clc_mask <- terra::rast("corine/CorineMask_2018_clundefined.tif")

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

# --- #
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

# --------------- #
#### Ensemble projections #####
# Get all projections for ensemble calculations
ll <- list.files(path_output, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("Projection", ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

# Make a container of all projections
df <- data.frame(ifname = ll, split = tools::file_path_sans_ext(basename(ll))) |> 
  tidyr::separate(split,sep = "__",into = c("Proj", "variable", "model", "timeperiod", "ssp", "gcm")) |> 
  dplyr::select(-Proj)

# Apply mask directly?
apply_mask <- TRUE
if(apply_mask) assertthat::assert_that(exists("bb"))

# Create ensembles
for(vars in unique(df$variable)){
  for(g in unique(df$gcm)){
    for(s in unique(df$ssp)){
      for(p in unique(df$timeperiod)){
        
        sub <- df |> dplyr::filter(variable == {{vars}},
                                   gcm == {{g}}, ssp == {{s}}, 
                                   timeperiod == {{p}})
        # Get the weights too
        w <- vals |> dplyr::filter(variable == {{vars}}) |> 
          # Aggregate per model 
          dplyr::filter(metric == "f1") |> dplyr::group_by(model) |> 
          dplyr::summarise(value = mean(value))
        if(!all(sub$model == w$model)){ message("skipped"); next()} 
        assertthat::assert_that(all(sub$model == w$model))
        
        # Load both rasters and average
        ras <- terra::rast(sub$ifname)
        # If raster lacks time, add here.
        if(is.na(terra::time(ras)[1])) terra::time(ras) <- rep(as.Date(as.numeric(p)), terra::nlyr(ras))
        
        for(n in unique(names(ras)) ){
          if(n %in% c("suitability_q05","suitability_q50","suitability_q95","suitability_mean")){
            # Now create folder structure and save
            # Scenario / metric /entity / time
            # Build output folder
            ofg <- paste0(path_ensemble, "/", 
                          paste0(s,"__",g), # Scenario
                          "/",
                          paste0(n), # Metric
                          "/",
                          paste0(tolower(vars)) # Entity
            )
            dir.create(ofg, showWarnings = FALSE,recursive = TRUE)
            # Build output name
            ofname <- paste0(ofg, "/", lubridate::year(terra::time(ras)[1]), ".tif")
            if(file.exists(ofname)) { next() }
            
            message(vars, "  -  " , n)
            out <- terra::weighted.mean(ras[[which(names(ras)==n)]],w = w$value)
            names(out) <- n
            if(apply_mask){
              # Crop
              out <- terra::crop(out, bb)
              out <- terra::extend(out, bb)
              # Now mask
              out <- terra::mask(out, bb)
            }
            if(is.na(terra::time(out))){
              terra::time(out) <- terra::time(ras)[1]
            } 
            terra::writeRaster(out, ofname, datatype = "FLT4S",overwrite = TRUE)
          }
        }
      }
    }
  }
}

# --- #
# List all items in the directory and filter to include only folders
folders <- list.dirs(path_ensemble, full.names = TRUE, recursive = FALSE)

# Loop through each folder and create a zip file
for (folder in folders) {
  # Get the folder name
  folder_name <- basename(folder)
  
  # Define the output zip file path
  zip_file <- file.path(path_ensemble, paste0(folder_name, ".zip"))
  
  # Create the zip file
  zip(zipfile = zip_file, files = list.files(folder, full.names = TRUE),)
  
  cat("Zipped folder:", folder, "->", zip_file, "\n")
}

# ---------------- #
#### Calculate ensemble projections per SSP ####
# Here we simply collate all layers and calculate an ensemble per ssp
# averaging across GCMs

# Get all projections for ensemble calculations
ll <- list.files(path_output, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("Projection", ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

# Make a container of all projections
df <- data.frame(ifname = ll, split = tools::file_path_sans_ext(basename(ll))) |> 
  tidyr::separate(split,sep = "__",into = c("Proj", "variable", "model", "timeperiod", "ssp", "gcm")) |> 
  dplyr::select(-Proj)

# Apply mask directly?
apply_mask <- TRUE
if(apply_mask) assertthat::assert_that(exists("bb"))

# Create ensembles
for(vars in unique(df$variable)){
    for(s in unique(df$ssp)){
      for(p in unique(df$timeperiod)){
        
        sub <- df |> dplyr::filter(variable == {{vars}},ssp == {{s}}, 
                                   timeperiod == {{p}})
        # Get the weights too
        w <- vals |> dplyr::filter(variable == {{vars}}) |> 
          # Aggregate per model 
          dplyr::filter(metric == "f1") |> dplyr::group_by(model) |> 
          dplyr::summarise(value = mean(value))
        sub <- dplyr::left_join(sub, w)

        # Load both rasters and average
        ras <- terra::rast(sub$ifname)
        # If raster lacks time, add here.
        if(is.na(terra::time(ras)[1])) terra::time(ras) <- rep(as.Date(as.numeric(p)), terra::nlyr(ras))
        
        for(n in unique(names(ras)) ){
          if(n %in% c("suitability_q05","suitability_q50","suitability_q95","suitability_mean")){
            # Now create folder structure and save
            # Scenario / metric /entity / time
            # Build output folder
            ofg <- paste0(path_ensemble, "/", 
                          paste0(s), # Scenario
                          "/",
                          paste0(n), # Metric
                          "/",
                          paste0(tolower(vars)) # Entity
            )
            dir.create(ofg, showWarnings = FALSE,recursive = TRUE)
            # Build output name
            ofname <- paste0(ofg, "/", lubridate::year(terra::time(ras)[1]), ".tif")
            if(file.exists(ofname)) { next() }
            
            message(vars, "  -  " , n)
            out <- terra::weighted.mean(ras[[which(names(ras)==n)]],w = sub$value)
            names(out) <- n
            if(apply_mask){
              # Crop
              out <- terra::crop(out, bb)
              out <- terra::extend(out, bb)
              # Now mask
              out <- terra::mask(out, bb)
            }
            if(is.na(terra::time(out))){
              terra::time(out) <- terra::time(ras)[1]
            } 
            terra::writeRaster(out, ofname, datatype = "FLT4S",overwrite = TRUE)
          }
        }
      }
    }
}

# --- #
# Also create/copy the starting year in there calculated as ensemble.
# Check results
ll <- list.files(path_output, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("Prediction", ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

ll_ex <- list.dirs(path_ensemble, full.names = TRUE, recursive = TRUE)
assertthat::assert_that(length(ll_ex)>0)

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

# Create ensembles
for(s in c("ssp126", "ssp370", "ssp585")){
  for(vars in unique(vals$variable)){
    ras <- terra::rast( grep(vars,ll,ignore.case = TRUE,value = TRUE) )
    w <- vals |> dplyr::filter(variable == vars) |> 
      # Aggregate per model 
      dplyr::filter(metric == "f1") |> dplyr::group_by(model) |> 
      dplyr::summarise(value = mean(value)) |> dplyr::pull(value)
    
    # Subset to relevant entries and rename
    ras <- ras[[grep(pattern = "mean|q05|q50|q95", names(ras))]]
    for(n in unique(names(ras))){
      message(vars, "-", n)
      out <- terra::weighted.mean(ras[[which(names(ras)==n)]],w = w)
      names(out) <- switch (n,
        "mean" = "suitability_mean",
        "q05" = "suitability_q05",
        "q50" = "suitability_q50",
        "q95" = "suitability_q95"
      )
      # Now write to each folder a copy
      ofg <- grep(names(out),ll_ex, value = TRUE)
      ofg <- grep(tolower(make.names(vars)), ofg, value = TRUE)
      assertthat::assert_that(length(ofg)==3)
      
      for(k in ofg){
        writeRaster(out, paste0(k, "/2020.tif"), datatype = "FLT4S", overwrite = TRUE)
      }
    }
  }
}
gc(verbose = TRUE)

