# Load packages and define paths
library(dataclima)
library(ibis.iSDM)
# devtools::load_all()
library(assertthat)
library(terra)
library(sf)
# library(blockCV)
library(gdalUtilities)
library(spatialsample)
library(readxl)
library(dplyr)

#' Helper function to get species from filename
#' @param input A `character()` for a species name.
#' @export
file_to_speciesname <- function(input){
  assertthat::assert_that(is.character(input))
  return(
    stringr::str_replace(
      tools::file_path_sans_ext( 
        basename(input)
      ),
      pattern = "_",replacement = " "
    )
  )
}

#' Calculate data volume in a list of lists
#' 
#' @param data A [`list`] of observed data. Usually containing presence and absence points
#' @export
count_dataobs <- function(data){
  assertthat::assert_that(length(data)>0,
                          is.list(data))
  out <- vector()
  for(n in names(data)){
    for(k in names(data[[n]])){
      v <- nrow(data[[n]][[k]])
      names(v) <- k
      out <- append(out, v)
    }
  }
  return(out)
}

#' Stars replace NA with 0
#' 
#' @description 
#' This function replaces all NA values in a given stars layer to replace with 0
#' @param obj A [`stars`] object.
#' @param tempalte Another [`stars`] object
#' @keywords internal, utils
#' @noRd
stars_replace_na <- function(obj, template){
  assertthat::assert_that(
    inherits(obj, "stars")
  )
  # Set all 
  obj[is.na(obj)] <- 0
  
  # Get template,
  template[template>=0] <- 1
  template[is.na(template)] <- 0
  
  # Mask the layer by template
  obj[template==0] <- NA
  
  return(obj)  
}

#' Replace NA with 0
#' 
#' @description 
#' This function replaces all NA values in a given Rasterlayer or Rasterstack with the respective
#' @param obj A [`SpatRaster`] object.
#' @keywords internal, utils
#' @noRd
raster_replace_na <- function(obj){
  assertthat::assert_that(
    is.Raster(obj)
  )
  
  for(n in 1:terra:nlyr(obj)){
    if(anyNA(obj[[n]][])){
      obj[[n]][is.na(obj[[n]])] <- 0 
    }
  }
  return(obj)  
}

#' Purpose of this script is to create blocked subsets of presence only data
#' Mainly relying on the blockCV package and doing a systematic division of sites
#' @param species A sf object with a Species column
#' @param speciesCol The column with the species name in the dataset
#' @param template A [`SpatRaster`] containing the background modelling extent.
#' @param k How many folds (default 3).
#' @param iteration Number of iterations for stability of classes
#' @param selection How blocks are selected (Default: Random).
#' @param return_sf A [`logical`] indication of whether the spatialfolds should be returned instead
#' @import blockCV
#' @import sf
spatialBlock <- function(species, speciesCol = 'Observed', k = 3, iteration = 10,
                         template = NULL,
                         selection = "random", return_sf = FALSE, ...){
  assertthat::assert_that(has_name(species,speciesCol),
                          is.Raster(template),
                          is.logical(return_sf),
                          k > 1)
  # For some reason the blocking does not work with projections other WGS84, thus reproject for this purpose
  if(sf::st_crs(species) != sf::st_crs(template)){
    species <- species %>%  sf::st_transform(crs = sf::st_crs(template))
  }
  
  # spatial blocking by specified range and random assignment
  sb <- try({
    blockCV::spatialBlock(
      speciesData = species, # sf or SpatialPoints
      species = speciesCol, # the response column (binomial or multi-class)
      rasterLayer = template,
      k = k, # number of folds
      selection = selection,
      iteration = iteration, # find evenly dispersed folds
      rows = 10, # For checkerboard
      cols = 10,
      verbose = FALSE,
      biomod2Format = TRUE,
      progress = TRUE,
      showBlocks = FALSE
    )
  })
  if(inherits(sb, "try-error")){
    sb <- try({
      blockCV::spatialBlock(
        speciesData = species, # sf or SpatialPoints
        species = speciesCol, # the response column (binomial or multi-class)
        k = k-1, # number of folds
        selection = selection,
        iteration = iteration, # find evenly dispersed folds
        rows = 10, # For checkerboard
        cols = 10,
        verbose = FALSE,
        biomod2Format = TRUE,
        progress = TRUE
      )
    })
    k <- k - 1
  }
  # Construct output data.frame
  out <- as.data.frame(sb$biomodTable) %>% 
    dplyr::rename_with(~ paste0('fold',1:k), paste0('RUN',1:k))
  if(nrow(species)!=nrow(out)){
    out2 <- as.data.frame(as.matrix(t( rep(sample(c('FALSE','TRUE'), ncol(out), replace = TRUE),nrow(species) - nrow(out)) ),ncol = 3))
    names(out2) <- names(out)
    out <- rbind(out, out2)
  }
  out[out==TRUE] <- 'Train';out[out==FALSE] <- 'Test'
  if(return_sf){
    o <- sf::st_as_sf(sb$blocks)["folds"]
    return( o )
  } else {
    out
  }
}


#' Spatial blocking using spatialsample
spatialBlock2 <- function(df, template, k = 3, method = "random"){
  assertthat::assert_that(is.Raster(template),
                          k > 1)
  # Spatial random blocking
  spcv <- spatialsample:::spatial_block_cv(df, method = method, v = k, relevant_only = TRUE)
  if(nrow(spcv)<k) stop("Something went wrong with the crossvalidation!")
  
  out <- data.frame()
  for(id in 1:nrow(spcv)){
    # Get points
    df_train <- df[spcv$splits[[id]]$in_id,]
    
    # Get Cells from template for those points
    copy <- emptyraster(template)
    copy[terra::cellFromXY(template, sf::st_coordinates(df_train))] <- 1
    
    new <- terra::as.polygons(copy) |> sf::st_as_sf() |>
      sf::st_union() |> sf::st_as_sf()
    new$folds <- spcv$id[[id]]

    out <- rbind(out,new)
  }
  
  return(out)
}

#' Get a mask of non-NA grid cells across the full stack
#' 
#' @description 
#' In order to add additional predictors to a GLOBIOM stack, we need to ensure that there
#' aren't values in areas that are missing in others. Since we are constrained to GLOBIOM
#' as common reference, we will use all non-data in the supplied stack to create a mask for stacking.
#' @param stk A [`SpatRaster`] object with multiple layers
#' @export
createNonNAmask <- function(stk, template){
  assertthat::assert_that(
    terra::nlyr(stk)>2
  )
 # Check ibis is there
 require(ibis.iSDM)
  
 mat <- as.data.frame(stk) 
 out <- emptyraster(stk) 
 out[] <- apply(mat, 1, function(x) !anyNA(x))
 # Mask template
 out <- terra::mask(out, template)
 out[out==0] <- NA
 # Return mask
 return(out)
}

#' Plot the NA data of a given rasterLayer
#' 
#' @description 
#' Small helper function to plot the NA data for a given raster
#' @param layer A [`SpatRaster`] object with a single layer
#' @export
plotNAdata <- function(layer, template){
  assertthat::assert_that(
    terra::nlyr(layer)==1,
    is.Raster(layer)
  )
  
  mat <- as.data.frame(layer) 
  out <- ibis.iSDM:::emptyraster(layer) 
  out[] <- apply(mat, 1, function(x) !anyNA(x))
  # Mask template
  out <- terra::mask(out, template)
  out <- terra::as.factor(out)
  plot(out, col = c("red", "lightgreen"))
}
