# --- --- #
# This script aggregates some historic land-use trajectories for a first figure in the manuscript.
# Idea is to take the HURTT data and assess for Europe the trajectory of natural land
# Script expects a folder with downloaded HYDE data
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")

# Set run cores for parallization
options('ibis.nthread' = cores)
options('ibis.runparallel' = ifelse(cores>1, FALSE, FALSE) ) # Set to FALSE for slower, but memory proof runs

# Temporary outputs saved for later
path_res <- "results/"
dir.create(path_res, showWarnings = FALSE)

# -------------------- #
# Get reference raster for the given grain (WGS84 projection here)
background <- terra::rast( paste0(path_background, "background_refWGS84.tif") )
bb <- terra::as.polygons(background) |> sf::st_as_sf()

# Path to HURTT data
path_hurt <- "/mnt/pdrive/bec_data/100_rawdata/LUH2_v2h/states.nc"
assert_that(file.exists(path_hurt))

# Load proxy
hu <- stars::read_ncdf(path_hurt, var = c("primf", "primn"),proxy = FALSE)

# Crop to background
hu <- st_crop(hu, bb)



