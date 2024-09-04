# --- --- #
# This script simply plots the uncertainty ranges (difference low to high) for
# the ensemble posterior
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
library(terra)
library(sf)
library(ibis.iSDM)
library(tidyterra)
library(ggplot2)

# Path to the create (depends on system)
path_output <- "corine/"
assertthat::assert_that(dir.exists(path_output))

# Colours
cols <- c("Woodland.and.forest"  = "#117733",
          "Heathland.and.shrub" = "#B65719",
          "Grassland" = "#DDCC77",
          "Sparsely.vegetated.areas" = "#AA4499",
          "Wetlands"  = "#332288",
          "Marine.inlets.and.transitional.waters" = "#88CCEE"
)

# Path to Ensemble
path_ensemble <- "ensemble/Masked/"

# -------------------- #
#### Load data ####
# Get reference raster for the given grain (WGS84 projection here)
background <- terra::rast( "data/background_ref.tif" )
bb <- terra::as.polygons(background) |> sf::st_as_sf()

# Load all the varying predictions
ll <- list.files(path_ensemble, full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "tif")]
assert_that(length(ll)>0)

# Load raster and get only median
ras <- terra::rast(ll)
ras <- ras[[grep("sd",names(ras))]]
names(ras) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Ensemble__")
names(ras) <- stringr::str_remove(names(ras),"__masked")

# Mask with background
ras <- terra::crop(ras, bb)
ras <- terra::mask(ras, bb)

# Now plot
gm <- ggplot() +
  geom_spatraster(data = ras) +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_viridis_c(
    na.value = NA,direction = -1
  ) +
  labs(
    fill = "",
    title = "Current potential natural vegetation (PNV)",
    subtitle = "Standard devation of posterior"
  )
 
ggsave(plot = gm, filename = "figures/SI_3_PNV_sd.png", width = 8, height = 10)
