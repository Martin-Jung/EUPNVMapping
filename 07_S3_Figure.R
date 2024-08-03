# --- --- #
# This script simply plots the aggregated corine layers for the supplementary
# materials
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

# -------------------- #
#### Load data ####
# Get reference raster for the given grain (WGS84 projection here)
background <- terra::rast( "data/background_ref.tif" )
bb <- terra::as.polygons(background) |> sf::st_as_sf()

# Load all the varying predictions
ll <- list.files(path_output, full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "tif")]
assert_that(length(ll)>0)

# Load raster and get only median
ras <- rast(ll)
names(ras) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")

# Mask with background
ras <- terra::crop(ras, bb)
ras <- terra::mask(ras, bb)

# Convert to fractions
ras <- ras / 10000

# TODO:

# Reorder by cols
# ras_med <- ras_med[[match(names(cols),names(ras_med))]]
# 
# # Now plot
# gm <- ggplot() +
#   geom_spatraster(data = ras_med) +
#   facet_wrap(~lyr, ncol = 2) +
#   scale_fill_whitebox_c(
#     palette = "atlas",
#     # n.breaks = 12,
#     breaks = seq(0,1,length.out = 11),
#     direction = -1,
#     guide = guide_legend(reverse = TRUE)
#   ) +
#   labs(
#     fill = "",
#     title = "Current potential natural vegetation (PNV)",
#     subtitle = "Median predictions"
#   )
# 
# ggsave(plot = gm, filename = "figures/SI_2__CurrentPNV.png", width = 8, height = 10)
