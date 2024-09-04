# --- --- #
# This script makes a simple plot of the predicted estimates for the
# supplementary materials
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
library(terra)
library(sf)
library(ibis.iSDM)
library(tidyterra)
library(ggplot2)

# Path to the create (depends on system)
# path_output <- "PNVHabitats__ClimateRun/"
path_output <- "ensemble/Masked/"
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
# ll <- ll[grep("Prediction_",ll)]
# ll <- ll[grep("Threshold_",ll)]
assert_that(length(ll)>0)

# Load raster and get only median
ras_med <- rast(ll)
ras_med <- ras_med[[which(names(ras_med)=="q50")]]
# ras_med <- ras_med[[which(names(ras_med)=="threshold_q50_percentile")]]
# names(ras_med) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")
names(ras_med) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Ensemble__")
names(ras_med) <- stringr::str_remove(names(ras_med),"__masked")

# Reorder by cols
ras_med <- ras_med[[match(names(cols),names(ras_med))]]

# Now plot
gm <- ggplot() +
  geom_spatraster(data = ras_med) +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(
    palette = "atlas",
    # n.breaks = 12,
    breaks = seq(0,1,length.out = 11),
    direction = -1,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    fill = "",
    title = "Current potential natural vegetation (PNV)",
    subtitle = "Median predictions"
  )

ggsave(plot = gm, filename = "figures/SI_2__CurrentPNV.png", width = 8, height = 10)
