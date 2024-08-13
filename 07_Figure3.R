# --- 
# The purpose of this figure is to present a map showing both 
# - a) the challenge of moving from current to potential vegetation
# - b) the most likely transition from Current land systems to potential vegetation.
# ---
library(terra)
library(ggplot2)
library(sf)
library(tidyterra)
library(readxl)
library(assertthat)

source("00_Functions.R")

# Temporary result folder
path_res <- "results/"
dir.create(path_res,showWarnings = FALSE)

# Load Dou et al. raster
dou <- terra::rast("DouEtAl/EU_landSystem.tif")

# Get reference raster for the given grain
background <- terra::rast( "data/background_ref.tif")

# The SI table 1 with transition challenges
dou_trans <- read_xlsx("C:/Users/tuete/OneDrive - IIASA/Manuscripts/EUPNV/SITable1.xlsx",col_types = "guess")
dou_trans$pnv_forest <- as.numeric(dou_trans$pnv_forest)
dou_trans$pnv_heath <- as.numeric(dou_trans$pnv_heath)
dou_trans$pnv_grass <- as.numeric(dou_trans$pnv_grass)
dou_trans$pnv_sparse <- as.numeric(dou_trans$pnv_sparse)
dou_trans$pnv_wetland <- as.numeric(dou_trans$pnv_wetland)
dou_trans$pnv_marine <- as.numeric(dou_trans$pnv_marine)

# Load stack of current potential natural vegetation
ll <- list.files("PNVHabitats__ClimateRun/",recursive = TRUE,full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
ll <- ll[grep("Prediction__", ll)]
pnv <- rast(ll)
pnv_sd <- pnv[[which(names(pnv)=="sd")]]
names(pnv_sd) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")
pnv <- pnv[[which(names(pnv)=="mean")]]
names(pnv) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")

# Checks
assertthat::assert_that(ibis.iSDM::is.Raster(pnv),
                        ibis.iSDM::is.Raster(pnv_sd),
                        nrow(dou_trans)>0,
                        all(dou_trans$code %in% unique(dou)[,1])
                        )
# ---- #
#### 1. Align all layers to the background ####
# Clip and mask
dou <- terra::crop(dou, terra::ext(background))
dou <- terra::mask(dou, background)
dou <- terra::project(dou, background)

pnv <- terra::extend(pnv, background)
pnv <- terra::crop(pnv, terra::ext(background))
pnv <- terra::mask(pnv, background)
pnv <- terra::project(pnv, background)

pnv_sd <- terra::extend(pnv_sd, background)
pnv_sd <- terra::crop(pnv_sd, terra::ext(background))
pnv_sd <- terra::mask(pnv_sd, background)
pnv_sd <- terra::project(pnv_sd, background)

#### 2. Calculate transition challenge to natural vegetation ####
# This is calculated as follows:
# > min( average score / (probability estimate))

out1 <- terra::rast() 
for(cl in dou_trans$code){
  message("Processing ", cl)
  sub <- dou_trans |> dplyr::filter(code == cl)
  
  # Get the target class
  ras <- ibis.iSDM::emptyraster(dou)
  ras[] <- NA
  ras[dou == cl] <- 1

  # Mask and multiply the sub
  pp <- terra::mask(pnv, ras)
  # pp <- (1-pp) # Invert since to align direction of score and probability
  pp <- terra::app(
    c(
      1/pp$Grassland * sub$pnv_grass,
      1/pp$Heathland.and.shrub * sub$pnv_heath,
      1/pp$Marine.inlets.and.transitional.waters * sub$pnv_marine,
      1/pp$Sparsely.vegetated.areas * sub$pnv_sparse,
      1/pp$Wetlands * sub$pnv_wetland,
      1/pp$Woodland.and.forest * sub$pnv_forest
    ),"min"
  )
  names(pp) <- paste0('code_',cl)
  # Now add to out
  suppressWarnings( out1 <- c(out1, pp))
  rm(pp, ras)
}
gc() # Cleanup

# Combine all (there are no overlaps anyway)
out1 <- mean(out1, na.rm = TRUE)
out1 <- terra::mask(out1, background)

# Save a copy of the output
writeRaster(out1, paste0(path_res, "CurrentPNV_challenge.tif"),overwrite=TRUE)

#### 2. Calculate most likely natural vegetation transition####
# This is calculated as follows:
# > which.min( average score x probability estimate)

out2 <- terra::rast() 
for(cl in dou_trans$code){
  message("Processing ", cl)
  sub <- dou_trans |> dplyr::filter(code == cl)
  
  # Get the target class
  ras <- ibis.iSDM::emptyraster(dou)
  ras[] <- NA
  ras[dou == cl] <- 1
  
  # Mask and multiply the sub
  pp <- terra::mask(pnv, ras)
  # pp <- (1-pp) # Invert since to align direction of score and probability
  pp <- terra::which.min(
    c(
      1/pp$Woodland.and.forest * sub$pnv_forest,
      1/pp$Heathland.and.shrub * sub$pnv_heath,
      1/pp$Grassland * sub$pnv_grass,
      1/pp$Sparsely.vegetated.areas * sub$pnv_sparse,
      1/pp$Wetlands * sub$pnv_wetland,
      1/pp$Marine.inlets.and.transitional.waters * sub$pnv_marine
    )
  )
  names(pp) <- paste0('code_',cl)
  # Now add to out
  suppressWarnings( out2 <- c(out2, pp))
  rm(pp, ras)
}
gc() # Cleanup

# Combine all (there are no overlaps anyway)
out2 <- terra::modal(out2, na.rm = TRUE)
out2 <- terra::mask(out2, background)
names(out2) <- "MostLikelyPNV"
# Classify
out2 <- as.factor(out2)
levels(out2) <- data.frame(
  ID = seq(1,6),
  MostLikelyPNV = c("Woodland.and.forest","Heathland.and.shrub","Grassland","Sparsely.vegetated.areas","Wetlands","Marine.inlets.and.transitional.waters")
)

# Save a copy of the output
writeRaster(out2, paste0(path_res, "CurrentPNV_mostlikely.tif"), overwrite=TRUE)

#### 3. Build plot ####
# Simply visualize both maps side by side

# Save a copy of the output
out1 <- terra::rast( paste0(path_res, "CurrentPNV_challenge.tif"))
out2 <- terra::rast( paste0(path_res, "CurrentPNV_mostlikely.tif"))

assertthat::assert_that(ibis.iSDM::is.Raster(out1),
                        ibis.iSDM::is.Raster(out2))

#TODO: Built a plot
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

ggsave(plot = gm, filename = "figures/Figure3.png", width = 15, height = 10)
