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
ll <- list.files("ensemble/Masked/",recursive = TRUE,full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
# ll <- ll[grep("Prediction__", ll)]
pnv <- rast(ll)
pnv_sd <- pnv[[which(names(pnv)=="sd")]]
names(pnv_sd) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Ensemble__")
names(pnv_sd) <- stringr::str_remove(names(pnv_sd),"__masked")

pnv <- pnv[[which(names(pnv)=="mean")]]
names(pnv) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Ensemble__")
names(pnv) <- stringr::str_remove(names(pnv),"__masked")

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

# Combine all (no overlaps)
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
library(tidyterra)
library(patchwork)

# Save a copy of the output
out1 <- terra::rast( paste0(path_res, "CurrentPNV_challenge.tif"))
out2 <- terra::rast( paste0(path_res, "CurrentPNV_mostlikely.tif"))
# Reclassify for visualization
out1 <- ibis.iSDM::predictor_transform(out1,option = "perc")

names(out1) <- "Challenge of transitioning to PNV"
names(out2) <- "Most likely PNV"

# Colours
cols <- c("Woodland.and.forest"  = "#117733",
          "Heathland.and.shrub" = "#B65719",
          "Grassland" = "#DDCC77",
          "Sparsely.vegetated.areas" = "#AA4499",
          "Wetlands"  = "#332288",
          "Marine.inlets.and.transitional.waters" = "#88CCEE"
)

assertthat::assert_that(ibis.iSDM::is.Raster(out1),
                        ibis.iSDM::is.Raster(out2))

# Now plot
gm1 <- ggplot() +
  theme_grey(base_size = 20) +
  geom_spatraster(data = out2,maxcell = 1e6) +
  facet_wrap(~lyr, ncol = 2) +
    theme(strip.text.x.top = element_text(size = 20),
          strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = cols,na.value = NA,
                    na.translate = FALSE) +
    guides(fill = guide_legend(title = "")) +
    theme(legend.position = "bottom",legend.text = element_text(22))+
  labs(tag = "a)")
gm1

gm2 <- ggplot() +
  theme_grey(base_size = 20) +
  geom_spatraster(data = out1) +
  facet_wrap(~lyr, ncol = 2) +
    theme(strip.text.x.top = element_text(size = 20),
          strip.background = element_rect(fill = "white")) +
  scale_fill_viridis_c(option = "F",na.value = NA,direction = -1,
                    guide = guide_colorbar(title = "Challenge")
  ) + 
    theme(legend.position = "bottom",legend.title.position = "left",
          legend.title = element_text(vjust = 1),
          legend.key.width  = unit(0.5, 'in') ) +
  # Transparent
  theme( panel.background = element_rect(fill = 'grey90'),
         plot.background = element_rect(fill = "transparent",colour = NA),
         legend.background = element_rect(fill = "transparent")) +
  labs(tag = "b)")
gm2

gm <- gm1 + gm2

ggsave(plot = gm, filename = "figures/Figure3.png", width = 20, height = 10)


# Make a summary for the paper on the most likely one
ex <- as.data.frame(out2)
names(ex) <- 'class'

# Group and summarize
ex |> dplyr::group_by(class) |> dplyr::summarise(n = n()) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(n = n/sum(n))
