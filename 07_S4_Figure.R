# --- --- #
# Here we simply plot the relative change (subtraction) in suitability between
# the current and future conditions per ssp.
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
library(terra)
library(sf)
library(ibis.iSDM)
library(tidyterra)
library(ggplot2)

# Colours
cols <- c("Woodland.and.forest"  = "#117733",
          "Heathland.and.shrub" = "#B65719",
          "Grassland" = "#DDCC77",
          "Sparsely.vegetated.areas" = "#AA4499",
          "Wetlands"  = "#332288",
          "Marine.inlets.and.transitional.waters" = "#88CCEE"
)

# Path to Ensemble
path_output <- "ensemble/Masked/"

# Get reference raster for the given grain
background <- terra::rast( "data/background_ref.tif")


#### Get all data ####

ll <- list.files(path_output, full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
ll <- ll[grep("EnsembleProjection",ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

# Make a container of all projections
df <- data.frame(ifname = ll, split = tools::file_path_sans_ext(basename(ll))) |> 
  tidyr::separate(split,sep = "__",into = c("Proj", "variable", "timeperiod", "ssp", "gcm")) |> 
  dplyr::select(-Proj)

# Convert date
df$timeperiod <- as.Date(as.numeric(df$timeperiod))

ll_cont <- list.files(path_output, full.names = TRUE)
ll_cont <- ll_cont[has_extension(ll_cont, "tif")]
ll_cont <- ll_cont[grep("Ensemble__",ll_cont)]
assertthat::assert_that(length(ll_cont)>0, msg = "No results found...")

df_cont <- data.frame(ifname = ll_cont, split = tools::file_path_sans_ext(basename(ll_cont)) ) |> 
  tidyr::separate(split,sep = "__",into = c("Proj", "variable")) |> 
  dplyr::select(-Proj)

# --- #
results <-list()

# Assuming single GCM at the moment
assertthat::assert_that(n_distinct(df$gcm)==1)
# Computed here only for a single GCM given their similarities

# Get mean suitability per ssp and variable
for(vars in unique(df$variable)){
  message(vars)
  for(s in unique(df$ssp)[1]){
    sub <- df |> dplyr::filter(ssp == {{s}}, variable == {{vars}})
    
    # Get future
    ras <- rast(sub$ifname)
    ras <- ras[[grep("suitability_mean", names(ras))]]
    terra::time(ras) <- sub$timeperiod
    # Add current
    cur <- terra::rast(df_cont$ifname[df_cont$variable==vars])
    cur <- subset(cur, "mean")
    names(cur) <- "suitability_mean"
    terra::time(cur) <- as.Date("2010-01-01")
    ras <- c(cur,ras); rm(cur)
    
    out <- terra::regress(ras, 1:nlyr(ras))
    out <- out[["x"]]
    names(out) <- vars
    
    # Save
    writeRaster(out, paste0("results/Regression_",vars,"_",s,".tif"))
    rm(out,ras)
  }
}

#### Plot ####

# Load layers
ll <- list.files("results/",full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
ll <- ll[grep("Regression", ll)]

ras <- terra::rast(ll)
ras <- ras[[match(names(cols),names(ras))]]

# Now plot
gm <- ggplot() +
  theme_gray(base_size = 20) +
  geom_spatraster(data = ras) +
  facet_wrap(~lyr, ncol = 3) +
  scale_fill_gradient2(
    na.value = NA
  ) + 
  guides(fill = guide_legend("Trend"))

ggsave(plot = gm, filename = "figures/SI_4_RegFuturePNV.png", width = 14, height = 10)


