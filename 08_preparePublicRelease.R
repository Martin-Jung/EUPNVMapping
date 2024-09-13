# Script to prepare a public release of the created layers to zenodo
# This script assembles the layers and prepares for upload
# To be uploaded:
# - Most likely transition layer
# - A png layer of such
# - The ensemble current potential natural vegetation layers
# - The ensemble future potential natural vegetation layers
# Formats:
# - Cloud-optimized for all current PNV layers without time dimension
# - netcdf and EBVcube for all future ones

library(terra)
library(ncmeta)
library(ebvcube)
library(gdalUtilities)
library(assertthat)
library(ibis.iSDM)
library(tidyterra)
library(stars)

source("00_Functions.R")

# Path to masked current ensemble
path_cPNV <- "ensemble/Masked/"
path_fPNV <- ""
path_mostlikely <- "results/CurrentPNV_mostlikely.tif"

# Get reference raster for the given grain
background <- terra::rast( "data/background_ref.tif")

# Writing path
path_export <- "export/"
dir.create(path_export,showWarnings = FALSE)

assertthat::assert_that(
  is.dir(path_cPNV),
  file.exists(path_mostlikely)
)

#### Export most likely transition and screenshot ####

# Copy layer
# file.copy(from = path_mostlikely,
          # to = paste0(path_export,"pnv_mostlikely_current_laea_1km.tif"),
          # overwrite = TRUE)

# Convert to cloud-optimized geoTIFF
gdalUtilities::gdal_translate(
  src_dataset = path_mostlikely,
  dst_dataset = paste0(path_export,"pnv_mostlikely_current_laea_1km.tif"),
  of = "COG",
  co = matrix(c("TILED=YES", 
           "COPY_SRC_OVERVIEWS=YES", 
           "COMPRESS=DEFLATE"), 
         ncol = 1)
)
assertthat::assert_that(file.exists(paste0(path_export,"pnv_mostlikely_current_laea_1km.tif")))

# Colours
cols <- c("Woodland.and.forest"  = "#117733",
          "Heathland.and.shrub" = "#B65719",
          "Grassland" = "#DDCC77",
          "Sparsely.vegetated.areas" = "#AA4499",
          "Wetlands"  = "#332288",
          "Marine.inlets.and.transitional.waters" = "#88CCEE"
)

ras <- rast(paste0(path_export,"pnv_mostlikely_current_laea_1km.tif"))
assertthat::assert_that(ibis.iSDM::is.Raster(ras))

if(!is.factor(ras)){
  ras <- as.factor(ras)
  levels(ras) <- data.frame(
    ID = seq(1,6),
    MostLikelyPNV = c("Woodland.and.forest","Heathland.and.shrub","Grassland","Sparsely.vegetated.areas","Wetlands","Marine.inlets.and.transitional.waters")
  )
}

# Now plot
gm1 <- ggplot() +
  theme_void() +
  geom_spatraster(data = ras, maxcell = 1e6) +
  theme(strip.text.x.top = element_text(size = 20),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = cols,na.value = NA,
                    na.translate = FALSE) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.position = "bottom",legend.text = element_text(22))
gm1

ggsave(plot = gm1, filename = paste0(path_export, "pnv_mostlikely_current_laea_1km.png"),
       width = 13, height = 13)

#### Copy current PNV as cloud-optimized geotiffs ####
# We only save the median posterior estimates for fast extraction

# Get list of ensemble data
ll <- list.files(path_cPNV, full.names = TRUE)
ll <- ll[grep("Ensemble__", ll)]
ll <- ll[has_extension(ll, "tif")]

# Loop through each and export
for(lyr in ll){
  # Get ofname
  ofname <- stringr::str_remove(tools::file_path_sans_ext(basename(lyr)),"Ensemble__")
  ofname <- stringr::str_remove(ofname,"__masked")
  message(ofname)
  
  # Convert to cloud-optimized geoTIFF
  gdalUtilities::gdal_translate(
    src_dataset = lyr,
    dst_dataset = paste0(path_export,"pnv_",ofname,"_current_laea_1km.tif"),
    of = "COG",
    co = matrix(c("TILED=YES", 
                  "COPY_SRC_OVERVIEWS=YES", 
                  "COMPRESS=DEFLATE"), 
                ncol = 1)
  )
}

# Test
rast(lyr)

#### Create netcdfs for export ####
# Export netcdf files

toint <- F

# Get all proejctions
ll <- list.files(path_cPNV, full.names = TRUE, recursive = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("Projection", ll)]
assertthat::assert_that(length(ll)>0, msg = "No results found...")

# Make a container of all projections
df <- data.frame(ifname = ll, split = tools::file_path_sans_ext(basename(ll))) |> 
  tidyr::separate(split,sep = "__",into = c("Proj", "variable","timeperiod", "ssp", "gcm")) |> 
  dplyr::select(-Proj)

# Load in for each ssp and gcm
for(s in unique(df$ssp)){ # s = df$ssp[1]
  for(g in unique(df$gcm)){ # g = df$gcm[1]
    sub <- df |> dplyr::filter(ssp == {{s}}, gcm == {{g}})
    # Format date
    sub$timeperiod <- as.Date(as.numeric(sub$timeperiod))

    for(vars in unique(sub$variable)){ # vars = "Grassland"
      message(vars, '-', s)
      ofname <- paste0(path_export, "pnv_projection__",vars,"__",s,"__",g,".nc")
      if(file.exists(ofname)) next()
      out <- NULL
      for(metric in c("suitability_q05","suitability_q50",
                      "suitability_q95","suitability_mean")){ # metric = "suitability_mean"
        ras <- terra::rast(sub$ifname[sub$variable==vars])
        ras2 <- ras[[grep(metric, names(ras))]]
        terra::time(ras2) <- sub$timeperiod[sub$variable==vars]
        # Convert to integer
        if(toint) ras2 <- terra::as.int(round(ras2 * 10000))
        z <- lapply(ras2, function(z) ibis.iSDM:::raster_to_stars(z) )
        z <- Reduce("c", z)
        names(z) <- metric
        if(is.null(out)) out <- z else {
          if( all(!stars::st_get_dimension_values(out, 3) != stars::st_get_dimension_values(z, 3))){
            stars::st_dimensions(z) <- stars::st_dimensions(out)
          }
          out <- c(out, z)
        }
        rm(z)
      }
      # Write mdim area
      stars::write_mdim(out, ofname, NA_value = NA)
      rm(out)
      gc()
    }
  }
}
