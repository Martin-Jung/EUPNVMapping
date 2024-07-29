# Small script to generate Figure S1
# Purpose:
# - Showing the sampling distribution of each of the different datasets

# Load plotting packages
library(ggplot2)
library(scales)
library(forcats)
library(ggthemes)
library(ggsci)
library(patchwork)
library(sf)

# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")

# --------------------- #

# Get reference raster for the given grain
background <- terra::rast( paste0(path_background, "background_ref.tif") )
# Aggregate to speed up and ease visual interpretation of coverage
background <- terra::aggregate(background, fact = 10, fun= terra::modal)

# Get all habitat databases
habitatdbs <- list.files(path = path_habitat, full.names = TRUE, recursive = TRUE)
habitatdbs <- habitatdbs[has_extension(habitatdbs, "gpkg")]
habitatdbs <- habitatdbs[grep("Cropland|Urban|River|Marine.gpkg", habitatdbs,invert = TRUE)]
assertthat::assert_that(length(habitatdbs)>0)

#### Load and simulate all data ####

out <- NULL
for(hab in sort(habitatdbs)){ # hab <- habitatdbs[5] 
  
  # Habitat name
  sname <- tools::file_path_sans_ext(basename(hab))
  
  message("Processing ", sname)  
  
  # Loading data for the target habitat
  lyrs <- sf::st_layers(hab)
  
  # --- #
  # Now load all cleaned data as list
  data <- list()
  for(name in lyrs$name) {
    o <- sf::st_read(hab, layer = name, quiet = TRUE) |> dataclima:::laea_transform() |> 
      ibis.iSDM:::rename_geometry("geometry") 
    if(!has_name(o, "observed")){
      # For polygon:
      if(any(unique(sf::st_geometry_type(o)) == "POLYGON")){
        o <- sf::st_sample(o, nrow(o), type = "regular") %>% sf::st_as_sf() %>% 
          ibis.iSDM:::rename_geometry("geometry") %>% 
          dplyr::mutate(observed = 1, dataset = "Natura2000")
      }
      if(has_name(o, "population_low")){
        o <- dplyr::rename(observed = population_low) %>% dplyr::mutate(observed = if_else(observed > 0, 1, 0))
      } else { o <- o %>% dplyr::mutate(observed = 1) }
    }
    data[[name]] <- o
  }
  assert_that(length(data) > 0)
  # --- #
  # For PNV we simply take the all point information. These will be concatenated later
  df <- data.frame()
  for(l in 1:length(data)){
    df <- rbind(df,
                data[[l]] %>% dplyr::select(observed) %>% 
                  dplyr::mutate(type = lyrs$name[l]) |> 
                  sf::st_transform(crs = sf::st_crs(background))
    )
  }
  if(nrow(df)==0) next() # return(NULL) # No point data found...
  df <- ibis.iSDM:::guess_sf(df) |> ibis.iSDM:::rename_geometry("geometry")
  
  # Aggregate observations to grid as to remove duplicates
  df <- ibis.iSDM:::aggregate_observations2grid(df, template = background, field_occurrence = "observed")
  if(is.na(sf::st_crs(df))) df <- sf::st_set_crs(df, value = sf::st_crs(background))
  df <- df |> sf::st_transform(df, crs = sf::st_crs(background))
  
  # Add name
  df$name <- sname
  
  message("Collated ",sname," (", sum(df$observed), " points)")
  
  if(is.null(out)) out <- df else out <- dplyr::bind_rows(out, df)
  rm(df)
}

#### Make a coverage figure ####

# Intersect with background again
ex <- terra::extract(background, out)
out <- out[-which(is.na(ex[,2])),]
rm(ex)

# Recategorizing among into 7 percentile
# out$cuts <- cut_interval(out$observed, 6)
out$cuts <- cut(out$observed, c(0,1, 10, 50, 100, max(out$observed)) ,
                labels = c("1", "1-10","10-50", "50-100", ">100"))
assertthat::assert_that(!anyNA(out$cuts))

# Simply plot all per type, showing a gradient per amount
g <- ggplot() + 
  # Add background
  geom_sf(data = sf::st_as_sf(as.polygons(background))) +
  # Add points
  geom_sf(data = out, aes(colour = cuts), alpha = 0.6) +
  # geom_point() +
  facet_wrap(~name) + 
  scale_color_viridis_d(option = "F",direction = -1)
g
ofname <- paste0(path_figures, "SI_1__DataCoverage.png")
ggsave(plot = g, filename = ofname, width = 14, height = 10)
