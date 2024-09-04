# --- --- #
# This script creates the map overview figure based on current potential vegetation
# The overview is plotted through a hexgrid visualized across europe.
# Inspiration and credit to  https://labs.ala.org.au/posts/2024-01-25_hex_point_maps
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)
library(terra)
library(exactextractr)
normalize <- function(x) return( (x - min(x)) / (max(x) - min(x)) )

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
#### Make and load hexgrid ####
# Get reference raster for the given grain (WGS84 projection here)
background <- terra::rast( "data/background_refWGS84.tif" )
# Crop put canaries for now
bbox <- ext(background) |> as.vector()
bbox['xmin'] <- -17; bbox['xmax'] <- 36; bbox['ymin'] <- 34
background <- terra::crop(background, bbox )
bb <- terra::as.polygons(background) |> sf::st_as_sf()
background <- terra::rast( "data/background_ref.tif" )
bb <- bb |> sf::st_transform(background, crs = sf::st_crs(background))

# Make a global grid 
hex_grid <- st_make_grid(bb,
                         # cellsize = 1, # Half a degree
                         cellsize = 200000, # Resoluton of hex grid
                         what = "polygons",
                         square = FALSE,
                         flat_topped = TRUE) |> 
  st_as_sf() |> 
  st_set_geometry("hex_geometry") |> 
  tibble::rowid_to_column(var = "hex_id")

# Remove empty grid
hex_grid <- st_join(x = hex_grid,
               y = bb,
               join = st_intersects,
               left = FALSE) |> 
  dplyr::select(-background)

# Short plot
ggplot() +
  geom_sf(data = bb,
          colour = "darkgrey",
          fill = NA,
          linewidth = 0.3) +
  geom_sf(data = hex_grid, 
          fill = NA, 
          col = "deepskyblue4", 
          linewidth = 0.2) +
  theme_void()

#### Get and extract predictions ####
# Load all the varying predictions
ll <- list.files(path_output, full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "tif")]
# ll <- ll[grep("Threshold_",ll)]
assert_that(length(ll)>0)

# Load raster and get only median
ras_med <- rast(ll)
ras_med <- ras_med[[which(names(ras_med)=="q50")]]
# ras_med <- ras_med[[which(names(ras_med)=="threshold_q50_percentile")]]
names(ras_med) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Ensemble__")
names(ras_med) <- stringr::str_remove(names(ras_med),"__masked")

# Make smaller hexagon within
vertex_coords <- hex_grid |> 
  dplyr::mutate(vertices = purrr::pmap(
    .l = list(x = hex_geometry),
    .f = function(x) {
      x |>
        # st_buffer(dist = -25000) |>       # STEP 2: set size of smaller hex
        st_buffer(dist = -50000) |>       # STEP 2: set size of smaller hex
        st_coordinates() |>               # STEP 3: get vertex coordinates of smaller hex        
        as_tibble() |>                    # convert matrix to tibble  
        st_as_sf(coords = c("X", "Y")) |> # convert tibble to simple features
        dplyr::select(-L1, -L2) |>        # remove unnecessary columns (created from st_coordinates)
        mutate(vertex_position = 1:7)     # STEP 4: number vertices 
    })) |> 
  tidyr::unnest(cols = vertices)

# Get vertex coordinates
vertex_centroid_coords <- vertex_coords |> 
  mutate(geometry = ifelse(vertex_position == 7,      
                           sf::st_centroid(hex_geometry), 
                           geometry)) |> 
  sf::st_drop_geometry()

# --- #
# Now extract values per grid
z1 <- exactextractr::exact_extract(ras_med, hex_grid, "mean")
assert_that(nrow(z1) == nrow(hex_grid))
names(z1) <- stringr::str_remove(names(z1), "Threshold__")
# --- #

zz <- bind_rows(
  z1 |> 
    # dplyr::rename(hex_id = ID) |> 
    dplyr::mutate(hex_id = hex_grid$hex_id) |>
    reshape2::melt(id.vars = "hex_id") |> 
    # Drop missing data
    tidyr::drop_na() |> 
    dplyr::mutate(variable = stringr::str_remove(variable, "mean.")) |>
    # Assign corners
    dplyr::mutate(vertex_position = case_when(
      variable == "Woodland.and.forest" ~ 1,
      variable == "Heathland.and.shrub" ~ 2,
      variable == "Grassland" ~ 3,
      variable == "Sparsely.vegetated.areas" ~ 4,
      variable == "Wetlands" ~ 5,
      variable == "Marine.inlets.and.transitional.waters" ~ 6
    ))
  # z2 |> 
  #   # dplyr::rename(hex_id = ID) |> 
  #   dplyr::mutate(hex_id = 1:nrow(z)) |>
  #   reshape2::melt(id.vars = "hex_id") |> 
  #   # Drop missing data
  #   tidyr::drop_na() |> 
  #   dplyr::mutate(variable = str_remove(variable, "mean.")) |>
  #   # Assign corners
  #   dplyr::mutate(vertex_position = case_when(
  #     variable == "Biodiversity" ~ 1,
  #     variable == "Food" ~ 3,
  #     variable == "Water" ~ 5,
  #     variable == "Climate" ~ 6
  #   )) |> 
  #   dplyr::mutate(ambition = "high")
)

# Join in extracted indicator values
ind_points <- zz |>
  left_join(vertex_centroid_coords |> sf::st_as_sf(crs = sf::st_crs(bb)),
            by = join_by(vertex_position, hex_id)) |> 
  sf::st_as_sf() 

ind_points$variable <- factor(ind_points$variable, levels = names(cols))

# Normalize within group
ind_points <- ind_points |> dplyr::group_by(hex_id) |> dplyr::mutate(value = normalize(value))
# Small filter, considering only points with values above 0.5
# ind_points <- ind_points |> dplyr::filter(value > 0.05)

# Build the plot
gm <- ggplot() +
  theme_gray(base_size = 20) +
  # cowplot::theme_map(font_size = 20) +
  geom_sf(data = bb, alpha = 0.8, fill = "white", 
          colour = "black") +
  geom_sf(data = hex_grid, fill = NA) +
  # Add points
  geom_sf(data = ind_points,
          aes(colour = variable, size = value)) +
  scale_color_manual(values = cols) +
  # scale_size_binned(range = c(0.5,1.5),n.breaks = 4,nice.breaks = TRUE,transform = "asn") +
  scale_size_binned_area(max_size = 4) +
  # Guides
  guides(colour = guide_legend(title = "", override.aes = list(size=7)),
         size = "none") +
    theme(legend.position = 'bottom',legend.text = element_text(size = 20)) +
  labs(title = "Probability of current PNV") +
  theme(plot.title = element_text(hjust = 0.5,size = 30))
#gm
ggsave(plot = gm, paste0("figures/Figure2_full.png"), width = 14,height = 13,dpi = 400)

#### Basic histograms ####
# Here simply make an area estimate relative to the total land area
ar_total <- cellSize(background, unit = "km") * background

# Also get uncertainty as posterior estimate
ll <- list.files(path_output, full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "tif")]
# ll <- ll[grep("Prediction_",ll)]
#ll <- ll[grep("Threshold_",ll)]
ras <- rast(ll)
ras_low <- ras[[which(names(ras)=="q05")]]
# ras_low <- ras[[which(names(ras)=="threshold_q05_percentile")]]
# names(ras_low) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")
ras_med <- ras[[which(names(ras)=="q50")]]
# ras_med <- ras[[which(names(ras)=="threshold_q50_percentile")]]
# names(ras_med) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")
ras_high <- ras[[which(names(ras)=="q95")]]
# ras_high <- ras[[which(names(ras)=="threshold_q95_percentile")]]
# names(ras_high) <- stringr::str_remove(tools::file_path_sans_ext(basename(ll)),"Prediction__")

ar1 <- ras_low * cellSize(ras_low,unit = "km")
ar2 <- ras_med * cellSize(ras_med,unit = "km")
ar3 <- ras_high * cellSize(ras_high,unit = "km")

names(ar1) <- stringr::str_remove( tools::file_path_sans_ext(basename(ll)), "Ensemble__")
names(ar1) <- stringr::str_remove( names(ar1), "__masked")
names(ar2) <- stringr::str_remove( tools::file_path_sans_ext(basename(ll)), "Ensemble__")
names(ar2) <- stringr::str_remove( names(ar2), "__masked")
names(ar3) <- stringr::str_remove( tools::file_path_sans_ext(basename(ll)), "Ensemble__")
names(ar3) <- stringr::str_remove( names(ar3), "__masked")

# Summarize
ss1 <- terra::global(ar1, "sum", na.rm = TRUE) |> as.data.frame() |> 
  tibble::rownames_to_column(var = "variable") |> rename(low = sum)
ss2 <- terra::global(ar2, "sum", na.rm = TRUE) |> as.data.frame() |> 
  tibble::rownames_to_column(var = "variable") |> rename(med = sum)
ss3 <- terra::global(ar3, "sum", na.rm = TRUE) |> as.data.frame() |> 
  tibble::rownames_to_column(var = "variable") |> rename(high = sum)

ss <- left_join(ss1,ss2) |> left_join(ss3)

ss$variable <- stringr::str_remove(ss$variable, "Threshold__")

ss$variable <- factor(ss$variable, levels = names(cols))

ss <- ss |> mutate(landarea = terra::global(ar_total, "sum", na.rm = TRUE)[,1])
ss[,c(2:5)] <- ss[,c(2:5)] * 1e-6
# Save result
saveRDS(ss, "results/TotalAreaStats.rds")

# Convert to proportion for visualization
ss[,c(2:4)] <- ss[,c(2:4)] / terra::global(ar_total, "sum", na.rm = TRUE)[,1]

# Build histogram
g <- ggplot(ss, aes(x = variable, y = med, fill = variable)) +
  theme_light(base_size = 22) +
  # Bar and uncertainty
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low, ymax = high),linewidth = 1.25,width = 0.5) +
  scale_y_continuous(limits = c(0, 0.55),expand = c(0,Inf),breaks = pretty_breaks(5)) +
  scale_fill_manual(values = cols) +
    # guides(fill = guide_legend(title = "",nrow = 3)) +
    guides(fill = "none") +
    theme(legend.position = "bottom",legend.text = element_text(size = 12)) +
  theme(axis.text.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) +
  labs(x = "", y = "Share of land area occupied (%)") +
  # Transparent
  theme( panel.background = element_rect(fill = 'transparent'),
         plot.background = element_rect(fill = "transparent",colour = NA),
         legend.background = element_rect(fill = "transparent"))
g
ggsave(plot = g, filename = "figures/Figure2_bar.png", width = 4, height = 6,bg = "transparent")
