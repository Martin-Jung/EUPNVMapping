# Purpose:
# Show the transition of PNV per SSP and type
# Visualized as https://corybrunson.github.io/ggalluvial/index.html

library(terra)
library(ggplot2)
library(sf)
library(tidyterra)
library(assertthat)
library(ggalluvial)
library(scales)

source("00_Functions.R")

# Get reference raster for the given grain
background <- terra::rast( "data/background_ref.tif")

# Get masked ensemble projections
path_output <- "ensemble/Masked/"

# Get them all
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

assertthat::assert_that(nrow(df)>0)

# --------------------- #
#### Extract per time step      ####
# Idea is to summarize ("Total sum")

# Results container
results <- data.frame()
  
pb <- progress::progress_bar$new(total = nrow(df))
for(i in 1:nrow(df)){
  pb$tick()
  ras <- terra::rast( df$ifname[i] )
  ss <- terra::global(ras, "sum", na.rm = TRUE) |> 
    tibble::rownames_to_column("metric")
  results <- dplyr::bind_rows(
    results,
    cbind(df[i,], ss)
  )
}

# Add contemporary conditions
ll_cont <- list.files(path_output, full.names = TRUE)
ll_cont <- ll_cont[has_extension(ll_cont, "tif")]
ll_cont <- ll_cont[grep("Ensemble__",ll_cont)]
assertthat::assert_that(length(ll_cont)>0, msg = "No results found...")

df_cont <- data.frame(ifname = ll_cont, split = tools::file_path_sans_ext(basename(ll_cont)) ) |> 
  tidyr::separate(split,sep = "__",into = c("Proj", "variable")) |> 
  dplyr::select(-Proj)

results_con <- data.frame()
# Extract
for(i in 1:nrow(df_cont)){
  print(i)
  ras <- terra::rast(df_cont$ifname[i])
  names(ras) <- paste0("suitability_", names(ras))
  ras <- subset(ras, c("suitability_q05","suitability_q50","suitability_q95","suitability_mean"))
  ss <- terra::global(ras, "sum", na.rm = TRUE) |> 
    tibble::rownames_to_column("metric")
  results_con <- bind_rows(
    results_con,
    data.frame(cbind(df_cont[i,], ss), timeperiod = as.Date("2010-01-01"))
  )
}

# Join in
results <- rbind(
  results,
  results_con |> tidyr::crossing(ssp = c("ssp126", "ssp370", "ssp585"), gcm = unique(df$gcm) )
)
rm(results_con)

# --- #
# Also add the overall area from the masked background
# Background as basis
bb <- background

#Load latest CORINE layer 
clc_mask <- terra::rast("corine/CorineMask_2018_clundefined.tif")

# Align with background
if(!terra::compareGeom(bb, clc_mask, stopOnError = FALSE)){
  if(sf::st_crs(bb) != sf::st_crs(clc_mask)) clc_mask <- terra::project(clc_mask, bb)
  # Crop and mask
  clc_mask <- terra::crop(clc_mask, bb)
  clc_mask <- terra::mask(clc_mask, bb)
}

# Any grid cells with more than 50% water bodies or urban surface -> mask out from bb
m <- clc_mask
clc_mask[clc_mask <=5000] <- NA
clc_mask[clc_mask>0] <- 1

bb <- terra::mask(bb, clc_mask, inverse = TRUE)
rm(clc_mask)

# Add
results$totalarea <- terra::global(bb, 'sum', na.rm = TRUE)[,1]

# Save the result
saveRDS(results, "results/FutureAreaStats.rds")

# --------------------- #
#### Build the alluvial plot ####
# Load data and plot as alluvial figure
results <- readRDS("results/FutureAreaStats.rds")

# Colours
cols <- c("Woodland.and.forest"  = "#117733",
          "Heathland.and.shrub" = "#B65719",
          "Grassland" = "#DDCC77",
          "Sparsely.vegetated.areas" = "#AA4499",
          "Wetlands"  = "#332288",
          "Marine.inlets.and.transitional.waters" = "#88CCEE"
)
# Variable factor loadings
results$variable <- factor(results$variable, levels = names(cols))

# area formatting to million ha
results$sum <- (results$sum * 100) / 1e6

# Focus on mean only and for a specific ssp 
oo <- results |> dplyr::filter(metric == "suitability_mean") |> 
  arrange(timeperiod, ssp)

# Convert to Date
oo$timeperiod <- factor(oo$timeperiod, labels = c("2010", "2040", "2070", "2100"))

# Load
oo$count <- rep(1:6,( dplyr::n_distinct(oo$timeperiod))*(dplyr::n_distinct(oo$ssp)))

# Format SSPs
oo$ssp <- factor(oo$ssp, levels = c("ssp126", "ssp370", "ssp585"),
                      labels = c("SSP1-2.6", "SSP3-7.0", "SSP5-8.5"))

g <- ggplot(oo,
       aes(x = timeperiod, y = sum, stratum = variable, alluvium = count,
           fill = variable)) +
  theme_grey(base_size = 20) +
  scale_x_discrete(expand = c(.2, 0)) +
  geom_flow(width = 2/4) +
  geom_stratum(alpha = .5, width = 2/4) +
  scale_linetype_manual(values = c("blank", "solid")) +
  scale_fill_manual(values = cols) + 
    guides(fill = guide_legend(title = "")) +
    theme(legend.position = "bottom") +
  # Facets
  facet_wrap(~ssp) + 
  # Scale y axis
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  geom_text(stat = "stratum",
            aes(label = round(sum, 1) )) +
  # ggrepel::geom_text_repel(
  #   aes(label = ifelse(as.numeric(survey) == 1, as.character(response), NA)),
  #   stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  # ) +
  # ggrepel::geom_text_repel(
  #   aes(label = ifelse(as.numeric(survey) == 3, as.character(response), NA)),
  #   stat = "stratum", size = 4, direction = "y", nudge_x = .5
  # ) +
  theme(legend.position = "none") +
  labs(x = "", y = "Total area (million ha)")
g

ggsave(plot = g, filename = "figures/Figure4.png", width = 14,height = 8,dpi = 400)

