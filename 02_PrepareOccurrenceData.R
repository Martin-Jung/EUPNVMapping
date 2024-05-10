# Prepare habitat data
# We would want to create pnv projections of habitats sensu Hengl et al. 
# To do so we need training data for these habitats which we can obtain from LUCAS, the EEA habitat predictions and also the shapefiles 
# as part of the EU habtitat assessment. 

# We will use the MAES system for all classes and not distinguish further between the secondary habitats. 
# We will also ignore marine habitats for now.

# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")
# NOTE: The code below is not well reproducible owing to the range of different data types and storage locations.

# ---- #

# Load reference grids as list
background <- terra::rast( paste0(path_background, "ReferenceGrid_Europe_bin_",grain,"m.tif") )

# Get habitat data crosswalks from dataclima
data("crosswalk_eunisclc")

# 20 - 2Heathland&Shrub
# 30 - 3Grassland
# 40 - 5Cropland
# 50 - 6Urban
# 60 - 4Sparsely vegetated areas
# 70 - 4Sparsely vegetated areas
# 80 - 8Rivers and lakes
# 90 - 7Wetlands
# 100 - 4Sparsely vegetated areas
# 111-116 - 1Woodland and forest
# 121- 126 - 1Woodland and forest
# 200 - 10Marine
## Format EU Habitat data

# Load these data
eu_habs <- sf::st_read(paste0(path_rawdata, "EU_HabitatsStatusTrends/2020/art17_2013_2018_public.gpkg"), "Art17_habitats_distribution_2013_2018_EU")
eu_habs_desc <- sf::st_read(paste0(path_rawdata, "EU_HabitatsStatusTrends/2020/art17_2013_2018_public.gpkg"), "habitatsEUassessment")

# Get the crosswalk
data("crosswalk_annexIeunis")
data("crosswalk_eunisclc")

# Join them to map
# Remapped cr
remcr <- dplyr::left_join(
  crosswalk_annexIeunis %>% dplyr::select(Natura.2000.code, EUNIS.code) %>% distinct(), 
  crosswalk_eunisclc %>% dplyr::select(EUNIS.Code, CLC.Code, CLC.Name) %>% distinct(),
  by = c("EUNIS.code" = "EUNIS.Code")
) 
# Get the dominant class name for unmapped ones
for(i in 1:nrow(remcr)){
  if(is.na(remcr$CLC.Name[i])){
    # Get all codes for the N2k code
    sub <- subset(remcr, Natura.2000.code == remcr$Natura.2000.code[i])
    val <- terra::modal(sub$CLC.Name, na.rm = TRUE)
    if(is.na(val)){
      sub <- subset(remcr, stringr::str_sub(Natura.2000.code, 0, 2) == stringr::str_sub(remcr$Natura.2000.code[i], 0, 2) )
      val <- terra::modal(sub$CLC.Name, na.rm = TRUE)
    }
    remcr$CLC.Name[i] <- val
  }
}

# Reclassify all habitats using the corine to MAES crosswalk 
# Sadly not all habitats are covered by the crosswalk, so we will do this exercise manually
# https://biodiversity.europa.eu/ecosystems/mapping-and-assessment-of-ecosystems-and-their-services-maes-1/correspondence-between-corine-land-cover-classes-and-ecosystem-types
remcr$MAES <- forcats::fct_collapse(remcr$CLC.Name,
                                    Marine = c("Sea and ocean"),
                                    "Marine inlets and transitional waters" = c("Intertidal flats", "Salt marshes", "Estuaries", "Coastal lagoons"),
                                    "Rivers and lakes" = c("Water bodies","Water courses"),
                                    "Cropland" = c("Permanently irrigated land", "Agro-forestry area"),
                                    "Heathland and shrub" = c("Moors and heathland", "Sclerophyllous vegetation"),
                                    "Grassland" = c("Pasture", "Natural grassland"),
                                    "Woodland and forest" = c("Coniferous fores", "Broad-leaved fores",
                                                              "Coniferous forest","Transitional woodland shrub", "Broad-leaved forest"),
                                    "Urban" = c(),
                                    "Sparsely vegetated areas" = c("Bare rock","Beaches, dunes, and sand plains", "Sparsely vegetated areas", "Glaciers and perpetual snow"),
                                    "Wetlands" = c("Peatbogs", "Inland marshes")
)

# Join in 
eu_habs  <- dplyr::left_join(eu_habs, remcr, by = c("habitatcode" = "Natura.2000.code"))
# Remove any missing ones
eu_habs <- eu_habs %>% dplyr::filter(!is.na(MAES))

# Grassland, Urban, Cropland, Grassland, Woodland and forest, Heathland and shrub, Sparsely vegetated areas, Wetlands, Marine inlets and transitional waters, Rivers and lakes, Marine
n_distinct(eu_habs$category)

# Now process
for(h in unique(eu_habs$MAES)){
  print(h)
  sub <- subset(eu_habs, MAES == h) %>% laea_transform()
  sub$record_type <- "po"
  
  hname <- paste0(make.names(h), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eu_n2k_habitats", append = FALSE, quiet = TRUE)
  
}


## Format LUCAS data -------
# The LUCAS records cover a large number of land cover training points for Europe

path_lucas <- "/mnt/pdrive/bec_data/100_rawdata/EU_LUCAS/harmonized/"

# Load the LUCAS geogemtries
lucas <- sf::st_read(paste0(path_lucas, "/2_geometry/LUCAS_gps_geom/LUCAS_gps_geom.shp"))

# Get the latest survey data
lucas_srv2018 <- readr::read_csv(paste0(path_lucas, "/1_table/lucas_harmo_uf_2018/lucas_harmo_uf_2018.csv")) %>%
  dplyr::select(id, point_id, year, survey_date, nuts0:nuts3, th_lat, th_long, lc1:lc2_perc)

# Get only those where the primary habitat is larger than 50%?

lucas_srv2018$MAES_cl1 <- forcats::fct_collapse(lucas_srv2018$lc1_label,
                                                "Marine inlets and transitional waters" = c("Salt marshes", 
                                                                                            "Inland salty water bodies", "Intertidal flats", "Transitional water bodies","Inland salty running water"),
                                                "Cropland" = c("Common wheat", "Clovers", "Durum wheat",
                                                               "Rye", "Permanent crops (only pi)",
                                                               "Arable land (only pi)",
                                                               "Potatoes", "Barley", "Olive groves","Maize",
                                                               "Vineyards", "Oats", "Apple fruit",
                                                               "Sunflower", "Mixed cereals for fodder", "Oranges",
                                                               "Cherry fruit", "Sunflower", "Rape and turnip rape",
                                                               "Dry pulses", "Other root crops", "Rice"
                                                ),
                                                "Heathland and shrub" = c("Shrubland without tree cover"),
                                                "Grassland" = c("Grassland with sparse tree/shrub cover",
                                                                "Grassland without tree/shrub cover",
                                                                "Temporary grasslands"),
                                                "Woodland and forest" = c("Broadleaved woodland",
                                                                          "Pine dominated coniferous woodland",
                                                                          "Spruce dominated mixed woodland",
                                                                          "Shrubland with sparse tree cover",
                                                                          "Spruce dominated coniferous woodland",
                                                                          "Pine dominated mixed woodland",
                                                                          "Other mixed woodland",
                                                                          "Other coniferous woodland"
                                                ),
                                                "Urban" = c("Other artificial areas", 
                                                            "Buildings with more than 3 floors",
                                                            "Buildings with 1 to 3 floors",
                                                            "Non built-up linear features"),
                                                "Sparsely vegetated areas" = c("Other bare soil",
                                                                               "Spontaneously vegetated surfaces", "Sand",
                                                                               "Lichens and moss", 
                                                                               "Glaciers, permanent snow"),
                                                "Rivers and lakes" = c("Inland fresh water bodies", 
                                                                       "Inland fresh running water"),
                                                "Wetlands" = c("Peatbogs", "Inland marshes"),
                                                other_level = "other"
)
lucas_srv2018$MAES_cl1[(lucas_srv2018$MAES_cl1 == "other")] <- "Cropland"
lucas_srv2018$MAES_cl1 <- droplevels(lucas_srv2018$MAES_cl1)

# Join the survey information with the point data
lucas <- dplyr::left_join(lucas, lucas_srv2018, by = c("POINT_ID" = "point_id"))
lucas <- lucas %>% dplyr::filter(!is.na(id))
# Rename lucas ids
lucas <- lucas %>% dplyr::rename(lucas_id = id) %>% dplyr::select(-ID)
lucas <- lucas %>% dplyr::select(-YEAR)

# Manual writing hack (since writing to network gpkg fails)
# outpath <- "/media/martin/data/habitat_occurrence/"
rename_geometry <- function(g, name){
  current = attr(g, "sf_column")
  names(g)[names(g)==current] = name
  st_geometry(g)=name
  g
}

# Now process
for(h in unique(lucas$MAES_cl1)){
  if(is.na(h)) next()
  print(h)
  sub <- subset(lucas, MAES_cl1 == h) %>% laea_transform()
  sub$record_type <- "po"
  # sub <- rename_geometry(sub, "geom")
  
  hname <- paste0(make.names(h), ".gpkg")
  
  sf::st_write(obj = sub, dsn = paste0(outpath, hname),
               layer = "lucas", append = FALSE, quiet = F)
}

## EU Ecoregions -------
# Load the ecoregion layer

## Format EUNIS habitat data

# Temporary copy for adding
outpath <- "/media/martin/data/habitat_occurrence/"

# Load the crosswalk
eunis_path <- paste0("/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/")

# Load Coastal habitat data
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(N) COASTAL HABITAT TYPES/n_distribution_gpkg/N_distribution.gpkg"
lyrs <- sf::st_layers(wp)
cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"),sheet = "Coastal") %>% 
  dplyr::select(Level, Code, Name) %>% distinct() %>% 
  dplyr::filter(Code %in% lyrs$name)

# Manually add MAES classes
cr$MAES <- NA
cr$MAES[grep("beach|cliff|slack|dune",cr$Name)] <- "Sparsely vegetated areas"
cr$MAES[grep("forest",cr$Name)] <- "Woodland and forest"
cr$MAES[grep("heath|scrub",cr$Name)] <- "Heathland and shrub"
cr$MAES[grep("grassland",cr$Name)] <- "Grassland"
assertthat::assert_that(!anyNA(cr$MAES))

# Now loop through each class and add to the datafiles
for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("fme_feature_type" = "Code")) %>% 
    dplyr::rename(Code = "fme_feature_type")
  
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# --- #
# Grasslands
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(R) GRASSLAND AND LANDS DOMINATED BY FORBS, MOSSES OR LICHENS HABITAT TYPES/r_distribution_gpkg/R_distribution.gpkg"
lyrs <- sf::st_layers(wp)

cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"),sheet = "Grassland") %>% 
  dplyr::select(Level, Code, Name) %>% distinct() %>% 
  dplyr::filter(Code %in% lyrs$name)

# Manually add MAES classes
cr$MAES <- NA
cr$MAES[grep("steppe|grassland|meadow|pasture|vegetation|fern",cr$Name)] <- "Grassland"
cr$MAES[grep("Snow-bed vegetation|desert",cr$Name)] <- "Sparsely vegetated areas"
cr$MAES[grep("heath|scrub",cr$Name)] <- "Heathland and shrub"
cr$MAES[grep("marsh",cr$Name)] <- "Wetlands"
cr$MAES[grep("Forest|forest",cr$Name)] <- "Woodland and forest"
# cr[which(is.na(cr$MAES)),]
assertthat::assert_that(!anyNA(cr$MAES))

for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("fme_feature_type" = "Code")) %>% 
    dplyr::rename(Code = "fme_feature_type")
  
  if(hasName(sub, "UNCERTAINT")) {
    sub <- dplyr::rename(sub, "PRECISION" = "UNCERTAINT")
    if(nrow(sub)==0) next()
  }
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# ---- #
# Heathland and shrub datasets
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(S) HEATHLAND, SCRUB AND TUNDRA HABITAT TYPES/s_distribution_gpkg/S_distribution.gpkg"
lyrs <- sf::st_layers(wp)

cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"),sheet = "Heathland") %>% 
  dplyr::select(Level, Code, Name) %>% distinct() %>% 
  dplyr::filter(Code %in% lyrs$name)

# Manually add MAES classes
cr$MAES <- NA
cr$MAES[grep("Shrub|shrub|heath|garrigue|scrub|maquis",cr$Name)] <- "Heathland and shrub"
cr$MAES[grep("desert|Moss and lichen",cr$Name)] <- "Sparsely vegetated areas"
# cr[which(is.na(cr$MAES)),]
assertthat::assert_that(!anyNA(cr$MAES))

for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("fme_feature_type" = "Code")) %>% 
    dplyr::rename(Code = "fme_feature_type")
  
  if(hasName(sub, "UNCERTAINT")) {
    sub <- dplyr::rename(sub, "PRECISION" = "UNCERTAINT")
    if(nrow(sub)==0) next()
  }
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# ---- #
# Forests
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(T) FOREST AND OTHER WOODED LAND HABITAT TYPE/T_distribution.gpkg/T_distribution.gpkg"
lyrs <- sf::st_layers(wp)

cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"),sheet = "Forest") %>% 
  dplyr::select(Level, Code, Name) %>% distinct() %>% 
  dplyr::filter(Code %in% lyrs$name)

# Manually add MAES classes
cr$MAES <- NA
cr$MAES[grep("forest|trees|vegetation|taiga",cr$Name)] <- "Woodland and forest"
# cr[which(is.na(cr$MAES)),]
assertthat::assert_that(!anyNA(cr$MAES))

for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("dest_featureclass" = "Code")) %>% 
    dplyr::rename(Code = "dest_featureclass")
  
  if(hasName(sub, "UNCERTAINT")) {
    sub <- dplyr::rename(sub, "PRECISION" = "UNCERTAINT")
    if(nrow(sub)==0) next()
  }
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# Sparsely vegetated stuff
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(U) INLAND HABITATS WITH NO OR LITTLE SOIL AND MOSTLY WITH SPARSE VEGETATION (Sparsely vegetated habitats)/u_distribution_gpkg/U_distribution.gpkg"
lyrs <- sf::st_layers(wp)

cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"),sheet = "Sparsely vegetated") %>% 
  dplyr::select(Level, Code, Name) %>% distinct() %>% 
  dplyr::filter(Code %in% lyrs$name)

# Manually add MAES classes
cr$MAES <- NA
cr$MAES[grep("scree|cliff|desert|volcanic",cr$Name)] <- "Sparsely vegetated areas"
# cr[which(is.na(cr$MAES)),]
assertthat::assert_that(!anyNA(cr$MAES))

for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("dest_featureclass" = "Code")) %>% 
    dplyr::rename(Code = "dest_featureclass")
  
  if(hasName(sub, "UNCERTAINT")) {
    sub <- dplyr::rename(sub, "PRECISION" = "UNCERTAINT")
    if(nrow(sub)==0) next()
  }
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# ---- #
# Marine stuff / wetlands
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(MA2) LITTORAL BIOGENIC HABITAT TYPES (SALT MARCHES)/m_distribution_gpkg/M_distribution.gpkg"
lyrs <- sf::st_layers(wp)

# cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"), sheet = "Grassland") %>% 
#   dplyr::select(Level, Code, Name) %>% distinct() %>% 
#   dplyr::filter(Code %in% lyrs$name)
# Not existing, thus map manually
cr <- data.frame(Level = 1, Code = lyrs$name, Name = NA)

# Manually add MAES classes
# Here looked up from the documents!
cr$MAES <- NA
cr$MAES[cr$Code == "MA211"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA221"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA222"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA223"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA224"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA225"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA232"] <- "Sparsely vegetated areas"
cr$MAES[cr$Code == "MA241"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA251"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA252"] <- "Marine inlets and transitional waters"
cr$MAES[cr$Code == "MA253"] <- "Marine inlets and transitional waters"
# cr[which(is.na(cr$MAES)),]
assertthat::assert_that(!anyNA(cr$MAES))

for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("dest_featureclass" = "Code")) %>% 
    dplyr::rename(Code = "dest_featureclass")
  
  if(hasName(sub, "UNCERTAINT")) {
    sub <- dplyr::rename(sub, "PRECISION" = "UNCERTAINT")
    if(nrow(sub)==0) next()
  }
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# Man-made stuff
wp <- "/mnt/pdrive/bec_data/100_rawdata/EU_EUNIS2021_PredictedRasterVector/(V) VEGETATED MAN-MADE HABITATS/v_distribution_gpkg/V_distribution.gpkg"
lyrs <- sf::st_layers(wp)

cr <- read_xlsx(paste0(eunis_path, "EUNIS terrestrial habitat classification 2021 including crosswalks.xlsx"),sheet = "Man-made") %>% 
  dplyr::select(Level, `Code 2018`, `Name 2018`) %>% distinct() %>% 
  dplyr::rename(Code = `Code 2018`, Name = `Name 2018`) %>% 
  dplyr::filter(Code %in% lyrs$name)

# Manually add MAES classes
cr$MAES <- NA
cr$MAES[grep("Intensive|crop|fallow|vegetation",cr$Name)] <- "Cropland"
cr$MAES[grep("grassland",cr$Name)] <- "Grassland"
# cr[which(is.na(cr$MAES)),]
assertthat::assert_that(!anyNA(cr$MAES))

for(m in lyrs$name){
  print(m)
  sub <- sf::st_read(wp, m, quiet = T) %>% laea_transform()
  if(!(m %in% cr$Code)) next()
  sub <- dplyr::left_join(sub, cr, by = c("dest_featureclass" = "Code")) %>% 
    dplyr::rename(Code = "dest_featureclass")
  
  if(hasName(sub, "UNCERTAINT")) {
    sub <- dplyr::rename(sub, "PRECISION" = "UNCERTAINT")
    if(nrow(sub)==0) next()
  }
  if(hasName(sub, "PRECISION")) {
    sub <- sub %>% dplyr::filter(PRECISION <= 2000) # Accuracy
    if(nrow(sub)==0) next()
  }
  sub$record_type <- "po"
  
  sub <- sub %>% dplyr::select(YEAR, PRECISION, Code, Name, MAES, record_type)
  
  hname <- paste0(make.names(unique(sub$MAES)), ".gpkg")
  
  sf::st_write(obj = sub, dsn = file.path(outpath, hname),
               layer = "eunis_habitats", append = TRUE, quiet = T)
}

# ------------------------------- #
#### Assign European habitat indicator species here ####
library(readxl)
library(ggplot2)

# Get the vegetation list from Chrity et al and filter.
charspecies <- readxl::read_xlsx(paste0(
  path_rawdata, "EU_EUNIS_VegetationPlotTaxonomicCrosswalks/Characteristic-species-combinations-EUNIS-habitats-2021-06-01.xlsx"
))
names(charspecies) <- make.names(names(charspecies))
charspecies <- charspecies |> 
  # Add habitat group
  dplyr::mutate(Habitat.group = stringr::str_sub(Habitat.code, 0, 1)) |> 
  # Ignore Man-made habitats
  dplyr::filter(Habitat.group != "V") |> 
  # Split
  dplyr::mutate(Speciesname = Species) |> 
  tidyr::separate(Species, into = c("Genus", "Species", "Epithet"),sep = " ") 
# Filter to only Diagnostic species
charspecies <- charspecies |> dplyr::filter(Species.type == "Diagnostic")

# Assess coverage across habitat.groups
ggplot2::ggplot(charspecies, aes(x = Habitat.code, y = Value)) +
  theme_classic(base_size = 18) +
  geom_hline(yintercept = c(25,50, 75), color = "red") +
  geom_boxplot() + 
  labs(y = "Phi-coefficient") +
  # stat_summary(fun = "mean_se", color = "red", size = 2, geom = "point") +
  facet_wrap(~Habitat.group,scales = "free_x") +
  theme(axis.text.x.bottom = element_text(angle = 90))

# Now take for each Habitat.type the top 2 ranked species to check for.
splist <- charspecies |> dplyr::arrange(desc(Value)) |>
  dplyr::group_by(Habitat.code) |> 
  dplyr::slice_head(n = 2) |> dplyr::ungroup() |> 
  dplyr::select(Habitat.group, Habitat.code:Habitat.name, Genus:Epithet, Speciesname) |> 
  distinct()
# Assign habitat codes
splist$Habitat.group <- factor(splist$Habitat.group, levels = c("M","N","Q","R","S","T","U"),
                               labels = c("Marine inlets and transitional waters", "Marine inlets and transitional waters",
                                 "Wetlands", "Grassland", "Heathland and shrub", "Woodland and forest", "Sparsely vegetated areas") )
# ---- #
# Now search for occurrence records respectively
library(rgbif)
library(fst)
library(data.table)
# Create a bounding box
bb <- project(background, terra::crs("WGS84"))
# ⁠min-longitude, min-latitude, max-longitude, max-latitude
wkt <- sf::st_bbox(bb) |> sf::st_as_sfc() |> sf::st_as_text()
spnames <- splist$Speciesname

# prepare download
pwd <- read.csv("gbifpasswd.txt")
options("gbif_user" = pwd$username,
        "gbif_pwd" = pwd$pw,
        "gbif_email" = pwd$email)

# Query key names
keys <- vector()
pb <- progress::progress_bar$new(total = length(spnames))
for(s in spnames){
  pb$tick()
  q <- rgbif::name_backbone(name = s)
  if(hasName(q, "usageKey")){
    qq <- q["usageKey"]
    names(qq) <- q$verbatim_name[1]
    keys <- append(keys, qq)
  }
  rm(q)
}
keys <- do.call(c, keys)

gbif_download <- rgbif::occ_download(
  rgbif::pred_in("taxonKey", keys),
  # rgbif::pred_in("scientificName", spnames),
  rgbif::pred("hasCoordinate", TRUE),
  rgbif::pred("hasGeospatialIssue", FALSE),
  # rgbif::pred_gte("coordinateUncertaintyInMeters", 0.0),
  rgbif::pred_lte("coordinateUncertaintyInMeters", 1000),
  rgbif::pred("occurrenceStatus","PRESENT"),
  rgbif::pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION")),
  rgbif::pred_gte("year", 2000), rgbif::pred_lte("year", 2024),
  rgbif::pred_within( wkt ),
  format = "SIMPLE_CSV")

# Reference
# GBIF Occurrence Download https://doi.org/10.15468/dl.6z2f6w  Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-01-30
# Check download status
occ_download_wait('0178396-240321170329656')

# Split externally into smaller batches
# cat 0178396-240321170329656.csv | parallel --header : --pipe -N999999 'cat >split_file_{#}.csv'

assertthat::assert_that(exists("splist"))
ll <- list.files("data/", pattern = "split_file", full.names = TRUE)
ll <- ll[has_extension(ll, "csv")]
library(fst)

# Now process overall to number of observations per habitat type
for(f in ll){
  print(basename(f))
  a <- data.table::fread(f)
  assertthat::assert_that(utils::hasName(a, "scientificName")) # Check that names have correctly parsed
  # Anything less or equal than 1km
  sps <- a |> dplyr::filter(coordinateUncertaintyInMeters <= 1000) |> 
    # After 2000
    dplyr::filter(year >= 2000) |> 
    dplyr::select(kingdom:species, scientificName,verbatimScientificName,countryCode,
                  decimalLatitude,decimalLongitude, coordinateUncertaintyInMeters,
                  eventDate, day:year, taxonKey, speciesKey) |> dplyr::distinct()
  rm(a)
  
  # Join with habitat codes
  sps <- dplyr::inner_join(splist |> dplyr::select(Habitat.code, Speciesname) |> distinct(), sps, by = c("Speciesname"="species")) 
  if(nrow(sps)==0) next()
  # Now save per Habitat. code
  pb <- progress::progress_bar$new(total = length(sps$Habitat.code))
  for(code in unique(sps$Habitat.code)){
    pb$tick()
    dir.create("data/Habitatsplits",showWarnings = FALSE)
    sub <- subset(sps, Habitat.code == code)
    ofname <- paste0("data/Habitatsplits/", code, ".fst")
    if(file.exists(ofname)){
      d <- fst::read_fst(ofname)
      fst::write_fst( bind_rows(d, sub), ofname ) # overwrite
    } else {
      fst::write_fst(sub, ofname) # Create output
    }
  }
}

# Now per habitatcode rasterize at 1km
ll <- list.files("data/Habitatsplits", full.names = TRUE)
ll <- ll[has_extension(ll, "fst")]
pb <- progress::progress_bar$new(total = length(ll))
for(f in ll){
  pb$tick()
  sub <- read_fst(f) |> sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = sf::st_crs(4326)) |> 
    sf::st_transform(crs = sf::st_crs(background)) |> 
    dplyr::select(Speciesname) |> distinct()
  # Count the number of species per grid cell and remove duplicates
  sub$cellid <- terra::cellFromXY(background, sf::st_coordinates(sub))
  sub <- sub |> dplyr::filter(!duplicated(cellid))
  if(nrow(sub)==0) next()
  ras <- terra::rasterize(sub, background, field = "Speciesname", fun = "count")
  ras <- terra::mask(ras, background)
  
  # Save the result as new tif under habitats split
  writeRaster(ras, paste0("data/Habitatsplits/",basename(tools::file_path_sans_ext(f)), ".tif"), overwrite = TRUE)
  rm(ras)
}

# Now finally decide which grid cells to use for reference conditions of aggregated habitats.
# Here we aggregate per tif in a given group and consider the cell if at least 5 descriptive species have been 
# found in the respective grid cell. This is assumed to be indicative of community representation.
# Assess per PU how many occur. If 3 or more, take as evidence of appropriate site (community)
ll <- list.files("data/Habitatsplits", full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
# Make processing groups
gr <- splist |> dplyr::select(Habitat.group,Habitat.code) |> distinct()
gr$ofname <- paste0("/media/martin/AAB4A0AFB4A08005/habitat_occurrence/", make.names(gr$Habitat.group), ".gpkg")

# Join with file names
gr <- left_join( data.frame(ifname = ll, Habitat.code = tools::file_path_sans_ext(basename(ll))),
           gr)

assertthat::assert_that(
  nrow(gr)>0,
  all(file.exists(unique(gr$ofname)))
)

for(e in unique(gr$Habitat.group)){
  print(e)
  sub <- subset(gr, Habitat.group == e)
  
  # Load all and 
  o <- rast(sub$ifname) |> sum(na.rm = TRUE)
  
  # Cutoff of 5
  o[o<5] <- NA
  
  # Convert to point
  o <- as.data.frame(o, xy = TRUE)
  o <- o |> sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(background)) |> 
    dplyr::mutate(Habitat.group = e)
  
  # Get output file
  sf::write_sf(o, unique(sub$ofname), layer = "GBIF_CharacteristicVeg")
}
