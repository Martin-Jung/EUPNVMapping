# Define parameters for the runs
# it is expected that species gpkg are in the input folder and

# -------------------------------------------------- #
# -------------------------------------------------- #

# Runname
runname <- "ClimateRun"

# Scale/grain of analysis in m
# It is expected that covariates are available in this grain size
grain <- c("1000", "5000", "10000")[1]

user <- c('martin_local', 'martin_hpg901')[1]


# Cross-validation strategy
strategy_cv <- c("blocks")
strategy_cvnr <- c("blocks" = 3, "metric" = "F1")

# Projections?
doproj <- TRUE
# GCMs
gcms <- c("GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")[c(1)]

# Save models?
modelsave <- TRUE

# Cores and parallel options
cores <- 7

# Be chatty
verbose <- TRUE

# Define home folder. Gets overwritten if not used
path_home <- here::here()

if(user == 'martin_local') td <- "/media/martin/AAB4A0AFB4A08005/tmp/"
if(user == 'martin_hpg901') td <- "~/tmp/"
  
if(dir.exists(td))
terra::terraOptions(tempdir = td) # Overwrite temporary directory in raster

# Path output
if(user == 'martin_local') path_output = "/media/martin/AAB4A0AFB4A08005/"
if(user == 'martin_hpg901') path_output = "~/"
assertthat::assert_that(dir.exists(path_output))

# Path background
path_background <- "data/"
assertthat::assert_that(dir.exists(path_background))

# Figure path
path_figures <- "figures/"
dir.create(path_figures, showWarnings = FALSE)

# Path processed
if(user == 'martin_local') {
  path_rawdata = "/mnt/pdrive/bec_data/100_rawdata/"
  path_processed = "/mnt/pdrive/bec_data/200_processeddata/"
} else if(user == 'martin_hpg901') {
  path_rawdata = "~/100_rawdata/"
  path_processed = "~/200_processeddata/"
}

# Finally create an analysis folder and subfolder in the output path
outdirname <- paste0("PNVHabitats__", runname)
dir.create(file.path(path_output, outdirname), showWarnings = FALSE)
path_output <- file.path(path_output, outdirname)
  
# Default Projection to be used for all inputs and outputs
proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

# Path to NUTS shapefile
if(user == 'martin_local') {
  path_nuts <- "/mnt/hdrive/NaturaConnect/GridsAndBoundaries/Outputs/Vector/PG_gadm_octsaccession_name0_mollweide.gpkg"
} else if(user == 'martin_hpg901') {
  path_nuts <- "~/NaturaConnect/GridsAndBoundaries/Outputs/Vector/PG_gadm_octsaccession_name0_mollweide.gpkg"
}
# path_nuts <- paste0(path_rawdata, "/EU_NUTS2021/NUTS_RG_var_2021_4326/NUTS_RG_01M_2021_4326.shp")

# Clip to NUTS (default, no)
clip_nuts <- FALSE

# Create temporary folder
dir.create(paste0(path_home,"/", "resSaves"), showWarnings = FALSE)

# --------- #
# Path variables to input data. This one is expected to point to a directory with gpkg files
path_habitat <- paste0(path_processed, "/dataclima/habitat_occurrence/")

# Predictors present
# path_presentcovs <- paste0(path_processed, "dataclima/predictors/", grain)
if(user == 'martin_local') path_presentcovs <- paste0("/media/martin/AAB4A0AFB4A08005/PNV_covs")
if(user == 'martin_hpg901') path_presentcovs <- paste0(path_processed, "/PNV_covs")

dir.create(path_presentcovs, showWarnings = FALSE)

# Other predictors related to pnv
# path_pnvcovs <- paste0(path_processed, "dataclima/pnv_baselines/", grain)
path_pnvcovs <- paste0(path_presentcovs, "/OtherVariables/")

# Predictors future
path_future <- c(
  "ssp126" = paste0(path_presentcovs, "/ssp126"),
  "ssp370" = paste0(path_presentcovs, "/ssp370"),
  "ssp585" = paste0(path_presentcovs, "/ssp585")
)

# -------------------------------------------------- #
# -------------------------------------------------- #
# Make an security check that all input files and folders and parameters are correctly set
require(assertthat)
assert_that(
  is.character(runname),
  is.character(grain) && is.numeric(as.numeric(grain)),
  is.dir(path_output),
  is.dir(td),
  is.dir(path_processed),
  is.numeric(cores),
  is.dir(path_presentcovs), length(list.files(path_presentcovs)) > 0,
  is.dir(path_pnvcovs),
  all(sapply(path_future, dir.exists))
)
