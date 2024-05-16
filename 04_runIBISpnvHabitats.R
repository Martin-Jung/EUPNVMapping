# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")

path_habitat <- "/media/martin/AAB4A0AFB4A08005/habitat_occurrence"
# Path variables to input data. This one is expected to point to a 

# Set run cores for parallization
options('ibis.nthread' = cores)
options('ibis.runparallel' = ifelse(cores>1, FALSE, FALSE) ) # Set to FALSE for slower, but memory proof runs

# -------------------- #
#### Prepare files for modelling ####

# Get reference raster for the given grain
background <- terra::rast( paste0(path_background, "background_ref.tif") )

# Get all habitat databases
habitatdbs <- list.files(path = path_habitat, full.names = TRUE, recursive = TRUE)
habitatdbs <- habitatdbs[has_extension(habitatdbs, "gpkg")]
habitatdbs <- habitatdbs[grep("Cropland|Urban|River", habitatdbs,invert = TRUE)]
assertthat::assert_that(length(habitatdbs)>0)

# Get selected static covariates
ll <- list.files(paste0(path_presentcovs),"tif", full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
# Get ERA BIOCLIM variables
ll <- ll[grep("BIO|era5|potential|DistanceCoast|volumetric|growing|frost-days|EUDEM|cloud-cover|aridity|GroundWater|Lithology",ll)]
ll <- ll[grep("aspect",ll,invert = TRUE)] # Remove aspect as we have sinsuiodal transformed eastness/northness
prescovs <- terra::rast(ll)

# Load full stack of all PNV covariates
ll <- list.files(paste0(path_pnvcovs),"tif", full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
pnvcovs <- terra::rast(ll)

# Clip to nuts if set
if(clip_nuts){
  # Clip with NUTS (excluding Turkey in the process)
  nuts <- sf::st_read(path_nuts) %>% dataclima:::laea_transform() %>%
    dplyr::filter(CNTR_CODE != "TR") %>%
    dplyr::filter(CNTR_CODE != "NO", CNTR_CODE != "IS")
  background <- dataclima:::clip_nuts(background, nuts)
  background <- terra::mask(background, nuts)
}

# --- #
varspnv <- c("EuroVegMap_closedForest",
             "EuroVegMap_Grassland",
             "EuroVegMap_Shrub", "EuroVegMap_Wetlands", "TootchiFloodedWetlands",
             "EuroVegMap_openForest", "EuroVegMap_SparseVeg")

# Now fill every layer with 0 unless set
fullcovs <- c(prescovs, pnvcovs[[varspnv]])
fullcovs <- terra::crop(fullcovs, background)
fullcovs <- terra::mask(fullcovs, background)
fullcovs <- ibis.iSDM:::predictor_homogenize_na(fullcovs)

# Scale all variables except the factor ones for unit consistency
new <- predictor_transform(
  fullcovs[[grep("koeppen|Lithology", names(fullcovs),invert = T)]],option = "scale"
)
# Rasterize the lithology and koeppen class
new <- c(new,
         fullcovs$Lithology_classes,
         fullcovs$koeppen.geiger.class_era5.to.1km_1979.2018.mean_v1.0)
fullcovs <- new
rm(new)
assertthat::assert_that(compareRaster(fullcovs, background),
                        terra::global(fullcovs$aridity_annual.mean_era5.to.1km_1979.2018.mean_v1.0, "max")[,1] !=
                          terra::global(fullcovs$BIO05_era5.to.1km_1979.2018.mean_v1.0, "max")[,1]
                        )

# Match PNV layers to IUCN codes
pnv_legends <- data.frame(layer = names(pnvcovs)[1:16], IUCN.code = as.character(c(1, 4, 3, 6, 1, 4, 1, 3, 6, 5,4, 3, 5, 6, 5, 1)) )
pnv_legends <- dplyr::left_join(pnv_legends, 
                                crosswalk_iucnmaes %>% dplyr::filter(IUCN.level==1) %>%  dplyr::select(IUCN.code, MAES.code, MAES.name)
)
assertthat::assert_that(pnv_legends$layer[8] == "EuroVegMap_Shrub",
                        pnv_legends$MAES.name[8] == "Heathland and shrub") # Security checks

# Default pseudo-absence information
abss <- pseudoabs_settings(background = NULL, nrpoints = 0, min_ratio = 1, method = "random",
                           buffer_distance = 10000)

assertthat::assert_that(
  length(habitatdbs)>0
)
# Fire up a cluster for parallel processing
#library(doMC); library(doParallel)
#cores <- 3
# cores <- min(
#  floor(parallel::detectCores()*2/3),
#  floor((memory.limit()*0.8/1000)/10)
# )
# cl <- parallel::makeCluster(cores, outfile = "")
# invisible(parallel::clusterEvalQ(cl, sink(file.path(getwd(), paste0(Sys.getpid(), ".txt")))))
#doParallel::registerDoParallel(cl)
#registerDoMC(cores = cores)

# -------------------- #
#### Run sequential ####
#proc <- foreach(species = iter(sort(speciesdbs)),
#        .errorhandling = 'stop',
#        .inorder = TRUE,
#        .packages = c('sf','dataclima', "ibis.iSDM", "assertthat", "raster", "modEvA", "sf"),
#        .export = c('path_output','fullcovs')
#        ) %dopar% {
for(hab in sort(habitatdbs)){ # hab <- habitatdbs[1] #speciesdbs[50]->species # Rana_macrocnemis
  
  # Process per directory of files
  d <- basename(dirname(hab))
  outdir <- paste0(path_output, outdirname, "/", d)
  dir.create(outdir, showWarnings = FALSE)
  assertthat::assert_that(dir.exists(outdir))
  # Also create subfolders for individual runs and ensembles in there
  dir.create(paste0(outdir, "/modelresults"),showWarnings = FALSE)
  
  # Group name
  grp <- basename(dirname(hab))
  # Species name
  sname <- tools::file_path_sans_ext(basename(hab))
  # Variables
  vars <- names(fullcovs)
  assertthat::assert_that(length(vars)>0)
  
  # If file already existing, skip
  if(file.exists(
    paste0(outdir, "/", "Threshold__",tools::file_path_sans_ext(basename(hab)), ".tif")
  )) next() # return(NULL)
  
  myLog("[Starting]", "green", "Now processing ", sname)
  # Loading data for the target species
  lyrs <- sf::st_layers(hab)
  
  # --- #
  # Now load all cleaned data as list
  data <- list()
  # stop("Grep population_lower to remove absence points")
  for(name in lyrs$name) {
    o <- sf::st_read(hab, layer = name, quiet = TRUE)
    if(!has_name(o, "observed")){
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
  for(l in grep("Point", unlist(lyrs$geomtype),value = F)){
    df <- rbind(df,
                data[[l]] %>% dplyr::select(observed) %>% 
                  dplyr::mutate(type = lyrs$name[l])
    )
  }
  if(nrow(df)==0) next() # return(NULL) # No point data found...
  df <- ibis.iSDM:::guess_sf(df)
  # Also if we have Natura 2000 layers, add simulated points (nr presence * 2) to this

  # Aggregate observations to grid as to remove duplicates
  df <- ibis.iSDM:::aggregate_observations2grid(df, template = background, field_occurrence = "observed") %>% 
    dplyr::mutate(observed = 1)
  if(is.na(sf::st_crs(df))) df <- sf::st_set_crs(df, value = sf::st_crs(background))

  # Identify spatial cross-validation blocks
  # if(nrow(df) < 50){
  #   # Simply select at random
  #   spc <- data.frame(fold1 = rep("Train", nrow(df)))
  #   spc$fold1[sample(1:nrow(df), size = ceiling(nrow(df) *.25))] <- "Test"
  #   
  # } else {
  #   spc <- spatialBlock(df, speciesCol = "observed", k = 3, iteration = 10)
  # }
  
  # Results containers for saving
  results <- list()
  results_validation <- list()
  # for(mm in c("bart")){ 
  mm = "bart"
  print(mm)
    
    # Now build ibis modelling object for all iterations. Using the zones layer for projections
    basemodel <- distribution(background) %>% 
      add_predictors(env = fullcovs[[vars]], transform = "none")
    
    if(mm == "xgboost") basemodel <- basemodel %>%  engine_xgboost(nrounds = 10000, gamma = 4)
    if(mm == "gdb") basemodel <- basemodel %>% engine_gdb(boosting_iterations = 2500, learning_rate = 0.001)
    if(mm == "bart") basemodel <- basemodel %>% engine_bart(nburn = 50, ntree = 100, chains = 4)
    
    # IUCN range data present for breeding season?
    if(any(c("eu_n2k_habitats") %in% lyrs$name) && grain == 10000){
      ind <- grep("eu_n2k_habitats", names(data))
      for(i in ind){
        # basemodel <- basemodel %>% add_offset_range(layer = data$polygon_iucn, type = "binomial")
        basemodel <- basemodel %>%
          add_predictor_range(layer = data[[i]],
                              method = "distance")
      }
    }
    # Check if any distance ranges are added, then provide prior
    if(length( grep("distance_range", basemodel$get_predictor_names()) )>0){
      vv <- grep("distance_range", basemodel$get_predictor_names(),value = TRUE)
      # Add prior to the range
      if(basemodel$engine$name == "<XGBOOST>"){
        if(length(vv)<2) pp <- XGBPrior(variable = vv, hyper = "increasing") else pp <- XGBPriors(variable = vv, hyper = "increasing")
      } else if(basemodel$engine$name == "<GDB>"){
        if(length(vv)<2) pp <- GDBPrior(variable = vv, hyper = "increasing") else pp <- GDBPriors(variable = vv, hyper = "increasing")
      } else if(basemodel$engine$name == "<BART>"){
        if(length(vv)<2) pp <- BARTPrior(variable = vv, hyper = 0.75) else pp <- BARTPriors(variable = vv, hyper = 0.75)
      }
      basemodel <- basemodel %>% add_priors( priors(pp) )
    }
    
    # Now compute for each fold
    # for(fold in names(spc) ){ # fold = "fold1"
      
      # Subset dataset to training and testing
      # sub_train <- df %>% dplyr::slice(which(spc[,fold]=="Train")) %>% 
      #   add_pseudoabsence(field_occurrence = "observed",template = background,
      #                     settings = abss) %>% unique()
      # sub_test <- df %>% dplyr::slice(which(spc[,fold]=="Test")) %>% 
      #   add_pseudoabsence(field_occurrence = "observed",template = background,
      #                     settings = abss) %>% unique()
      train <- df %>% add_pseudoabsence(field_occurrence = "observed",template = background,
                                                            settings = abss) %>% unique()
      
      # For this exercise we will estimate binomial distributed data
      x <- basemodel %>% 
        add_biodiversity_poipa(train, name = sname, field_occurrence = "observed", docheck = F)
      
      # Add nearest neighbour distance?
      # x <- x %>% add_latent_spatial(method = "nnd")
      
      mod1 <- try({
        train(x, runname = paste0("pnv","_",sname), varsel = "none", rm_corPred = TRUE,
              # bias_variable = paste0("OccurrenceBias__",grain), bias_value = cellStats(x$predictors$data[[paste0("OccurrenceBias__",grain)]],"min"),
              only_linear = ifelse(mm=="xgboost", TRUE, FALSE), verbose = verbose)
      })
      if(inherits(mod1, "try-error")) next()
      
      # Save model summary for later
      suppressWarnings(
        write_summary(mod = mod1,
                      fname = paste0(outdir, "/modelresults/", paste0(sname, "__", paste0(mm)), ".rds"))
      )
      ofname <- paste0(outdir, "/modelresults/", paste0(sname, "__", mm), ".tif")
      write_output(mod1, ofname)

    # Calculate thresholds
    # Use the F1 score for thresholding 
    o <- try({threshold(mod1, method = "F1")}, silent = TRUE)
    o_tru <- try({threshold(mod1, method = "F1", truncate = TRUE)}, silent = TRUE)
    if(inherits(o, "try-error")) {
      message("Thresholding went wrong..")
      # Use quantiles from 0.05 to 0.95 and take the one that maximizes the f1
      o1 <- threshold(mod1, method = "percent",value = .05)
      o2 <- threshold(mod1, method = "percent",value = .50)
      o3 <- threshold(mod1, method = "percent",value = .95)
      m1 <- validate(o1, method = "discrete", point_column = "observed")
      m2 <- validate(o2, method = "discrete", point_column = "observed")
      m3 <- validate(o3, method = "discrete", point_column = "observed")
      # Now pick whichever has highest TSS
      i <- which.max(c(m1$value[m1$metric=="f1"],m2$value[m2$metric=="f1"],m3$value[m3$metric=="f1"]))
      if(i==1) {
        o <- threshold(mod1, method = "percent",value = .05,truncate = FALSE)
        o_tru <- threshold(mod1, method = "percent",value = .05,truncate = TRUE)} else if(i==2) {
          o <- threshold(mod1, method = "percent",value = .5,truncate = FALSE)
          o_tru <- threshold(mod1, method = "percent",value = .5,truncate = TRUE)
        } else {
          o <- threshold(mod1, method = "percent",value = .95,truncate = FALSE)
          o_tru <- threshold(mod1, method = "percent",value = .95,truncate = TRUE)
          }
      }
    
    
    # Create output layers
    ofname <- paste0(outdir, "/", "Threshold__lower_",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(deratify(o$fits[[grep("threshold", names(o$fits))]][[grep("q05", names(o$fits[[grep("threshold", names(o$fits))]]))]],complete = TRUE), ofname, "INT2S")
    ofname <- paste0(outdir, "/", "Threshold__median_",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(deratify(o$fits[[grep("threshold", names(o$fits))]][[grep("q50", names(o$fits[[grep("threshold", names(o$fits))]]))]],complete = TRUE), ofname, "INT2S")
    ofname <- paste0(outdir, "/", "Threshold__upper_",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(deratify(o$fits[[grep("threshold", names(o$fits))]][[grep("q95", names(o$fits[[grep("threshold", names(o$fits))]]))]],complete = TRUE), ofname, "INT2S")
    
    # Truncated thresholds
    ofname <- paste0(outdir, "/", "Threshold__lowerTruncated_",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(o$fits[[grep("threshold", names(o$fits))]][[grep("q05", names(o$fits[[grep("threshold", names(o$fits))]]))]], ofname, "FLT4S")
    ofname <- paste0(outdir, "/", "Threshold__medianTruncated_",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(o$fits[[grep("threshold", names(o$fits))]][[grep("q50", names(o$fits[[grep("threshold", names(o$fits))]]))]], ofname, "FLT4S")
    ofname <- paste0(outdir, "/", "Threshold__upperTruncated_",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(o$fits[[grep("threshold", names(o$fits))]][[grep("q95", names(o$fits[[grep("threshold", names(o$fits))]]))]], ofname, "FLT4S")
    
    ofname <- paste0(outdir, "/", "Prediction__",tools::file_path_sans_ext(basename(hab)), ".tif")
    write_output(o$fits$prediction, ofname, "FLT4S")
    
    # Create csv
    val <- validate(o, method = "disc")
    ofname <- paste0(outdir, "/", "Validation__",tools::file_path_sans_ext(basename(hab)), ".rds")
    saveRDS(val, ofname)
    
    # Delete and clean up
    try({rm(results, results_validation, basemodel, df, data, tr, o, o_tru, mod1)}, silent = TRUE)
    invisible(gc())
}

# Finally calculate the potential restorable land 
message("Finally, calculate potential restorable land.")

# Get the ensemble outputs
ll_pred <- list.files(paste0(path_output, outdirname),recursive = TRUE, full.names = T,ignore.case = T)
ll_pred <- ll_pred[has_extension(ll_pred, "tif")]
ll_pred <- ll_pred[grep("Prediction", ll_pred)]

# Then load the prepared actual land cover
# ll_curr <- list.files(paste0(path_processed, "dataclima/predictors/",grain),full.names = TRUE)
# ll_curr <- ll_curr[has_extension(ll_curr, "tif")]
# ll_curr <- ll_curr[grep("Corine", ll_curr)]

assertthat::assert_that(length(ll_pred)>0
                        # length(ll_curr)>= length(ll_pred)
                        )
# Loop through the predictions and process them
for(f in ll_pred){
  print(basename(f))
  fname <- stringr::str_replace(tools::file_path_sans_ext(basename(f)), "Prediction__", "")
  
  # Get corine layer and convert to fractional
  # ras_cor <- raster(grep(fname, ll_curr,ignore.case = T,value = T)) / 10000
  
  # Load and normalize
  ras <- raster::stack(f)
  ras <- subset(ras, 3:5)
  names(ras) <- c("lower", "median", "upper")
  
  if(all(cellStats(ras,"max")!=1)){
    ras <- predictor_transform(ras, option = "norm")
  }
  # Now get the minimum consensus with current
  # ras <- raster::resample(ras, ras_cor)
  # assertthat::assert_that(compareRaster(ras,ras_cor))
  
  # new_lower <- max(ras$lower, ras_cor)
  # new_median <- max(ras$median, ras_cor)
  # new_upper <- max(ras$upper, ras_cor)
  # 
  # # Difference between lower and upper
  # new_uncert <- raster::clamp(new_upper - new_lower,0,Inf)
  new_lower <- ras$lower
  new_median <- ras$median
  new_upper <- ras$upper
  new_uncert <- new_upper - new_lower
  
  # Write new outputs
  ibis.iSDM:::writeGeoTiff(new_lower, paste0(path_output, outdirname, "/", "PotentialMAES__lower__",fname,".tif"),
               dt = "FLT4S",varNA = -1)
  ibis.iSDM:::writeGeoTiff(new_median, paste0(path_output, outdirname, "/", "PotentialMAES__median__",fname,".tif"),
               dt = "FLT4S",varNA = -1)
  ibis.iSDM:::writeGeoTiff(new_upper, paste0(path_output, outdirname, "/", "PotentialMAES__upper__",fname,".tif"),
               dt = "FLT4S",varNA = -1)
  ibis.iSDM:::writeGeoTiff(new_uncert, paste0(path_output, outdirname, "/", "PotentialMAES__uncert__",fname,".tif"),
               dt = "FLT4S",varNA = -1)
  rm(new_lower, new_median, new_upper, new_uncert)
}
stop("Done!")

# -------------------------- #
#### Recalculate thresholds ####
# Parameter to recalculate thresholds 

if(TRUE == FALSE){
  for(species in rev(speciesdbs)){ 
    
    # Process per directory of files
    d <- basename(dirname(species))
    outdir <- paste0(path_output, outdirname, "/", d)
    
    # Group name
    grp <- basename(dirname(species))
    # Species name
    sname <- tools::file_path_sans_ext(basename(species))
    
    # If file already existing, then recalculate threshold
    if(file.exists(
      paste0(outdir, "/", "EnsembleThreshold__",tools::file_path_sans_ext(basename(species)), ".tif")
    )) {
      # Load raster stack
      o <- raster::stack( paste0(outdir, "/", "Ensemble__",tools::file_path_sans_ext(basename(species)), ".tif") )
      names(o)[1] <- "ensemble_mean"
      myLog("[Starting]", "green", "Now processing ", sname)
      
      # -- #
      # Loading data for the target species
      lyrs <- sf::st_layers(species)
      
      # --- #
      # Now load all cleaned data as list
      data <- list()
      for(name in lyrs$name) data[[name]] <- sf::st_read(species, layer = name, quiet = TRUE) %>% dplyr::mutate(observed = 1)
      assert_that(length(data) > 0)
      # --- #
      # For PNV we simply take the all point information. These will be concatenated later
      df <- data.frame()
      for(l in grep("point", lyrs$name,value = T)){
        df <- rbind(df,
                    data[[l]] %>% dplyr::select(observed) %>% dplyr::mutate(type = l)
        )
      }
      if(nrow(df)==0) next() # No point data found...
      df <- ibis.iSDM:::guess_sf(df)
      # Also if we have Natura 2000 layers, add simulated points (nr presence * 2) to this
      if("polygon_natura2000" %in% lyrs$name){
        a <- sf::st_sample(data$polygon_natura2000,
                           nrow(data$polygon_natura2000) * 2,
                           type = "regular") %>% sf::st_as_sf() %>% 
          ibis.iSDM:::rename_geometry("geom") %>% 
          dplyr::mutate(observed = 1, type = "Natura2000"
                        # x = sf::st_coordinates(.)[,1],
                        # y = sf::st_coordinates(.)[,2]
          )
        df <- rbind(df,a)
      }
      # Same if  BIA are present
      if("polygon_iba" %in% lyrs$name){
        o <- data$polygon_iba #%>% dplyr::filter(season == "breeding")
        if(nrow(o)>1){
          a <- sf::st_sample(o,
                             nrow(o) * 2,
                             type = "regular") %>% sf::st_as_sf() %>% 
            ibis.iSDM:::rename_geometry("geom") %>% 
            dplyr::mutate(observed = 1, type = "IBA"
                          # x = sf::st_coordinates(.)[,1],
                          # y = sf::st_coordinates(.)[,2]
            )
          df <- rbind(df,a)
        }
      } 
      
      # Aggregate observations to grid as to remove duplicates
      df <- ibis.iSDM:::aggregate_observations2grid(df, template = background, field_occurrence = "observed") %>% 
        dplyr::mutate(observed = 1)
      
      # -------------- #
      # New Threshold calculation
      # Load raster stack
      tr <- threshold(o$ensemble_mean, method = "percentile", value = 0.05, poi = df, truncate = F)
      tr <- raster::deratify(tr, complete = TRUE)
      # plot(tr, col = ibis.iSDM:::ibis_colours$distinct_random[5:6])
      
      # Overwrite the existing threshold
      # Create output layer
      ofname <- paste0(outdir, "/", "EnsembleThreshold__",tools::file_path_sans_ext(basename(species)), ".tif")
      write_output(tr, ofname, "INT2S")
      
    }
  }
}
