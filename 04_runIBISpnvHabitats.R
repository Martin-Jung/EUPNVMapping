# --- --- #
# This script runs the EUPNV predictions by loading all prepared variables in
# then computing per class and using BART as predictive framework.
# - The predicted probability (low, median, upper) and threshold for the current covariates
# - Predict per SSP and GCM
# - Also make a mean ensemble for each
# - Also run a cross-validation to evaluate model predictive performance on current variables (spatial block, 3). Save results as csv
# --- --- #

# Load in functions and packages 
source("00_Functions.R")
# Source in parameter settings
source("01_ParameterSettings.R")

# Set run cores for parallization
options('ibis.nthread' = cores)
options('ibis.runparallel' = ifelse(cores>1, FALSE, FALSE) ) # Set to FALSE for slower, but memory proof runs

# -------------------- #
#### Prepare files for modelling ####

# Get reference raster for the given grain
background <- terra::rast( paste0(path_background, "background_ref.tif") )

# For future check variables and buld a dataframe to load them
ll <- Reduce("c", sapply(path_future, function(z) list.files(z, recursive = TRUE, full.names = TRUE)) )
ll <- ll[has_extension(ll, "tif")]
df_future <- data.frame(ifname = ll, gcm = basename(dirname(ll)), ssp = basename(dirname(dirname(ll)))) |> 
  dplyr::mutate(period = stringr::str_split(basename(ifname), "_",simplify = TRUE)[,3],
                varname = paste0(stringr::str_split(basename(ll), "_",simplify = TRUE)[,1],
                                 "_",
                                 stringr::str_split(basename(ll), "_",simplify = TRUE)[,2]
                )
  )
assertthat::assert_that(all( c("ssp126","ssp370","ssp585") %in% df_future$ssp ))
# Create a time variable per period
df_future$time <- ifelse(df_future$period == "2011-2040", "2040-01-01",
                         ifelse(df_future$period == "2041-2070", "2070-01-01",
                                ifelse(df_future$period == "2071-2100", "2100-01-01", NA))
                         ) |> as.Date()

# Filter to target gcms only
df_future <- df_future |> dplyr::filter(gcm %in% gcms)

# 46 per GCM|SSP|PERIOD
df_future |> dplyr::group_by(gcm, ssp, period) |> summarise(n = n()) |> as.data.frame()

# Get selected static covariates
ll <- list.files(paste0(path_presentcovs),"tif", full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
# Get Copernicus BIOCLIM variables
# ll <- ll[grep("aridity|BIO|era5|potential|DistanceCoast|volumetric|growing|frost-days|EUDEM|cloud-cover|GroundWater|Lithology",ll)]
# Get only projected variables
ll <- ll[grep("min|max|range", ll,invert = TRUE)]
ll <- ll[grep(paste0(unique(df_future$varname), collapse= "|"), ll)]
assertthat::assert_that(length(ll) == dplyr::n_distinct(df_future$varname))

prescovs <- terra::rast(ll)
# Aggregate/resample to background
# prescovs <- terra::resample(prescovs, background)
assertthat::assert_that(terra::compareGeom(prescovs, background))
# Rename the variables
names(prescovs) <- stringr::str_remove(names(prescovs), "_1981-2010")
names(prescovs) <- stringr::str_remove(names(prescovs), "_V.2.1")
assertthat::assert_that(anyDuplicated(names(prescovs))==0)
# For snow (SWE) correct by setting all NA to 0
prescovs$CHELSA_swe[is.na(prescovs$CHELSA_swe)] <- 0
# Same for gst
prescovs$CHELSA_gst[is.na(prescovs$CHELSA_gst)] <- 0
# Check variables are all in future too
assertthat::assert_that( all( names(prescovs) %in% df_future$varname))

# Load full stack of all PNV covariates
ll <- list.files(paste0(path_pnvcovs),"tif", full.names = TRUE)
ll <- ll[has_extension(ll, "tif")]
# For this exercise we rely specifically on 
ll <- ll[grep("DistanceCoast|EUDEM|EuroVegMap|Lithology|Tootchi",ll)]
pnvcovs <- terra::rast(ll)
# Align with prescovs
# pnvcovs <- terra::resample(pnvcovs, prescovs)
assertthat::assert_that( terra::compareGeom(pnvcovs, prescovs, stopOnError = FALSE) )

# Clip to nuts if set
if(clip_nuts){
  # Clip with NUTS (excluding Turkey in the process)
  nuts <- sf::st_read(path_nuts) %>% dataclima:::laea_transform() %>%
    dplyr::filter(CNTR_CODE != "TR") %>%
    dplyr::filter(CNTR_CODE != "NO", CNTR_CODE != "IS")
  background <- dataclima:::clip_nuts(background, nuts)
  background <- terra::mask(background, nuts)
}

# Get all habitat databases
habitatdbs <- list.files(path = path_habitat, full.names = TRUE, recursive = TRUE)
habitatdbs <- habitatdbs[has_extension(habitatdbs, "gpkg")]
habitatdbs <- habitatdbs[grep("Cropland|Urban|River|Marine.gpkg", habitatdbs,invert = TRUE)]
assertthat::assert_that(length(habitatdbs)>0)

# Now align all covariates to be sure
if(file.exists(paste0(dirname(path_presentcovs),"/pnvcovariates.tif"))){
  fullcovs <- terra::rast(paste0(dirname(path_presentcovs),"/pnvcovariates.tif"))
  state <- readRDS(paste0(dirname(path_presentcovs),"/pnvcovariates.rds"))
} else {
  fullcovs <- c(prescovs, pnvcovs)
  fullcovs[is.na(fullcovs)] <- 0 # Replace all with
  fullcovs <- terra::crop(fullcovs, background)
  fullcovs <- terra::mask(fullcovs, background)
  fullcovs <- ibis.iSDM:::predictor_homogenize_na(fullcovs)
  
  # Scale all variables except the factor ones for unit consistency
  plot(fullcovs[[1]]) # Plot for visual cross-check
  # names(fullcovs)
  new <- predictor_transform(
    fullcovs[[!is.factor(fullcovs)]], option = "scale"
  )
  
  # Get the state variable here from the transformation
  state <- attr(new, "transform_params")
  
  # Rasterize the lithology and koeppen class
  new <- c(new, fullcovs[[is.factor(fullcovs)]] )
  fullcovs <- new
  writeRaster(fullcovs, paste0(dirname(path_presentcovs),"/pnvcovariates.tif"))
  saveRDS(state, paste0(dirname(path_presentcovs),"/pnvcovariates.rds"))
  rm(new); gc()
}
assertthat::assert_that(terra::compareGeom(fullcovs, background),
                        nrow(state)>1)

# Match PNV Maes codes layers to IUCN codes
pnv_legends <- data.frame(layer = grep("VegMap|Wetland", names(pnvcovs),value = TRUE),
                          IUCN.code = as.character(c(1, 4, 4, 3, 6, 5, 5)) )
pnv_legends <- dplyr::left_join(pnv_legends, 
                                crosswalk_iucnmaes %>% dplyr::filter(IUCN.level==1) %>%  dplyr::select(IUCN.code, MAES.code, MAES.name)
)

# Default pseudo-absence information if not available elsewhere.
# Sample absence points randomly equal to the number of presence points
abss <- pseudoabs_settings(background = background, nrpoints = 10000, min_ratio = 1, method = "random",
                           buffer_distance = 10000)

assertthat::assert_that(
  length(habitatdbs)>0
)

# -------------------- #
#### Execute prediction sequential ####

mm = "bart" # Stays the same. We use BART for estimation throughout
doFuture <- T # Future projections
doValidation <- F # Validation run to assess model performance
testing <- F # Use aggregated 10km data instead
doPrediction <- FALSE
## Parallel code
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
#proc <- foreach(species = iter(sort(speciesdbs)),
#        .errorhandling = 'stop',
#        .inorder = TRUE,
#        .packages = c('sf','dataclima', "ibis.iSDM", "assertthat", "raster", "modEvA", "sf"),
#        .export = c('path_output','fullcovs')
#        ) %dopar% {
for(hab in sort(habitatdbs)){ # hab <- habitatdbs[3] 
  
  # Habitat name
  sname <- tools::file_path_sans_ext(basename(hab))

  message("Processing ", sname)  
  
  # TODO:
  # Sate variable only for transformed, thus merge of different predictors.
  # Make test prediction 
  # Future projection per SSP and class

  # Process per directory of files
  outdir <- paste0(path_output, "/", sname )
  dir.create(outdir, showWarnings = FALSE)
  assertthat::assert_that(dir.exists(outdir))
  # Also create subfolders for individual runs and ensembles in there
  dir.create(paste0(outdir, "/modelresults"),showWarnings = FALSE)
  
  # Variables to use
  vars <- names(fullcovs)
  assertthat::assert_that(length(vars)>0)
  
  # Use different ones depending on the type
  # For forests only use treeclim projected variables (fgd, gsl, gsp, gst, lgd, scd), https://doi.org/10.1007/s00035-014-0124-0
  if(sname == "Wetlands"){
    vars <- grep("fgd|gsl|gsp|gst|lgd|scd", vars, value = TRUE,invert = TRUE)
    vars <- grep("Forest|Grassland|Shrub|Sparse", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  } else if( sname == "Woodland.and.forest"){
    # vars <- grep("kg1|kg2|kg3|kg4|kg5", vars, value = TRUE,invert = TRUE)
    vars <- grep("Wetland|Grassland|Shrub|Sparse", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  } else if( sname == "Heathland.and.shrub") {
    vars <- grep("fgd|gsl|gsp|gst|lgd|scd", vars, value = TRUE,invert = TRUE)
    vars <- grep("Wetland|Grassland|Forest|Sparse", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  } else if( sname == "Marine.inlets.and.transitional.waters"){
    vars <- grep("fgd|gsl|gsp|gst|lgd|scd", vars, value = TRUE,invert = TRUE)
    vars <- grep("Wetland|Grassland|Forest|Sparse|Shrub", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  } else if( sname == "Sparsely.vegetated.areas"){
    vars <- grep("fgd|gsl|gsp|gst|lgd|scd", vars, value = TRUE,invert = TRUE)
    vars <- grep("Wetland|Grassland|Shrub|Forest", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  } else if( sname == "Grassland"){
    vars <- grep("fgd|gsl|gsp|gst|lgd|scd", vars, value = TRUE,invert = TRUE)
    vars <- grep("Wetland|Shrub|Forest|Sparse", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  } else {
    vars <- grep("fgd|gsl|gsp|gst|lgd|scd", vars, value = TRUE,invert = TRUE)
    vars <- grep("Wetland|Grassland|Shrub|Forest|Sparse", vars, value = TRUE, invert = TRUE, ignore.case = TRUE)
  }
  # Subset variables
  covs <- fullcovs[[vars]]
  sstate <- state[,colnames(state) %in% vars] # Same for state
  
  # ---
  # Check if existing, otherwise skip
  tex <- c( paste0(outdir, "/", "Model__",sname, ".tif"),
            paste0(outdir, "/", "Prediction__",sname, ".tif"),
            paste0(outdir, "/", "Threshold__",sname, ".tif")
  )
  if(all(file.exists(tex))) next()
  # ---

  if(testing){
    covs <- terra::aggregate(covs, fact = 10)
    background <- terra::aggregate(background, fact = 10, fun = terra::modal)
  }
  
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
  df <- ibis.iSDM:::aggregate_observations2grid(df, template = background, field_occurrence = "observed") %>% 
    dplyr::mutate(observed = 1)
  if(is.na(sf::st_crs(df))) df <- sf::st_set_crs(df, value = sf::st_crs(background))
  df <- df |> sf::st_transform(df, crs = sf::st_crs(background))

  message("Running ",mm, " model for ", sname, " (", sum(df$observed), " points)")
    
  # Now build ibis modelling object for all iterations. Using the zones layer for projections
  basemodel <- distribution(background) %>% 
    add_predictors(env = covs, transform = "none",explode_factors = TRUE)
  
  # if(mm == "xgboost") basemodel <- basemodel %>%  engine_xgboost(nrounds = 10000, gamma = 4)
  # if(mm == "gdb") basemodel <- basemodel %>% engine_gdb(boosting_iterations = 2500, learning_rate = 0.001)
  if(mm == "bart") basemodel <- basemodel %>% engine_bart(iter = 1000,nburn = 250)
  
  # N2k data present, add simulated points
  # Now done further above
  # if(any(c("eu_n2k_habitats") %in% lyrs$name)){
  #   ind <- grep("eu_n2k_habitats", names(data))
  #   for(i in ind){
  #     # basemodel <- basemodel %>% add_offset_range(layer = data$polygon_iucn, type = "binomial")
  #     # basemodel <- basemodel %>% add_predictor_range(layer = data[[i]], method = "distance")
  #     # FIX: Removed owing to biases in certain areas
  #   }
  # }
  
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
  
  # Add other priors related to pnv_legends
  if(any(pnv_legends$layer %in% basemodel$get_predictor_names())){
    vv <- grep(paste0(pnv_legends$layer,collapse = "|"), basemodel$get_predictor_names(),value = TRUE)
    if(basemodel$engine$name == "<BART>"){
      if(length(vv)<2) pp <- BARTPrior(variable = vv, hyper = 0.75) else pp <- BARTPriors(variable = vv, hyper = 0.75)
    }
    basemodel <- basemodel %>% add_priors( priors(pp) )
  }
  
  # Set pseudo-absence data to full model
  train <- df %>% add_pseudoabsence(field_occurrence = "observed",template = background,
                                                        settings = abss) %>% unique()
  
  # For this exercise we will estimate binomial distributed data
  x <- basemodel %>% 
    add_biodiversity_poipa(train, name = sname, field_occurrence = "observed", docheck = F)
  
  ofname <- paste0(outdir, "/", "Model__",sname, ".rds")
  # if(file.exists(ofname)){
  #   mod1 <- load_model(ofname)
  # } else {
    mod1 <- try({
      train(x, runname = paste0("pnv","_",sname), filter_predictors = "pear",
            only_linear = FALSE, verbose = verbose)
    })
  # }
  if(inherits(mod1, "try-error")) next() # Overall model fitting failed
  
  # Resave model
  if(doPrediction){
    # Save model summary for later
    suppressWarnings(
      write_summary(mod = mod1,
                    fname = paste0(outdir, "/modelresults/", paste0(sname, "__", paste0(mm)), ".rds"))
    )
  
    # Calculate thresholds based on a simple 5% percentile
    mod1 <- try({threshold(mod1, method = "percentile",value = 0.05)}, silent = TRUE)
    if(inherits(mod1, "try-error")) {
      message("Thresholding went wrong..")
    }
    ofname <- paste0(outdir, "/", "Model__",sname, ".rds")
    write_model(mod1, ofname)
    
    # Create output layers #
    # Prediction
    ofname <- paste0(outdir, "/", "Prediction__",sname, ".tif")
    write_output(mod1$get_data(), ofname, "FLT4S")
    
    # Thresholds
    ofname <- paste0(outdir, "/", "Threshold__",sname, ".tif")
    write_output(mod1$get_data("threshold_percentile")[[c(1,3,4,5)]], ofname, "INT2S")
  } else {
    # But check threshold and set if necessary
    if(is.Waiver(mod1$get_thresholdvalue())){
      mod1 <- try({threshold(mod1, method = "percentile",value = 0.05)}, silent = TRUE)
    }
  }
  
  #### Scenario projection ----
  if(doFuture){
    # Load future covariates 
    sub <- df_future |> dplyr::filter(varname %in% basemodel$get_predictor_names())
    
    # Now for each SSP and gcm combo
    for(s in unique(sub$ssp)){ # s = unique(sub$ssp)[1]
      subs <- subset(sub, ssp == s)
      for(g in unique(subs$gcm)){ # g = unique(sub$gcm)[1]
        subss <- subset(subs, gcm == g)
        message(s, " - ", g)
        ofname <- paste0(outdir, "/", "Projection__",sname, "__", s, "__", g ,".nc")
        if(file.exists(ofname)) next()
        
        # Now load the stacks per timeframe
        ras <- terra::rast(subss$ifname)
        if(testing){
          ras <- terra::aggregate(ras, fact = 10)
        }
        if(!terra::compareGeom(ras, basemodel$predictors$get_data(),stopOnError = FALSE)){
          ras <- terra::resample(ras, basemodel$predictors$get_data())
        }
        terra::time(ras) <- subss$time
        names(ras) <- subss$varname
        ras <- terra::subst(ras, NA, 0) # Replace all NA with 0
        # Mask again
        # if(!terra::compareGeom(ras, background, stopOnError = FALSE)){
        #   ras <- terra::resample(ras, background)
        # }
        ras <- terra::mask(ras, basemodel$background)
        # # For snow (SWE) correct by setting all NA to 0
        # ras[[grep("CHELSA_swe", names(ras))]] <- terra::subst(ras[[grep("CHELSA_swe", names(ras))]],NA,0)
        # # Same for gst
        # if(length(grep("CHELSA_gst", names(ras)))>0){
        #   ras[[grep("CHELSA_gst", names(ras))]] <- terra::subst(ras[[grep("CHELSA_gst", names(ras))]],NA,0)
        # }

        # Convert to stars
        # fut <- list()
        # for(n in unique(names(ras))){
        #   fut[[n]] <- stars::st_as_stars(ras[[which(names(ras)==n)]]) 
        # }
        # fut <- Reduce('c', fut)
        # # Transform the variables by scale depending on state
        # fut2 <- predictor_transform(fut[which(names(fut) %in% colnames(state))],
        #                     option = "scale",state = state[,which(colnames(state) %in% names(fut))])
        # fut <- c(fut2, fut[which(!names(fut) %in% colnames(state))] )
        # # --- #
        # try({ rm(ras, fut2) }, silent = TRUE)
        # 
        # # Now add all other vars that are missing from the model object
        # miss <- basemodel$get_predictor_names()[which(!(basemodel$get_predictor_names() %in% names(fut)))]
        # if(length(miss)>0){
        #   fut <- ibis.iSDM:::st_add_raster(fut, basemodel$predictors$get_data()[[miss]])
        # }
        ## Alternative implementation per timeslot
        for(p in unique(terra::time(ras))){ # p = unique(terra::time(ras))[1]
          message("Modelling time period:", p)
          
          fut <- ras[[which(terra::time(ras)==p)]]
          # Transform
          fut2 <- predictor_transform(fut[[which(names(fut) %in% colnames(state))]],
                              option = "scale",state = state[,which(colnames(state) %in% names(fut))])
          fut <- c(fut2, fut[[which(!names(fut) %in% colnames(state))]] )
          rm(fut2)
          # Now add all other vars that are missing from the model object
          miss <- basemodel$get_predictor_names()[which(!(basemodel$get_predictor_names() %in% names(fut)))]
          if(length(miss)>0) fut <- c(fut, basemodel$predictors$get_data()[[miss]])
          # Add time dimension again
          terra::time(fut) <- rep(as.Date(p), terra::nlyr(fut))
          # Security check
          assertthat::assert_that(all(names(fut) %in% basemodel$get_predictor_names() ))
          
          message("Starting with scenario...")
          # ---- #
          # Now create the projection model
          sc <- scenario(mod1) |> 
            threshold() |> 
            add_predictors(fut,transform = "none",explode_factors = TRUE)
          proj1 <- try({ project(sc, layer = "q05") })
          if(inherits(proj1, "try-error")){ warning("Projection failed..."); next() }
          proj2 <- try({ project(sc, layer = "q50") })
          proj3 <- try({ project(sc, layer = "q95") })
          proj4 <- try({ project(sc, layer = "mean") })
          # ---- #
          # Combine all
          p1 <- proj1$get_data()
          names(p1) <- c("suitability_q05", "threshold_q05")
          p2 <- proj2$get_data()
          names(p2) <- c("suitability_q50", "threshold_q50")
          p3 <- proj3$get_data()
          names(p3) <- c("suitability_q95", "threshold_q95")
          p4 <- proj4$get_data()
          names(p4) <- c("suitability_mean", "threshold_mean")
          
          # Save the outputs
          message("Saving the outputs...")
          ofname <- paste0(outdir, "/", "Projection__",sname, "__", p ,"__", s, "__", g ,".nc")
          write_output(c(p1,p2,p3,p4), ofname)
          try({ rm(proj1, sc, proj2, proj3, proj4, p1,p2,p3,p4) },silent = TRUE)
          gc()
        } # End of period loop
      }
    }
    # --- #
    message("Projection done for ", sname)
  }
  
  #### Validation ----
  if(doValidation){
    
    # Get training data from object
    spc <- spatialsample::spatial_block_cv(data = train, repeats = 2, method = "random", v = 3)

    # Now compute for each fold
    for(fold in unique(spc$id2) ){ # fold = "fold1"
      message("Computing ", fold)
     
      # Get testing data
      for(r in unique(spc$id[spc$id2==fold])){
        sub <- spc |> dplyr::filter(id2 == fold, id == r)
        
        train <- spatialsample::analysis(sub$splits[[1]])
        test <- spatialsample::assessment(sub$splits[[1]])
        
        x <- basemodel %>% 
          add_biodiversity_poipa(train, name = sname, field_occurrence = "observed", docheck = F)
        
        # Train with slighly fewer iterations
        mod1 <- try({
          train(x |> engine_bart(iter = 100,nburn = 25), runname = paste0("cv","_",sname), filter_predictors = "pear",
                only_linear = FALSE, verbose = verbose)
        })
        
        # Threshold with independent data 
        mod1 <- threshold(mod1, method = strategy_cvnr[['metric']], point = test,field_occurrence = 'observed')
        # Validate as well
        val <- validate(mod1, method = "discrete", point = test, field_occurrence = 'observed')
        # Also append the threshold data
        val <- dplyr::bind_rows(val, data.frame(metric = "threshold", value = mod1$get_thresholdvalue() ))
        
        ofname <- paste0(outdir, "/", "Validation__",sname, "__", paste0(fold,"-",r), ".rds")
        saveRDS(val, ofname)
        # Clean up
        try({ rm(sub, mod1, x) })
      }
      
    }
    
    message("Cross-validation done for ", sname)
  }
  
  # Delete and clean up
  try({rm(mod1, ofname, basemodel)})
  invisible(gc())
  
}

stop("DONE!")
