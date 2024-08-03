# --- --- #
# This script creates a formatted table with the prediction accuracy estimates
# for each class (and model)
# --- --- #

# Load in functions and packages 
source("00_Functions.R")

library(purrr)
library(gt)
library(gtExtras)

# Path to the create (depends on system)
path_output <- "PNVHabitats__ClimateRun/"
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
#### Load prediction accuracy estimates ####
# Load all the varying predictions
ll <- list.files(path_output, full.names = TRUE,recursive = TRUE)
ll <- ll[has_extension(ll, "rds")]
ll <- ll[grep("Fold",ll)]
assert_that(length(ll)>0)

# Load them in
vals <- ll |> map_dfr( function(z) readRDS(z) |> dplyr::mutate(filename = tools::file_path_sans_ext(basename(z))) ) |> 
  separate(filename, c("name", "variable", "Fold"), sep = "__")

# Format factor
vals$variable <- factor(vals$variable, levels = names(cols))

# Filter to relevant metric and summarize across folds
vals <- vals |> dplyr::filter(metric %in% c("f1","tss")) |> 
  dplyr::group_by(variable, metric) |> 
  dplyr::summarise(y = ggplot2::mean_se(value)[,'y'],
                   ymin = ggplot2::mean_se(value)[,'ymin'],
                   ymax = ggplot2::mean_se(value)[,'ymax'],
                   value_sd = sd(value))

# Build plot
gg <- vals |> dplyr::select(-ymin,-ymax) |> 
  gt(groupname_col = 'variable') |> 
  tab_header(
    title = 'Prediction performance',
    subtitle = 'Spatial blocks crossvalidation \n(3 blocks and 2 repeats each)'
  ) |> 
  # cols_label(.list = metric) |> 
  fmt_number(columns = where(is.numeric), decimals = 2) |> 
  gt_theme_pff() #|> 
  # gt_color_rows(
  #   columns = c(y,value_sd), 
  #   domain = c(0, 1)
  # )
gg
gt::gtsave(as_word(gg), "figures/SI_3_AccuracyTable.html")
