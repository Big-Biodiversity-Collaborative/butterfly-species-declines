# Modeling
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(tidyverse)
require(ENMeval)

#' Modeling using Maxent
#' 
#' Function that takes training data, environmental rasters and information 
#' about parallelization and builds and tunes SDMs. Will try maxnet algorithm #' first, and if that fails, the old maxent.jar
#' 
#' @param data training data generated from the \code{\link{train_test_split}} function. 
#' @param env_raster the cropped env_raster for a given species and time period. 
#' @param num_cores the number of cores to use in parallelization
#'  
#' @return an ENMeval object - contains information on all of the modeling tuning and testing (internal cross validation)
#' 
#' @examples

model_func = function(data = NULL, env_raster, num_cores = to_use) {
  data_occ = data %>%  #Generating occurence lat long
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude) %>%
    drop_na()
  
  bg_data = data %>% #Generating background lat long
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude) %>%
    drop_na()
  
  #Running the model
  eval = try(ENMevaluate(occ = data_occ, 
                         bg.coords = bg_data,
                         env = env_raster,
                         method = 'randomkfold', 
                         kfolds = 5,
                         parallel = TRUE,
                         numCores = num_cores,
                         algorithm = 'maxnet'))
  
  if (class(eval) == "try-error") {
    cat("Caught an error running maxnet, trying maxent")
    
    eval = try(ENMevaluate(occ = data_occ, 
                           bg.coords = bg_data,
                           env = env_raster,
                           method = 'randomkfold', 
                           kfolds = 5,
                           # parallel = TRUE,
                           # numCores = num_cores,
                           algorithm = 'maxent.jar'))
  }
  return(eval)
}