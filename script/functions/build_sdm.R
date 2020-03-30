# build_sdm 2.0
# Keaton Wilson
# keatonwilson@me.com
# 2020-01-16  

# packages
require(raster)
require(tidyverse)
require(parallel)

# List of operations and associated functions
# 1. Prepping data - prep_data
source("./script/functions/prep_data.R")
# 2. Running blockCV - run_block_cv
source("./script/functions/run_block_cv.R")
# 3. Prepping Data more - prep_data_2
source("./script/functions/prep_data_2.R")
# 4. Training and testing split - train_test_split
source("./script/functions/train_test_split.R")
# 5. Modeling - model_func
source("./script/functions/model_func.R")
# 6. Evaluation plots - eval_plots
source("./script/functions/eval_plots.R")
# 7. Choosing the best model - best_mod
source("./script/functions/best_mod.R")
# 8. Evaluating the best model - evaluate_models
source("./script/functions/evaluate_models.R")
# 9. Extracting arguments from the best model - make_args
source("./script/functions/make_args.R")
# 10. Building the full model on all data - full_model
source("./script/functions/full_model.R")

# importing env rasters into the workspace
bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")

#' build_sdm - master function for building temporally explicit models
#' Goes through data manipulation and year-splitting, blockCV, training and 
#' testing splits, model development and evaluation.
#'
#' @param filename the path of the single-species dataframe as a .rds file
#' @param env_raster_t1 the environmental raster associated with the first 
#' time period
#' @param env_raster_t2 the environmental raster associated with the second 
#' time period
#' @param year_split what year you want to split the data by (e.g. if you 
#' pick 2000, your two groups would years up to 1999, and 2000 to current)
#' @param full_or_minimal parameter for selecting whether you want the full 
#' complement of objects returned, or just the raw data, evaluation metrics, 
#' full models for each time period. 
#'
#' @return a list (either 9 units long (full), or 3 units long (minimal))
#' @export
#'
#' @examples
build_sdm = function(filename, 
                     env_raster_t1, 
                     env_raster_t2,
                     year_split, 
                     full_or_minimal = "full",
                     cores = NULL){
  
  # Setting seed for reproducibility
  set.seed(42)
  
  # Reading in the data
  raw_data = readRDS(filename)
  
  # Prepping Data
    prepped_data = prep_data(data = raw_data, year_split = year_split,
                             env_raster_t1 = env_raster_t1, 
                             env_raster_t2 = env_raster_t2)
    
  # Creating the master list to feed stuff into - nice to have it in one place
    master_list = list("raw_data" = raw_data, "prepped_data" = prepped_data)
    rm(raw_data)
    rm(prepped_data)
    
    
  # Block CV for each time piece 
    block_t1 = try(run_block_cv(prepped_data = master_list[[2]][[1]][[1]], 
                                  bv_raster = master_list[[2]][[2]][[1]]))
    block_t2 = try(run_block_cv(prepped_data = master_list[[2]][[1]][[2]], 
                                  bv_raster = master_list[[2]][[2]][[2]]))
  # writing block objects to data list
    master_list$block_objs = list("block_t1" = block_t1, 
                                  "block_t2" = block_t2)
    rm(block_t1)
    rm(block_t2)
    
  # Second round of data prep - 
    prepped_2_t1 = try(prep_data_2(data = master_list[[2]][[1]][[1]], 
                                    env_raster = master_list[[2]][[2]][[1]]))
    prepped_2_t2 = try(prep_data_2(data = master_list[[2]][[1]][[2]], 
                                   env_raster = master_list[[2]][[2]][[2]]))
    master_list$extra_prepped = list("extra_prepped_t1" = prepped_2_t1, 
                                     "extra_prepped_t2" = prepped_2_t2)
    rm(prepped_2_t1)
    rm(prepped_2_t2)
    
  # Training and testing split  
    training_list_t1 = try(train_test_split(master_list$extra_prepped[[1]],
                                            blocked_obj = master_list$block_objs$block_t1))
    training_list_t2 = try(train_test_split(master_list$extra_prepped[[2]],
                                            blocked_obj = master_list$block_objs$block_t2))
    
  # writing training and test data to master list  
    master_list$train_test = list("train_test_t1" = training_list_t1,
                                  "train_test_t2" = training_list_t2)
    rm(training_list_t1)
    rm(training_list_t2)
    
  # Setting up internal parallelization within a single species for model 
  # building
  if(is.null(cores)){
    total_cores = parallel::detectCores()
    to_use = total_cores - 2
    doParallel::registerDoParallel(to_use)
  } else {
    to_use = cores
    doParallel::registerDoParallel(to_use)
  }
    
  # Modeling
    models_t1 = try(model_func(data = master_list$train_test$train_test_t1$training_data, 
                                 env_raster = master_list$prepped_data$env_data[[1]],
                                 num_cores = to_use))
    models_t2 = try(model_func(data = master_list$train_test$train_test_t2$training_data, 
                                 env_raster = master_list$prepped_data$env_data[[2]],
                                 num_cores = to_use))
  # Writing to master list
    master_list$model_objs = list("models_t1" = models_t1, 
                                    "models_t2" = models_t2)
    rm(models_t1)
    rm(models_t2)
    
  # Model selection
    best_mod_t1 = try(best_mod(model_obj = master_list$model_objs$models_t1))
    best_mod_t2 = try(best_mod(model_obj = master_list$model_objs$models_t2))
    
    master_list$best_mods = list("best_mod_t1" = best_mod_t1, 
                                 "best_mod_t2" = best_mod_t2)
    rm(best_mod_t1)
    rm(best_mod_t2)
    
  # evaluating models on test data
    ev_t1 = try(evaluate_models(test_data = master_list$train_test$train_test_t1$test_data,
                                model = master_list$best_mods$best_mod_t1[[1]],
                                env_raster = master_list$prepped_data$env_data[[1]]))
    ev_t2 = try(evaluate_models(test_data = master_list$train_test$train_test_t2$test_data,
                                model = master_list$best_mods$best_mod_t2[[1]],
                                env_raster = master_list$prepped_data$env_data[[2]]))
  # Writing evaluate objects to master list  
    master_list$eval_objs = list("eval_t1" = ev_t1, 
                                 "eval_t2" = ev_t2)
    rm(ev_t1)
    rm(ev_t2)
    
    
  # Building full models on all data
    full_mod_t1 = try(full_model(models_obj = master_list$model_objs$models_t1,
                                 best_model_index = master_list$best_mods$best_mod_t1[[2]],
                                 full_data = master_list$extra_prepped$extra_prepped_t1,
                                 env_raster = master_list$prepped_data$env_data[[1]]))
    
    full_mod_t2 = try(full_model(models_obj = master_list$model_objs$models_t2,
                                 best_model_index = master_list$best_mods$best_mod_t2[[2]],
                                 full_data = master_list$extra_prepped$extra_prepped_t2,
                                 env_raster = master_list$prepped_data$env_data[[2]]))
  # Writing full model objs to master list  
    master_list$full_mods = list("full_mod_t1" = full_mod_t1, 
                                 "full_mod_t2" = full_mod_t2)
    rm(full_mod_t1)
    rm(full_mod_t2)
    
  # slimming down list if minimal is selected in function argument
    if(full_or_minimal == "minimal"){
      master_list = master_list[[c(1,8,9)]]
    }
    
    return(master_list)
}


