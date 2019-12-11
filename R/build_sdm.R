# Master function to iterate through species and time periods and build SDMs
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(raster)
require(tidyverse)
require(parallel)

# List of operations and associated functions
# 1. Prepping data - prep_data
source("./R/prep_data.R")
# 2. Running blockCV - run_block_cv
source("./R/run_block_cv.R")
# 3. Prepping Data more - prep_data_2
source("./R/prep_data_2.R.R")
# 4. Training and testing split - train_test_split
source("./R/train_test_split.R")
# 5. Modeling - model_func
source("./R/model_func.R")
# 6. Evaluation plots - eval_plots
source("./R/eval_plots.R")
# 7. Choosing the best model - best_mod
source("./R/best_mod.R")
# 8. Evaluating the best model - evaluate_models
source("./R/evaluate_models.R")
# 9. Extracting arguments from the best model - make_args
source("./R/make_args.R")
# 10. Building the full model on all data - full_model
source("./R/full_model.R")

# Creating a new folder to house individual lists/stuff in
dir.create("./output/")

# importing env rasters into the workspace
bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")

# Full master function - calls other functions in the pipieline and steps through
# the whole analysis of multiple species over the two time periods. 

build_sdm = function(multi_species_df, year_split, env_raster_t1, env_raster_t2, num_cores = NULL){
  # Setting seed for reproducibility
    set.seed(42)
  
  # QC to make sure we have enough records for each species
  data_summary = multi_species_df %>%
    mutate(year = lubridate::year(date)) %>%
    filter(!is.na(year)) %>%
    mutate(time_frame = ifelse(year > 1999, "T1", "T2")) %>%
    group_by(name, time_frame) %>%
    summarize(n = n()) %>%
    print(n = 50)
  
  invisible(readline(prompt="Do this data summary look reasonable?
Check to make sure you have enough data for each species to continue.
Press [enter] to continue"))
  
  
  # split multi-species dataframe into a list
  butt_list = split(multi_species_df, f = multi_species_df$name)
  
  #Writing each component of the list out to file
  files = c()
  for(i in 1:length(names(butt_list))) {
    files[i] = paste0("./output/", names(butt_list)[i], ".rds")
  }
  
  for(j in 1:length(names(butt_list))){
    saveRDS(butt_list[[j]], files[j])
  }
  
  #Paralell stuff
  num_cores = detectCores() - 2
  
  # Iterating the prep_data function over the list of species dataframes
  # This overwrites original data to a list that incldues original data and new
  # objects
  file_list = as.list(list.files("./output", full.names=TRUE))
  
  prep_over_list = function(file_name){
    data = readRDS(file_name)
    prepped_data = prep_data(data = data, year_split = year_split,
                             env_raster_t1 = env_raster_t1, 
                             env_raster_t2 = env_raster_t2)
    data_list = list(data, prepped_data)
    saveRDS(data_list, file_name)
  }
  
  #using multi-core lapply 
  mclapply(file_list, prep_over_list, mc.cores = num_cores)  
  
  
  # Generating blockCV objects for each time period for each species and attaching to master list
  
  block_over_list = function(filename){
    data = readRDS(file_name)
    block_t1 = try(run_block_cv(prepped_data = data[[2]][[1]][[1]], 
                                bv_raster = data[[2]][[2]][[1]]))
    block_t2 = try(run_block_cv(prepped_data = data[[2]][[1]][[2]], 
                                bv_raster = data[[2]][[2]][[2]]))
    data[[2]]$blocks = list(block_t1, block_t2)
    saveRDS(data, file_name)
  }
  
  #using multi-core lapply
  mclapply(file_list, block_over_list, mc.cores = num_cores)
  
  # Running second prepping data function
  
  for(j in 1:length(prepped_data_list)){
    prepped_df_t1 = try(prep_data_2(data = prepped_data_list[[j]]$data[[1]], 
                                    env_raster = prepped_data_list[[j]]$env_data[[1]]))
    prepped_df_t2 = try(prep_data_2(data = prepped_data_list[[j]]$data[[2]], 
                                    env_raster = prepped_data_list[[j]]$env_data[[2]]))
    data_df = list(prepped_df_t1, prepped_df_t2)
    prepped_data_list[[j]]$prepped_dfs = data_df
  }
  
  # Training and test split
  for(j in 1:length(prepped_data_list)){
    training_list_t1 = try(train_test_split(prepped_data_list[[j]]$prepped_dfs[[1]],
                                            blocked_obj = prepped_data_list[[j]]$blocks[[1]]))
    training_list_t2 = try(train_test_split(prepped_data_list[[j]]$prepped_dfs[[2]],
                                            blocked_obj = prepped_data_list[[j]]$blocks[[2]]))
    prepped_data_list[[j]]$train_test_split[[1]] = training_list_t1
    prepped_data_list[[j]]$train_test_split[[2]] = training_list_t2
    
  }
  
  # parallelization
  total_cores = parallel::detectCores()
  to_use = total_cores - 2
  doParallel::registerDoParallel(to_use)
  
  # Modeling
  for(k in 1:length(prepped_data_list)){
    models_t1 = try(model_func(data = prepped_data_list[[k]]$train_test_split[[1]]$training_data, 
                               env_raster = prepped_data_list[[k]]$env_data[[1]],
                               num_cores = to_use))
    models_t2 = try(model_func(data = prepped_data_list[[k]]$train_test_split[[2]]$training_data, 
                               env_raster = prepped_data_list[[k]]$env_data[[2]],
                               num_cores = to_use))
    prepped_data_list[[k]]$models = list(models_t1, models_t2)
  }
  
  # Choosing best model
  for(l in 1:length(prepped_data_list)){
    best_mod_t1 = try(best_mod(model_obj = prepped_data_list[[l]]$models[[1]]))
    best_mod_t2 = try(best_mod(model_obj = prepped_data_list[[l]]$models[[2]]))
    prepped_data_list[[l]]$best_mod = list(best_mod_t1, best_mod_t2)
  }
  
  # evaluating best model on blockCV test data
  for(m in 1:length(prepped_data_list)){
    ev_t1 = try(evaluate_models(test_data = prepped_data_list[[m]]$train_test_split[[1]]$test_data,
                                model = prepped_data_list[[m]]$best_mod[[1]][[1]],
                                env_raster = prepped_data_list[[m]]$env_data[[1]]))
    ev_t2 = try(evaluate_models(test_data = prepped_data_list[[m]]$train_test_split[[2]]$test_data,
                                model = prepped_data_list[[m]]$best_mod[[2]][[1]],
                                env_raster = prepped_data_list[[m]]$env_data[[2]]))
    
    prepped_data_list[[m]]$evaluations = list(ev_t1, ev_t2)
  }
  
  # full mods on all data
  for(n in 1:length(prepped_data_list)){
    full_mod_t1 = try(full_model(models_obj = prepped_data_list[[n]]$models[[1]],
                                 best_model_index = prepped_data_list[[n]]$best_mod[[1]][[2]],
                                 full_data = prepped_data_list[[n]]$prepped_dfs[[1]],
                                 env_raster = prepped_data_list[[n]]$env_data[[1]]
    ))
    
    full_mod_t2 = try(full_model(models_obj = prepped_data_list[[n]]$models[[2]],
                                 best_model_index = prepped_data_list[[n]]$best_mod[[2]][[2]],
                                 full_data = prepped_data_list[[n]]$prepped_dfs[[2]],
                                 env_raster = prepped_data_list[[n]]$env_data[[2]]))
    
    prepped_data_list[[n]]$full_mods = list(full_mod_t1, full_mod_t2)
    
  }
  return(prepped_data_list)
  saveRDS(prepped_data_list, "./data/model_data_list.rds")
  doParallel::stopImplicitCluster()
}


# Testing

#loading in multi-species

test_df = read_csv("./data/candidate_occurences.csv")

data_summary = multi_species_df %>%
  mutate(year = lubridate::year(date)) %>%
  filter(!is.na(year)) %>%
  mutate(time_frame = ifelse(year > 1999, "T1", "T2")) %>%
  group_by(name, time_frame) %>%
  summarize(n = n()) %>%
  print(n = 50)

invisible(readline(prompt="Do this data summary look reasonable?
Check to make sure you have enough data for each species to continue.
Press [enter] to continue"))


# split multi-species dataframe into a list
butt_list = split(multi_species_df, f = multi_species_df$name)

#Writing each component of the list out to file
files = c()
for(i in 1:length(names(butt_list))) {
  files[i] = paste0("./output/", names(butt_list)[i], ".rds")
}

for(j in 1:length(names(butt_list))){
  saveRDS(butt_list[[j]], files[j])
}

#Paralell stuff
num_cores = detectCores() - 2

# Iterating the prep_data function over the list of species dataframes
file_list = as.list(list.files("./output", full.names=TRUE))

prep_over_list = function(file_name){
  data = readRDS(file_name)
  prepped_data = prep_data(data = data, year_split = 2000,
                           env_raster_t1 = bv_t1, 
                           env_raster_t2 = bv_t2)
  data_list = list(data, prepped_data)
  saveRDS(data_list, file_name)
}

mclapply(file_list, prep_over_list, mc.cores = 2)  


block_over_list = function(file_name){
  data = readRDS(file_name)
  block_t1 = try(run_block_cv(prepped_data = data[[2]][[1]][[1]], 
                              bv_raster = data[[2]][[2]][[1]]))
  block_t2 = try(run_block_cv(prepped_data = data[[2]][[1]][[2]], 
                              bv_raster = data[[2]][[2]][[2]]))
  data[[2]]$blocks = list(block_t1, block_t2)
  saveRDS(data, file_name)
}

#using multi-core lapply
mclapply(file_list, block_over_list, mc.cores = num_cores)

### EDITS START BELOW HERE
