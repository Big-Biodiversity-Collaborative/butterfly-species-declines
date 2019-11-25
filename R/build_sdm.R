# Master function to iterate through species and time periods and build SDMs
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(raster)
require(tidyverse)
require(parallel)


# Full master function - calls other functions in the pipieline and steps through
# the whole analysis of multiple species over the two time periods. 

build_sdm = function(multi_species_df, year_split, env_raster_t1, env_raster_t2, num_cores = NULL){
  # Setting seed for reproducibility
  
  set.seed(42)
  
  # QC to make sure we have enough records for each species
  data_summary = multi_species_df %>%
    filter(!is.na(year)) %>%
    mutate(time_frame = ifelse(year > 1999, "T1", "T2")) %>%
    group_by(true_name, time_frame) %>%
    summarize(n = n()) 
  
  if(any(data_summary$n < 10)) {
    
    offender_list = data_summary %>%
      filter(n <= 10) %>%
      select(true_name) %>%
      pull()
    
    print(paste("Removing species (", offender_list,") as these species have less than 10 records for a given period"))
    
    # Defining not-in funcion
    '%!in%' <- function(x,y)!('%in%'(x,y))
    
    multi_species_df = multi_species_df %>% 
      filter(true_name %!in% offender_list)
    
    
  }
  
  # split multi-species dataframe into a list
  butt_list = split(multi_species_df, f = multi_species_df$true_name)
  
  # Iterating the prep_data function over the list of species dataframes
  prepped_data_list = lapply(butt_list, 
                             try(prep_data), 
                             env_raster_t1 = env_raster_t1, 
                             env_raster_t2 = env_raster_t2, 
                             year_split = year_split)
  
  # Generating blockCV objects for each time period for each species and attaching to master list
  
  for(i in 1:length(prepped_data_list)){
    # initializing the list
    block_list = list()
    
    # Calculating stuff
    block_t1 = try(run_block_cv(prepped_data = prepped_data_list[[i]][[1]][[1]], 
                                bv_raster = prepped_data_list[[i]][[2]][[1]]))
    block_t2 = try(run_block_cv(prepped_data = prepped_data_list[[i]][[1]][[2]], 
                                bv_raster = prepped_data_list[[i]][[2]][[2]]))
    
    block_list = list(t1_block = block_t1, t2_block = block_t2)
    prepped_data_list[[i]]$blocks = block_list
  }
  
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
