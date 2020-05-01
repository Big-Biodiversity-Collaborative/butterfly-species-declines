# Monte carlo methods master function
# Keaton Wilson
# keatonwilson@me.com
# 2020-04-20

# packages
require(raster)
require(tidyverse)
require(parallel)

# loading in required environmental data
bv_full = readRDS("./data/bioclim_full.rds")
bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")
names(bv_t1) = paste0("Bio", seq(1:19))
names(bv_t2) = paste0("Bio", seq(1:19))


# sourcing other functions
# List of operations and associated functions
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

# function
# Lots going on here - but essentially will take species data and output either
# a. A comparison of the t2 and t1 time periods where t2 is downsampled to 
# the number of occurences in t1 multiple times (you can set this parameter), 
# each downsampling results in a separate model, and model summaries are 
# combined and a distribution of two summary statistics (average area predicted, 
# and average latitude) are generated. This is this compared to the single 
# statistic from T1. 
# 
# Alternatively, the downsampled data is drawn from ALL occurences (t1+t2). 
#

monte_carlo_sdm = function(path_to_rds_data, post_or_all = "all", 
                           number_of_iterations = 100, num_cores = 10){

# Setup -------------------------------------------------------------------

  # unpacking the data
  species = readRDS(path_to_rds_data)
  
  # calculating extent of occurences
  max_lat = ceiling(max(species$latitude))
  min_lat = floor(min(species$latitude))
  max_lon = ceiling(max(species$longitude))
  min_lon = floor(min(species$longitude))
  
  # added a 1º buffer in every direction
  geographic_extent <- extent(x = c(min_lon-1, max_lon+1, min_lat-1, max_lat+1))
  
  # Crop bioclim data to geographic extent of species
  bioclim_full_cropped = crop(x = bv_full, y = geographic_extent)
  bv_t1_cropped = crop(x = bv_t1, y = geographic_extent)
  bv_t2_cropped = crop(x = bv_t2, y = geographic_extent)
  
# T1 Building -------------------------------------------------------------
  #t1 work - we'll do this now because we always do t1, no matter the choices 
  # in the function
  
  small_data = species %>%
  mutate(year = lubridate::year(date)) %>%
    dplyr::select(name, longitude, latitude, date, year) %>%
    mutate(time_frame = ifelse(year < 2000, "t1", "t2"))
  
  df_t1 = small_data %>% 
    dplyr::filter(time_frame == "t1")
  df_t2 = small_data %>%
    dplyr::filter(time_frame == "t2")
  
  #function to generate background points and add to an existing df
  bg_add = function(df){
    bg_1 = dismo::randomPoints(bv_t1_cropped, 10000)
    colnames(bg_1) = c("longitude", "latitude")
    
    df_comb = data.frame(df) %>%
      mutate(pb = 1) %>%
      dplyr::select(pb, longitude, latitude) %>%
      bind_rows(data.frame(bg_1) %>% 
                  mutate(pb = 0))  %>%
      mutate(Species = as.integer(pb)) %>%
      dplyr::select(-pb)
    
    return(df_comb)
  }
  
  # doing the same for t1
  species_t1_with_bg = bg_add(df_t1)
  
  df_sp_t1 = SpatialPointsDataFrame(species_t1_with_bg[,c("longitude","latitude")], 
                                    species_t1_with_bg, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  df_sp_t1$time_frame = "t1"
  
  # run blockCV
  block_obj_t1 = run_block_cv(prepped_data = df_sp_t1, 
                              bv_raster =bv_t1_cropped)
  
  # more prep
  names(bv_t1) = paste0("Bio", seq(1:19))
  prepped_t1 = prep_data_2(data = df_sp_t1, env_raster = bv_t1_cropped)
  
  # train test split
  training_test_split_t1 = train_test_split(extra_prepped_data = prepped_t1, 
                                            blocked_obj = block_obj_t1)
  # modeling
  t1_mod = model_func(data = training_test_split_t1$training_data, 
                      env_raster = bv_t1_cropped, 
                      num_cores = num_cores)
  # best mod
  t1_best_mod = best_mod(t1_mod)
  
  # evaluating models
  evaluation_t1 = evaluate_models(test_data = training_test_split_t1$test_data, 
                                  model = t1_best_mod[[1]], 
                                  env_raster = bv_t1_cropped)

  # full model
  auc_mod = t1_mod@results[t1_best_mod[[2]],]                         
  FC_best = tolower(as.character(auc_mod$fc[1]))  
  rm_best = as.numeric(auc_mod$rm)
  best_mod_t1 = maxnet(p = prepped_t1$Species, data = prepped_t1[,1:19],
                       maxnet.formula(prepped_t1$Species, 
                                      prepped_t1[,1:19], 
                                      classes = FC_best),
                       regmult = rm_best)
  
  # bundling data, full model, and evaluations and saving
  t1_list = list(df_sp_t1, best_mod_t1, evaluation_t1)
  
  species_name = str_remove(str_remove(path_to_rds_data, "./data/split_data/"),
                                           ".rds")
  
  saveRDS(t1_list, paste0("./output/", species_name, "t1_output.rds"))
  

# T2 or All Iterating and Building ----------------------------------------

  if(post_or_all == "all"){
    # Building on ALL ---------------------------------------------------------
    # downsample data n times
    species_all_downsample_list = list()
    
    downsample_num = nrow(df_t1)
    
    for(i in 1:number_of_iterations){
      downsample = sample_n(species, downsample_num)
      species_all_downsample_list[[i]] = downsample
    }
    
    # lapply bg_add function over the list
    #function to generate background points and add to an existing df
    bg_add = function(df){
      bg_1 = dismo::randomPoints(bioclim_full_cropped, 10000)
      colnames(bg_1) = c("longitude", "latitude")
      
      df_comb = data.frame(df) %>%
        mutate(pb = 1) %>%
        dplyr::select(pb, longitude, latitude) %>%
        bind_rows(data.frame(bg_1) %>% 
                    mutate(pb = 0))  %>%
        mutate(Species = as.integer(pb)) %>%
        dplyr::select(-pb)
      
      return(df_comb)
    }
    species_all_list_with_bg = lapply(species_all_downsample_list, bg_add)
    
    # changing to a spatial points data frame
    # function
    spdf_func = function(df){
      df_sp = SpatialPointsDataFrame(df[,c("longitude","latitude")], 
                                        df, 
                                        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      df_sp$time_frame = "all"
      return(df_sp)
    }
    
    # applying
    df_sp_all_list = lapply(species_t2_list_with_bg, 
                            spdf_func)
    
    # block CV
    block_obj_all_list = lapply(df_sp_all_list, 
                                run_block_cv, 
                                bv_raster = bioclim_full_cropped)
    
    # more prep
    prepped_all_list = lapply(df_sp_all_list, 
                             prep_data_2, 
                             env_raster = bioclim_full_cropped)
    
    # training and test data split
    training_test_split_all_list = list()
    for(i in 1:length(prepped_all_list)){
      training_test_split_all_list[[i]] = train_test_split(extra_prepped_data = prepped_all_list[[i]], 
                                                          blocked_obj = block_obj_all_list[[i]])
    }
    
    all_training_data_only_list = lapply(training_test_split_all_list, "[[", 1)
    all_test_data_only_list = lapply(training_test_split_all_list, "[[", 2)
    
    # modeling
    all_mod_list = list()
    for(i in 1:length(all_training_data_only_list)){
      all_mod_list[[i]] = model_func(data = all_training_data_only_list[[i]], 
                                    env_raster = bioclim_full_cropped, 
                                    num_cores = num_cores)
    }
    
    # model selection
    all_best_mod_list = (lapply(all_mod_list, best_mod))
    
    # evaluating
    evaluation_all_list = list()
    for(i in 1:length(all_best_mod_list)){
      evaluation_all_list[[i]] = evaluate_models(test_data = all_test_data_only_list[[i]], 
                                                model = all_best_mod_list[[i]][[1]], 
                                                env_raster = bioclim_full_cropped)
    }
    
    # best models full
    best_mods_all_list = list()
    for(i in 1:length(evaluation_all_list)){
      auc_mod = all_mod_list[[i]]@results[all_best_mod_list[[i]][[2]],]                         
      FC_best = tolower(as.character(auc_mod$fc[1]))  
      rm_best = as.numeric(auc_mod$rm)
      best_mods_all_list[[i]] = maxnet(p = prepped_all_list[[i]]$Species, data = prepped_all_list[[i]][,1:19],
                                      maxnet.formula(prepped_all_list[[i]]$Species, 
                                                     prepped_all_list[[i]][,1:19], 
                                                     classes = FC_best),
                                      regmult = rm_best)
    }
    
    # bundling data, full model, and evaluations and saving
    all_list = list(df_sp_all_list, best_mods_all_list, evaluation_all_list)
    
    species_name = str_remove(str_remove(path_to_rds_data, "./data/split_data/"),
                              ".rds")
    
    saveRDS(all_list, paste0("./output/", species_name, "all_output.rds"))
    
    
  } else {
    # Building on T2 ---------------------------------------------------------
    # downsample data n times
    species_t2_downsample_list = list()
    
    downsample_num = nrow(df_t1)
    
    for(i in 1:number_of_iterations){
      downsample = sample_n(df_t2, downsample_num)
      species_t2_downsample_list[[i]] = downsample
    }
    
    # lapply bg_add function over the list
    #function to generate background points and add to an existing df
    bg_add = function(df){
      bg_1 = dismo::randomPoints(bv_t2_cropped, 10000)
      colnames(bg_1) = c("longitude", "latitude")
      
      df_comb = data.frame(df) %>%
        mutate(pb = 1) %>%
        dplyr::select(pb, longitude, latitude) %>%
        bind_rows(data.frame(bg_1) %>% 
                    mutate(pb = 0))  %>%
        mutate(Species = as.integer(pb)) %>%
        dplyr::select(-pb)
      
      return(df_comb)
    }
    species_t2_list_with_bg = lapply(species_t2_downsample_list, bg_add)
    
    # changing to a spatial points data frame
    # function
    spdf_func = function(df){
      df_sp = SpatialPointsDataFrame(df[,c("longitude","latitude")], 
                                     df, 
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      df_sp$time_frame = "t2"
      return(df_sp)
    }
    
    # applying
    df_sp_t2_list = lapply(species_t2_list_with_bg, 
                            spdf_func)
    
    # block CV
    block_obj_t2_list = lapply(df_sp_t2_list, 
                                run_block_cv, 
                                bv_raster = bv_t2_cropped)
    
    # more prep
    prepped_t2_list = lapply(df_sp_t2_list, 
                              prep_data_2, 
                              env_raster = bv_t2_cropped)
    
    # training and test data split
    training_test_split_t2_list = list()
    for(i in 1:length(prepped_t2_list)){
      training_test_split_t2_list[[i]] = train_test_split(extra_prepped_data = prepped_t2_list[[i]], 
                                                           blocked_obj = block_obj_t2_list[[i]])
    }
    
    t2_training_data_only_list = lapply(training_test_split_t2_list, "[[", 1)
    t2_test_data_only_list = lapply(training_test_split_t2_list, "[[", 2)
    
    # modeling
    t2_mod_list = list()
    for(i in 1:length(t2_training_data_only_list)){
      t2_mod_list[[i]] = model_func(data = t2_training_data_only_list[[i]], 
                                     env_raster = bv_t2_cropped, 
                                     num_cores = num_cores)
    }
    
    # model selection
    t2_best_mod_list = (lapply(t2_mod_list, best_mod))
    
    # evaluating
    evaluation_t2_list = list()
    for(i in 1:length(t2_best_mod_list)){
      evaluation_t2_list[[i]] = evaluate_models(test_data = t2_test_data_only_list[[i]], 
                                                 model = t2_best_mod_list[[i]][[1]], 
                                                 env_raster = bv_t2_cropped)
    }
    
    # best models full
    best_mods_t2_list = list()
    for(i in 1:length(evaluation_t2_list)){
      auc_mod = t2_mod_list[[i]]@results[t2_best_mod_list[[i]][[2]],]                         
      FC_best = tolower(as.character(auc_mod$fc[1]))  
      rm_best = as.numeric(auc_mod$rm)
      best_mods_t2_list[[i]] = maxnet(p = prepped_t2_list[[i]]$Species, data = prepped_t2_list[[i]][,1:19],
                                       maxnet.formula(prepped_t2_list[[i]]$Species, 
                                                      prepped_t2_list[[i]][,1:19], 
                                                      classes = FC_best),
                                       regmult = rm_best)
    }
    
    # bundling data, full model, and evaluations and saving
    t2_list = list(df_sp_t2_list, best_mods_t2_list, evaluation_t2_list)
    
    species_name = str_remove(str_remove(path_to_rds_data, "./data/split_data/"),
                              ".rds")
    
    saveRDS(t2_list, paste0("./output/", species_name, "t2_output.rds"))
  }
  
}