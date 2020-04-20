# Function that builds a single sdm with all data and average (1970-2000 bioclim)
# Keaton Wilson
# keatonwilson@email.arizona.edu
# 2020-04-17

# packages
require(raster)
require(tidyverse)
require(parallel)

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

# importing complete environmental raster
complete_bv = readRDS("./data/bioclim_full.rds")

build_complete_sdm = function(path_to_rds_data){
  
  # unpacking the data
  species = readRDS(path_to_rds_data)
  bioclim = readRDS("./data/bioclim_full.rds")
  
  # calculating extent of occurences
  max_lat = ceiling(max(pieris$latitude))
  min_lat = floor(min(pieris$latitude))
  max_lon = ceiling(max(pieris$longitude))
  min_lon = floor(min(pieris$longitude))
  
  # added a 1º buffer in every direction
  geographic_extent <- extent(x = c(min_lon-1, max_lon+1, min_lat-1, max_lat+1))
  
  # Crop bioclim data to geographic extent of species
  bioclim_cropped = crop(x = bioclim, y = geographic_extent)
  
  # checking
  plot(bioclim_cropped)
  
  # background points
  bg_1 = dismo::randomPoints(bioclim_cropped, 10000)
  colnames(bg_1) = c("longitude", "latitude")
  
  species_w_bg = data.frame(species) %>%
    mutate(pb = 1) %>%
    dplyr::select(pb, longitude, latitude) %>%
    bind_rows(data.frame(bg_1) %>% 
                mutate(pb = 0))  %>%
    mutate(Species = as.integer(pb)) %>%
    dplyr::select(-pb)
  
  # Turning into a spatial points data frame
  species_sp = SpatialPointsDataFrame(species_w_bg[,c("longitude","latitude")], 
                                     species_w_bg, 
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  # block cv
  source("./script/functions/run_block_cv.R")
  
  species_block_obj = run_block_cv(species_sp, bioclim_cropped)
  
  # data prep
  source("./script/functions/prep_data_2.R")
  
  species_extra_prepped = prep_data_2(species_sp, bioclim_cropped)
  
  # train test split
  source("./script/functions/train_test_split.R")
  training_test_split = train_test_split(extra_prepped_data = species_extra_prepped, 
                                         blocked_obj = species_block_obj)
  
  # modeling
  source("./script/functions/model_func.R")
  species_full_mod = model_func(data = training_test_split$training_data, 
                               env_raster = bioclim_cropped, 
                               num_cores = 10)
  
  # model selection
  species_best_mod = best_mod(species_full_mod)
  
  # evaluation
  species_mod_eval = evaluate_models(training_test_split$test_data, 
                                     model = species_best_mod, 
                                     env_raster = bioclim_cropped)
  # full model
  species_full_mod = full_model(species_full_mod, 
                                best_model_index = species_best_mod[[2]], 
                                full_data = species_extra_prepped, 
                                env_raster = bioclim_cropped)
  
  # The two things we really need are the full model and the evaluation obj
  species_name = str_remove(str_remove(path_to_rds_data, "./data/split_data/"),
                            ".rds")
  
  saveRDS(species_mod_eval, paste0("./output/", species_name, "_eval.rds"))
  saveRDS(species_full_mod, paste0("./output/", species_name, "_full_mod.rds"))
}