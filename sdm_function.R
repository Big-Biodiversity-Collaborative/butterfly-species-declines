# Function for building butterfly SDMs
# Keaton Wilson
# keatonwilson@me.com
# 2019-09-20
# 
# packages ----------------------------------------------------------------

library(dismo)
library(raster)
library(tidyverse)
library(blockCV)
library(tidyverse)
library(maxnet)
library(ENMeval)
if (Sys.getenv("JAVA_HOME")!="")
  Sys.setenv(JAVA_HOME="")
library(rJava)


# Importing big bioclim data ----------------------------------------------

bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")

# Prepping Occurrence Data ------------------------------------------------

#' Initial prepping of occurence data to create SDMs
#'
#' @param data A dataframe outputted from the butt_obs script. Generated from 
#' \code{\link[spocc]{occ}}
#' @param year_split The year to split the data by. Non inclusive (e.g. 2000 will
#'  split the  data into everything through 1999, and 2000-everything after).
#'  Default is the year 2000. 
#'
#' @return A list four elements long: the first two are the occurence data split
#'  by the year_split argument with 10k background points added. Additionally, 
#'  they are converted to a spatial points dataframe. 
#'  The second two items are the env rasters cropped 
#'  to the area of the occurences for each subset split by year_split.
#'
#' @examples
prep_data = function(data, year_split = 2000, env_raster_t1, env_raster_t2) {
  
  # selecting the pieces we want and separating by time
  small_data = data %>%
    select(name = true_name, longitude, latitude, date = eventDate, year) %>%
    mutate(time_frame = ifelse(year < year_split, "t1", "t2"))
  
  # calculating extent of occurences
  max_lat = ceiling(max(small_data$latitude))
  min_lat = floor(min(small_data$latitude))
  max_lon = ceiling(max(small_data$longitude))
  min_lon = floor(min(small_data$longitude))
  
  # added a 1º buffer in every direction
  geographic_extent <- extent(x = c(min_lon-1, max_lon+1, min_lat-1, max_lat+1))
  
  # Crop bioclim data to geographic extent of species
  bv_t1_cropped <- crop(x = env_raster_t1, y = geographic_extent)
  bv_t2_cropped <- crop(x = env_raster_t2, y = geographic_extent)
  
  # Split by time period into two data frames
  df_t1 = small_data %>% 
    filter(time_frame == "t1")
  df_t2 = small_data %>%
    filter(time_frame == "t2")
  
  # print each to make sure it looks ok
  print(glimpse(df_t1))
  print(glimpse(df_t2))
  
  # Generate 10k background points for each one. 
  bg_t1 = dismo::randomPoints(bv_t1_cropped, 10000)
  colnames(bg_t1) = c("longitude", "latitude")
  
  bg_t2 = randomPoints(bv_t2_cropped, 10000)
  colnames(bg_t2) = c("longitude", "latitude")
  
  # Merging background data and occurence data
  df_comb_t1 = data.frame(df_t1) %>%
    mutate(pb = 1) %>%
    dplyr::select(pb, longitude, latitude) %>%
    bind_rows(data.frame(bg_t1) %>% 
                mutate(pb = 0))  %>%
    mutate(Species = as.integer(pb)) %>%
    dplyr::select(-pb)
  
  df_comb_t2 = data.frame(df_t2) %>%
    mutate(pb = 1) %>%
    dplyr::select(pb, longitude, latitude) %>%
    bind_rows(data.frame(bg_t2) %>% 
                mutate(pb = 0)) %>%
    mutate(Species = as.integer(pb)) %>%
    dplyr::select(-pb)
  
  # Changing to a spatial points data frame
  df_sp_t1 = SpatialPointsDataFrame(df_comb_t1[,c("longitude","latitude")], 
                                    df_comb_t1, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) 
  df_sp_t1$time_frame = "t1"
  df_sp_t2 = SpatialPointsDataFrame(df_comb_t2[,c("longitude","latitude")], 
                                    df_comb_t2, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  df_sp_t2$time_frame = "t2"
  
  #Converting to a list with the two dataframes
  prepared_data_list = list(data = list(t1 = df_sp_t1, t2 = df_sp_t2),
                            env_data = list(bv_t1_cropped, bv_t2_cropped))
  #Names
  bio_names = c()
  for(i in 1:19){
    bio_names[i] = paste0("Bio", i)
  }
  
  names(prepared_data_list[[2]][[1]]) = bio_names
  names(prepared_data_list[[2]][[2]]) = bio_names
  
  return(prepared_data_list)
}


# Block CV ----------------------------------------------------------------
#' Running blockCV with a preset config for this project
#'
#' @param prepped_data The prepped spatial points dataframe created by 
#' \code{link{prep_data}}
#' @param bv_raster the cropped raster associated with the same time period 
#' as the prepped_data above. 
#'
#' @return a blockCV object that we will use in later analysis
#'
#' @examples
run_block_cv = function(prepped_data, bv_raster){
  
  blocked = spatialBlock(speciesData = prepped_data,
                                 species = "Species",
                                 rasterLayer = bv_raster,
                                 theRange = 400000,
                                 k = 5, 
                                 selection = "random", 
                                 iteration = 250, 
                                 biomod2Format = TRUE, 
                                 xOffset = 0, 
                                 yOffset = 0, 
                                 progress = T, 
                                 showBlocks = F
                         )
  return(blocked)
}


# Preparing data 2 ----------------------------------------------------------
#' More data preparation prior to SDM building
#'
#' @param data Prepped spatial points datframe created by \code{link{prep_data}}
#' @param env_raster the cropped raster associated with the same time period 
#' as the prepped_data above.
#'
#' @return a dataframe with extracted environmental variables along with presence 
#' for all of the occurence and background data
#'
#' @examples
prep_data_2 = function(data, env_raster){
extra_prepped = raster::extract(env_raster, data, df = TRUE) %>%
  bind_cols(as.data.frame(data)) %>%
  drop_na() %>%
  dplyr::select(-ID, Species, longitude, latitude, Bio1:Bio19)
return(extra_prepped)
}


# Train and test split ---------------------------------------

train_test_split = function(extra_prepped_data, blocked_obj){

extract_index = function(list_of_folds = NULL) {
  for(k in 1:length(list_of_folds)){
    train_index <- unlist(list_of_folds[[k]][1]) # extract the training set indices
    test_index <- unlist(list_of_folds[[k]][2])# extract the test set indices
  }
  mini_list = list(train_index, test_index)
  return(mini_list)
}

indices = extract_index(blocked_obj$folds)
print(length(indices[[1]]))
print(length(indices[[2]]))

#applying indexes and splitting data
training_data = extra_prepped_data[indices[[1]],]
test_data = extra_prepped_data[-indices[[2]],]

return(list(training_data = training_data, 
            test_data = test_data))
}


# Modeling ----------------------------------------------------------------

model_func = function(data = NULL, env_raster) {
  data_occ = data %>%  #Generating occurence lat long
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude) %>%
    drop_na()
  
  data_bg = data %>% #Generating background lat long
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude) %>%
    drop_na()
  
  #Running the model
  eval = ENMevaluate(occ = data_occ, 
                     bg.coords = data_bg,
                     env = env_raster,
                     method = 'randomkfold', 
                     kfolds = 5, 
                     algorithm = 'maxent.jar')
  return(eval)
}


# Evaluation plots ---------------------------------------------------------

eval_plots = function(eval_object = NULL) {
  par(mfrow=c(2,3))
  eval.plot(eval_object@results)
  eval.plot(eval_object@results, 'avg.test.AUC', legend = F)
  eval.plot(eval_object@results, 'avg.diff.AUC', legend = F)
  eval.plot(eval_object@results, 'avg.test.or10pct', legend = F)
  eval.plot(eval_object@results, 'avg.test.orMTP', legend = F)
  plot(eval_object@results$avg.test.AUC, eval_object@results$delta.AICc, bg=eval_object@results$features, pch=21, cex= eval_object@results$rm/2, xlab = "avg.test.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
  legend("topright", legend=unique(eval_object@results$features), pt.bg=eval_object@results$features, pch=21)
  mtext("Circle size proportional to regularization multiplier value", cex = 0.6)
  
}


# Model Selection ---------------------------------------------------------

best_mod = function(model_obj){
  best_index = as.numeric(row.names(model_obj@results[which(model_obj@results$avg.test.AUC== max(model_obj@results$avg.test.AUC)),]))[1]

  best_mod = model_obj@models[[best_index]]
  return(list(best_mod, best_index))
}


# Evaluating on test data -------------------------------------------------

evaluate_models = function(test_data, model, env_raster) {
  test_data_occ = test_data %>%
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude)
  
  bg_data = test_data %>%
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude)
  
  ev = evaluate(test_data_occ, a = bg_data, model = model, x = env_raster)
  return(ev)
}


# Building full models on all data ----------------------------------------

full_model = function(models_obj, best_model_index, full_data = NULL, env_raster) {
  auc_mod = models_obj@results[best_model_index,]
  FC_best = as.character(auc_mod$features[1])
  rm_best = auc_mod$rm
  
  
  maxent.args = ENMeval::make.args(RMvalues = rm_best, fc = FC_best)
  
  # re calculating environmental raster
  

  full_mod = maxent(env_raster, as.matrix(full_data[,1:2]), args = maxent.args[[1]])
  return(full_mod)
}


# Master Function - build_sdm() -------------------------------------------

build_sdm = function(multi_species_df, year_split, env_raster_t1, env_raster_t2){
  # split multi-species dataframe into a list
  butt_list = split(multi_species_df, f = multi_species_df$true_name)
  
  # Iterating the prep_data function over the list of species dataframes
  prepped_data_list = lapply(butt_list, 
                             prep_data, 
                             env_raster_t1 = env_raster_t1, 
                             env_raster_t2 = env_raster_t2, 
                             year_split = year_split)
  
  # Generating blockCV objects for each time period for each species and attaching to master list

  for(i in 1:length(prepped_data_list)){
    # initializing the list
    block_list = list()
    
    # Calculating stuff
    block_t1 = run_block_cv(prepped_data = prepped_data_list[[i]][[1]][[1]], 
                            bv_raster = prepped_data_list[[i]][[2]][[1]])
    block_t2 = run_block_cv(prepped_data = prepped_data_list[[i]][[1]][[2]], 
                            bv_raster = prepped_data_list[[i]][[2]][[2]])
    
    block_list = list(t1_block = block_t1, t2_block = block_t2)
    prepped_data_list[[i]]$blocks = block_list
  }
  
  return(prepped_data_list)
  }

# Testing Sandbox ---------------------------------------------------------
# Loading saved temp objects
best_model = readRDS("./tmp/best_mod.rds")
block_test = readRDS("./tmp/block_test.rds")
ev = readRDS("./tmp/evaluation_obj.rds")
eval = readRDS("./tmp/models.rds")


test_data = read_csv("./data/candidate_occurences.csv") %>%
  filter(true_name == "Leptotes marina")

test_prepped = prep_data(test_data)
#block_test = run_block_cv(test_prepped[[1]][[1]], test_prepped[[2]][[1]])

prep_2_test = prep_data_2(data = test_prepped[[1]][[1]], env_raster = test_prepped[[2]][[1]])

split_test = train_test_split(prep_2_test, block_test)

# Takes forever to run
#eval = model_func(data = split_test[[1]], env_raster = test_prepped[[2]][[1]])

# plots
eval_plots(eval)

# best mod
#best_model = best_mod(eval)

# ev obj

# ev = evaluate_models(test_data = prep_data_2(test_prepped[[1]][[2]], 
#                                              env_raster = test_prepped[[2]][[2]]),
#                      model = best_model, 
#                      env_raster = test_prepped[[2]][[2]])


full_mod = full_model(models_obj = eval, 
                      best_model_index = best_model[[2]], 
                      full_data = test_data[,2:3],
                      env_raster = test_prepped[[2]][[1]])

#Testing master function
full_data = read_csv("./data/candidate_occurences.csv")
unique(full_data$true_name)

small_test_data = full_data %>%
  filter(true_name == "Leptotes marina" | true_name == "Vanessa atalanta")

t = build_sdm(multi_species_df = small_test_data, year_split = 2000, env_raster_t1 = bv_t1, env_raster_t2 = bv_t2)
