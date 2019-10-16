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
    select(name = true_name, longitude, latitude, date, year) %>%
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
    dplyr::select(-ID, Species, longitude, latitude, Bio1:Bio19) %>%
    filter_all(all_vars(. != -Inf))
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
  training_data = extra_prepped_data[indices[[1]],] %>%
    drop_na()
  test_data = extra_prepped_data[-indices[[2]],] %>%
    drop_na()
  
  return(list(training_data = training_data, 
              test_data = test_data))
}


# Modeling ----------------------------------------------------------------

model_func = function(data = NULL, env_raster, num_cores = NULL) {
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
                         # parallel = TRUE,
                         # numCores = num_cores,
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


# Evaluation plots ---------------------------------------------------------

eval_plots = function(eval_object = NULL) {
  par(mfrow=c(2,3))
  eval.plot(eval_object@results)
  eval.plot(eval_object@results, 'auc.test.avg', legend = F)
  eval.plot(eval_object@results, 'auc.diff.avg', legend = F)
  eval.plot(eval_object@results, 'or.mtp.avg', legend = F)
  eval.plot(eval_object@results, 'or.10p.avg', legend = F)
  plot(eval_object@results$auc.test.avg, eval_object@results$delta.AICc, bg=as.factor(eval_object@results$fc), pch= 21, cex = eval_object@results$rm/2, xlab = "avg.test.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
  legend("topright", legend=unique(eval_object@results$fc), pt.bg=as.factor(eval_object@results$fc), pch=21)
  mtext("Circle size proportional to regularization multiplier value", cex = 0.6)
  
}

# Model Selection ---------------------------------------------------------

best_mod = function(model_obj){
  best_index = as.numeric(row.names(model_obj@results[which(model_obj@results$auc.test.avg== max(model_obj@results$auc.test.avg)),]))[1]
  
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
  
  ev = evaluate(test_data_occ, a = bg_data, model = model, x = env_raster, type = 'cloglog')
  return(ev)
}


# make args ---------------------------------------------------------------
make.args <- function(RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), labels=FALSE) {
  
  other.args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
  args.list <- list()
  
  for (i in 1:length(fc)) {
    args.list[[i]] <- other.args
    if(!grepl("L", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "nolinear")
    if(!grepl("Q", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "noquadratic")
    if(!grepl("H", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "nohinge")
    if(!grepl("P", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "noproduct")
    if(!grepl("T", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "nothreshold")
  }
  
  RM.lab <- rep(RMvalues, each=length(fc))
  RM.arg <- paste("betamultiplier=", RM.lab, sep="")
  fc.lab <- rep(fc, times=length(RMvalues))
  fc.arg <- rep(args.list, times=length(RMvalues))
  
  args <- list()
  feats.lab <- c()
  rms.lab <- c()
  for (i in 1:length(fc.lab)) {
    args[[i]] <- c(RM.arg[i], fc.arg[[i]])
    feats.lab <- c(feats.lab, fc.lab[[i]])
    rms.lab <- c(rms.lab, RM.lab[i])
  }
  args.lab <- list(feats.lab, rms.lab)
  
  if(labels==FALSE) {
    return(args)
  } else {
    return(args.lab)
  }
}
# Building full models on all data ----------------------------------------

full_model = function(models_obj, best_model_index, full_data = NULL, env_raster) {
  auc_mod = models_obj@results[best_model_index,]
  FC_best = as.character(auc_mod$fc[1])
  rm_best = auc_mod$rm
  
  
  maxent.args = make.args(RMvalues = rm_best, fc = FC_best)
  
  # re calculating environmental raster
  
  
  full_mod = maxent(env_raster, as.matrix(full_data[,1:2]), args = maxent.args[[1]])
  return(full_mod)
}


# Master Function - build_sdm() -------------------------------------------

build_sdm = function(multi_species_df, year_split, env_raster_t1, env_raster_t2){
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
  
  # Modeling
  for(k in 1:length(prepped_data_list)){
    models_t1 = try(model_func(data = prepped_data_list[[k]]$train_test_split[[1]]$training_data, 
                               env_raster = prepped_data_list[[k]]$env_data[[1]]))
    models_t2 = try(model_func(data = prepped_data_list[[k]]$train_test_split[[2]]$training_data, 
                               env_raster = prepped_data_list[[k]]$env_data[[2]]))
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
                                 full_data = butt_list[[n]] %>%
                                   filter(year < year_split), 
                                 env_raster = prepped_data_list[[n]]$env_data[[1]]
    ))
    
    full_mod_t2 = try(full_model(models_obj = prepped_data_list[[n]]$models[[2]],
                                 best_model_index = prepped_data_list[[n]]$best_mod[[2]][[2]],
                                 full_data = butt_list[[n]] %>%
                                   filter(year >= year_split), 
                                 env_raster = prepped_data_list[[n]]$env_data[[2]]))
    
    prepped_data_list[[n]]$full_mods = list(full_mod_t1, full_mod_t2)
    
  }
  return(prepped_data_list)
  saveRDS(prepped_data_list, "./data/model_data_list.rds")
}

# Testing Sandbox ---------------------------------------------------------
# # Loading saved temp objects
# best_model = readRDS("./tmp/best_mod.rds")
# block_test = readRDS("./tmp/block_test.rds")
# ev = readRDS("./tmp/evaluation_obj.rds")
# eval = readRDS("./tmp/models.rds")

# Need to do a bit of data wrangling before feeding into the model
# #
#  test_data = read_csv("./data/candidate_occurences.csv") %>%
#    filter(name == "Leptotes marina") %>%
#    mutate(true_name = name,
#           year = lubridate::year(date)) %>%
#      select(-name)
# # #
# test_prepped = prep_data(test_data, env_raster_t1 = bv_t1, env_raster_t2 = bv_t2)
# # 
# block_test = run_block_cv(test_prepped[[1]][[1]], test_prepped[[2]][[1]])
# # 
# prep_2_test = prep_data_2(data = test_prepped[[1]][[1]], env_raster = test_prepped[[2]][[1]])
# # #
# split_test = train_test_split(prep_2_test, block_test)
# # #
# # # # Takes forever to run
# # #
# doParallel::registerDoParallel(cores = 2)
# eval = model_func(data = split_test[[1]], env_raster = test_prepped[[2]][[1]], num_cores = 4)
# # 
# # # #
# # # # # plots
# eval_plots(eval)
# # #
# # # # best mod
# best_model = best_mod(eval)
# # #
# # # # ev obj
# # #
# ev = evaluate_models(test_data = split_test$test_data,
#                      env_raster = test_prepped[[2]][[1]],
#                      model = best_model[[1]])
# # #
# # #
# full_mod = full_model(models_obj = eval,
#                       best_model_index = best_model[[2]],
#                       full_data = test_data[,1:2],
#                       env_raster = test_prepped[[2]][[1]])
# #

# #Testing master function
# full_data = read_csv("./data/candidate_occurences.csv")
#
#
# sum(is.na(full_data$year))
# sum(is.na(full_data$eventDate))
#
#
# full_data %>%
#   group_by(true_name, time_frame) %>%
#   summarize(n = n())
#
# # write_csv(full_data, "./data/candidate_occurences.csv")
#
# small_test_data = full_data %>%
#   filter(true_name == "Danaus plexippus" | true_name == "Colias eurytheme")
#


# Run ---------------------------------------------------------------------

'%!in%' <- function(x,y)!('%in%'(x,y))

read_csv("./data/candidate_occurences.csv") %>%
  mutate(true_name = name,
         year = lubridate::year(date), 
         time_frame = ifelse(year < 2000, "t1", "t2")) %>%
  drop_na() %>%
  group_by(true_name, time_frame) %>%
  summarize(n = n()) %>%
  print(n = 45)

data = read_csv("./data/candidate_occurences.csv") %>%
  mutate(true_name = name,
         year = lubridate::year(date)) %>%
  drop_na() %>%
  filter(true_name == "Parnassius clodius") %>%
  select(-name)

test_cases_2 = build_sdm(multi_species_df = data, year_split = 2000, env_raster_t1 = bv_t1, env_raster_t2 = bv_t2)
saveRDS(test_cases_2, "./output/clodius_new_ev.rds")



test_ev_2 = evaluate_models(test_data = test_cases_2$`Parnassius clodius`$prepped_dfs[[2]], 
                          model = test_cases_2$`Parnassius clodius`$best_mod[[2]][[1]], 
                          env_raster = test_cases_2$`Parnassius clodius`$env_data[[2]])
