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

prepped = prep_data(data = readRDS("./data/split_data/colias_eurytheme.rds"), 
                    year_split = 2000, 
                    env_raster_t1 = bv_t1, 
                    env_raster_t2 = bv_t2)

block_1 = run_block_cv(prepped$data$t1, bv_t1)

prepped_2 = prep_data_2(data = prepped$data$t1, prepped$env_data[[1]])

train_test = train_test_split(extra_prepped_data = prepped_2, 
                              blocked_obj = block_1)

total_cores = parallel::detectCores()
to_use = total_cores - 2
doParallel::registerDoParallel(to_use)
doSNOW::registerDoSNOW(cl = makeCluster(spec = ))

t1_mod = model_func(data = train_test$training_data, 
                    env_raster = prepped$env_data[[1]])

java_test = ENMevaluate(occ = train_test$training_data %>%
                          filter(Species == 1) %>%
                          select(longitude, latitude) %>%
                          drop_na(), 
                       bg.coords = train_test$training_data %>%
                         filter(Species == 0) %>%
                         select(longitude, latitude) %>%
                         drop_na(),
                       env = prepped$env_data[[1]],
                       method = 'randomkfold', 
                       kfolds = 5,
                       parallel = TRUE,
                       numCores = to_use,
                       algorithm = 'maxent.jar')

