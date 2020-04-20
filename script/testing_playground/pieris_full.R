# Pieris full bioclim model build
# Keaton Wilson
# keatonwilson@me.com
# 2020-04-17

# packages
library(ENMeval)
library(blockCV)
library(raster)
library(sf)
library(dismo)
library(sp)
library(tidyverse)

# loading in data
pieris = readRDS("./data/split_data/pieris_rapae.rds")
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

pieris_w_bg = data.frame(pieris) %>%
  mutate(pb = 1) %>%
  dplyr::select(pb, longitude, latitude) %>%
  bind_rows(data.frame(bg_1) %>% 
              mutate(pb = 0))  %>%
  mutate(Species = as.integer(pb)) %>%
  dplyr::select(-pb)

# Turning into a spatial points data frame
pieris_sp = SpatialPointsDataFrame(pieris_w_bg[,c("longitude","latitude")], 
                                  pieris_w_bg, 
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# block cv
source("./script/functions/run_block_cv.R")

pieris_block_obj = run_block_cv(pieris_sp, bioclim_cropped)

# data prep
source("./script/functions/prep_data_2.R")

pieris_extra_prepped = prep_data_2(pieris_sp, bioclim_cropped)

# train test split
source("./script/functions/train_test_split.R")
training_test_split = train_test_split(extra_prepped_data = pieris_extra_prepped, 
                                       blocked_obj = pieris_block_obj)

# modeling
source("./script/functions/model_func.R")
pieris_full_mod = model_func(data = training_test_split$training_data, 
                             env_raster = bioclim_cropped, 
                             num_cores = 10)
