# Function for building butterfly SDMs
# Keaton Wilson
# keatonwilson@me.com
# 2019-09-20


# packages ----------------------------------------------------------------

library(dismo)
library(raster)
library(tidyverse)
library(blockCV)
library(tidyverse)
library(maxnet)
library(ENMeval)


# Importing big bioclim data ----------------------------------------------

bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")

# Prepping Occurrence Data ------------------------------------------------

prep_data = function(data) {
  
  # selecting the pieces we want and separating by time
  small_data = data %>%
    select(name = true_name, longitude, latitude, date = eventDate, year) %>%
    mutate(time_frame = ifelse(year < 2000, "t1", "t2"))
  
  # calculating extent of occurences
  max_lat = ceiling(max(small_data$latitude))
  min_lat = floor(min(small_data$latitude))
  max_lon = ceiling(max(small_data$longitude))
  min_lon = floor(min(small_data$longitude))
  
  # added a 1º buffer in every direction
  geographic_extent <- extent(x = c(min_lon-1, max_lon+1, min_lat-1, max_lat+1))
  
  # Crop bioclim data to geographic extent of species
  bv_t1_cropped <- crop(x = bv_t1, y = geographic_extent)
  bv_t2_cropped <- crop(x = bv_t2, y = geographic_extent)
  
  # Split by time period into two data frames
  df_t1 = small_data %>% 
    filter(time_frame == "t1")
  df_t2 = small_data %>%
    filter(time_frame == "t2")
  
  # print each to make sure it looks ok
  print(glimpse(df_t1))
  print(glimpse(df_t2))
  
  # Generate 10k background points for each one. 
  bg_t1 = dismo::randomPoints(bv_t1, 10000)
  colnames(bg_t1) = c("longitude", "latitude")
  
  bg_t2 = randomPoints(bv_t2, 10000)
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
  return(prepared_data_list)
}


# Block CV ----------------------------------------------------------------
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
                                 progress = T)
  return(blocked)
  }

# Testing Sandbox ---------------------------------------------------------

test_data = read_csv("./data/candidate_occurences.csv") %>%
  filter(true_name == "Leptotes marina")

test_prepped = prep_data(test_data)
block_test = run_block_cv(test_prepped[[1]][[1]], test_prepped[[2]][[1]])

