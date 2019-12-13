# Prepping Occurence Data for Modeling
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(tidyverse)
require(raster)
require(dismo)
require(sp)

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
    mutate(year = lubridate::year(date)) %>%
    dplyr::select(name, longitude, latitude, date, year) %>%
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
