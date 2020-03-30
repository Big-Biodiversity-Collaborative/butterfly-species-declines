# Making threshold maps
# Keaton Wilson
# keatonwilson@me.com
# 2020-03-23

# packages
library(tidyverse)
library(lubridate)
library(sp)
library(raster)
library(maptools)
library(ggmap)
library(viridis)
library(ggthemes)
library(rgeos)
library(maps)
library(ggpubr)
library(ENMeval)
library(maxnet)
library(stringr)
library(dismo)



# Geographic Mapping Data ---------------------------------------

#Pulling in polygons for states and provinces
#Getting map data
usa = getData(country = 'USA', level = 1, path = "./data/")
#extract states (need to uppercase everything)
to_remove = c("Alaska", "Hawaii")

#filtering
usa = usa[-match(toupper(to_remove), toupper(usa$NAME_1)),]

#simplying polygons
simple_map_US = gSimplify(usa, tol = 0.01, topologyPreserve = TRUE)

#Pulling Canada Province data
can = getData(country = 'CAN', level = 1, path = "./data/")
simple_map_can = gSimplify(can, tol = 0.01, topologyPreserve = TRUE)

#Pulling Mexico data
mex = getData(country = 'MEX', level = 1, path = "./data/")
simple_map_mex = gSimplify(mex, tol = 0.01, topologyPreserve = TRUE)

# Function
make_map_time_compare_threshold = function(species_rds_out, 
                                           env_raster_t1, 
                                           env_raster_t2){
  
  # reading in list of stuff
  spec_list = readRDS(species_rds_out)
  
  #Cropping data to actual occurences
  max_lat = ceiling(max(spec_list$raw_data$latitude))
  min_lat = floor(min(spec_list$raw_data$latitude))
  max_lon = ceiling(max(spec_list$raw_data$longitude))
  min_lon = floor(min(spec_list$raw_data$longitude))
  
  # added a 1º buffer in every direction
  geographic_extent <- extent(x = c(min_lon-1, max_lon+1, min_lat-1, max_lat+1))
  
  # Crop bioclim data to geographic extent of species
  bv_t1_cropped <- crop(x = bv_t1, y = geographic_extent)
  bv_t2_cropped <- crop(x = bv_t2, y = geographic_extent)
  
  # Create new data to predict on
  newdata_t1 = as(bv_t1_cropped, "SpatialPixelsDataFrame")
  newdata_t1 = as.data.frame(newdata_t1) %>%
    drop_na()
  
  newdata_t2 = as(bv_t2_cropped, "SpatialPixelsDataFrame")
  newdata_t2 = as.data.frame(newdata_t2) %>%
    drop_na()
  
  names = c(names(spec_list$prepped_data$env_data[[1]]), "x", "y")
  names(newdata_t1) = names
  names(newdata_t2) = names
  
  # predictions
  pred_t1 = dismo::predict(object = spec_list$full_mods$full_mod_t1,
                           newdata = newdata_t1,
                           x = spec_list$prepped_data$env_data[[1]],
                           ext = spec_list$prepped_data$env_data[[1]]@extent,
                           type = "cloglog")
  
  if(class(pred_t1) == "RasterLayer"){
    pred_t1_df = as(pred_t1, "SpatialPixelsDataFrame")
    pred_t1_df = as.data.frame(pred_t1_df) %>%
      drop_na()
    
    colnames(pred_t1_df) = c("value", "x", "y")
    
  } else {
    pred_t1_df = newdata_t1 %>%
      dplyr::select(x, y) %>%
      cbind(as.data.frame(pred_t1)) %>%
      as_tibble() 
    
    colnames(pred_t1_df) = c("x", "y", "value")
  }
  
  #T2 Predictions
  pred_t2 = dismo::predict(object = spec_list$full_mods$full_mod_t1,
                           newdata = newdata_t2,
                           x = spec_list$prepped_data$env_data[[1]],
                           ext = extent(spec_list$prepped_data$env_data[[1]]),
                           type = "cloglog")
  
  if(class(pred_t2) == "RasterLayer"){
    pred_t2_df = as(pred_t2, "SpatialPixelsDataFrame")
    pred_t2_df = as.data.frame(pred_t2_df) %>%
      drop_na()
    
    colnames(pred_t2_df) = c("value", "x", "y")
    
  } else {
    pred_t2_df = newdata_t2 %>%
      dplyr::select(x, y) %>%
      cbind(as.data.frame(pred_t2)) %>%
      as_tibble() 
    
    colnames(pred_t2_df) = c("x", "y", "value")
  }
  
  # Extracting the species name
  species = str_remove(str_to_sentence(str_replace(str_remove(str_remove(species_rds_out, 
                                                                         "_output.rds"), 
                                                              "./output/"), 
                                                   "_", " ")
  ), ".rds")
  
  # loading evaluation objs and thresholds
  eval_t1 = spec_list$eval_objs$eval_t1
  eval_t2 = spec_list$eval_objs$eval_t2
  
  thresh_t1 = threshold(eval_t1, 'spec_sens')
  thresh_t2 = threshold(eval_t2, 'spec_sens')
  
  # thresholded filter
  t1_threshold = pred_t1_df %>%
    filter(value > thresh_t1)
  
  t2_threshold = pred_t2_df %>%
    filter(value > thresh_t2)
  
  
  # Plotting
  # 
  # # Plotting
  g1 = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=t1_threshold, aes(x=x, y=y), fill = "lightgrey") + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.20, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.20, fill = NA) +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.20, fill = NA) +
    scale_fill_viridis(name = "Probability of Occurence") +
    theme(legend.position="right") +
    theme(legend.key.width=unit(2, "cm"),
          plot.title = element_text(hjust = 0.5, size = 24)) +
    theme_nothing(legend = TRUE) +
    coord_quickmap() +
    ggtitle("pre-2000")
  
  #T2
  g2 = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=t2_threshold, aes(x=x, y=y), fill = "lightgrey") + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.20, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.20, fill = NA) +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.20, fill = NA) +
    scale_fill_viridis(name = "Probability of Occurence") +
    theme(legend.position="right") +
    theme(legend.key.width=unit(2, "cm"),
          plot.title = element_text(hjust = 0.5, size = 24)) +
    theme_nothing(legend = TRUE) +
    coord_quickmap() +
    ggtitle("post-2000")
  
  gfull = ggarrange(g1, g2, common.legend = TRUE)
  gfull = annotate_figure(gfull,
                          top = text_grob(species, face = "italic", size = 22))
  plotname = paste0("./output/thresh_maps/", species,"_thresh", ".png")
  print(species)
  print(plotname)
  ggsave(plotname, gfull)
  
}


# Testing
# make_map_time_compare_threshold(species_rds_out = "./output/strymon_melinus.rds",
#                       bv_t1, bv_t2)