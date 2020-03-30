# make map comparisons (of T1 and T2) of a species range from model outputs
# Keaton Wilson
# keatonwilson@me.com
# 2020-02-19

# packages
library(tidyverse)
library(lubridate)
#library(rgdal)
library(sp)
library(raster)
library(maptools)
library(ggmap)
library(viridis)
library(ggthemes)
library(rgeos)
library(maps)
library(ggpubr)
#library(blockCV)
library(ENMeval)
library(maxnet)
library(stringr)

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

# making maps function
make_map_time_compare = function(species_rds_out, 
                                 env_raster_t1, 
                                 env_raster_t2){
  # unpacking rds
  spec_list = readRDS(species_rds_out)
  
  # building newdata to feed into predict
  newdata_t1 = as(spec_list$prepped_data$env_data[[1]], "SpatialPixelsDataFrame")
  newdata_t1 = as.data.frame(newdata_t1) %>%
    drop_na()
  
  newdata_t2 = as(spec_list$prepped_data$env_data[[2]], "SpatialPixelsDataFrame")
  newdata_t2 = as.data.frame(newdata_t2) %>%
    drop_na()
  
  # T1 Predictions
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
  species = str_to_sentence(str_replace(str_remove(str_remove(species_rds_out, ".rds"), "./output/"), "_", " ")) 
  
  # Plot building
  # Panel 1
  g1 = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=pred_t1_df, aes(x=x, y=y, fill=value)) + 
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
  
  g2 = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=pred_t2_df, aes(x=x, y=y, fill=value)) + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.20, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.25, fill = NA) +
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
  plotname = str_remove(paste0("./output/prob_maps/", species, ".png"), ".rds")
  ggsave(plotname, gfull)
}

# Testing
# make_map_time_compare(species_rds_out = "./output/strymon_melinus.rds",
#                       bv_t1, bv_t2)
