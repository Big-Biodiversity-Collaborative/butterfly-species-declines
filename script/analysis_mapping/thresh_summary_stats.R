# Threshold summary stats
# Keaton Wilson
# keatonwilson@me.com
# 2020-03-30

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

# env rasters
bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")


# threshold summary stats
thresh_summary_stats = function(input_folder) {
  
  # getting list of things to run through
  file_list = str_subset(list.files(input_folder, full.names = TRUE), ".rds")
  
  # initializing dataframe to thresholded occurences into
  thresh_df = data.frame()
  
  # Looping through
  for(i in 1:length(file_list)) {
    
    spec_list = readRDS(file_list[i])
    
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
    species = str_remove(str_to_sentence(str_replace(str_remove(str_remove(file_list[i], 
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
    
    # combining
    all_data = t1_threshold %>% bind_rows(t2_threshold, .id = "time_frame") %>%
      mutate(species_name = species)
    
    # binding to master data
    thresh_df = bind_rows(thresh_df, all_data)
    
    
    # plot
    g1 = ggplot(all_data, aes(x = y, fill = as.factor(time_frame))) +
      geom_density(alpha = 0.7) +
      theme_classic() +
      xlab("Latitude (º)") +
      ylab("Density") +
      scale_fill_discrete(name = "Time Frame", 
                          labels = c("1959-1999", "2000-2020")) +
      coord_flip()
    
    ggsave(filename = paste0("./output/density_plots/", species, ".png"), 
           plot = g1)
    
  }
  
  # summary table
  thresh_df_summary = thresh_df %>%
    group_by(species_name, time_frame) %>%
    summarize(total_sq_km = n()*5) %>%
    ungroup() %>%
    group_by(species_name) %>%
    summarize(diff = last(total_sq_km) - first(total_sq_km))
  
  return(thresh_df_summary)
  
}

# testing
# test = thresh_summary_stats("./output")
# test
