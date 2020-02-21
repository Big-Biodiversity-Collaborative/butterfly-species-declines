# make map comparisons (of T1 and T2) of a species range from model outputs
# Keaton Wilson
# keatonwilson@me.com
# 2020-02-19

# packages
library(tidyverse)
library(lubridate)
library(rgdal)
library(sp)
library(raster)
library(maptools)
library(ggmap)
library(viridis)
library(ggthemes)
library(rgeos)
library(maps)
library(ggpubr)
library(blockCV)
library(ENMeval)
library(maxnet)
library(stringr)

# Geographic Mapping Data ---------------------------------------

#Pulling in polygons for states and provinces
#Getting map data
usa = getData(country = 'USA', level = 1)
#extract states (need to uppercase everything)
to_remove = c("Alaska", "Hawaii")

#filtering
usa = usa[-match(toupper(to_remove), toupper(usa$NAME_1)),]

#simplying polygons
simple_map_US = gSimplify(usa, tol = 0.01, topologyPreserve = TRUE)

#Pulling Canada Province data
can = getData(country = 'CAN', level = 1)
simple_map_can = gSimplify(can, tol = 0.01, topologyPreserve = TRUE)

# making maps function
make_map_time_compare = function(species_rds_out, 
                                 env_raster_t1, 
                                 env_raster_t2){
  # unpacking rds
  spec_list = readRDS(species_rds_out)
  
  #T1 Predictions
  pred_t1 = dismo::predict(object = spec_list$full_mods$full_mod_t1,
                           newdata = spec_list$extra_prepped$extra_prepped_t1,
                           x = spec_list$prepped_data$env_data[[1]],
                           ext = extent(spec_list$prepped_data$env_data[[1]]),
                           args = 'outputformat=cloglog')
  
  pred_t1_df = spec_list$extra_prepped$extra_prepped_t1 %>%
    dplyr::select(longitude, latitude) %>%
    cbind(pred_t1) %>%
    as_tibble() 
  
  colnames(pred_t1_df) = c("value", "x", "y")
  
  #T2 Predictions
  pred_t2 = dismo::predict(object = spec_list$full_mods$full_mod_t2,
                           newdata = spec_list$extra_prepped$extra_prepped_t2,
                           x = spec_list$prepped_data$env_data[[2]],
                           ext = extent(spec_list$prepped_data$env_data[[2]]),
                           args = 'outputformat=cloglog')
  
  pred_t2_df = spec_list$extra_prepped$extra_prepped_t2 %>%
    dplyr::select(longitude, latitude) %>%
    cbind(pred_t2) %>%
    as_tibble() 
  
  colnames(pred_t2_df) = c("value", "x", "y")
  
  # Extracting the species name
  species = str_to_sentence(str_replace(str_remove(str_remove(species_rds_out, 
                                                              "_output.rds"), 
                                                   "./output/"), 
                                        "_", " ")
  ) 
  
  # Plot building
  # Panel 1
  g1 = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=pred_t1, aes(x=x, y=y, fill=value)) + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.25, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.25, fill = NA) +
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
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=pred_t2, aes(x=x, y=y, fill=value)) + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.25, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
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
  plotname = paste0("./output/", species, ".png")
  ggsave(plotname, gfull)
}

# Testing
spec_list_test = readRDS("./output/leptoties_marina_output.rds")

pred_t1 = dismo::predict(object = spec_list_test$full_mods$full_mod_t1,
                          newdata = spec_list_test$extra_prepped$extra_prepped_t1,
                          x = spec_list_test$prepped_data$env_data[[1]],
                          ext = spec_list_test$prepped_data$env_data[[1]]@extent,
                          type = "cloglog")


pred_t1_df = spec_list_test$extra_prepped$extra_prepped_t1 %>%
  dplyr::select(longitude, latitude) %>%
  cbind(pred_t1) %>%
  as_tibble() 

colnames(pred_t1_df) = c("x", "y", "value")


species = str_to_sentence(str_replace(str_remove(str_remove("./output/leptoties_marina_output.rds", 
                                "_output.rds"), 
                     "./output/"), 
                     "_", " "))
# calculating grid size
pred_t1_df$x[[3]] - pred_t1_df$x[[2]]
# Plot building
# Panel 1
g1 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
               color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_t1_df, aes(x=x, y=y, fill=value), 
            height = 0.5, width = 0.5) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
               color = "grey50", size = 0.25, fill = NA) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="right") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  theme_nothing(legend = TRUE) +
  coord_quickmap() +
  ggtitle("pre-2000")