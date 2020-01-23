# Initial Maps
# Keaton Wilson
# keatonwilson@me.com
# 2019-10-14
options(java.parameters = "-Xmx8000m")
# Packages ----------------------------------------------------------------
#packages
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

# Full env rasters
bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")
new_names = paste0("Bio", rep(1:19))
names(bv_t1) = new_names
names(bv_t2) = new_names

# loading in the big objects and binding together into a big ol' list
files = list.files("./output/", full.names = TRUE)

master_list = list()
for(x in 1:length(files)) {
  data = readRDS(files[x])
  master_list = c(master_list, data)
}


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

# #Great lakes issues
# lakes <- rgdal::readOGR("./data/raw_data/ne_10m_lakes/ne_10m_lakes.shp")
# lakes = lakes[lakes$scalerank==0,]
# lakes = crop(lakes, geographic.extent)

# Predictions and maps-------------------------------------------------------------

# iterating through and adding predictions
for(x in 1:length(master_list)){
  pred_t1 = dismo::predict(object = master_list[[x]]$full_mods[[1]],
                           x = bv_t1,
                           ext = master_list[[x]]$env_data[[1]]@extent,
                           args = 'outputformat=cloglog')
  pred_t1 = as(pred_t1, "SpatialPixelsDataFrame")
  pred_t1 = as.data.frame(pred_t1)
  colnames(pred_t1) = c("value", "x", "y")
  
  pred_t2 = dismo::predict(object = master_list[[x]]$full_mods[[2]],
                           x = bv_t2,
                           ext = master_list[[x]]$env_data[[2]]@extent,
                           args = 'outputformat=cloglog')
  pred_t2 = as(pred_t2, "SpatialPixelsDataFrame")
  pred_t2 = as.data.frame(pred_t2)
  colnames(pred_t2) = c("value", "x", "y")
  
  species = names(master_list)[x]
  
  g1 = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=pred_t1, aes(x=x, y=y, fill=value)) + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.25, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
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

# Threshold and point maps
evaluate_models = function(test_data, model, env_raster) {
  test_data_occ = test_data %>%
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude)
  
  bg_data = test_data %>%
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude)
  
  ev = evaluate(test_data_occ, a = bg_data, model = model, x = env_raster, type = 'response')
  return(ev)
}
ev_test = evaluate_models(test_data = master_list$`Parnassius clodius`$train_test_split[[1]][[2]], 
                          model = master_list$`Parnassius clodius`$best_mod[[1]][[1]], 
                          env_raster = master_list$`Parnassius clodius`$env_data[[1]])

ev_test = evaluate(p = master_list$`Pieris rapae`$train_test_split[[1]][[2]] %>%
                     filter(Species == 1), 
                   a = master_list$`Pieris rapae`$train_test_split[[1]][[2]] %>%
                     filter(Species == 0), 
                   model = master_list$`Pieris rapae`$best_mod[[1]][[1]], 
                   x = master_list$`Pieris rapae`$env_data[[1]], type = 'response')

thresh = threshold(master_list$`Parnassius clodius`$evaluations[[2]])
test_pred = dismo::predict(object = master_list[[x]]$full_mods[[1]],
                           x = bv_t1,
                           ext = master_list[[x]]$env_data[[1]]@extent)

test_pred = as(test_pred, "SpatialPixelsDataFrame")
test_pred = as.data.frame(test_pred)
colnames(test_pred) = c("value", "x", "y")

filtered = test_pred %>%
  filter(value < thresh)

ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "grey10") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "grey10") +
  geom_tile(data=filtered, aes(x=x, y=y), fill = "lightgrey") + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_point(data = master_list$`Parnassius clodius`$prepped_dfs[[1]] %>%
               filter(Species == 1), 
             aes(x = longitude, y = latitude), 
             alpha = 0.2, 
             color = "yellow", 
             shape = 3, 
             size = 0.5) +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("2000-2019") +
  coord_quickmap()