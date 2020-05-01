# downsampling and looping test script
# Keaton Wilson
# keatonwilson@me.com
# 2020-04-07

# Packages
library(tidyverse)
library(dismo)

# data
pieris = readRDS("./data/split_data/pieris_rapae.rds")
bv_t1 = readRDS("./data/bioclim_t1.rds")
bv_t2 = readRDS("./data/bioclim_t2.rds")

# will have to modify or manually go through the steps in prep_data.R as it 
# adds background points... we want to downsample before adding bg points. 

# selecting the pieces we want and separating by time
small_data = pieris %>%
  mutate(year = lubridate::year(date)) %>%
  dplyr::select(name, longitude, latitude, date, year) %>%
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

# downsample 10 times
pieris_t2_list = list()
downsample_num = nrow(df_t1)
for(i in 1:10){
  downsample = sample_n(df_t2, downsample_num)
  pieris_t2_list[[i]] = downsample
}

#function to generate background points and add to an existing df
bg_add = function(df){
  bg_1 = dismo::randomPoints(bv_t1_cropped, 10000)
  colnames(bg_1) = c("longitude", "latitude")
  
  df_comb = data.frame(df) %>%
    mutate(pb = 1) %>%
    dplyr::select(pb, longitude, latitude) %>%
    bind_rows(data.frame(bg_1) %>% 
                mutate(pb = 0))  %>%
    mutate(Species = as.integer(pb)) %>%
    dplyr::select(-pb)
  
  return(df_comb)
}

#lapply over the list
pieris_t2_list_with_bg = lapply(pieris_t2_list, bg_add)

# doing the same for t1
pieris_t1_with_bg = bg_add(df_t1)

# changing both to spatial points data frames
# Changing to a spatial points data frame
df_sp_t1 = SpatialPointsDataFrame(pieris_t1_with_bg[,c("longitude","latitude")], 
                                  pieris_t1_with_bg, 
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
df_sp_t1$time_frame = "t1"

spdf_func = function(df){
  df_sp_t2 = SpatialPointsDataFrame(df[,c("longitude","latitude")], 
                                    df, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  df_sp_t2$time_frame = "t2"
  return(df_sp_t2)
}

df_sp_t2_list = lapply(pieris_t2_list_with_bg, spdf_func)

glimpse(df_sp_t2_list)

df_sp_t2_list[[1]]

## Three important objects going forward - 
# 1. t1 data
df_sp_t1
# 2. list of downsampled t2 data
df_sp_t2_list
# 3. cropped environmental rasters
bv_t1_cropped
bv_t2_cropped

# run blockCV
source("./script/functions/run_block_cv.R")
block_obj_t1 = run_block_cv(prepped_data = df_sp_t1, 
                            bv_raster =bv_t1)

block_obj_t2_list = lapply(df_sp_t2_list, run_block_cv, bv_raster = bv_t2)

# prep data round 2
source("./script/functions/prep_data_2.R")
names(bv_t1) = paste0("Bio", seq(1:19))
names(bv_t2) = paste0("Bio", seq(1:19))

prepped_t1 = prep_data_2(data = df_sp_t1, env_raster = bv_t1)

prepped_t2_list = lapply(df_sp_t2_list, prep_data_2, env_raster = bv_t2)

# train test split
source("./script/functions/train_test_split.R")
training_test_split_t1 = train_test_split(extra_prepped_data = prepped_t1, 
                                          blocked_obj = block_obj_t1)

training_test_split_t2_list = list()
for(i in 1:length(prepped_t2_list)){
  training_test_split_t2_list[[i]] = train_test_split(extra_prepped_data = prepped_t2_list[[i]], 
                                                      blocked_obj = block_obj_t2_list[[i]])
}
t2_training_data_only_list = lapply(training_test_split_t2_list, "[[", 1)
t2_test_data_only_list = lapply(training_test_split_t2_list, "[[", 2)
# Modeling
source("./script/functions/model_func.R")
# pulling out just training data


t1_mod = model_func(data = training_test_split_t1$training_data, 
                    env_raster = bv_t1_cropped, 
                    num_cores = 14)


t2_mod_list = list()
for(i in 1:length(t2_training_data_only_list)){
  t2_mod_list[[i]] = model_func(data = t2_training_data_only_list[[i]], 
                                env_raster = bv_t2_cropped, 
                                num_cores = 14)
}

# model selection
source("./script/functions/best_mod.R")
t1_best_mod = best_mod(t1_mod)

t2_best_mod_list = (lapply(t2_mod_list, best_mod))

# evaluate models
source("./script/functions/evaluate_models.R")
evaluation_t1 = evaluate_models(test_data = training_test_split_t1$test_data, 
                                model = t1_best_mod[[1]], 
                                env_raster = bv_t1_cropped)

evaluation_t2_list = list()
for(i in 1:length(t2_best_mod_list)){
  evaluation_t2_list[[i]] = evaluate_models(test_data = t2_test_data_only_list[[i]], 
                                            model = t2_best_mod_list[[i]][[1]], 
                                            env_raster = bv_t2_cropped)
}


# full model
source("./script/functions/full_model.R")
names(bv_t1_cropped) = paste0("Bio", seq(1:19))
names(bv_t2_cropped) = paste0("Bio", seq(1:19))

auc_mod = t1_mod@results[t1_best_mod[[2]],]                         
FC_best = tolower(as.character(auc_mod$fc[1]))  
rm_best = as.numeric(auc_mod$rm)
best_mod_t1 = maxnet(p = prepped_t1$Species, data = prepped_t1[,1:19],
                     maxnet.formula(prepped_t1$Species, 
                                    prepped_t1[,1:19], 
                                    classes = FC_best),
                     regmult = rm_best)

best_mods_t2_list = list()
for(i in 1:length(evaluation_t2_list)){
  auc_mod = t2_mod_list[[i]]@results[t2_best_mod_list[[i]][[2]],]                         
  FC_best = tolower(as.character(auc_mod$fc[1]))  
  rm_best = as.numeric(auc_mod$rm)
  best_mods_t2_list[[i]] = maxnet(p = prepped_t2_list[[i]]$Species, data = prepped_t2_list[[i]][,1:19],
                                  maxnet.formula(prepped_t2_list[[i]]$Species, 
                                                 prepped_t2_list[[i]][,1:19], 
                                                 classes = FC_best),
                                  regmult = rm_best)
}

full_mod_t1 = readRDS("./output/best_mod_t1_pieris_looping.rds")
full_mods_t2 = readRDS("./output/best_mods_t2_pieris_looping.rds")



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

# Create new data to predict on
newdata_t1 = as(bv_t1_cropped, "SpatialPixelsDataFrame")
newdata_t1 = as.data.frame(newdata_t1) %>%
  drop_na()

newdata_t2 = as(bv_t2_cropped, "SpatialPixelsDataFrame")
newdata_t2 = as.data.frame(newdata_t2) %>%
  drop_na()

names = c(paste0("Bio", seq(1:19)), "x", "y")
names(newdata_t1) = names
names(newdata_t2) = names

# predictions
pred_t1 = dismo::predict(object = full_mod_t1,
                         newdata = newdata_t1,
                         x = bv_t1_cropped,
                         ext = extent(newdata_t1),
                         type = "cloglog")

pred_t2_list = list()
for(i in 1:length(full_mods_t2)){
  pred_t2_list[[i]] = dismo::predict(object = full_mods_t2[[i]], 
                                     newdata = newdata_t2, 
                                     x = bv_t2_cropped,
                                     ext = extent(newdata_t2), 
                                     type = "cloglog")
}


pred_t1_df = newdata_t1 %>%
  dplyr::select(x, y) %>%
  cbind(as.data.frame(pred_t1)) %>%
  as_tibble() 

colnames(pred_t1_df) = c("x", "y", "value")

pred_t2_df_list = list()
for(i in 1:length(pred_t2_list)){
  df = newdata_t2 %>%
    dplyr::select(x,y) %>%
    cbind(as.data.frame(pred_t2_list[[i]])) %>%
    as_tibble()
  colnames(df) = c("x", "y", "value")
  
  pred_t2_df_list[[i]] = df
}

# evaluations and filtering
eval_t1 = readRDS("./output/evaluation_t1_pieris_looping.rds")
eval_t2_list = readRDS("./output/evaluations_t2_pieris_looping.rds")

thresh_t1 = threshold(eval_t1, 'spec_sens')

thresh_t2_list = lapply(eval_t2_list, threshold, 'spec_sens')

# filtering
t1_threshold = pred_t1_df %>%
  filter(value > thresh_t1)

t2_thresholded_list = list()
for(i in 1:length(thresh_t2_list)){
  t2_thresholded_list[[i]] = pred_t2_df_list[[i]] %>% filter(value > thresh_t2_list[[i]])
}

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
  geom_point(data = df_t1, aes(x = longitude, y = latitude), size = 0.3, 
             color = "white", 
             alpha = 0.7) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="right") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  theme_nothing(legend = TRUE) +
  coord_quickmap() +
  ggtitle("pre-2000")


post_2000_list = list()
for(i in 1:length(t2_thresholded_list)){
  post_2000_list[[i]] = ggplot() +  
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color=NA, size=0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = NA, size = 0.25, fill = "#440154FF") +
    geom_tile(data=t2_thresholded_list[[i]], aes(x=x, y=y), fill = "lightgrey") + 
    geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
                 color="grey50", size=0.20, fill = NA) +
    geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.20, fill = NA) +
    geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
                 color = "grey50", size = 0.20, fill = NA) +
    geom_point(data = pieris_t2_list[[i]], aes(x = longitude, y = latitude), size = 0.3, 
               color = "white", 
               alpha = 0.7) +
    scale_fill_viridis(name = "Probability of Occurence") +
    theme(legend.position="right") +
    theme(legend.key.width=unit(2, "cm"),
          plot.title = element_text(hjust = 0.5, size = 16)) +
    theme_nothing(legend = TRUE) +
    coord_quickmap() +
    ggtitle(paste0("post-2000 iteration: ", i))
}


test = ggarrange(g1, plotlist = post_2000_list)
ggsave(plot = test, filename = "./output/pieris_test.png", width = 11, height = 8.5, units = "in")

# trying to summarize into 1 model
dim(pred_t1_df)
dim(pred_t2_df_list[[1]])

# generating a column that changes to a binary prediction based on threshold
thresh_t1

pred_t1_df_binary = pred_t1_df %>%
  mutate(pred = ifelse(value >= thresh_t1, 1, 0))

pred_t2_df_binary_list = list()
for(i in 1:length(pred_t2_df_list)){
  pred_t2_df_binary_list[[i]] = pred_t2_df_list[[i]] %>%
    mutate(pred = ifelse(value >= thresh_t2_list[[i]], 1, 0))
}

all_pred_t2 = bind_rows(pred_t2_df_binary_list, .id = "iter")
summary_pred_t2 = all_pred_t2 %>%
  group_by(x,y) %>%
  summarize(pred_stack = sum(pred)/10) %>%
  ungroup()
  
# summarizing predictions
g2 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
               color = NA, size = 0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
               color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=summary_pred_t2, aes(x=x, y=y, fill = pred_stack)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.20, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), 
               color = "grey50", size = 0.20, fill = NA) +
  geom_polygon(data = simple_map_mex, aes(x = long, y = lat, group = group), 
               color = "grey50", size = 0.20, fill = NA) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="right") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  theme_nothing(legend = TRUE) +
  coord_quickmap() +
  ggtitle("Summarized Stack of T2")
  
ggsave("./output/summarized_threshold_stack_pieris.png", plot = g2)




