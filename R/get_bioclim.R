#New take on get_bioclim
#Keaton Wilson
#keatonwilson@me.com
#2019-09-08

#Packages
library(tidyverse)
library(RNetCDF)
library(dismo)
library(raster)
library(progress)

# lat-lon range splitting function ----------------------------------------

#Function to split lat-lon range into 8x8 grid to reduce size for memory
split_lat_lon = function(lat_range, lon_range) {
  #Breaking the lat range into 8 chunks
  lat_diff = (max(lat_range)-min(lat_range))/8
  
  #initializing a list to put the chunks into
  
  lat_list = list()
  min = min(lat_range)
  
  #looping through 8 chunks
  for (i in 1:8) {
    lat_range_temp = c(min, min+lat_diff)
    lat_list[[i]] = lat_range_temp
    min = min+lat_diff
  }
  
  #breaking the lon range into 8 chunks
  lon_diff = (max(lon_range)-min(lon_range))/8
  
  #initializing a list to put the chunks into
  
  lon_list = list()
  min = min(lon_range)
  
  
  #looping through 8 chunks
  for (i in 1:8) {
    lon_range_temp = c(min, min+lon_diff)
    lon_list[[i]] = lon_range_temp
    min = min+lon_diff
  }
  
  return(master_lat_lon_list = list(lat_list, lon_list))
}


# Terraclim Pull ----------------------------------------------------------
#Function to get all three environmental variables from terraclim based on 
#lat-lon range

pull_terraclim = function(lat_range, lon_range){

    env_vars = c("ppt", "tmax", "tmin")
    temp_list = list()
    
    for(j in 1:3){
      
      #iterating through environmental variables    
      env_variable = env_vars[j]
      
      lon.range = lon_range
      lat.range = lat_range
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",env_variable),"_1958_CurrentYear_GLOBE.nc")
      
      nc <- open.nc(baseurlagg)
      lon <- var.get.nc(nc, "lon")
      lat <- var.get.nc(nc, "lat")
      lat.range <- sort(lat.range)                              
      lon.range <-sort(lon.range)
      lat.index <- which(lat>=lat.range[1]&lat<=lat.range[2])    
      lon.index <- which(lon>=lon.range[1]&lon<=lon.range[2])    
      lat.n <- length(lat.index)                                
      lon.n <- length(lon.index)
      start <- c(lon.index[1], lat.index[1], 1)
      count <- c(lon.n, lat.n, NA)                            
      
      # read in the full period of record using aggregated files
      
      temp_list[[j]] <- var.get.nc(nc, variable = env_variable,start = start, count,unpack=TRUE)    #! argument change: 'variable' instead of 'varid'  # Output is now a matrix
      
    }
    return(temp_list)
}


# Separating years -------------------------------------------

env_var_year_split = function(terraclim_data, year_split) {
  #Getting rid of weirdly high values in prcp
  terraclim_data = ifelse(terraclim_data > 10000, NA, terraclim_data)
  
  #
  months = (year_split-1958)*12
  
  terraclim_data_t1 = terraclim_data[,,1:months]
  terraclim_data_t2 = terraclim_data[,,(months+1):732]
  
  terraclim_data_t1_t2_list = list(terraclim_data_t1, terraclim_data_t2)
  return(terraclim_data_t1_t2_list)
}



# Bioclim Var Calculation -------------------------------------------------

bioclim_calc = function(prcp, tmax, tmin, lat_range, lon_range) {
  #rasterizing
  prcp_raster = (brick(prcp, 
                             crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", 
                             xmn = min(lon_range),
                             xmx = max(lon_range),
                             ymn = min(lat_range),
                             ymx = max(lat_range), transpose = TRUE))
  tmax_raster = (brick(tmax, 
                        crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", 
                        xmn = min(lon_range),
                        xmx = max(lon_range),
                        ymn = min(lat_range),
                        ymx = max(lat_range), transpose = TRUE))
  tmin_raster = (brick(prcp, 
                        crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0", 
                        xmn = min(lon_range),
                        xmx = max(lon_range),
                        ymn = min(lat_range),
                        ymx = max(lat_range), transpose = TRUE))
  
  biovar_list = list() # Generating empty list to feed into
  length = dim(prcp_raster)[3]/12 #How many times are we going to cycle through?
  seq = 1:12 #initalizing sequence
  
  for(i in 1:length) {
    precip_sub = prcp_raster[[seq]]
    tmin_sub = tmin_raster[[seq]]
    tmax_sub = tmax_raster[[seq]]
    
    biovar_list[[i]] = biovars(prec = precip_sub,
                                  tmin = tmin_sub,
                                  tmax = tmax_sub)
    seq = seq + 12
  }
                       
return(biovar_list)       
  
}


# Averaging over entire time period ---------------------------------------

bioclim_averaging = function(biovar_list, nrows, ncols, lat_range, lon_range){
  biovar_avg_combined = raster::brick(nl = 19, nrows = nrows, ncols = ncols, 
                                       xmn = min(lon_range),
                                       xmx = max(lon_range),
                                       ymn = min(lat_range),
                                       ymx = max(lat_range))
    for(i in 1:19) {
      biovar_sublist = lapply(biovar_list, '[[', i) #pulls out each bioclim variable iteratively
      biovar_substack = stack(biovar_sublist) #combines all years into a raster stack
      biovar_avg = calc(biovar_substack, fun = mean) #Calculates the average for each var
      biovar_avg_combined[[i]] = biovar_avg
    }
  return(biovar_avg_combined)
}

# Master ------------------------------------------------------------------
pb = progress_bar$new(total = 64, format = "    working [:bar] :percent eta: :eta", 
                      clear = FALSE, width = 60)

get_bioclim = function(lat_range, lon_range, year_split) {
  pb$tick(0)
  # breaking up into chunks
  lat_lon_chunks = split_lat_lon(lat_range, lon_range)
  
  # initializing list to put stuff into
  bioclim_final = list()
  
  
  
  for(i in 1:8){
    sub_lat_range = lat_lon_chunks[[1]][[i]]
    
    sub_final = list()
    
    for(g in 1:8){
    sub_lon_range = lat_lon_chunks[[2]][[g]]
    
    temp_terraclim = pull_terraclim(lat_range = sub_lat_range,
                                    lon_range = sub_lon_range)
    
    
    year_split_terraclim = lapply(temp_terraclim, env_var_year_split, year_split = year_split)
    
    
    
    bioclims_t1 = bioclim_calc(prcp = year_split_terraclim[[1]][[1]],
                               tmax = year_split_terraclim[[2]][[1]],
                               tmin = year_split_terraclim[[3]][[1]], 
                               lat_range = sub_lat_range, 
                               lon_range = sub_lon_range)
    
    
    
    bioclims_t2 = bioclim_calc(prcp = year_split_terraclim[[1]][[2]],
                               tmax = year_split_terraclim[[2]][[2]],
                               tmin = year_split_terraclim[[3]][[2]], 
                               lat_range = sub_lat_range, 
                               lon_range = sub_lon_range)
    
   
    
    dims_1 = dim(bioclims_t1[[1]])
    dims_2 = dim(bioclims_t2[[1]])
    
    avg_t1 = bioclim_averaging(biovar_list = bioclims_t1, 
                               nrows = dims_1[1], 
                               ncols = dims_1[2], 
                              lat_range = sub_lat_range, 
                              lon_range = sub_lon_range)
    
    avg_t2 = bioclim_averaging(biovar_list = bioclims_t2, 
                               nrows = dims_2[1], 
                               ncols = dims_2[2], 
                               lat_range = sub_lat_range, 
                               lon_range = sub_lon_range)
    
    avg_list = list(avg_t1, avg_t2)
    
    sub_final[[g]] = avg_list
    pb$tick(1)
    }
    
    bioclim_final[[i]] = sub_final
  }
  #Splitting out t1 and t2 - two nested lists each with 64 individual chunks
  t1 = lapply(bioclim_final, function(x) lapply(x, '[[', 1))
  t2 = lapply(bioclim_final, function(x) lapply(x, '[[', 2))

  #unlisting
  t1_full = unlist(t1)
  t2_full = unlist(t2)

  #mergning
  t1_apogee = do.call(merge, t1_full)
  t2_apogee = do.call(merge, t2_full)

  final_list = list(t1_apogee, t2_apogee)
  return(final_list)
}

# Testing -----------------------------------------------------------------

#Testing on small lat/lon range by the coast to make sure the transpose issue is functioning
# florida_test = get_bioclim(lat_range = c(25, 32), lon_range = c(-90, -80))
# 
# florida_t1 = lapply(florida_test, function(x) lapply(x, '[[', 1))
# florida_t1_full = unlist(florida_t1)
# florida_t1_merged = do.call(merge, florida_t1_full)
# plot(florida_t1_merged)
# 
# florida_transposed = lapply(florida_t1_full, t)
# florida_flipped = lapply(florida_transposed, flip, 2)
# florida_t1_merged = do.call(merge, florida_flipped)
# plot(florida_t1_merged)
# 
# 
# # Bigger problems testing
# terraclim_dat = pull_terraclim(lat_range = c(25, 26), lon_range = c(-81, -80))
# year_split_terraclim_dat = lapply(terraclim_dat, env_var_year_split)
# bioclims_1 = bioclim_calc(prcp = year_split_terraclim_dat[[1]][[2]],
#                           tmax = year_split_terraclim_dat[[2]][[2]],
#                           tmin = year_split_terraclim_dat[[3]][[2]], 
#                           lat_range = c(25,26), 
#                           lon_range = c(-81, -80))
# dims_1 = dim(bioclims_1[[1]])
# avg_1 = bioclim_averaging(biovar_list = bioclims_1, 
#                         nrows = dims_1[1], 
#                         ncols = dims_1[2], 
#                         lat_range = c(25, 26), 
#                         lon_range = c(-81, -80))
# library(maps)
# library(mapdata)
# library(rgdal)
# library(fields)
# states = map_data("state")
# 
# florida = states %>%
#   filter(region == "florida")
# 
# ggplot(data = florida) +
#   geom_polygon(aes(x = long, y = lat)) +
#   geom_raster(data = avg_1$layer.1.1.1.1)
# 
# plot((avg_1$layer.1.1.1.1))
# US(add=TRUE)
# plot(avg_2$layer.1.1.1.1, add=TRUE)
# 
# 
# #Part 2
# terraclim_dat = pull_terraclim(lat_range = c(25, 26), lon_range = c(-82, -81))
# year_split_terraclim_dat = lapply(terraclim_dat, env_var_year_split)
# bioclims_2 = bioclim_calc(prcp = year_split_terraclim_dat[[1]][[2]],
#                           tmax = year_split_terraclim_dat[[2]][[2]],
#                           tmin = year_split_terraclim_dat[[3]][[2]], 
#                           lat_range = c(25,26), 
#                           lon_range = c(-82, -81))
# dims_1 = dim(bioclims_2[[1]])
# avg_2 = bioclim_averaging(biovar_list = bioclims_2, 
#                           nrows = dims_1[1], 
#                           ncols = dims_1[2], 
#                           lat_range = c(25, 26), 
#                           lon_range = c(-82, -81))

#test
# t = get_bioclim(lat_range = c(65, 66), lon_range = c(-101, -100), year_split = 2000)

# big_bioclim_1 = get_bioclim(lat_range = c(15, 66), lon_range = c(-140, -100))
# saveRDS(big_bioclim_1, "./data/big_bioclim_1.RDS")
big_bioclim_1 = readRDS("./data/big_bioclim_1.RDS")
big_bioclim_2 = get_bioclim(lat_range = c(15, 66), lon_range = c(-100, -40))
saveRDS(big_bioclim_2, "./data/big_bioclim_2.RDS")

#merging
bioclim_t1 = raster::merge(big_bioclim_1[[1]], big_bioclim_2[[1]])
bioclim_t2 = raster::merge(big_bioclim_1[[2]], big_bioclim_2[[2]])

saveRDS(bioclim_t1, "./data/bioclim_t1.rds")
saveRDS(bioclim_t2, "./data/bioclim_t2.rds")

