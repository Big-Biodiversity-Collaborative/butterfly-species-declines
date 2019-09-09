# Reworking terraclim function - need to iterate and build biovar chunks 
# and then bind together. Will help reduce the memory requirements for these 
# giant arrays. 
# Keaton Wilson
# keatonwilson@me.com
# 2019-09-03

# packages
library(tidyverse)
library(RNetCDF)
library(abind)


get_bioclim = function(lat_range, lon_range) {
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
  
  #unit test to make sure above works
  print(lat_list)
  print(lon_list)
  
  
  # Iterating through blocks and pulling terraclim data
  
  # Initializing master list to feed everything into
  master_bioclim_list = list()
  
  for(n in 1:8){
    #iterating through each lat range in the outer loop
    lat.range = lat_list[[n]]
    
    # initializing the sub_list
    
    # iterating through each lon within each lat for the inner loop
    # initalizing list that contains three items: an array of data for each 
    # environmental variable
    for(l in 1:8) {
      
      env_vars = c("ppt", "tmax", "tmin")
      
      
      temp_list = list()
      for(j in 1:3){
      
        #iterating through environmental variables    
        env_variable = env_vars[j]
        
        lon.range = lon_list[[l]]
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
        #Splitting into time chunks
        #


      temp_list[[1]] = ifelse(temp_list[[1]] > 10000, NA, temp_list[[1]])

      prcp_t1 = temp_list[[1]][,,1:492]
      prcp_t2 = temp_list[[1]][,,493:732]

      print(prcp_t1)
      print(prcp_t2)
#       tmin_t1 = temp_list[[3]][,,1:492]
#       tmin_t2 = temp_list[[3]][,,493:732]
# 
#       tmax_t1 = temp_list[[2]][,,1:492]
#       tmax_t2 = temp_list[[2]][,,493:732]
# 
#       #test
#       biovar_t1 = biovars(prec = prcp_t1,
#                           tmin = tmin_t1,
#                           tmax = tmax_t1)
# 
#       print(biovar_t1)
  }
  }
}


#Tests
get_bioclim(lat_range = c(50,54), lon_range = c(-81, -78))
