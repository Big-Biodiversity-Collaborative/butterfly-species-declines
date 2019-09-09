#New take on get_bioclim
#Keaton Wilson
#keatonwilson@me.com
#2019-09-08

#Packages
library(tidyverse)
library(RNetCDF)
library(dismo)
library(raster)

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
}