#Script to get and combine big terraclim data
#Keaton Wilson
#keatonwilson@me.com
#2019-08-21

#RNetCDF
library(RNetCDF)
library(tidyverse)
library(abind)


get_terraclim = function(lat_range, lon_range, env_variable, path=NULL) {
  #Breaking the lat range into 4 chunks
  lat_diff = (max(lat_range)-min(lat_range))/8
  
  #initializing a list to put the chunks into
  
  lat_list = list()
  min = min(lat_range)
  
  #looping through 4 chunks
  for (i in 1:8) {
    lat_range_temp = c(min, min+lat_diff)
    lat_list[[i]] = lat_range_temp
    min = min+lat_diff
  }
  
  #breaking the lon range into 4 chunks
  lon_diff = (max(lon_range)-min(lon_range))/8
  
  #initializing a list to put the chunks into
  
  lon_list = list()
  min = min(lon_range)
  
  #looping through 4 chunks
  for (i in 1:8) {
    lon_range_temp = c(min, min+lon_diff)
    lon_list[[i]] = lon_range_temp
    min = min+lon_diff
  }
  
  #unit test to make sure above works
  print(lat_list)
  print(lon_list)
  
  #Terraclim
  
  #initializing a list to pipe into
  ppt_list = list()
  
  for(n in 1:8){
    #iterating through each lat range in the outer loop
    lat.range = lat_list[[n]]
    
    #initializing the sub_list
    ppt_sub_list = list()
    
    #iterating through each lon within each lat for the inner loop
    for(l in 1:8) {
      lon.range = lon_list[[l]]
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",env_variable),"_1958_CurrentYear_GLOBE.nc")
      
      nc <- open.nc(baseurlagg)
      lon <- var.get.nc(nc, "lon")
      lat <- var.get.nc(nc, "lat")
      lat.range <- sort(lat.range)                              #!sort user input values from low to high
      lon.range <-sort(lon.range)
      lat.index <- which(lat>=lat.range[1]&lat<=lat.range[2])    #! index values within specified range
      lon.index <- which(lon>=lon.range[1]&lon<=lon.range[2])    
      lat.n <- length(lat.index)                                #!value for count
      lon.n <- length(lon.index)
      start <- c(lon.index[1], lat.index[1], 1)
      count <- c(lon.n, lat.n, NA)                            #! parameter change: 'NA' instead of '-1' to signify entire dimension
      
      
      # read in the full period of record using aggregated files
      
      data <-var.get.nc(nc, variable = env_variable,start = start, count,unpack=TRUE)    #! argument change: 'variable' instead of 'varid'  # Output is now a matrix
  ppt_sub_list[[l]] = data
  
    }
    #binding each sub list into the main list
    ppt_list[[n]] = ppt_sub_list

  }
  #unit test to see if above is working
  print(str(ppt_list))
  # return(ppt_list)
  
  #Binding the matrices together
  # Binding each latitude range together
  # Col bind first, then r bind for one big matrix
  test_col_bound = list()
  for(i in 1:8) {
    test_col_bound[[i]] = abind(ppt_list[[i]], along=2)
  }
  #unit test
  print(str(test_col_bound))
  
  glued = abind(test_col_bound, along=1)
  print(str(glued))
  return(glued)
  
  if(is.null(path)) {
    
  } else{
  saveRDS(glued, file = path)
  }
}



#Testing
# lat.range=c(5, 66)        #! Ranges instead of point values. Order does not matter
# lon.range=c(-166, -120)
# 
# test = get_terraclim(env_variable = "ppt", lat_range = c(40, 66), lon_range = c(-166, -160))
# 
big_test = get_terraclim(env_variable = "ppt", lat_range = c(15, 66), lon_range = c(-140, -50), path = "./data/prcp.rds")
