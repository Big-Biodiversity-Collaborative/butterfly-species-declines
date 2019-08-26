#Script to get and combine big terraclim data
#Keaton Wilson
#keatonwilson@me.com
#2019-08-26

#RNetCDF
library(RNetCDF)
library(tidyverse)
library(abind)


#A note - don't put the file (in the file argument below in the terraclim folder - as this get's cleaned out once the function finishes buildin)

get_terraclim = function(lat_range, lon_range, env_variable, file=NULL) {
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
  
 #going to loop through the blocks and write each chunk to disk
 
  for(n in 1:8){
    #iterating through each lat range in the outer loop
    lat.range = lat_list[[n]]
    
    #initializing the sub_list
    
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
      save_name = paste0("./data/terraclim/", paste(paste("chunk", n, l, sep = "_"), "rds", sep = "."))
      saveRDS(data, file = save_name)
      rm(data)
    }
  }
  #unit test to see if above is working
  print(list.files("./data/terraclim/"))
  # return(ppt_list)

  #Binding the chunks together across longitude first and writing to disk
  file_list = list.files("./data/terraclim", full.names = TRUE)
  
  for(i in 1:8) {
    #Setting up filenames for each "row"
    str_search = file_list[str_detect(file_list, paste0("_", i, "_"))]
    
    #initializing list
    chunk_list = list()
    
    #iterating
      for(a in 1:8) {
        chunk_list[[a]] = readRDS(str_search[i])
      }
    #Binding columns together 
    big_wide_chunk = abind(chunk_list, along=2)
    
    
    saveRDS(big_wide_chunk, paste0("./data/terraclim/", paste(paste("wide_chunk", i, sep = "_"), "rds", sep = ".")))
    rm(big_wide_chunk)
  }
  #Binding all the wide chunks together and writing to disk
  
  #initializing list
  wide_chunk_list = list()
  
  #Making names
  wide_chunk_names = list.files("./data/terraclim", full.names = TRUE)
  wide_chunk_names = wide_chunk_names[str_detect(wide_chunk_names, "wide")]
  
  #reading in
  for(i in 1:8){
    wide_chunk_list[[i]] = readRDS(wide_chunk_names[i])
  }
  
  #binding
  glued = abind(wide_chunk_list, along=1)
  rm(wide_chunk_list)
  print(str(glued))
  saveRDS(glued, file)
  rm(glued)
  
  #Cleaning up
  file.remove(list.files("./data/terraclim", full.names = TRUE))
}

#Testing
# lat.range=c(5, 66)        #! Ranges instead of point values. Order does not matter
# lon.range=c(-166, -120)
# 
# test = get_terraclim(env_variable = "ppt", lat_range = c(40, 66), lon_range = c(-166, -160))
# 
big_test = get_terraclim(env_variable = "ppt", lat_range = c(15, 66), lon_range = c(-140, -60), file = "./data/prcp.rds")
