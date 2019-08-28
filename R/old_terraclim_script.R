#Creating biovars from terraclim for both time periods
#Keaton Wilson
#keatonwilson@me.com
#2019-04-26

#libraries
library(dismo)
library(tidyverse)
library(raster)

#Each of these files is an array. 4 km grid of each variable (prcp, tmin, tmax), with 720 layers - months between 1958 and 2018. 
#First thing is to slice them into t1 and t2 chunks. 

#41 years between 1958 and 1999. 
41*12

#Reading in the files
prcp = readRDS("./data/terraclim/prcp.rds")
#Getting rid of weirdly high values in prcp
prcp = ifelse(prcp > 10000, NA, prcp)

prcp_t1 = prcp[,,1:492]
prcp_t2 = prcp[,,493:720]
rm(prcp)

tmin = readRDS("./data/terraclim/tmin.rds")
tmin_t1 = tmin[,,1:492]
tmin_t2 = tmin[,,493:720]
rm(tmin)

tmax = readRDS("./data/terraclim/tmax.rds")
tmax_t1 = tmax[,,1:492]
tmax_t2 = tmax[,,493:720]
rm(tmax)


lat.range=c(25, 55)        #! Ranges instead of point values. Order does not matter
lon.range=c(-94, -65)
#Changing to rasterbricks
prcp_t1_raster = t(brick(prcp_t1, crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         xmn = 25,
                         xmx = 55,
                         ymn = -94,
                         ymx = -65))
tmin_t1_raster = t(brick(tmin_t1, crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         xmn = 25,
                         xmx = 55,
                         ymn = -94,
                         ymx = -65))
tmax_t1_raster = t(brick(tmax_t1, crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         xmn = 25,
                         xmx = 55,
                         ymn = -94,
                         ymx = -65))

##Great. This works - but it outputs a rasterbrick, not sure how that will work with stuff feeding into a raster stack below. 
biovar_list_t1 = list() # Generating empty list to feed into
length = dim(prcp_t1_raster)[3]/12 #How many times are we going to cycle through?
seq = 1:12 #initalizing sequence

for(i in 1:length) {
  precip_sub = prcp_t1_raster[[seq]]
  tmin_sub = tmin_t1_raster[[seq]]
  tmax_sub = tmax_t1_raster[[seq]]
  
  biovar_list_t1[[i]] = biovars(prec = precip_sub,
                             tmin = tmin_sub,
                             tmax = tmax_sub)
  seq = seq + 12
  print(seq)
}

#Workflow 
#pull out a list of raster layers for each bioclim variable
#turn those layers into a rasterStack
#Compute average
#Recombine

biovar_avg_combined_t1 = raster::brick(nl = 19, nrows = 720, ncols = 696, 
                                       xmn = -94, 
                                       xmx = -65,
                                       ymn = 25, 
                                       ymx = 55)
for(i in 1:19) {
  biovar_sublist = lapply(biovar_list_t1, '[[', i) #pulls out each bioclim variable iteratively
  biovar_substack = stack(biovar_sublist) #combines all years into a raster stack
  biovar_avg = calc(biovar_substack, fun = mean) #Calculates the average for each var
  biovar_avg_combined_t1[[i]] = biovar_avg
}


#Changing to rasterbricks
prcp_t2_raster = t(brick(prcp_t2, crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         xmn = 25,
                         xmx = 55,
                         ymn = -94,
                         ymx = -65))
tmin_t2_raster = t(brick(tmin_t2, crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         xmn = 25,
                         xmx = 55,
                         ymn = -94,
                         ymx = -65))
tmax_t2_raster = t(brick(tmax_t2, crs ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
                         xmn = 25,
                         xmx = 55,
                         ymn = -94,
                         ymx = -65))

##Great. This works - but it outputs a rasterbrick, not sure how that will work with stuff feeding into a raster stack below. 
biovar_list_t2 = list() # Generating empty list to feed into
length = dim(prcp_t2_raster)[3]/12 #How many times are we going to cycle through?
seq = 1:12 #initalizing sequence

for(i in 1:length) {
  precip_sub = prcp_t2_raster[[seq]]
  tmin_sub = tmin_t2_raster[[seq]]
  tmax_sub = tmax_t2_raster[[seq]]
  
  biovar_list_t2[[i]] = biovars(prec = precip_sub,
                             tmin = tmin_sub,
                             tmax = tmax_sub)
  seq = seq + 12
  print(seq)
}

#Workflow 
#pull out a list of raster layers for each bioclim variable
#turn those layers into a rasterStack
#Compute average
#Recombine

biovar_avg_combined_t2 = raster::brick(nl = 19, nrows = 720, ncols = 696, 
                                       xmn = -94, 
                                       xmx = -65,
                                       ymn = 25, 
                                       ymx = 55)
for(i in 1:19) {
  biovar_sublist = lapply(biovar_list_t2, '[[', i) #pulls out each bioclim variable iteratively
  biovar_substack = stack(biovar_sublist) #combines all years into a raster stack
  biovar_avg = calc(biovar_substack, fun = mean) #Calculates the average for each var
  biovar_avg_combined_t2[[i]] = biovar_avg
}

writeRaster(biovar_avg_combined_t1, "./data/terraclim/biovar_avg_t1.grd", overwrite = TRUE)
writeRaster(biovar_avg_combined_t2, "./data/terraclim/biovar_avg_t2.grd", overwrite = TRUE)

test = raster::brick("./data/terraclim/biovar_avg_t1.grd")
test2 = raster::brick("./data/terraclim/biovar_avg_t2.grd")
