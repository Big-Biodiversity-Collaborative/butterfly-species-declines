#Terraclim Script 
#Keaton Wilson
#keatonwilson@me.com
#2019-04-26
#
#packages
library(RNetCDF)

#Enter lat and lon ranges
lat.range=c(5, 66)        #! Ranges instead of point values. Order does not matter
lon.range=c(-166, -120)

# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html

var="ppt"

baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")

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

data <-var.get.nc(nc, variable = var,start = start, count,unpack=TRUE)    #! argument change: 'variable' instead of 'varid'  # Output is now a matrix
data_prcp = data

saveRDS(data_prcp, "./data/terraclim/prcp.rds")

# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html

var="tmin"

baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")

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

data_tmin <-var.get.nc(nc, variable = var,start = start, count,unpack=TRUE)    #! argument change: 'variable' instead of 'varid'  # Output is now a matrix
saveRDS(data_tmin, "./data/terraclim/tmin.rds")

# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html

var="tmax"

baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")

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

data_tmax <-var.get.nc(nc, variable = var,start = start, count,unpack=TRUE)    #! argument change: 'variable' instead of 'varid'  # Output is now a matrix

saveRDS(data_tmax, "./data/terraclim/tmax.rds")
