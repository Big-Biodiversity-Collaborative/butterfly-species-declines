# simple pull of worldclim/terraclim summary data
# Keaton Wilson
# keatonwilson@me.com
# 2020-04-16

# packages
library(raster)

# raster brick of all bioclim
files = list.files("./data/bioclim/", full.names = TRUE)
bioclim_full = stack(files)
names(bioclim_full) = paste0("Bio", seq(1:19))

saveRDS(bioclim_full, "./data/bioclim_full.rds")
