# R Packages for sdm work

# Normal CRAN installs
install.packages(c('dismo', 'raster', 'tidyverse', 'devtools', 'rJava', 'maxnet'))

# github installs
devtools::install_github('johnbaums/rmaxent')
devtools::install_github('rvalavi/blockCV')
devtools::install_github('bobmuscarella/ENMeval@dev')

# Installing maxent and putting it in the right place
library(rmaxent)
get_maxent()
