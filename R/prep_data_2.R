# Prep data 2: the preppening
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(tidyverse)
require(raster)


# Preparing data 2 ----------------------------------------------------------
#' More data preparation prior to SDM building
#'
#' @param data Prepped spatial points datframe created by \code{link{prep_data}}
#' @param env_raster the cropped raster associated with the same time period 
#' as the prepped_data above.
#'
#' @return a dataframe with extracted environmental variables along with presence 
#' for all of the occurence and background data
#'
#' @examples
prep_data_2 = function(data, env_raster){
  extra_prepped = raster::extract(env_raster, data, df = TRUE) %>%
    bind_cols(as.data.frame(data)) %>%
    drop_na() %>%
    dplyr::select(-ID, Species, longitude, latitude, Bio1:Bio19) %>%
    filter_all(all_vars(. != -Inf))
  return(extra_prepped)
}
