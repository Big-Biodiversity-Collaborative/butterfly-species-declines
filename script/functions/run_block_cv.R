# Running blockCV with a pre-defined grid for this project
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(blockCV)


#' Running blockCV with a preset config for this project
#'
#' @param prepped_data The prepped spatial points dataframe created by 
#' \code{link{prep_data}}
#' @param bv_raster the cropped raster associated with the same time period 
#' as the prepped_data above. 
#'
#' @return a blockCV object that we will use in later analysis
#'
#' @examples
run_block_cv = function(prepped_data, bv_raster){
  
  blocked = spatialBlock(speciesData = prepped_data,
                         species = "Species",
                         rasterLayer = bv_raster,
                         theRange = 400000,
                         k = 5, 
                         selection = "random", 
                         iteration = 250, 
                         biomod2Format = TRUE, 
                         xOffset = 0, 
                         yOffset = 0, 
                         progress = T, 
                         showBlocks = F
  )
  return(blocked)
}
