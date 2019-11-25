# Evaluating the best model on out-of-sample test data
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

require(raster)
require(tidyverse)
require(dismo)

# evalute function
evaluate_models = function(test_data, model, env_raster) {
  test_data_occ = test_data %>%
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude)
  
  bg_data = test_data %>%
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude)
  
  ev = evaluate(test_data_occ, a = bg_data, model = model, x = env_raster, type = 'cloglog')
  return(ev)
}
