# Building full final model on all of the data
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(raster)
require(ENMeval)
require(maxnet)

# function to build out full model on all data after hyperparameter tuning
full_model = function(models_obj, best_model_index, full_data = NULL, env_raster) {
  auc_mod = models_obj@results[best_model_index,]
  FC_best = tolower(as.character(auc_mod$fc[1]))
  rm_best = auc_mod$rm
  
  
  maxent.args = make.args(RMvalues = rm_best, fc = FC_best)
  
  # re calculating environmental raster
  
  
  # full_mod = maxent(env_raster, as.matrix(full_data[,1:2]), args = maxent.args[[1]])
  full_mod = maxnet(p = full_data$Species, data = full_data[,1:19],
                    maxnet.formula(full_data$Species, 
                                   full_data[,1:19], 
                                   classes = FC_best), 
                    regmult = rm_best)
  return(full_mod)
  
}