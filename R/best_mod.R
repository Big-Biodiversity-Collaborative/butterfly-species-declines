# Modeling selection - pulls out best model based on test AUC scores
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(ENMeval)

best_mod = function(model_obj){
  best_index = as.numeric(row.names(model_obj@results[which(model_obj@results$auc.test.avg== max(model_obj@results$auc.test.avg)),]))[1]
  
  best_mod = model_obj@models[[best_index]]
  return(list(best_mod, best_index))
}
