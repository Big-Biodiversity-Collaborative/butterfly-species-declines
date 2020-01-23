# Modeling selection - pulls out best model based on test AUC scores
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(ENMeval)

#' Model Selection
#' 
#' \code{best_mod} chooses the best model from the large ENMeval object of all 
#' models, based on the highest AUC value for test data. 
#'
#' @param model_obj A model object generated from the \code{\link{model_func}}.
#'
#' @return A list two elements long. The first item contains the full model 
#' model object of the best model and the second item is the index of the best 
#' model to reference back to all models in the full ENMeval object. 
#'
#' @examples

best_mod = function(model_obj){
  best_index = as.numeric(row.names(model_obj@results[which(model_obj@results$auc.test.avg== max(model_obj@results$auc.test.avg)),]))[1]
  
  best_model = model_obj@models[[best_index]]
  return(list(best_model, best_index))
}
