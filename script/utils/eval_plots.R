# Making evaluation plots that summarize hyperparameter tuning within Maxent
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(ENMeval)

#'  Generating evaluation plots for model tuning
#' 
#' \code{eval_plots} builds a multi-panel plot that highlights the tuning phase 
#' of model selection within ENMeval. It plots a variety of evaluation metrics 
#' across regularization and feature space.  
#' 
#' @param models_obj an evaluation object generated from the
#'  \code{\link{evaluate_models}} function.
#'
#' @return Generates a multi-panel plot using base R plotting, ready for save 
#' or export.  
#'
#' @examples

eval_plots = function(models_obj = NULL) {
  par(mfrow=c(2,3))
  eval.plot(models_obj@results)
  eval.plot(models_obj@results, 'auc.test.avg', legend = F)
  eval.plot(models_obj@results, 'auc.diff.avg', legend = F)
  eval.plot(models_obj@results, 'or.mtp.avg', legend = F)
  eval.plot(models_obj@results, 'or.10p.avg', legend = F)
  plot(models_obj@results$auc.test.avg, models_obj@results$delta.AICc, bg=as.factor(models_obj@results$fc), pch= 21, cex = models_obj@results$rm/2, xlab = "avg.test.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
  legend("topright", legend=unique(models_obj@results$fc), pt.bg=as.factor(models_obj@results$fc), pch=21)
  mtext("Circle size proportional to regularization multiplier value", cex = 0.6)
  
}