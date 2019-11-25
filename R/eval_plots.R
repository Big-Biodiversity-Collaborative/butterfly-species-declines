# Making evaluation plots that summarize hyperparameter tuning within Maxent
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25

# packages
require(ENMeval)

eval_plots = function(eval_object = NULL) {
  par(mfrow=c(2,3))
  eval.plot(eval_object@results)
  eval.plot(eval_object@results, 'auc.test.avg', legend = F)
  eval.plot(eval_object@results, 'auc.diff.avg', legend = F)
  eval.plot(eval_object@results, 'or.mtp.avg', legend = F)
  eval.plot(eval_object@results, 'or.10p.avg', legend = F)
  plot(eval_object@results$auc.test.avg, eval_object@results$delta.AICc, bg=as.factor(eval_object@results$fc), pch= 21, cex = eval_object@results$rm/2, xlab = "avg.test.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
  legend("topright", legend=unique(eval_object@results$fc), pt.bg=as.factor(eval_object@results$fc), pch=21)
  mtext("Circle size proportional to regularization multiplier value", cex = 0.6)
  
}