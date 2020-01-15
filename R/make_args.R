# Making arguments to feed into 'final' model
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-25


#' make_args - creates appropriate arguments that come from ENMeval to generate 
#' final model based on the parameters of the 'best model' before building the 
#' final model
#'
#' @param RMvalues regularization multiplier values
#' @param fc feature classes
#' @param labels labels as true or false
#'
#' @return
#' @export
#'
#' @examples
make_args <- function(RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), labels=FALSE) {
  
  other.args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
  args.list <- list()
  
  for (i in 1:length(fc)) {
    args.list[[i]] <- other.args
    if(!grepl("L", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "nolinear")
    if(!grepl("Q", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "noquadratic")
    if(!grepl("H", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "nohinge")
    if(!grepl("P", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "noproduct")
    if(!grepl("T", fc[[i]])) args.list[[i]] <- c(args.list[[i]], "nothreshold")
  }
  
  RM.lab <- rep(RMvalues, each=length(fc))
  RM.arg <- paste("betamultiplier=", RM.lab, sep="")
  fc.lab <- rep(fc, times=length(RMvalues))
  fc.arg <- rep(args.list, times=length(RMvalues))
  
  args <- list()
  feats.lab <- c()
  rms.lab <- c()
  for (i in 1:length(fc.lab)) {
    args[[i]] <- c(RM.arg[i], fc.arg[[i]])
    feats.lab <- c(feats.lab, fc.lab[[i]])
    rms.lab <- c(rms.lab, RM.lab[i])
  }
  args.lab <- list(feats.lab, rms.lab)
  
  if(labels==FALSE) {
    return(args)
  } else {
    return(args.lab)
  }
}