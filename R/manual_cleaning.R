# Examining output objects and cleaning/fixing
# Keaton Wilson
# keatonwilson@me.com
# 2019-10-14

# packages
library(dismo)
library(raster)
library(tidyverse)
library(blockCV)
library(tidyverse)
library(maxnet)
library(ENMeval)
if (Sys.getenv("JAVA_HOME")!="")
  Sys.setenv(JAVA_HOME="")
library(rJava)

# First object
t1 = readRDS("./output/eurytheme_philodice.rds")
glimpse(t1)
rm(t1)

# looks good no errors
t2 = readRDS("./output/marina_antiopa.rds")
glimpse(t2)

# Full mod not compelted for Nymphalis antiopa t1 - retry
full_model = function(models_obj, best_model_index, full_data = NULL, env_raster) {
  auc_mod = models_obj@results[best_model_index,]
  FC_best = as.character(auc_mod$fc[1])
  rm_best = auc_mod$rm
  
  
  maxent.args = make.args(RMvalues = rm_best, fc = FC_best)
  
  # re calculating environmental raster
  
  
  full_mod = maxnet(env_raster, as.matrix(full_data[,1:2]), args = maxent.args[[1]])
  return(full_mod)
}
make.args <- function(RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), labels=FALSE) {
  
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
glimpse(t2)

full_data = read_csv("./data/candidate_occurences.csv") %>%
  mutate(true_name = name,
         year = lubridate::year(date)) %>%
  drop_na() %>%
  filter(true_name == "Nymphalis antiopa") %>%
  filter(year < 2000) %>%
  select(-name)


new_full = full_model(models_obj = t2$`Nymphalis antiopa`$models[[1]], 
                      best_model_index = t2$`Nymphalis antiopa`$best_mod[[1]][[2]], 
                      full_data = full_data, 
                      env_raster = t2$`Nymphalis antiopa`$env_data[[1]])


auc_mod = t2$`Nymphalis antiopa`$models[[1]]@results[best_model_index,]
FC_best = as.character(auc_mod$fc[1])
rm_best = auc_mod$rm

maxent.args = make.args(RMvalues = rm_best, fc = FC_best)

maxent_test = maxent(x = t2$`Nymphalis antiopa`$prepped_dfs[[1]][,-c(22:25)],
                     p = t2$`Nymphalis antiopa`$prepped_dfs[[1]]$Species,
                     args = maxent.args[[1]])
t2$`Nymphalis antiopa`$full_mods[[1]] = maxent_test

saveRDS(t2, "./output/marina_antiopa.rds")

rm(t2)

t3 = readRDS("./output/phyleus_coenida.rds")
glimpse(t3)

best_model_index = t3$`Junonia coenia`$best_mod[[1]][[2]]
auc_mod = t3$`Junonia coenia`$models[[1]]@results[best_model_index,]
FC_best = as.character(auc_mod$fc[1])
rm_best = auc_mod$rm

maxent.args = make.args(RMvalues = rm_best, fc = FC_best)

maxent_test = maxent(x = t3$`Junonia coenia`$prepped_dfs[[1]][,-c(22:25)],
                     p = t3$`Junonia coenia`$prepped_dfs[[1]]$Species,
                     args = maxent.args[[1]])
t3$`Junonia coenia`$full_mods[[1]] = maxent_test

best_model_index = t3$`Junonia coenia`$best_mod[[2]][[2]]
auc_mod = t3$`Junonia coenia`$models[[2]]@results[best_model_index,]
FC_best = as.character(auc_mod$fc[1])
rm_best = auc_mod$rm

maxent.args = make.args(RMvalues = rm_best, fc = FC_best)

maxent_test = maxent(x = t3$`Junonia coenia`$prepped_dfs[[2]][,-c(22:25)],
                     p = t3$`Junonia coenia`$prepped_dfs[[2]]$Species,
                     args = maxent.args[[1]])

t3$`Junonia coenia`$full_mods[[2]] = maxent_test

saveRDS(t3, "./output/phyleus_coenida.rds")
rm(t)
rm(t3)


t4 = readRDS("./output/plexippus_clarus.rds")
glimpse(t4)
rm(t4)

t5 = readRDS("./output/vanillae_philenor.rds")
glimpse(t5)
