# testing build_sdm

source("./script/functions/build_sdm.R")

test = build_sdm(filename = "./data/split_data/colias_eurytheme.rds", 
                 bv_t1, 
                 bv_t2, 
                 year_split = 2000)

saveRDS(test, file = "./output/test_list.rds")

