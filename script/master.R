# Master function that runs build_sdm() for a single species

source("./script/functions/build_sdm.R")

spec = commandArgs()[[6]]
filename = paste0("./data/split_data/", spec)

# Running a specific species
test_obj = build_sdm(filename = filename, 
                     bv_t1, 
                     bv_t2, 
                     year_split = 2000,
		             cores = 14)

# writing list out
file_out = paste0("./output/", spec)
saveRDS(test_obj, file_out)
