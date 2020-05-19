#Identifying candidate species
#Keaton Wilson
#keatonwilson@me.com
#2019-08-06

#packages
library(tidyverse)
library(spocc)

#loading in data
declines = read_csv("./data/declines_tax_probs.csv")

glimpse(declines)

#Splitting species into decline quartiles and then taking the top 5 species 
# by occurences from each quartile
candidate_species = declines %>%
  mutate(quartile = ntile(mean, 4)) %>%
  group_by(quartile) %>%
  top_n(n = 5, wt = gbif_occs)
         
#Writing to csv
write_csv(candidate_species, "./data/candidate_species.csv")
