#Pulling data for canidate species
#Keaton Wilson
#keatonwilson@me.com
#2019-08-13

# packages
library(tidyverse)
library(spocc)

# loading in candidate species data
candidates = read_csv('./data/candidate_species.csv')

# Also loading in the full species list to check whether any of the candidate
# species are problems and we should rename

declines = read_csv("./data/declines.csv")

#joining
can_joins = candidates %>%
  left_join(declines, by = "species") %>%
  select(species, gbif_name, mean = mean.x, gbif_occs, quartile, tax_prob)

#fuction to feed in dataframe and build a dataframe with all of the gbif
#records
names = can_joins$gbif_name

butt_obs = function(names){
  df = data.frame()
  for (i in 1:length(names)){
    sub = occ(query = names[i], from = "gbif", limit = 100000, has_coords=TRUE, 
              gbifopts=list(continent='north_america'))
    bind_rows(df, sub$gbif$data[[1]])
  }
  return(df)
}

test = butt_obs(names)


