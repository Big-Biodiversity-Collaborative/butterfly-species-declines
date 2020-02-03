#Pulling data for canidate species
#Keaton Wilson
#keatonwilson@me.com
#2019-08-13

# packages
library(tidyverse)
library(spocc)
library(stringr)

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
names_small = names[c(1,2,4,6)]

butt_obs = function(names){
  df = data.frame()
  for (i in 1:length(names)){
    sub = occ(query = names[i], from = c("gbif", "inat"), limit = 10000, 
              has_coords=TRUE, 
              gbifopts=list(continent='north_america'), 
              geometry = c(-140, 20, -90, 60))
    
    sub_df = occ2df(sub)
    
    df = bind_rows(df, sub_df)
  }
  return(df)
}

#Running function above
butterfly_data = butt_obs(names_small)

# Will need to do some duplicate removal
butterfly_data_clean = butterfly_data %>%
  mutate(longitude = as.numeric(longitude), 
         latitude = as.numeric(latitude),
         name = word(name, 1, 2)) %>%
  distinct(longitude, latitude, date, name) 

write_csv(butterfly_data_clean, './data/candidate_occurences.csv')
