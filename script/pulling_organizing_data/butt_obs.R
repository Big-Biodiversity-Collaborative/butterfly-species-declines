#Pulling data for canidate species
#Keaton Wilson
#keatonwilson@me.com
#2019-08-13

# packages
library(tidyverse)
library(spocc)
library(stringr)
library(lubridate)

# loading in candidate species data
candidates = read_csv('./data/candidate_species.csv')

# Also loading in the full species list to check whether any of the candidate
# species are problems and we should rename

declines = read_csv("./data/declines.csv")

#joining
can_joins = candidates %>%
  left_join(declines, by = "species") %>%
  dplyr::select(species, gbif_name, mean = mean.x, gbif_occs, quartile, tax_prob)

#fuction to feed in dataframe and build a dataframe with all of the gbif
#records
names = can_joins$gbif_name
names = str_replace(names, " \\(eriphyle\\)", "")
names = str_replace(names, "Everes", "Cupido")

names_small = names[c(3,8)]

butt_obs = function(names){
  df = data.frame()
  for (i in 1:length(names)){
    sub = occ(query = names[i], from = c("gbif", "inat"), limit = 10000, 
              has_coords=TRUE, 
              gbifopts=list(continent='north_america'), 
              geometry = c(-140, 20, -90, 60))

    sub_gbif_df = sub$gbif$data[[1]] %>%
      dplyr::select(longitude, latitude, src_name = species, date = eventDate, key) %>%
      mutate(key = as.numeric(key)) %>%
      mutate(prov = "gbif") %>%
      mutate(name = names[i])
      
    
    sub_inat_df = sub$inat$data[[1]] %>%
      dplyr::select(longitude, latitude, src_name = name, date = observed_on, key = id) %>%
      mutate(longitude = as.numeric(longitude), 
             latitude = as.numeric(latitude), 
             key = as.numeric(key)) %>%
      mutate(prov = "inat") %>%
      mutate(name = names[i])

      dplyr::select(longitude, latitude, name = species, date = eventDate, key) %>%
      mutate(key = as.numeric(key)) %>%
      mutate(prov = "gbif")
    
    sub_inat_df = sub$inat$data[[1]] %>%
      dplyr::select(longitude, latitude, name, date = observed_on, key = id) %>%
      mutate(longitude = as.numeric(longitude), 
             latitude = as.numeric(latitude), 
             key = as.numeric(key)) %>%
      mutate(prov = "inat")

    
    df = bind_rows(sub_gbif_df, sub_inat_df, df)
    
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
  distinct(longitude, latitude, date, name, .keep_all = TRUE) %>%
  drop_na()

# Everything below was tailored cleaning, which is no longer relevant because of the 
# workflow implemented above

# 
# # dealing with more name issues
# butterfly_data_clean %>%
#   mutate(time_frame = ifelse(year(date) < 2000, "T1", "T2")) %>%
#   group_by(name, time_frame) %>%
#   count() %>%
#   print(n = 47)
# 
# # Issue 1 - Pyrgus and Burnsius are the same thing (and burnsius adepta is... something else) 
# butterfly_data_clean$name = str_replace(butterfly_data_clean$name, "Pyrgus communis", "Burnsius communis")
# butterfly_data_clean = butterfly_data_clean %>%
#   filter(name != "Burnsius adepta")
# 
# # Issue 2 - records for eklalyce comyntas should be merged with cupido comyntas
# butterfly_data_clean$name = str_replace(butterfly_data_clean$name, "Elkalyce comyntas", "Cupido comyntas")
# 
# # Issue 3 - Pieris oleracea shouldn't be in the west... primarily a northeastern species
# butterfly_data_clean = butterfly_data_clean %>%
#   filter(name != "Pieris oleracea")
# 
# # Checking again
# summary = butterfly_data_clean %>%
#   mutate(time_frame = ifelse(year(date) < 2000, "T1", "T2")) %>%
#   group_by(name, time_frame) %>%
#   count() %>%
#   print(n = 47)
# 
# # So we have one extra pair of species names compared to the candidate species list. 
# candidates = read_csv("./data/candidate_species.csv")
# 
# summary %>% anti_join(candidates, by = c("name" = "species"))
# 
# # Pieris marginalis is showing up in the records from the pull, but shouldn't be in there when we examine the candidate species
# butterfly_data_clean = butterfly_data_clean %>%
#   filter(name != "Pieris marginalis")

  distinct(longitude, latitude, date, name, .keep_all = TRUE) 


write_csv(butterfly_data_clean, './data/candidate_occurences.csv')
