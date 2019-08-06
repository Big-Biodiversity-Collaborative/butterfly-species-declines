#Identifying taxonomic problem-species
#Keaton Wilson
#keatonwilson@me.com
#2019-08-06

#packages
library(tidyverse)
library(spocc)

#Loading in decline data
declines = read_csv("./data/declines.csv")
glimpse(declines)

#Function that cycles through every name on the list and pulls out how many gbif
# occurence records there are
found = c()
for (i in 1:nrow(declines)) {
ocs = occ(query = declines[[i,1]], from = "gbif", limit = 1)
found[i] = ocs$gbif$meta$found
}

#binding onto dataframe
declines$gbif_occs = found

glimpse(declines)

#Tagging with another column that makes it easier to filter
declines = declines %>%
  mutate(tax_prob = ifelse(gbif_occs == 0, TRUE, FALSE))

write_csv(declines, path = "./data/declines_tax_probs.csv")
