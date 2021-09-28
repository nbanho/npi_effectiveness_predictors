library(tidyverse)

df <- read_csv("data/data_brauner.csv") %>%
  rename(date = Date, country = `Region Name`, cases = Confirmed) %>%
  dplyr::select(-`Mask Wearing`, -`Travel Screen/Quarantine`, -`Travel Bans`, -`Public Transport Limited`, 
                -`Internal Movement Limited`, -`Public Information Campaigns`, -`Symptomatic Testing`,
                -Active, -Deaths, -`Country Code`) %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(new_cases = dplyr::lead(cases) - cases) %>%
  mutate(new_cases = ifelse(new_cases < 0, floor(0.5 * (dplyr::lead(new_cases) + dplyr::lag(new_cases))), new_cases)) %>% # impute negs
  mutate(cases = dplyr::lag(cumsum(new_cases))) %>%
  slice(1:(n()-1)) %>%
  ungroup() %>%
  dplyr::filter(date <= '2020-05-30') %>%
  dplyr::filter(country != "Andorra") # Filter Andorra because many moderators are NA

country_id <- df %>%
  group_by(country) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(country_id = 1:nrow(.)) %>%
  dplyr::select(country, country_id)

df <- left_join(df, country_id) %>%
  dplyr::select(country_id, country, date, cases, new_cases, everything())

write.csv(df, file = "data/data_brauner_prep.csv", row.names = F)
