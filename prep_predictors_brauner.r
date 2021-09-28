library(tidyverse)
library(xlsx)


#### Brauner Data ####
brauner <- read.csv("data/data_brauner.csv") %>% 
  rename(country = `Region.Name`) %>%
  dplyr::select(country) %>%
  unlist() %>%
  unique() %>%
  data.frame() %>%
  set_names("country")


#### GDP per capita ####
gdp_per_capita <- xlsx::read.xlsx("data/brauner_predictors/gdp_per_capita_prep.xlsx", sheetIndex = 1,
                                  stringsAsFactors = F) %>%
  dplyr::select(`Country.Name`, `GDP.per.capita`) %>%
  set_names(c("country", "gdp_per_capita")) %>%
  mutate(gdp_per_capita = as.numeric(gdp_per_capita))


#### Age ####
age <- xlsx::read.xlsx("data/brauner_predictors/age.xlsx", sheetName = "prep", stringsAsFactors = F) %>%
  set_names(c("country", "share_age_young", "share_age_working", "share_age_elderly")) %>%
  mutate_at(vars(matches("share")), as.numeric)


#### Urban ####
urban <- xlsx::read.xlsx("data/brauner_predictors/urban_prep.xls", sheetIndex = 1, stringsAsFactors = F) %>%
  dplyr::select(`Country.Name`, `Share.Urban.2019`) %>%
  set_names(c("country", "share_urban"))


#### Services ####
services <- read_csv("data/brauner_predictors/service.csv") %>%
  rename(country = `ref_area.label`,
         value = obs_value) %>%
  dplyr::filter(`sex.label` == "Sex: Total") %>% 
  dplyr::filter(`classif1.label` %in% c("Economic activity (Detailed): Accommodation and food service activities ~ISIC rev.4 I" ,
                                        "Economic activity (Detailed): Total")) %>%
  dplyr::select(country, time, value, `classif1.label`) %>%
  spread(`classif1.label`, value) %>%
  set_names(c("country", "time", "service", "total")) %>%
  mutate(share_service = service / total) %>%
  dplyr::select(country, share_service)


#### Health expenditure ####
health_exp <- xlsx::read.xlsx("data/brauner_predictors/health_expenditure.xlsx", sheetName = "prep", stringsAsFactors = F) %>%
  mutate(share_health_exp_gdp = as.numeric(share_health_exp_gdp))


#### Informal employment ####
informal <- read_csv2("data/brauner_predictors/informal_employment.csv")


#### Household size ####
hhsize <- read.xlsx("data/brauner_predictors/household_size.xlsx", sheetIndex = 1, stringsAsFactors = F) %>%
  mutate(hh_size = `average.household.size` %>% as.character %>% as.numeric) %>%
  mutate(country = ifelse(country == "Bosnia-Herzegovina", "Bosnia and Herzegovina",
                          ifelse(country == "Czechia", "Czech Republic", country))) %>%
  dplyr::select(country, hh_size)


#### Global health security index ####
ghs <- xlsx::read.xlsx("data/brauner_predictors/jhu_ghs_index.xlsx", sheetIndex = 3, stringsAsFactors = F) %>%
  dplyr::select(`c1`, `ghs_overall`) %>% # Rapid Response Score
  set_names(c("country", "epi_prep_score"))


#### Governance Indicators ####
gov_ind <- read_csv2("data/brauner_predictors/governance_indicators_prep.csv") %>%
  mutate_at(vars(-country), as.numeric)
  

#### Population ####
pop <- xlsx::read.xlsx("data/brauner_predictors/pop_prep.xls", sheetIndex = 1, stringsAsFactors = F) %>%
  dplyr::select(`Country.Name`, `Pop`) %>%
  set_names(c("country", "pop"))


#### Population density ####
pop_density <- xlsx::read.xlsx("data/brauner_predictors/pop_density.xlsx", sheetName = "prep", stringsAsFactors = F) %>%
  mutate(pop_density = log(as.numeric(pop_density)))


#### Join data ####
predictors <- brauner %>%
  left_join(gdp_per_capita %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  left_join(age %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  left_join(informal) %>%
  left_join(hhsize) %>%
  left_join(health_exp %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  left_join(urban %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  left_join(services %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", ifelse(country == "Czechia", "Czech Republic", country)))) %>%
  left_join(gov_ind %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  left_join(ghs) %>%
  left_join(pop %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  left_join(pop_density %>% mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country))) %>%
  dplyr::filter(country != "Andorra") %>% # most predictors are missing for Andorra
  dplyr::select(country, pop, everything()) %>%
  dplyr::select(-reg_qual, -rule_of_law)  # government indicators are strongly correlated, only consider government effectiveness

#### Standardization ####
predictors <- predictors %>%
  mutate_at(vars(-country, -pop), scale) 

#### Missing values ####
# note that the imputation is not used in the main models
# predictors missing is an identifier for missing values that are given parameters to estimate them
predictors_missing <- apply(predictors, c(1,2), function(x) ifelse(is.na(x), 1, 0))
imp <- mice::mice(data = dplyr::select(predictors, share_informal_employment, gdp_per_capita, share_health_exp_gdp) %>% as.matrix, method = "norm.predict", m = 1, maxit = 1)
predictors$share_informal_employment <- complete(imp)[["share_informal_employment"]]

#### Save ####
write.csv(predictors, "data/predictors_brauner.csv", row.names = F)
write.csv(predictors_missing, "data/predictors_missing_brauner.csv", row.names = F)
