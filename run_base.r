run_base <- function() {
  
  # Settings
  args <- commandArgs(trailingOnly = TRUE)
  model_name <- args[1] # the name of the model
  sens <- ifelse(args[2] == "true", TRUE, FALSE) # whether to the leave-one-country-out analysis

  # Files
  prep_file <- "utils/covid.R"
  data_file <- "data/data_brauner_prep.csv"
  pred_file <- "data/predictors_brauner.csv"
  if (sens) { pred_file <- "fitted-models/brauner/theta_j_pred_factor/latent_factors.rds" }
  pred_missing_file <- "data/predictors_missing_brauner.csv"
  stan_file <- paste0("models/", model_name, ".stan")

  # Save models
  save_dir <- paste0(getwd(), "/fitted-models/brauner/", model_name)
  if (sens) { save_dir <- paste0(save_dir, "_factor/sens/") }
  
  # Libraries
  library(cmdstanr)
  library(tidyverse)
  source(prep_file)
  
  # Load data
  df <- read_csv(data_file)
  
  # NPI names
  npi_names <- colnames(dplyr::select(df, -country_id, -country, -date, -cases, -new_cases))
  
  # if leave-one-country-out analysis
  if (sens) { countries <- unique(df$country) } 
  else { countries <- "none" }
  for (loo_country in countries) {
    # Filter country
    df_loo <- dplyr::filter(df, country != loo_country) 
    
    # Create stan data file
    stanDF <- make_stan_dat(dat = df_loo , npis = npi_names)
    
    # Add predictors
    if (sens) { predictors <- readRDS(pred_file) } 
    else { predictors <- read_csv(pred_file) }
    pred_missing <- read_csv(pred_missing_file)
    if (model_name == "theta_j_pred") {
      if (sens) {
        stanDF$predictors <- predictors %>%
          dplyr::filter(country != loo_country) %>%
          dplyr::select(country, mean, d) %>% 
          reshape2::dcast(country ~ d, value.var = "mean") %>%
          dplyr::select(-country) %>%
          as.matrix()
        stanDF$predictors_missing <- matrix(0, nrow = nrow(stanDF$predictors), ncol = ncol(stanDF$predictors))
      } else {
        stanDF$predictors <- predictors %>% 
          dplyr::filter(country != loo_country) %>%
          dplyr::select(-country, -pop) %>% 
          as.matrix()
        stanDF$predictors_missing <- pred_missing %>% 
          dplyr::filter(country != loo_country) %>%
          dplyr::select(-country, -pop) %>% 
          as.matrix()
      }
      stanDF$K_missing <- sum(stanDF$predictors_missing)
      stanDF$K <- ncol(stanDF$predictors)
    } else if (model_name == "theta_j_pred_factor") {
      stanDF$predictors <- predictors %>% 
        dplyr::filter(country != loo_country) %>%
        dplyr::select(-country, -pop, -share_age_young) %>% 
        as.matrix()
      stanDF$predictors_missing <- pred_missing %>% 
        dplyr::filter(country != loo_country) %>%
        dplyr::select(-country, -pop, -share_age_young) %>% 
        as.matrix()
      stanDF$K <- ncol(stanDF$predictors)
      stanDF$D <- stanDF$K -1
      stanDF$K_missing <- sum(stanDF$predictors_missing)
    } 
    
    
    # Run model
    init_fun <- function() {
      list(mu_p_in = runif(1, 8.0, 14.0))
    }
    main <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = T), include_paths = "models/main_parts/")
    main_fit <- main$sample(
      data =  stan_dat_omit(stanDF),
      seed = 12345,
      chains = 4,
      iter_warmup = 1000, 
      iter_sampling = 1000,
      parallel_chains = 4,
      threads_per_chain = 2, 
      refresh = 200,
      max_treedepth = 15,
      adapt_delta = 0.9,
      save_warmup = F,
      init = init_fun
    )
    
    # Remove old files
    file.remove(list.files(save_dir, full.names = T, pattern = model_name))
    
    # Save model
    if (sens) { main_fit$save_object(file = paste0(save_dir, loo_country, ".rds")) } 
    else { main_fit$save_output_files(dir = save_dir, basename = model_name) }
  }
}


run_base()