run_mod <- function() {
  
  # Files
  prep_file <- "utils/covid.R"
  data_file <- "data/data_brauner_prep.csv"
  pred_file <- "data/predictors_brauner.csv"
  pred_missing_file <- "data/predictors_missing_brauner.csv"
  stan_file <- "models/theta_j_mod.stan"
  
  # Save models
  save_dir <- paste0(getwd(), "/fitted-models/brauner/theta_j_mod/single")
  
  # Remove old files
  file.remove(list.files(save_dir, full.names = T))
  
  # Libraries
  library(cmdstanr)
  library(tidyverse)
  library(tidybayes)
  library(posterior)
  source(prep_file)
  
  # Model init
  init_fun <- function() {
    list(mu_p_in = runif(1, 8.0, 14.0))
  }
  main <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = T), include_paths = "models/main_parts/")
  
  # Load data
  df <- read_csv(data_file)
  
  # NPI names
  npi_names <- colnames(dplyr::select(df, -country_id, -country, -date, -cases, -new_cases))
  
  # Load predictors
  predictors <- read_csv(pred_file)
  pred_vars <- colnames(predictors)
  pred_vars <- pred_vars[!grepl("country|pop$", pred_vars)]
  pred_missing <- read_csv(pred_missing_file)
  for (v in pred_vars) {
    # Filter moderator
    single_pred <- predictors %>% 
      dplyr::filter(!is.na(!! sym(v))) 
    single_pred_missing <- pred_missing %>% 
      dplyr::filter(!is.na(!! sym(v))) 
    
    # Subset of countries for which moderator is available
    missing_countries <- unique(df$country[!(df$country %in% single_pred$country)])
    df_sub <- dplyr::filter(df, !(country %in% missing_countries))
    
   
    # Prepare data for stan
    stanDF <- make_stan_dat(dat = df_sub , npis = npi_names)
    
    # Add predictors
    stanDF$predictors <- single_pred %>%
      dplyr::select(v) %>%
      as.matrix()
    stanDF$pop <- single_pred$pop
    stanDF$K <- 1
    
    # Save file
    save_model_name <- paste0(model_name, "_", v)
    
    # Run model
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
    
    # Save model
    single_betas <- main_fit$draws("beta[1]") %>% as_draws_df() %>% as.data.frame()
    saveRDS(single_betas, paste0(save_dir, "/", save_model_name, ".rds"))
  }
}


run_mod()