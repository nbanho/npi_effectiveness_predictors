run_simulation <- function() {
  
  # Settings
  args <- commandArgs(trailingOnly = TRUE)
  scenario <- as.numeric(args[1]) # 1
  model_spec <- args[2] # theta_j
  start_sim <- as.numeric(args[3]) # 1
  end_sim <- as.numeric(args[4]) # 50
  
  
  prep_file <- "utils/covid.R"
  data_file <- "data/data_brauner_prep.csv"
  scenario_file <- paste0("simulations/", model_spec, "_sim/scenario_", scenario, "/setup.rds")
  sim_file <- "utils/simulation.R"
  stan_file <- paste0("models/", model_spec, ".stan")
  model_save_file <- paste0("simulations/", model_spec, "_sim/scenario_", scenario, "/")
  data_save_file <- paste0("simulations/", model_spec, "_sim/scenario_", scenario, "/")
  
  library(cmdstanr)
  library(tidyverse)
  source(prep_file)
  source(sim_file)
  
  # Load data
  df <- read_csv(data_file)
  
  # NPI names
  npi_names <- colnames(dplyr::select(df, -country_id, -country, -date, -cases, -new_cases))
  
  # Create stan data file
  stanDF <- make_stan_dat(dat = df , npis = npi_names)
  
  # Scenario
  setup <- readRDS(scenario_file)
  
  model_save_file <- paste0(model_save_file, "fit_s")
  data_save_file <- paste0(data_save_file, "data_s")
  
  # Remove old files
  file.remove(list.files(model_save_file, full.names = T))
  file.remove(list.files(data_save_file, full.names = T))
  
  
  for (i in start_sim:end_sim) {
    
    print(sprintf("Simulation %i of %i", i, end_sim))
    
    sim_new_cases <- simulate_sim(seed = i,
                                  X = stanDF$X,
                                  p_g = stanDF$p_g,
                                  max_N2 = stanDF$max_N2, J = stanDF$J, N1 = stanDF$N1, cum_N1 = stanDF$cum_N1, beforeEpiStart = stanDF$beforeEpiStart,
                                  theta_m = setup$true_m$effect, theta_j = setup$true_j$effect, 
                                  priors = setup$priors)
    
    fake_dat <- stanDF
    fake_dat$new_cases <- sim_new_cases$nc
    
    saveRDS(fake_dat, paste0(data_save_file, i, ".rds"))
    
    init_fun <- function() {
      list(mu_p_in = runif(1, 8.0, 14.0))
    }
    main <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = T), include_paths = "models/main_parts/")
    main_fit <- main$sample(
      data =  stan_dat_omit(fake_dat),
      seed = 12345,
      chains = 4,
      iter_warmup = 1000, 
      iter_sampling = 1000,
      parallel_chains = 4,
      threads_per_chain = 2, 
      refresh = 200,
      max_treedepth = 15,
      adapt_delta = 0.9,
      init = init_fun
    )
    
    main_fit$save_object(file = paste0(model_save_file, i, ".rds"))
  }
}

run_simulation()