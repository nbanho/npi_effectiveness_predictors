require(posterior)
require(cmdstanr)
require(tidybayes)
require(bayesplot)
require(LaplacesDemon)
require(loo)


# Gather draws for posterior draws from read_cmdstan_csv
gather_draws_csv <- function(model, par, mi, ma) {
  model$post_warmup_draws[,,paste0(par, "[", mi:ma, "]")]%>% 
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(variable = stringi::stri_extract(variable, regex = "\\d+")) 
}

# To support tidybayes gather_draws for cmdstanr
tidy_draws.CmdStanMCMC <- function(model, ...) { return(as_draws_df(model$draws())) }

# LOO for ragged matrices (padded matrices that are originally of unequal length0)
loo.ragged_mat <- function(stan_data, model) {
  ll_pars_num <- expand.grid(1:stan_data$max_N1, 1:stan_data$J)
  ll_pars <- paste0("log_lik[", ll_pars_num$Var1, ",", ll_pars_num$Var2, "]")
  model$post_warmup_draws[,,ll_pars] %>% 
    apply(3, na.omit) %>%
    Filter(f = length) %>%
    abind::abind(along = 3) %>%
    loo
}

# Posterior summary
post_summary <- function(fit, pars) {
  fit$post_warmup_draws[,,pars] %>% 
    summarize_draws(mean, median, function(x) quantile2(x, probs = c(0.025, 0.975)), ess_bulk, Rhat) %>%
    mutate_if(is.numeric, round, 3)
}

# Sampler diagnostics
diagnostics <- function(fit) {
  np <- fit$post_warmup_sampler_diagnostics %>% 
    reshape2::melt() %>% 
    set_names(c("Iteration", "Chain", "Parameter", "Value")) %>%
    dplyr::select(Iteration, Parameter, Value, Chain)
  lp <- fit$post_warmup_draws[,,"lp__"] %>% reshape2::melt() %>% 
    dplyr::select(iteration, value, chain) %>% 
    set_names(c("Iteration", "Value", "Chain"))
  return(list(np = np, lp = lp))
}
