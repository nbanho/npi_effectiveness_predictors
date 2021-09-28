require(GGally)
require(corrplot)

plot_theta_m <- function(stan_data, npi_names, model, file) {
  npi_var <- as.character(1:stan_data$M)
  names(npi_names) <- npi_var
  
  thetas <- model$post_warmup_draws[,,paste0("theta_m[", 1:stan_data$M, "]")] %>%
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(variable = stringi::stri_extract(variable, regex = "\\d")) 
  
  thetas_pl <- thetas %>%
    mutate(npi = dplyr::recode(variable, !!! npi_names)) %>%
    mutate(reduction = 1 - exp(-value)) %>%
    group_by(variable) %>%
    ggplot(aes(x = reduction, y = reorder(npi, reduction))) +
    stat_pointinterval(point_interval = ci_interval, shape = 21, point_size = 2, fill = "white") +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "red") +
    coord_cartesian(xlim = c(-.25, .75)) +
    scale_x_continuous(labels = function(x) x * 100) +
    labs(x = "Estimated reduction (%)", y = "") +
    theme_bw2() +
    theme(axis.title.y = element_blank())
  
  save_plot(thetas_pl, file, w = 12, h = 6)
  return(thetas_pl)
}


plot_theta_j <- function(stan_data, cid_country, model, file) {
  cs_thetas <- model$post_warmup_draws[,,paste0("theta_j[", 1:stan_data$J, "]")] %>%
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(country_id = stringi::stri_extract(variable, regex = "\\d{1,2}")) %>%
    left_join(cid_country %>% mutate(country_id = as.character(country_id))) %>%
    dplyr::select(country, value, .chain, .iteration, .draw) %>%
    mutate(value = 100*(1-exp(-value))) %>%
    group_by(country) %>%
    mutate(mv = mean(value))
  
  cs_thetas_pl <-  ggplot() +
    stat_pointinterval(data = cs_thetas, 
                       mapping = aes(y = reorder(country, value), x = value, color = mv),
                       shape = 21, point_size = 2, point_fill = "white", .width = c(0.8, 0.95)) +
    scale_color_gradient2(low = "indianred4", mid = "lightgoldenrod2", high = "steelblue4", midpoint = 0,
                          breaks = c(-20, 0, 20), labels = c("Lower", "Average", "Higher"), limits = c(-20, 20)) +
    geom_vline(aes(xintercept = 0), linetype = "dotted", color = "grey") +
    scale_x_continuous(limits = c(-20, 20), breaks = seq(-20,20,5)) +
    labs(x = "Estimated change (%)", y = "", color = "NPI effect") +
    theme_bw2() +
    theme(axis.title.y = element_blank(), legend.position = "bottom")
  
  panel_width = unit(1,"npc") - sum(ggplotGrob(cs_thetas_pl)[["widths"]][-3]) - unit(1,"line")
  cs_thetas_pl <- cs_thetas_pl + guides(color = guide_colorbar(barwidth = panel_width, title.position = "top", title.hjust = 0.5))
  
  save_plot(cs_thetas_pl, file, w = 14, h = 16)
  return(cs_thetas_pl)
}


plot_model_fit <- function(stan_data, cid_country, model, file, ...) {
  print("Get parameters...")
  EN <- list()
  EI <- list()
  for (j in 1:stanDF$J) {
    EN[[j]] <- model$post_warmup_draws[,,paste0("mu_new_cases[", (stan_data$cum_N1[j]+1):(stan_data$cum_N1[j+1]), "]")] %>%
      as_draws_df() %>%
      data.frame() %>%
      reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
      mutate(day = as.numeric(stringi::stri_extract(variable, regex = "\\d{1,4}"))) %>%
      mutate(day = day - stan_data$cum_N1[j]) %>%
      mutate(country = j) %>%
      rename(.value = value) %>%
      dplyr::select(-variable)
    EI[[j]] <- model$post_warmup_draws[,,paste0("mu_new_infections[", (stan_data$cum_N1[j]+1):(stan_data$cum_N1[j+1]), "]")] %>%
      as_draws_df() %>%
      data.frame() %>%
      reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
      mutate(day = as.numeric(stringi::stri_extract(variable, regex = "\\d{1,4}"))) %>%
      mutate(day = day - stan_data$cum_N1[j]) %>%
      mutate(country = j) %>%
      rename(.value = value) %>%
      dplyr::select(-variable)
  }
  EN <- do.call(rbind, EN) %>% group_by(day, country)
  EI <- do.call(rbind, EI) %>% group_by(day, country)
  
  pdf(file, width = 12 / cm(1), height = 10 / cm(1))
  for (j in cid_country$country_id) {
    print(sprintf("Plot for country %i of %i", j, stan_data$J))
    print(plot_compartments(cid_country$country[j], stan_data, EN, EI, ...))
  }
  dev.off()
}


plot_post_pairs <- function(stan_data, model, npi_names, file) {
  npi_var <- as.character(1:stan_data$M)
  names(npi_names) <- npi_var
  
  thetas <- model$post_warmup_draws[,,paste0("theta_m[", 1:stan_data$M, "]")] %>%
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(variable = stringi::stri_extract(variable, regex = "\\d")) %>%
    mutate(variable = as.factor(recode(as.character(variable), !!! npi_names))) %>%
    dplyr::select(variable, value, .draw) %>%
    spread(variable, value) %>%
    dplyr::select(-.draw)
  
  corr_pl <- ggpairs(thetas, 
                     lower = list(continuous = contours),
                     diag = list(continuous = wrap(hist_fun)),
                     upper = list(continuous = wrap(cor_fun, sz = text_size*5/14, stars = FALSE))) +
    theme_bw2()
  
  save_plot(corr_pl, pdf_file = file, w = 20, h = 20)
  return(corr_pl)
}


plot_influential <- function(stan_data, cid_country, model, file) {
  ll_pars_num <- expand.grid(1:stan_data$max_N1, 1:stan_data$J)
  ll_pars <- paste0("log_lik[", ll_pars_num$Var1, ",", ll_pars_num$Var2, "]")
  
  ll_draws <- model$post_warmup_draws[,,ll_pars] %>%
    apply(3, na.omit) %>%
    Filter(f = length) 
  cid <- sapply(stringi::stri_extract_all(names(ll_draws), regex = "\\d{1,2}"), function(x) x[[2]])
  loos <- ll_draws %>% abind::abind(along = 3) %>% loo
  
  loo_dat <- data.frame(country_id =  as.numeric(cid),
                        date = as.Date(unlist(sapply(1:stan_data$J, function(j)
                          as.character(stan_data$date[[j]][(stan_data$EpiStart[j]+1):length(stan_data$date[[j]])])))),
                        k = loos$diagnostics$pareto_k) %>%
    left_join(cid_country)
  
  loo_pl <- loo_dat %>%
    ggplot(aes(x = date, y = k)) +
    geom_point(shape = 3) +
    facet_wrap(~ country, ncol = floor(sqrt(nrow(cid_country))), scales = "free_x") +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 0.7), linetype = "dashed", color = "red") +
    scale_y_continuous(limits = c(min(loo_dat$k)-0.05, ifelse(max(loo_dat$k) > 1, max(loo_dat$k), 1)), breaks = c(0, 0.5, 0.7)) +
    labs(y = "Pareto shape k", x = "Date") +
    theme_bw2() 
  
  save_plot(loo_pl, file, w = 30, h = 30)
  return(loo_pl)
}

plot_vi_vs <- function(model, vi, vs, J, J_lab, vi_lab = "intercept", vs_lab = "slope") {
  vi_sum <- gather_draws_csv(model, vi, 1, J) %>%
    group_by(variable) %>%
    ci_interval() 
  vs_sum <- gather_draws_csv(model, vs, 1, J) %>%
    group_by(variable) %>%
    ci_interval() 
  VD <- data.frame(
    group_name = J_lab,
    intercept = vi_sum$value,
    slope = vs_sum$value
  )
  ggplot(VD, aes(x = intercept, y = slope, label = J_lab)) +
    geom_point() +
    geom_density2d(color = "grey") + 
    geom_text(size = 6 / cm(1), nudge_y = diff(range(VD$slope)) / 25) +
    geom_hline(aes(yintercept = 0), color = "grey", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), color = "grey", linetype = "dashed") +
    labs(x = vi_lab, y = vs_lab) +
    theme_bw2()
}


plot_beta <- function(K, mod_names, par, model, file, xl = c(-.3, .3)) {
  mod_var <- as.character(1:K)
  names(mod_names) <- mod_var
  
  betas <- model$post_warmup_draws[,,paste0(par, "[", 1:K, "]")] %>%
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(variable = stringi::stri_extract(variable, regex = "\\d")) 
  
  betas_pl <- betas %>%
    mutate(mod = dplyr::recode(variable, !!! mod_names)) %>%
    group_by(variable) %>%
    ggplot(aes(x = value, y = reorder(mod, value))) +
    stat_pointinterval(point_interval = ci_interval, shape = 21, point_size = 2, fill = "white") +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "red") +
    coord_cartesian(xlim = xl) +
    labs(x = "Estimated reduction (%)", y = "") +
    theme_bw2() +
    theme(axis.title.y = element_blank())
  
  save_plot(betas_pl, file, w = 12, h = 6)
  return(betas_pl)
}


plot_p_in <- function(stan_data, model, file) {
  
  post_p_in <- model$post_warmup_draws[,,paste0("p_in[", 1:stan_data$max_N2, "]")] %>%
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(x = stringi::stri_extract(variable, regex = "\\d{1,3}")) %>%
    mutate(x = as.numeric(x))
  post_mu_p_in <- model$post_warmup_draws[,,"mu_p_in"] %>%
    as_draws_df() %>%
    data.frame() %>%
    rename(value = mu_p_in)
  post_sigma_p_in <- model$post_warmup_draws[,,"sigma_p_in"] %>%
    as_draws_df() %>%
    data.frame() %>%
    rename(value = sigma_p_in)
  
  post_mu_p_in_pl <- ggplot() + 
    geom_line(data = data.frame(x = seq(0, 22, .001), y = dnorm(seq(0, 22, .001), 10.92, 0.94)), 
              mapping = aes(x = x, y = y), color = "skyblue3") +
    geom_density(data = post_mu_p_in, 
                 mapping = aes(x = value), col = "tomato3") +
    geom_vline(aes(xintercept = mean(post_mu_p_in$value)), col = "tomato3", linetype = "dashed") +
    geom_vline(aes(xintercept = 10.92), col = "skyblue3", linetype = "dashed") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 22), breaks = seq(0, 22, 5)) +
    labs(x = expression(mu), y = "Density") +
    theme_bw2()
  
  post_sigma_p_in_pl <- ggplot() + 
    geom_line(data = data.frame(x = seq(0, 11, .001), y = dnorm(seq(0, 11, .001), 5.41, 0.27)), 
              mapping = aes(x = x, y = y), color = "skyblue3") +
    geom_density(data = post_sigma_p_in, mapping = aes(x = value), col = "tomato3") +
    geom_vline(aes(xintercept = mean(post_sigma_p_in$value)), col = "tomato3", linetype = "dashed") +
    geom_vline(aes(xintercept = 5.41), col = "skyblue3", linetype = "dashed") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 11), breaks = seq(0, 11, 2)) +
    labs(x = expression(sigma), y = "Density") +
    theme_bw2()
  
  post_p_in <- post_p_in %>% 
    mutate(x = max(x) - x) %>%
    dplyr::filter(x <= 30) %>%
    group_by(x) %>%
    summarize(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              mean = mean(value)) %>%
    ungroup() 
  prior_p_in <- data_delay(xT = 30, from0 = T, p_in, m = rnorm(4000, 10.92, 0.94), v = rnorm(4000, 5.41, 0.27)^2) %>%
    group_by(x) %>%
    summarize(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              mean = mean(value)) %>%
    ungroup() 
  
  posterior_p_in_pl <- ggplot(data = rbind(post_p_in, prior_p_in) %>% 
                                mutate(type = rep(c("Posterior", "Prior"), each = 31)), 
                              mapping = aes(x = x, y = mean, ymin = lower, ymax = upper, color = type, fill = type)) +
    geom_line(linetype = "dashed") +
    geom_ribbon(alpha = 0.2, color = "white") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("tomato3", "skyblue3"), name = "") +
    scale_color_manual(values = c("tomato3", "skyblue3"), name = "") +
    labs(y = "Prob. to be reported", x = "Days post infection") +
    theme_bw2()
  
  p_in_grid_pl <- gridExtra::grid.arrange(post_mu_p_in_pl, post_sigma_p_in_pl, posterior_p_in_pl, nrow = 1)
  ggsave(plot = p_in_grid_pl, filename = file, width = 24 / cm(1), height = 4 / cm(1))
  return(p_in_grid_pl)
}
