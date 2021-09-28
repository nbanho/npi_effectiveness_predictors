simulate_sim <- function(seed,
                         X,
                         p_g,
                         max_N2, N1, cum_N1, J, beforeEpiStart,
                         theta_m = NULL, theta_j = NULL, theta_mj = NULL, 
                         priors = list()) {
  
  set.seed(seed)
  
  # Initialization
  p_in <- numeric(max_N2)
  initially_infected <- matrix(NA, nrow = beforeEpiStart, ncol = J)
  long_n <- sum(N1)
  being_contagious <- numeric(long_n)
  new_infections <- numeric(long_n)
  new_cases <- numeric(long_n)
  mu_new_infections <- numeric(long_n)
  mu_new_cases <- numeric(long_n)
  sigma2_new_infections <- numeric(long_n)
  
  # Priors / Assumptions based on base.stan 
  alpha_0 <- ifelse(priors$alpha_0 %>% is.null, 0.97, priors$alpha_0)
  tau <- ifelse(priors$tau %>% is.null, 0.20, priors$tau)
  if (is.null(priors$alpha_j)) {
    alpha_j <- numeric(J)
    for (j in 1:J) {
      alpha_j[j] <- rnorm(1, 0, tau)
    }
  } else {
    alpha_j <- priors$alpha_j
  }
  phi_new_cases <- ifelse(priors$phi_new_cases %>% is.null, 3.5, priors$phi_new_cases)
  phi_new_infections <- ifelse(priors$phi_new_infections %>% is.null, 8.5, priors$phi_new_infections)
  mu_p_in <- ifelse(priors$mu_p_in %>% is.null, 14.76, priors$mu_p_in)
  sigma_p_in <- ifelse(priors$sigma_p_in %>% is.null, 4.76, priors$sigma_p_in)
  lambda <- ifelse(priors$lambda %>% is.null, 2.3, priors$lambda)
  
  half_distr <- function(distr, ...) {
    y <- -1
    while (y < 0) {
      y <- distr(1, ...)
    }
    return(y)
  }
  
  
  # Delay distribution p_IN
  for (s in 1:max_N2) {
    if (max_N2 - s == 0) {
      p_in[s] <- plnorm(0.5, log( mu_p_in ^ 2 / sqrt( sigma_p_in ^ 2 + mu_p_in ^ 2 ) ), sqrt( log( sigma_p_in ^ 2 / mu_p_in ^ 2 + 1 ) ) )
    } else {
      p_in[s] <- plnorm(max_N2 - s + 0.5, log( mu_p_in ^ 2 / sqrt( sigma_p_in ^ 2 + mu_p_in ^ 2 ) ), sqrt( log( sigma_p_in ^ 2 / mu_p_in ^ 2 + 1 ) )) - 
        plnorm(max_N2 - s - 0.5, log( mu_p_in ^ 2 / sqrt( sigma_p_in ^ 2 + mu_p_in ^ 2 ) ), sqrt( log( sigma_p_in ^ 2 / mu_p_in ^ 2 + 1 ) ))
    }
  }
  
  # Model
  for (j in 1:J) {
    # initially infected
    if (is.null(priors$new_infections_EpiStart)) {
      initially_infected[beforeEpiStart,j] <- rexp(1, 1 / lambda)
    } else {
      initially_infected[beforeEpiStart,j] <- priors$new_infections_EpiStart[j]
    }
    for (i in 1:(beforeEpiStart-1)) {
      initially_infected[i,j] <- initially_infected[beforeEpiStart,j] / exp(3.28 / 5.5 * (beforeEpiStart - i))
    }
    
    # initially contagious
    being_contagious[cum_N1[j]+1] <- initially_infected[1:beforeEpiStart,j] %*% tail(p_g, beforeEpiStart)
    mu_new_infections[cum_N1[j]+1] <- being_contagious[cum_N1[j]+1] * exp(alpha_0 + alpha_j[j])
    sigma2_new_infections[cum_N1[j]+1] <- mu_new_infections[cum_N1[j]+1] * (1 + mu_new_infections[cum_N1[j]+1] / phi_new_infections)
    new_infections[cum_N1[j]+1] <- half_distr(rnorm, mu_new_infections[cum_N1[j]+1], sqrt(sigma2_new_infections[cum_N1[j]+1]))
    mu_new_cases[cum_N1[j]+1] <- c(initially_infected[1:beforeEpiStart,j], new_infections[cum_N1[j]+1]) %*% tail(p_in, beforeEpiStart+1)
    new_cases[cum_N1[j]+1] <- rnbinom(1, mu = mu_new_cases[cum_N1[j]+1], size = phi_new_cases)
    
    # new cases over time
    for (i in (cum_N1[j]+2):cum_N1[j+1]) {
      being_contagious[i] <- c(initially_infected[1:beforeEpiStart,j], new_infections[(cum_N1[j]+1):(i-1)]) %*% tail(p_g, i-cum_N1[j]+beforeEpiStart-1)
      if (!is.null(theta_mj)) {
        mu_new_infections[i] <- being_contagious[i] * exp(alpha_0 + alpha_j[j] - X[i, ] %*% (theta_m + theta_j[j] + theta_mj[ ,j]))
      } else if (!is.null(theta_j)) {
        mu_new_infections[i] <- being_contagious[i] * exp(alpha_0 + alpha_j[j] - X[i, ] %*% (theta_m + theta_j[j]))
      } else if (!is.null(theta_m)) {
        mu_new_infections[i] <- being_contagious[i] * exp(alpha_0 + alpha_j[j] - X[i, ] %*% theta_m)
      } else {
        mu_new_infections[i] <- being_contagious[i] * exp(alpha_0 + alpha_j[j])
      }
      sigma2_new_infections[i] <- mu_new_infections[i] * (1 + mu_new_infections[i] / phi_new_infections)
      new_infections[i] <- half_distr(rnorm, mu_new_infections[i], sqrt(sigma2_new_infections[i]))
      mu_new_cases[i] <- c(initially_infected[1:beforeEpiStart,j], new_infections[(cum_N1[j]+1):(i)]) %*% tail(p_in, i-cum_N1[j]+beforeEpiStart)
      new_cases[i] <- rnbinom(1, mu = mu_new_cases[i], size = phi_new_cases)
    }
  }
  
  return(list(bc = being_contagious, mu_ni = mu_new_infections, ni = new_infections, nc = new_cases))
}


theta_boxplot <- function(theta_draws, true_theta, file, levs = true_theta$npi, ordered_levs = true_theta$npi, par = "theta_m", plt_w = 8, plt_h = 4) {
  theta_means <- sapply(theta_draws, function(X) X %>% mean_qi %>% dplyr::select(par) %>% unlist) %>%
    t %>%
    data.frame() %>%
    set_names(levs) %>%
    gather() %>%
    set_names(c("npi", "effect")) 
  
  if (ordered_levs == "by_post_var") {
    # theta_var <- sapply(theta_draws, function(X) X %>% group_by(country) %>% summarize(v = var(!!sym(par))) %>% dplyr::select(v) %>% unlist) %>%
    #   t %>%
    #   data.frame() %>%
    #   set_names(levs) %>%
    #   gather() %>%
    #   set_names(c("npi", "v"))  %>%
    #   group_by(npi) %>%
    #   summarize(mv = mean(v)) %>%
    #   ungroup()
    var_means <- theta_means %>%
      group_by(npi) %>%
      summarize(sd_effect = boxplot(effect)$stats[5] - boxplot(effect)$stats[1]) %>%
      ungroup()
    theta_means <- mutate(theta_means, npi = factor(npi, levels = arrange(var_means, sd_effect)$npi))
  } else {
    theta_means <- mutate(theta_means, npi = factor(npi, levels = ordered_levs))
  }
  
  plt <- ggplot() +
    geom_boxplot(data = theta_means, mapping = aes(x = npi, y = effect)) +
    geom_point(data = true_theta, mapping = aes(x = npi, y = effect), color = "red", size = 2, shape = 4) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey") +
    labs(y = "Estimated effect", x = "") +
    coord_flip() +
    theme_bw2() 
  
  ggsave(file, width = plt_w / cm(1), height = plt_h / cm(1))
  return(plt)
}


theta_barchart <- function(theta_draws, true_theta, file, levs = true_theta$npi, ordered_levs = true_theta$npi, thresh = length(stan_fits), plt_w = 8, plt_h = 4) {
  theta_in <- sapply(theta_draws, function(X, tt) {
    X <- X %>% mean_qi %>% mutate(true_value = tt, is_in = ifelse(true_value >= .lower & true_value <= .upper, 1, 0))
    return(X$is_in)
  }, tt = true_theta$effect) %>%
    t %>%
    data.frame() %>%
    set_names(levs) %>%
    summarize_all(function(x) round(sum(x) / n() * 100)) %>%
    gather() %>%
    set_names(c("npi", "prob")) %>%
    mutate(npi = factor(npi, levels = ordered_levs))
  
  plt <- ggplot() +
    geom_bar(data = theta_in, mapping = aes(x = npi, y = prob), stat = "identity") +
    geom_text(data = theta_in, mapping = aes(x = npi, y = prob, label = prob), hjust = -0.2, color = "black",
              size = 8 / cm(1)) +
    labs(y = "Prop. of CrIs containing true effect (%)", x = "") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 120), breaks = seq(0, 100, 20)) +
    coord_flip() +
    theme_bw2() 
  
  ggsave(file, width = plt_w / cm(1), height = plt_h / cm(1))
  return(plt)
}


theta_heatmap <- function(theta_draws, true_theta, file, npi_labels = NULL) {
  theta_cor <- lapply(thetas_all, function(X) {
    X %>% spread(m, theta_m) %>% dplyr::select(-.chain, -.iteration, -.draw) %>% cor
  })
  
  theta_cor_mean <- Reduce("+", theta_cor) / length(theta_cor)
  colnames(theta_cor_mean) <- rownames(theta_cor_mean) <- true_theta$npi
  
  theta_cor_mean <- theta_cor_mean %>%
    data.frame() %>%
    add_rownames() %>%
    reshape2::melt(id.vars = "rowname") %>%
    mutate(variable = gsub(".", " ", variable, fixed = T)) %>%
    mutate(variable = ifelse(variable == "Stay at home order", "Stay-at-home order", variable)) %>%
    mutate(variable = ifelse(variable == "Work from home order", "Work-from-home order", variable)) 
  
  if (!is.null(npi_labels)) {
    theta_cor_mean <- theta_cor_mean %>%
      mutate(variable = dplyr::recode(variable, !!! npi_labels), 
             rowname = dplyr::recode(rowname, !!! npi_labels)) %>%
      mutate(rowname = factor(rowname, levels = npi_labels),
             variable = factor(variable, levels = npi_labels))
  } else {
    theta_cor_mean <- theta_cor_mean %>%
      mutate(rowname = factor(rowname, levels = true_theta$npi),
             variable = factor(variable, levels = true_theta$npi))
  }
  
  plt <- ggplot(data = theta_cor_mean) +
    geom_tile(aes(x = rowname, y = variable, fill = value)) +
    geom_text(aes(rowname, variable, label = round(value, 2)), color = "black", size = 4) +
    labs(y = "", x = "") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="") +
    theme_bw2() 
  
  ggsave(file, width = 18 / cm(1), height = 15 / cm(1))
  return(plt)
}
