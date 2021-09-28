for (j in 1:J) {
  new_infections[(cum_N1[j]+1):cum_N1[j+1]] ~ normal(mu_new_infections[(cum_N1[j]+1):cum_N1[j+1]], sqrt(sigma2_new_infections[(cum_N1[j]+1):cum_N1[j+1]]));
  new_cases[(cum_N1[j]+1+EpiStart[j]):cum_N1[j+1]] ~ neg_binomial_2(mu_new_cases[(cum_N1[j]+1+EpiStart[j]):cum_N1[j+1]], phi_new_cases);
}