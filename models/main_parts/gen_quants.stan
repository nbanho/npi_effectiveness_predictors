// log-likelihood
matrix[max_N1,J] log_lik;
for (j in 1:J) {
  for (i in 1:N[j]) {
    log_lik[i,j] = neg_binomial_2_log_lpmf(new_cases[cum_N1[j]+i+EpiStart[j]] | log(mu_new_cases[cum_N1[j]+i+EpiStart[j]]), phi_new_cases);
  }
}