#include funs.stan

data {
#include data.stan
}

parameters {
#include pars.stan
  vector[J] alpha_j; // exp(alpha_0 + alpha_j) is the country-specific transmission rate in the absence of NPIs
  vector[J] theta_j; // country-specific NPI effect
  vector<lower=0>[2] tau; // between-country variation for alpha_j and theta_j
  real psi; // the effect of alpha_j on theta_j
}


transformed parameters {
#include trans_pars.stan
  // Compute discretized p_IN distribution
#include p_in.stan

  for (j in 1:J) {
#include initial.stan
    // number of contagious subjects (up to a normalizing constant), expected number of new infections and cases from day t=-EpiStart+2 onwards
    for (i in (cum_N1[j]+2):cum_N1[j+1]) {
      being_contagious[i] = dot_product(append_row(initially_infected[1:beforeEpiStart,j], new_infections[(cum_N1[j]+1):(i-1)]), tail(p_g, i-cum_N1[j]+beforeEpiStart-1));
      // add theta_j to each theta_m
      mu_new_infections[i] = being_contagious[i] * exp(alpha_0 + alpha_j[j] - X[i] * (theta_m + theta_j[j])); 
      sigma2_new_infections[i] = mu_new_infections[i] * (1 + mu_new_infections[i] / phi_new_infections);
      mu_new_cases[i] = dot_product(append_row(initially_infected[1:beforeEpiStart,j], new_infections[(cum_N1[j]+1):(i)]), tail(p_in, i-cum_N1[j]+beforeEpiStart));
    }
  }
}


model {
  // priors
#include priors.stan
  tau[1] ~ student_t(4., 0., 1.);
  tau[2] ~ student_t(4., 0., 0.1);
  alpha_j ~ normal(0., tau[1]);
  psi ~ student_t(4., 0., .625);
  theta_j ~ normal(psi * alpha_j, tau[2]);

  // likelihood
#include likelihood.stan
}

generated quantities {
  #include gen_quants.stan
}