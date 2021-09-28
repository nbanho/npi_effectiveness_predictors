#include funs.stan

data {
#include data.stan
  // Bayesian PCA:
  int<lower=1> K; // number of predictors
  int<lower=1> D; // number of latent dimensions
  matrix[J,K] predictors; // predictors
  matrix[J,K] predictors_missing; // binary matrix indicating whether value of moderator variable is missing
  int<lower=0> K_missing; // number of missing predictors 
}

parameters {
#include pars.stan
  vector[J] alpha_j; // exp(alpha_0 + alpha_j) is the country-specific transmission rate in the absence of NPIs
  vector<lower=0>[2] tau; // between-country variation for alpha_j and theta_j
  vector[J] theta_j; // country-specific effect
  real psi; // the effect of alpha_j on theta_j
  vector[D] beta; // the effect of the factor
  // Bayesian PCA:
  matrix[J, D] Z; // latent factors
	real<lower=0> eps; // additive noise
	cholesky_factor_cov[K, D] W; // the weight matrix
	vector<lower=0>[D] omega; // hierarchical prior for inverse variance of W_d
	vector[K_missing] omicron; // parameters for missing data
}


transformed parameters {
#include trans_pars.stan
  matrix[J,K] predictors_imp; // imputed predictors
  // Bayesian PCA:
  vector<lower=0>[D] t_omega = inv_sqrt(omega);

  // Compute discretized p_IN distribution
#include p_in.stan

  // imputation of missing predictors
  {
    int m_i = 0; 
    for (j in 1:J) {
      for (k in 1:K) {
        if (predictors_missing[j,k] > 0) {
          m_i = m_i + 1;
          predictors_imp[j,k] = omicron[m_i];
        } else {
          predictors_imp[j,k] = predictors[j,k];
        }
      }
    }
  }
  

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
  theta_j ~ normal(psi * alpha_j + Z * beta, tau[2]);
  psi ~ student_t(4., 0., .625);
  beta ~ student_t(4., 0., .125);
  
  // priors for imputation of missing predictors
  omicron ~ normal(0., 1.);
 
  // priors for Bayesian PCA
  eps ~ student_t(4., 0., 1.);
	to_vector(Z) ~ normal(0,1);
	omega ~ gamma(1e-3, 1e-3);
	for (d in 1:D) W[ ,d] ~ normal(0., t_omega[d]);
	
	//likelihood for Bayesian PCA
	to_vector(predictors_imp) ~ normal(to_vector(Z*W'), eps);

  // likelihood
#include likelihood.stan
}

generated quantities {
  #include gen_quants.stan
}
