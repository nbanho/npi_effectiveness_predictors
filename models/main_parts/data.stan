int<lower=0> max_N1; // maximum number of observations (days in seeding + modelling phase) across countries
int<lower=0> max_N2; // maximum number of observations (days before seeding phase + in seeding phase + modeling phase) across countries
int<lower=1> J; // number of countries
int<lower=1> EpiStart[J]; // -EpiStart is the start of the seeding phase (default is 33)
int<lower=1> beforeEpiStart; // 1 + the number of days before the seeding phase (default is 1+7)
int<lower=1> M; // number of NPIs
int<lower=1> N[J]; // number of days in the modeling phase for each country j
int<lower=1> N1[J]; // number of days in the seeding + modeling phase for each country j
int<lower=0> cum_N1[J+1]; // observations are concatenated to one vector, this is the index vector for the time series (including seeding phase)
int<lower=0> new_cases[sum(N1)]; // documented number of new cases
matrix[sum(N1),M] X; // feature matrix with dummies indicating whether NPI m was implemented in country j at time t 
real p_in_mu_mu0; // prior: location hyperparameter mu0 in mu^p_IN ~ Normal(mu0, sigma0)
real p_in_mu_sigma0; // prior: scale hyperparameter sigma0 in mu^p_IN ~ Normal(mu0, sigma0)
real p_in_sigma_mu0; // prior: shape hyperparameter mu0 in sigma^p_IN ~ Normal+(mu0, sigma0)
real p_in_sigma_sigma0; /// prior: inverse scale hyperparameter sigma0 in sigma^p_IN ~ Normal+(mu0, sigma0)
vector<lower=0,upper=1>[max_N2] p_g; // assumed generation time distribution in reversed order (rev(p_G))