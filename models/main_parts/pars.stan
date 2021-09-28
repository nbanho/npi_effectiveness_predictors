real alpha_0; // exp(alpha_0) is the cross-country average transmission rate in the absence of NPIs 
vector[M] theta_m; // parameters for the cross-country average NPI effects
real<lower=0> xi_new_cases; // over-dispersion on the parameter 1 / sqrt(phi^N) for the number of new cases
real<lower=0> phi_new_infections; // over-dispersion on the parameter 1 / sqrt(phi^I) for the number of new infections
real<lower=0> new_infections_EpiStart[J]; // number of new infections in country j at the start of the seeding phase, i.e., at day t = -EpiStart
real<lower=0> lambda; // expected number of new infections in country j at the start of the seeding phase, i.e., at day t = -EpiStart
vector<lower=0>[sum(N1)] new_infections; // number of new infections (latent / unobserved)
real<lower=0> mu_p_in; // log mean in p_IN ~ Lognormal(mu, sigma)
real<lower=0> sigma_p_in; // log standard deviation in p_IN ~ Lognormal(mu, sigma)