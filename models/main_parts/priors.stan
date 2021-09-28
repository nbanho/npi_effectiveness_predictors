alpha_0 ~ student_t(7., 0., 10.);
for (m in 1:M) {
  theta_m[m] ~ asymmetric_laplace(0, 10., 0.5);
}
xi_new_cases ~ normal(0., 1.);
phi_new_infections ~ normal(8.7, 1.0);
new_infections_EpiStart ~ exponential(1. / lambda);
lambda ~ exponential(1.);
mu_p_in ~ normal(p_in_mu_mu0, p_in_mu_sigma0);
sigma_p_in ~ normal(p_in_sigma_mu0, p_in_sigma_sigma0);