initially_infected[beforeEpiStart,j] = new_infections_EpiStart[j];
    
// initially infected before EpiStart
for (i in 1:(beforeEpiStart-1)) {
  initially_infected[i,j] = initially_infected[beforeEpiStart,j] / exp(3.28 / 5.5 * (beforeEpiStart-i));
}
    
// number of contagious subjects, expected number of new infections and cases at day t=-EpiStart+1
being_contagious[cum_N1[j]+1] = dot_product(initially_infected[1:beforeEpiStart,j], tail(p_g, beforeEpiStart));
mu_new_infections[cum_N1[j]+1] = being_contagious[cum_N1[j]+1] * exp(alpha_0 + alpha_j[j]); // estimate number of new infections at EpiStart
sigma2_new_infections[cum_N1[j]+1] = mu_new_infections[cum_N1[j]+1] * (1 + mu_new_infections[cum_N1[j]+1] / phi_new_infections);
mu_new_cases[cum_N1[j]+1] = dot_product(append_row(initially_infected[1:beforeEpiStart,j], mu_new_infections[cum_N1[j]+1]),  tail(p_in, beforeEpiStart+1));
    