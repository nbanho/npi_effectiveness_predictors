real phi_new_cases = inv_square(xi_new_cases); // over-dispersion parameter for new cases N
vector[sum(N1)] mu_new_cases; // expected number of new cases
vector[sum(N1)] being_contagious; // number of contagious subjects (up to a normalizing constant)
vector[sum(N1)] mu_new_infections; // expected number of new infections
vector[sum(N1)] sigma2_new_infections; // variance of the expected number of new infections
vector[max_N2] p_in; // probability distribution for the time from infection to reporting of a new case p_IN
matrix[beforeEpiStart,J] initially_infected; // number of initially infected before the start of the seeding phase (including new_infections_EpiStart)