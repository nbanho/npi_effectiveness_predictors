for (s in 1:max_N2) { 
  p_in[s] = diff_lnorm(max_N2-s, mu_p_in, sigma_p_in);
}