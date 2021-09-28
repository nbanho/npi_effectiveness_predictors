functions {
  // discretization of the cumulative lognormal distribution function for p_in
  real diff_lnorm(real x, real mu, real sigma) {
    if (x == 0) {
      return lognormal_cdf( 0.5, log( mu ^ 2 / sqrt( sigma ^ 2 + mu ^ 2 ) ), sqrt( log( sigma ^ 2 / mu ^ 2 + 1 ) ) );
    } else {
      return lognormal_cdf( x + 0.5, log( mu ^ 2 / sqrt( sigma ^ 2 + mu ^ 2 ) ), sqrt( log( sigma ^ 2 / mu ^ 2 + 1 ) ) ) - 
      lognormal_cdf( x - 0.5, log( mu ^ 2 / sqrt( sigma ^ 2 + mu ^ 2 ) ), sqrt( log( sigma ^ 2 / mu ^ 2 + 1 ) ) );
    }
  }
  // asymmetric laplace prior from Brauner et al for theta_m
  real asymmetric_laplace_lpdf(real y, real mu, real scale, real symmetry) {
    return log(scale) - log(symmetry + 1 / symmetry) + ((y < mu) ? (scale / symmetry * (y - mu)) : (-scale * symmetry * (y - mu))); 
  }
}