#include "BlobConditionalPrior.h"

#include <cmath>

#include "DNest4/code/DNest4.h"
#include "Data.h"

// TODO: DNest4 requires Triangle distribution
// TODO: DNest4 requires LogGaussian distribution
// TODO: Once above are performed just use a for-loop in from/to uniform

/*
  Public
*/
BlobConditionalPrior::BlobConditionalPrior(
  double fluxlim_min, double fluxlim_max,
  double flux_std_min, double flux_std_max,
  double radiuslim_min, double radiuslim_max,
  double wd_min, double wd_max,
  double qlim_min
  ) :fluxlim_min(fluxlim_min)
    ,fluxlim_max(fluxlim_max)
    ,flux_std_min(flux_std_min)
    ,flux_std_max(flux_std_max)
    ,radiuslim_min(radiuslim_min)
    ,radiuslim_max(radiuslim_max)
    ,wd_min(wd_min)
    ,wd_max(wd_max)
    ,qlim_min(qlim_min) {
  // Hyperprior distributions
  hyperprior_fluxmu = DNest4::Uniform(fluxlim_min, fluxlim_max);
  hyperprior_fluxstd = DNest4::LogUniform(flux_std_min, flux_std_max);
  hyperprior_radiusmax = DNest4::LogUniform(radiuslim_min, radiuslim_max);
  hyperprior_wd = DNest4::LogUniform(wd_min, wd_max);
  hyperprior_qmin = DNest4::Uniform(qlim_min, 1.0);
}

void BlobConditionalPrior::from_prior(DNest4::RNG& rng) {
  // Hyperparameters
  wd = hyperprior_wd.generate(rng);
  flux_mu = hyperprior_fluxmu.generate(rng);
  flux_std = hyperprior_fluxstd.generate(rng);
  radiusmax = hyperprior_radiusmax.generate(rng);
  q_min = hyperprior_qmin.generate(rng);

  // Prior Distributions
  prior_rc = DNest4::Exponential(wd);
  prior_theta = DNest4::Uniform(0.0, 2.0*M_PI);
  prior_logflux = DNest4::Gaussian(flux_mu, flux_std);
  prior_width = DNest4::Uniform(radiuslim_min, radiusmax);
  prior_phi = DNest4::Uniform(0.0, M_PI);
}

double BlobConditionalPrior::perturb_hyperparameters(DNest4::RNG& rng) {
  double logH = 0.0;
  int which = rng.rand_int(5);
  switch (which) {
    case 0:
      logH += hyperprior_wd.perturb(wd, rng);
      prior_rc = DNest4::Exponential(wd);
      break;
    case 1:
      logH += hyperprior_fluxmu.perturb(flux_mu, rng);
      prior_logflux = DNest4::Gaussian(flux_mu, flux_std);
      break;
    case 2:
      logH += hyperprior_fluxstd.perturb(flux_std, rng);
      prior_logflux = DNest4::Gaussian(flux_mu, flux_std);
      break;
    case 3:
      logH += hyperprior_radiusmax.perturb(radiusmax, rng);
      prior_width = DNest4::Uniform(radiuslim_min, radiusmax);
      break;
    case 4:
      logH += hyperprior_qmin.perturb(q_min, rng);
      break;
    }

  return logH;
}

double BlobConditionalPrior::log_pdf(const std::vector<double>& vec) const {
  if (
    vec[0] < 0.0 ||
    vec[1] < 0.0 || vec[1] > 2.0*M_PI ||
    vec[3] < radiuslim_min || vec[3] > radiusmax ||
    vec[4] < q_min || vec[4] > 1.0 ||
    vec[5] < 0.0 || vec[5] > M_PI
    )
    return -1E300;

  double logp = 0.0;

  // Exponential for radius
  logp += prior_rc.log_pdf(vec[0]);

  // Lognormal for flux
  logp += prior_logflux.log_pdf(log(vec[2]));

  // Uniform for width with changing boundaries
  logp += prior_width.log_pdf(vec[3]);

  // Triangular distribution for q
  logp += 2.0*(vec[4] - q_min)/pow(1.0 - q_min, 2);

  return logp;
}

void BlobConditionalPrior::from_uniform(std::vector<double>& vec) const {
  vec[0] = prior_rc.cdf_inverse(vec[0]);
  vec[1] = prior_theta.cdf_inverse(vec[1]);
  vec[2] = exp(prior_logflux.cdf_inverse(vec[2]));
  vec[3] = prior_width.cdf_inverse(vec[3]);
  vec[4] = (1.0 - q_min)*sqrt(vec[4]) + q_min;
  vec[5] = prior_phi.cdf_inverse(vec[5]);
}

void BlobConditionalPrior::to_uniform(std::vector<double>& vec) const {
  vec[0] = prior_rc.cdf(vec[0]);
  vec[1] = prior_theta.cdf(vec[1]);
  vec[2] = prior_logflux.cdf(log(vec[2]));
  vec[3] = prior_width.cdf(vec[3]);
  vec[4] = pow(vec[4] - q_min, 2)/pow(1.0 - q_min, 2);
  vec[5] = prior_phi.cdf(vec[5]);
}

void BlobConditionalPrior::print(std::ostream& out) const {
  // Print hyperparameters
  out<<wd<<' '
     <<exp(flux_mu)<<' '<<flux_std<<' '
     <<radiuslim_min<<' '<<radiusmax<<' '
     <<q_min<<' ';
}
