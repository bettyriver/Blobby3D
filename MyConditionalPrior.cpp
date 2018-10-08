#include "MyConditionalPrior.h"

#include <cmath>
#include "boost/math/special_functions/erf.hpp"

#include "DNest4/code/DNest4.h"
#include "Data.h"

/*
  Public
*/
MyConditionalPrior::MyConditionalPrior(
  double fluxlim_min, double fluxlim_max,
  double flux_std_min, double flux_std_max,
  double radiuslim_min, double radiuslim_max,
  double wd_min, double wd_max,
  double qlim_min,
  double hp_step
  ) :fluxlim_min(fluxlim_min)
    ,fluxlim_max(fluxlim_max)
    ,flux_std_min(flux_std_min)
    ,flux_std_max(flux_std_max)
    ,radiuslim_min(radiuslim_min)
    ,radiuslim_max(radiuslim_max)
    ,wd_min(wd_min)
    ,wd_max(wd_max)
    ,qlim_min(qlim_min)
    ,hp_step(hp_step)
    ,fluxlim_width(fluxlim_max - fluxlim_min)
    ,flux_std_width(flux_std_max/flux_std_min)
    ,radiuslim_width(radiuslim_max/radiuslim_min)
    ,wd_width(wd_max/wd_min)
    ,qlim_width(1.0 - qlim_min) {}

void MyConditionalPrior::from_prior(DNest4::RNG& rng) {
  wd = exp(log(wd_min) + rng.rand()*log(wd_width));

  flux_mu = fluxlim_min + fluxlim_width*rng.rand();

  flux_std = exp(log(flux_std_min) + log(flux_std_width)*rng.rand());
  flux_var = pow(flux_std, 2);

  radiusmax = exp(log(radiuslim_min) + log(radiuslim_width)*rng.rand());

  q_min = qlim_min + qlim_width*rng.rand();
}

double MyConditionalPrior::perturb_hyperparameters(DNest4::RNG& rng) {
  double logH = 0.0;
  int which = rng.rand_int(5);

  switch (which) {
    case 0:
      wd = log(wd);
      wd += hp_step*log(wd_width)*rng.randh();
      wd = DNest4::mod(wd - log(wd_min), log(wd_width));
      wd += log(wd_min);
      wd = exp(wd);
      break;
    case 1:
      flux_mu += hp_step*fluxlim_width*rng.randh();
      flux_mu = DNest4::mod(flux_mu - fluxlim_min, fluxlim_width);
      flux_mu += fluxlim_min;
      break;
    case 2:
      flux_std = log(flux_std);
      flux_std += hp_step*log(flux_std_width)*rng.randh();
      flux_std = DNest4::mod(flux_std - log(flux_std_min), log(flux_std_width));
      flux_std += log(flux_std_min);
      flux_std = exp(flux_std);
      flux_var = pow(flux_std, 2);
      break;
    case 3:
      radiusmax = log(radiusmax);
      radiusmax += hp_step*log(radiuslim_width)*rng.randh();
      radiusmax = DNest4::mod(radiusmax - log(radiuslim_min), log(radiuslim_width));
      radiusmax += log(radiuslim_min);
      radiusmax = exp(radiusmax);
      break;
    case 4:
      q_min += hp_step*qlim_width*rng.randh();
      q_min = DNest4::mod(q_min - qlim_min, qlim_width);
      q_min += qlim_min;
      break;
    }

  return logH;
}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const {
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
  logp += -log(wd) - vec[0]/wd;

  // Lognormal for flux
  logp += -log(vec[2]*sqrt(2.0*M_PI*flux_var));
  logp += -0.5*pow(log(vec[2]) - flux_mu, 2)/flux_var;

  // Uniform for width with changing boundaries
  logp += -log(radiusmax - radiuslim_min);

  // Triangular distribution for q
  logp += 2.0*(vec[4] - q_min)/pow(1.0 - q_min, 2);

  return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const {
  vec[0] = -wd*log(1.0 - vec[0]);
  vec[1] = 2.0*M_PI*vec[1];
  vec[2] = exp(sqrt(2.0*flux_var)*boost::math::erf_inv((2.0*vec[2] - 1.0)*(1.0 - 1E-15)) + flux_mu);
  vec[3] = radiuslim_min + (radiusmax - radiuslim_min)*vec[3];
  vec[4] = (1.0 - q_min)*sqrt(vec[4]) + q_min;
  vec[5] = M_PI*vec[5];
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const {
  vec[0] = 1.0 - exp(-vec[0]/wd);
  vec[1] = 0.5*vec[1]/M_PI;
  vec[2] = 0.5 + 0.5*erf((log(vec[2]) - flux_mu)/sqrt(2.0*flux_var));
  vec[3] = (vec[3] - radiuslim_min)/(radiusmax - radiuslim_min);
  vec[4] = pow(vec[4] - q_min, 2)/pow(1.0 - q_min, 2);
  vec[5] = vec[5]/M_PI;
}

void MyConditionalPrior::print(std::ostream& out) const {
  out<<wd<<' '
     <<exp(flux_mu)<<' '<<flux_std<<' '
     <<radiuslim_min<<' '<<radiusmax<<' '
     <<q_min<<' ';
}
