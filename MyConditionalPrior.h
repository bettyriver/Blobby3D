#ifndef BLOBBY3D_MYCONDITIONALPRIOR_H_
#define BLOBBY3D_MYCONDITIONALPRIOR_H_

#include "DNest4/code/DNest4.h"

// Hyperparameters setting interim prior for galaxy properties
class MyConditionalPrior:public DNest4::ConditionalPrior {
 private:
  /*
    Hyperparameter priors
  */
  // Limits
  double fluxlim_min, fluxlim_max;
  double flux_std_min, flux_std_max;
  double radiuslim_min, radiuslim_max;
  double wd_min, wd_max;
  double qlim_min;
  double hp_step;

  //  Corresponding widths
  double fluxlim_width;
  double flux_std_width;
  double radiuslim_width;
  double wd_width;
  double qlim_width;

  /*
    Hyperparameters
  */
  // Lognormal for flux
  double flux_mu;
  double flux_std, flux_var;

  // Loguniform with changing limits for W
  double radiusmax;

  // Exponential for wd
  double wd;

  // Triangular with changing lower limit
  double q_min;

  double perturb_hyperparameters(DNest4::RNG& rng);

 public:
   MyConditionalPrior(
    double fluxlim_min, double fluxlim_max,
    double flux_std_min, double flux_std_max,
    double radiuslim_min, double radiuslim_max,
    double wd_min, double wd_max,
    double qlim_min,
    double hp_step
    );

   void from_prior(DNest4::RNG& rng);

   double log_pdf(const std::vector<double>& vec) const;
   void from_uniform(std::vector<double>& vec) const;
   void to_uniform(std::vector<double>& vec) const;

   void print(std::ostream& out) const;
};

#endif  // BLOBBY3D_MYCONDITIONALPRIOR_H_

