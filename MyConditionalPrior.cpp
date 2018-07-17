#include "MyConditionalPrior.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include <cmath>
#include "boost/math/special_functions/erf.hpp"

using namespace DNest4;

MyConditionalPrior::MyConditionalPrior(double x_min, double x_max,
				       double y_min, double y_max,
				       double r_min, double r_max,
				       double dx, double dy,
				       double dr,
				       double x_pad_dx, double y_pad_dy,
				       double fluxlim_min, double fluxlim_max,
				       double n2fluxlim_min, double n2fluxlim_max)
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,r_min(r_min)
,r_max(r_max)
,dx(dx)
,dy(dy)
,dr(dr)
,x_pad_dx(x_pad_dx)
,y_pad_dy(y_pad_dy)
,fluxlim_min(fluxlim_min)
,fluxlim_max(fluxlim_max)
,n2fluxlim_min(n2fluxlim_min)
,n2fluxlim_max(n2fluxlim_max)
,radiuslim_min(sqrt(dx*dy))
,radiuslim_max(2.0*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy)))
{



}

void MyConditionalPrior::from_prior(RNG& rng)
{

  /*
    Initialise
  */
  
  // Limits
  fluxlim_min = log(fluxlim_min);
  fluxlim_width = log(fluxlim_max) - fluxlim_min;
  
  radiusmin = sqrt(Data::get_instance().get_dx()*Data::get_instance().get_dy());
  radiuslim_width = radiuslim_max/radiuslim_min;

  qlim_min = Data::get_instance().get_qlim_min();

  flux_std_min = Data::get_instance().get_lnfluxsd_min();
  flux_std_width = Data::get_instance().get_lnfluxsd_max()/flux_std_min;
  
  // Initial hyperparameters
  flux_mu = fluxlim_min + fluxlim_width*rng.rand();
  flux_std = exp(log(flux_std_min) + log(flux_std_width)*rng.rand());
  flux_var = pow(flux_std, 2);

  radiusmax = exp(log(radiuslim_min) + log(radiuslim_width)*rng.rand());

  q_min = qlim_min + (1.0 - qlim_min)*rng.rand();

}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
  const double hp_step = Data::get_instance().get_hp_step();

  double logH = 0.0;

  int which;
  which = rng.rand_int(4);

  switch(which)
    {
    case 0:
      flux_mu += hp_step*fluxlim_width*rng.randh();
      flux_mu = mod(flux_mu - fluxlim_min, fluxlim_width);
      flux_mu += fluxlim_min;
      break;
    case 1:
      flux_std = log(flux_std);
      flux_std += hp_step*log(flux_std_width)*rng.randh();
      flux_std = mod(flux_std - log(flux_std_min), log(flux_std_width));
      flux_std += log(flux_std_min);
      flux_std = exp(flux_std);
      flux_var = pow(flux_std, 2);
      break;
    case 2:
      radiusmax = log(radiusmax);
      radiusmax += hp_step*log(radiuslim_width)*rng.randh();
      radiusmax = mod(radiusmax - log(radiuslim_min), log(radiuslim_width));
      radiusmax += log(radiuslim_min);
      radiusmax = exp(radiusmax);
      break;
    case 3:
      q_min += hp_step*(1.0 - qlim_min)*rng.randh();
      q_min = mod(q_min - qlim_min, 1.0 - qlim_min);
      q_min += qlim_min;
      break;
    }

  return logH;
}


double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{

  if(vec[0] < 0.0 ||
     vec[1] < 0.0 || vec[1] > 2.0*M_PI ||
     vec[3] < radiusmin ||
     vec[3] > radiusmax ||
     vec[4] < q_min || vec[4] > 1.0 ||
     vec[5] < 0.0 || vec[5] > M_PI)
    return -1E300;
  
  double logp = 0.0;
  
  // Exponential for radius
  logp += -vec[0];

  // Lognormal for flux
  logp += -log(vec[2]*sqrt(2.0*M_PI*flux_var)) - 0.5*pow(log(vec[2]) - flux_mu, 2)/flux_var;
  
  // LogUniform for width with changing boundaries
  logp += -log(radiusmax - radiusmin);

  // Triangular distribution for q
  logp += 2.0*(vec[4] - q_min)/pow(1.0 - q_min, 2);

  return logp;

}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{ 

  vec[0] = -log(1.0 - vec[0]);
  vec[1] = 2.0*M_PI*vec[1];
  vec[2] = exp(sqrt(2.0*flux_var)*boost::math::erf_inv((2.0*vec[2] - 1.0)*(1.0 - 1E-15)) + flux_mu); 
  vec[3] = radiusmin + (radiusmax - radiusmin)*vec[3];
  vec[4] = (1.0 - q_min)*sqrt(vec[4]) + q_min;
  vec[5] = M_PI*vec[5];
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{

  vec[0] = 1.0 - exp(-vec[0]);
  vec[1] = 0.5*vec[1]/M_PI;
  vec[2] = 0.5 + 0.5*erf((log(vec[2]) - flux_mu)/sqrt(2.0*flux_var));
  vec[3] = (vec[3] - radiusmin)/(radiusmax - radiusmin);
  vec[4] = pow(vec[4] - q_min, 2)/pow(1.0 - q_min, 2);
  vec[5] = vec[5]/M_PI;
}

void MyConditionalPrior::print(std::ostream& out) const
{
  
  out<<exp(flux_mu)<<' '<<flux_std<<' '
     <<radiusmin<<' '<<radiusmax<<' '
     <<q_min<<' ';

}
