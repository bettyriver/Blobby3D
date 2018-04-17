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
,radiuslim_max(3.0*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy)))
,dispersion_min(Data::get_instance().get_vdispmu_min())
,dispersion_max(Data::get_instance().get_vdispmu_max())
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

  n2fluxlim_min = log(n2fluxlim_min);
  n2fluxlim_width = log(n2fluxlim_max) - n2fluxlim_min;
  
  radiuslim_width = radiuslim_max/radiuslim_min;

  qlim_min = Data::get_instance().get_qlim_min();

  dispersion_mu_min = log(dispersion_min);
  dispersion_mu_width = log(dispersion_max) - dispersion_mu_min;

  flux_std_min = Data::get_instance().get_lnfluxsd_min();
  flux_std_width = Data::get_instance().get_lnfluxsd_max()/flux_std_min;
  
  n2flux_std_min = Data::get_instance().get_lnfluxsd_min();
  n2flux_std_width = Data::get_instance().get_lnfluxsd_max()/n2flux_std_min;
  
  dispersion_std_min = Data::get_instance().get_lnvdispsd_min();
  dispersion_std_width = Data::get_instance().get_lnvdispsd_max()/dispersion_std_min;

  wd_min = 0.1*sqrt(dx*dy);
  wd_width = 3.0*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy))/wd_min;
  
  // Initial hyperparameters
  flux_mu = fluxlim_min + fluxlim_width*rng.rand();
  flux_std = exp(log(flux_std_min) + log(flux_std_width)*rng.rand());
  flux_var = pow(flux_std, 2);

  if(Data::get_instance().get_model_n2lines())
    {
      n2flux_mu = n2fluxlim_min + n2fluxlim_width*rng.rand();
      n2flux_std = exp(log(n2flux_std_min) + log(n2flux_std_width)*rng.rand());
      n2flux_var = pow(n2flux_std, 2);
    }
  else
    {
      n2flux_mu = 0.0;
      n2flux_std = 0.0;
      n2flux_var = 0.0;
    }

  radiusmax = exp(log(radiuslim_min) + log(radiuslim_width)*rng.rand());
  radiusmin_ratio = rng.rand();

  q_min = qlim_min + (1.0 - qlim_min)*rng.rand();

  dispersion_mu = dispersion_mu_min + dispersion_mu_width*rng.rand();
  dispersion_std = exp(log(dispersion_std_min) + log(dispersion_std_width)*rng.rand());
  dispersion_var = pow(dispersion_std, 2);

  wd = exp(log(wd_min) + log(wd_width)*rng.rand());

}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
  double logH = 0.0;

  int which;
  if(Data::get_instance().get_model_n2lines() == 1)
    which = rng.rand_int(10);
  else
    which = rng.rand_int(8);

  switch(which)
    {
    case 0:
      flux_mu += fluxlim_width*rng.randh();
      flux_mu = mod(flux_mu - fluxlim_min, fluxlim_width);
      flux_mu += fluxlim_min;
      break;
    case 1:
      flux_std = log(flux_std);
      flux_std += log(flux_std_width)*rng.randh();
      flux_std = mod(flux_std - log(flux_std_min), log(flux_std_width));
      flux_std += log(flux_std_min);
      flux_std = exp(flux_std);
      flux_var = pow(flux_std, 2);
      break;
    case 2:
      dispersion_mu += dispersion_mu_width*rng.randh();
      dispersion_mu = mod(dispersion_mu - dispersion_mu_min, dispersion_mu_width);
      dispersion_mu += dispersion_mu_min;
      break;
    case 3:
      dispersion_std = log(dispersion_std);
      dispersion_std += log(dispersion_std_width)*rng.randh();
      dispersion_std = mod(dispersion_std - log(dispersion_std_min), log(dispersion_std_width));
      dispersion_std += log(dispersion_std_min);
      dispersion_std = exp(dispersion_std);
      dispersion_var = pow(dispersion_std, 2);
      break;
    case 4:
      radiusmax = log(radiusmax);
      radiusmax += log(radiuslim_width)*rng.randh();
      radiusmax = mod(radiusmax - log(radiuslim_min), log(radiuslim_width));
      radiusmax += log(radiuslim_min);
      radiusmax = exp(radiusmax);
      break;
    case 5:
      radiusmin_ratio += rng.randh();
      radiusmin_ratio = mod(radiusmin_ratio, 1.0);
      break;
    case 6:
      wd = log(wd);
      wd += log(wd_width)*rng.randh();
      wd = mod(wd - log(wd_min), log(wd_width));
      wd += log(wd_min);
      wd = exp(wd);
      break;
    case 7:
      q_min += (1.0 - qlim_min)*rng.randh();
      q_min = mod(q_min - qlim_min, 1.0 - qlim_min);
      q_min += qlim_min;
      break;
    case 8:
      n2flux_mu += n2fluxlim_width*rng.randh();
      n2flux_mu = mod(n2flux_mu - n2fluxlim_min, n2fluxlim_width);
      n2flux_mu += n2fluxlim_min;
      break;
    case 9:
      n2flux_std = log(n2flux_std);
      n2flux_std += log(n2flux_std_width)*rng.randh();
      n2flux_std = mod(n2flux_std - log(n2flux_std_min), log(n2flux_std_width));
      n2flux_std += log(n2flux_std_min);
      n2flux_std = exp(n2flux_std);
      n2flux_var = pow(n2flux_std, 2);
      break;
    }

  return logH;
}


double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
  
  double radiusmin = exp(log(radiuslim_min) + log(radiusmax/radiuslim_min)*radiusmin_ratio);

  if(vec[0] < 0.0 ||
     vec[1] < 0.0 || vec[1] > 2.0*M_PI ||
     vec[3] < radiusmin ||
     vec[3] > radiusmax ||
     vec[4] < q_min || vec[4] > 1.0 ||
     vec[5] < 0.0 || vec[5] > M_PI)
    return -1E300;
  
  double logp = 0.0;
  
  // Lognormal for flux and dispersion
  logp += -vec[0]/wd - log(wd);
  logp += -log(vec[2]*sqrt(2.0*M_PI*flux_var)) - 0.5*pow(log(vec[2]) - flux_mu, 2)/flux_var;
  logp += -log(vec[6]*sqrt(2.0*M_PI*dispersion_var)) - 0.5*pow(log(vec[6]) - dispersion_mu, 2)/dispersion_var; 
  
  // LogUniform for width with changing boundaries
  logp += -log(vec[3]*log(radiusmax/radiusmin));

  // Triangular distribution for q
  logp += 2.0*(vec[4] - q_min)/pow(1.0 - q_min, 2);

  if(Data::get_instance().get_model_n2lines() == 1)
    logp += -log(vec[7]*sqrt(2.0*M_PI*n2flux_var)) - 0.5*pow(log(vec[7]) - n2flux_mu, 2)/n2flux_var;

  return logp;

}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{ 

  double radiusmin = exp(log(radiuslim_min) + log(radiusmax/radiuslim_min)*radiusmin_ratio);

  vec[0] = -wd*log(1.0 - vec[0]);
  vec[1] = 2.0*M_PI*vec[1];
  vec[2] = exp(sqrt(2.0*flux_var)*boost::math::erf_inv((2.0*vec[2] - 1.0)*(1.0 - 1E-15)) + flux_mu); 
  vec[3] = exp(log(radiusmin) + log(radiusmax/radiusmin)*vec[3]);

  vec[4] = (1.0 - q_min)*sqrt(vec[4]) + q_min;

  vec[5] = M_PI*vec[5];
  vec[6] = exp(sqrt(2.0*dispersion_var)*boost::math::erf_inv((2.0*vec[6] - 1.0)*(1.0 - 1E-15)) + dispersion_mu);

  if(Data::get_instance().get_model_n2lines() == 1)
    vec[7] = exp(sqrt(2.0*n2flux_var)*boost::math::erf_inv((2.0*vec[7] - 1.0)*(1.0 - 1E-15)) + n2flux_mu);

}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
  double radiusmin = exp(log(radiuslim_min) + log(radiusmax/radiuslim_min)*radiusmin_ratio);

  vec[0] = 1.0 - exp(-vec[0]/wd);
  vec[1] = 0.5*vec[1]/M_PI;
  vec[2] = 0.5 + 0.5*erf((log(vec[2]) - flux_mu)/sqrt(2.0*flux_var));
  vec[3] = log(vec[3]/radiusmin)/log(radiusmax/radiusmin);
  vec[4] = pow(vec[4] - q_min, 2)/pow(1.0 - q_min, 2);
  vec[5] = vec[5]/M_PI;
  vec[6] = 0.5 + 0.5*erf((log(vec[6]) - dispersion_mu)/sqrt(2.0*dispersion_var));

  if(Data::get_instance().get_model_n2lines() == 1)
    vec[7] = 0.5 + 0.5*erf((log(vec[7]) - n2flux_mu)/sqrt(2.0*n2flux_var));

}

void MyConditionalPrior::print(std::ostream& out) const
{
  
  double radiusmin = exp(log(radiuslim_min) + log(radiusmax/radiuslim_min)*radiusmin_ratio);
  
  out<<exp(flux_mu)<<' '<<flux_std<<' '
     <<radiusmin<<' '<<radiusmax<<' '
     <<exp(n2flux_mu)<<' '<<n2flux_std<<' '
     <<exp(dispersion_mu)<<' '<<dispersion_std<<' '
     <<q_min<<' '<<wd<<' ';

}
