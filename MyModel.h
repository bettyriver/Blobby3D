#ifndef DNest4_GalaxyField_MyModel
#define DNest4_GalaxyField_MyModel

#include <vector>
#include "DNest4/code/DNest4.h"
#include "MyConditionalPrior.h"
#include "Conv.h"

class MyModel
{
 private:
  
  Conv conv; // setup convolution kernels
  DNest4::RJObject<MyConditionalPrior> objects;


  /*
    Arrays
  */
  std::vector< std::vector< std::vector<double> > > image;
  std::vector< std::vector< std::vector<double> > > imageos;
  std::vector< std::vector< std::vector<double> > > convolved;
  std::vector< std::vector< std::vector<double> > > convolvedos;
  std::vector< std::vector<double> > rel_lambda;
  std::vector< std::vector<double> > rad;

  void calculate_image();

  /*
    Global parameters
  */

  // Disk parameters
  double xcd;
  double x_imagecentre;

  double ycd;
  double y_imagecentre;

  double gamma_pos;

  double Md;
  double Md_min, Md_width;

  double Mdn2;
  double Mdn2_min, Mdn2_width;

  double wxd;
  double wxd_min, wxd_width;

  double wxd_n2;
  double wxd_n2_min, wxd_n2_width;

  double vsys;
  double gamma_vsys;
  double vsys_min, vsys_width;

  double vmax;
  double vmax_min, vmax_width;
		
  double vslope;
  double vslope_min, vslope_width;

  double vgamma;
  double vgamma_min, vgamma_width;

  double vbeta;
  double vbeta_min, vbeta_width;
  
  double sigmad_lambda;
  double sigmad_lambda_min, sigmad_lambda_width;

  double inc, pa;

  // Noise parameter
  double sigma0;
  double sigma0_min, sigma0_width;

  double sigma1;
  double sigma1_min, sigma1_width;

  // perturb booleans
  bool disk_perturb, rot_perturb;


 public:
  MyModel();

  // Generate the point from the prior
  void from_prior(DNest4::RNG& rng);

  // Metropolis-Hastings proposals
  double perturb(DNest4::RNG& rng);

  // Likelihood function
  double log_likelihood() const;

  // Print to stream
  void print(std::ostream& out) const;

  // Return string with column information
  std::string description() const;
};

#endif

