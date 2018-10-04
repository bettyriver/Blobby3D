#ifndef BLOBBY3D_MyModel_H_
#define BLOBBY3D_MyModel_H_

#include <vector>

#include "DNest4/code/DNest4.h"
#include "MyConditionalPrior.h"
#include "Conv.h"

class MyModel {
  private:
    Conv conv; // setup convolution kernels
    DNest4::RJObject<MyConditionalPrior> objects;

    /*
       Arrays
    */
    std::vector< std::vector< std::vector<double> > > image;
    std::vector< std::vector< std::vector<double> > > imageos;
    std::vector< std::vector< std::vector<double> > > convolved;

    std::vector< std::vector<double> > x_shft;
    std::vector< std::vector<double> > y_shft;
    std::vector< std::vector<double> > rad;
    std::vector< std::vector<double> > cos_angle;

    std::vector< std::vector<double> > flux;
    std::vector< std::vector<double> > rel_lambda;
    std::vector< std::vector<double> > vdisp;

    void calculate_image();

    // Construct cube from maps
    void calculate_shifted_arrays();

    void calculate_flux();
    void add_disc_flux();
    void add_blob_flux(std::vector< std::vector<double> >& components);

    void calculate_vdisp();
    void calculate_rel_lambda();
    void construct_cube();
    void clear_cube();
    void clear_flux_map();

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

    double wxd;
    double wxd_min, wxd_width;

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

    int vdisp_order;
    double vdisp0_min, vdisp0_width;
    std::vector<double> vdisp_param;

    double gama_inc;
    double inc, pa;

    // Noise parameter
    double sigma0;
    double sigma0_min, sigma0_width;

    double sigma1;
    double sigma1_min, sigma1_width;

    // perturb booleans
    bool array_perturb;
    bool vel_perturb;
    bool vdisp_perturb;
    bool disc_flux_perturb;
    bool blob_perturb;
    bool noise_perturb;

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

#endif  // BLOBBY3D_MyModel_H_

