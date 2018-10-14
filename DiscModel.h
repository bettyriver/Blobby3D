#ifndef BLOBBY3D_DISCMODEL_H_
#define BLOBBY3D_DISCMODEL_H_

#include <vector>

#include "DNest4/code/DNest4.h"
#include "BlobConditionalPrior.h"
#include "Conv.h"
#include "Data.h"

class DiscModel {
  private:
    DNest4::RJObject<BlobConditionalPrior> blobs;

    // TODO: Find elegant way to initialise constructor in DiscModel.cpp
    Conv conv = Conv(
      Data::get_instance().get_convolve(),
      Data::get_instance().get_psf_amp(),
      Data::get_instance().get_psf_fwhm(),
      Data::get_instance().get_psf_beta(),
      Data::get_instance().get_psf_sigma(),
      Data::get_instance().get_psf_sigma_overdx(),
      Data::get_instance().get_psf_sigma_overdy(),
      Data::get_instance().get_ni(),
      Data::get_instance().get_nj(),
      Data::get_instance().get_nr(),
      Data::get_instance().get_dx(),
      Data::get_instance().get_dy(),
      Data::get_instance().get_x_pad(),
      Data::get_instance().get_y_pad()
      ); // setup convolution kernels

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
      Parameters
    */
    int model;

    double xcd;
    double x_imagecentre;

    double ycd;
    double y_imagecentre;

    double gamma_pos;

    double Md;
    double Md_min, Md_max;

    double wxd;
    double wxd_min, wxd_max;

    double vsys;
    double gamma_vsys;
    double vsys_max;

    double vmax;
    double vmax_min, vmax_max;

    double vslope;
    double vslope_min, vslope_max;

    double vgamma;
    double vgamma_min, vgamma_max;

    double vbeta;
    double vbeta_min, vbeta_max;

    int vdisp_order;
    std::vector<double> vdisp_param;
    double vdisp0_min, vdisp0_max;

    double gama_inc;
    double inc, pa;

    double sigma0;
    double sigma0_min, sigma0_max;

    double sigma1;
    double sigma1_min, sigma1_max;

    // Prior distributions
    DNest4::Uniform prior_pa;
    DNest4::TruncatedCauchy prior_xc;
    DNest4::TruncatedCauchy prior_yc;

    DNest4::TruncatedCauchy prior_vsys;
    DNest4::LogUniform prior_vmax;
    DNest4::LogUniform prior_vslope;
    DNest4::LogUniform prior_vgamma;
    DNest4::Uniform prior_vbeta;

    DNest4::Uniform prior_vdisp0;
    DNest4::Gaussian prior_vdisp;

    DNest4::LogUniform prior_sigma0;
    DNest4::LogUniform prior_sigma1;

    DNest4::LogUniform prior_Md;
    DNest4::LogUniform prior_wxd;

    // Perturb booleans
    bool array_perturb;
    bool vel_perturb;
    bool vdisp_perturb;
    bool disc_flux_perturb;
    bool blob_perturb;
    bool noise_perturb;

  public:
    DiscModel();

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

#endif  // BLOBBY3D_DISCMODEL_H_

