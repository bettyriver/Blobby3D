#ifndef DNest4_GalaxyField_Data
#define DNest4_GalaxyField_Data

#include <vector>
#include <string>

class Data
{
 private:

  // constants
  double C, HA, N2UPP, N2LOW, INFTY;

  // files
  std::string metadata_file;
  std::string cube_file;
  std::string var_file;

  // model parameters
  int model;
  int nmax;
  bool nfixed;
  int model_n2lines;
  double prior_inc;
  int convolve;
  double vsys_gamma, vsys_max;
  double vmax_min, vmax_max;
  double fluxmu_min, fluxmu_max;
  double lnfluxsd_min, lnfluxsd_max;
  double vdispmu_min, vdispmu_max;
  double lnvdispsd_min, lnvdispsd_max;
  double qlim_min;
  double sigma_min, sigma_max;
  
  // sampling
  double sample;
  
  // Number of pixels
  int ni, nj;
  
  // Number of wavelength bins
  int nr;
  
  // Number of valid spaxels
  int nv;

  // Coordinates of image boundaries
  double x_min, x_max, y_min, y_max;

  // Wavelength boundaries
  double r_min, r_max;

  // Pixel widths
  double dx, dy;
		
  // Wavelength bin widths
  double dr;

  // Total bin volume
  double db;

  // PSF
  double psf_fwhm;
  double psf_beta;
  double psf_sigma;
  double sigma_cutoff;

  // LSF
  double lsf_fwhm;
  double lsf_sigma;

  // Padding due to convolution
  double x_pad, y_pad;
  double x_pad_dx, y_pad_dy;
		
  // Analytical paramaters
  double sum_flux;

  // Spatial oversampling
  double dxos, dyos;
  double x_pados, y_pados;
  double nios, njos;
  double x_pad_dxos, y_pad_dyos;

  // Coordinates of pixel centers
  std::vector< std::vector<double> > x_rays;
  std::vector< std::vector<double> > y_rays;
  std::vector<double> r_rays;

  // The pixels
  std::vector< std::vector< std::vector<double> > > image;

  // Sigma map
  std::vector< std::vector< std::vector<double> > > var_cube;

  // Valid spaxels
  std::vector< std::vector<int> > valid;

  // Private functions
  std::vector< std::vector< std::vector<double> > > arr_3d();
  void compute_ray_grid();

 public:
  Data();
  void load(const char* moptions_file);

  // Getters
  double get_C() const { return C; }
  double get_HA() const { return HA; }
  double get_N2UPP() const { return N2UPP; }
  double get_N2LOW() const { return N2LOW; }
  int get_model() const { return model; }
  int get_nmax() const { return nmax; }
  bool get_nfixed() const { return nfixed; }
  double get_prior_inc() const { return prior_inc; }
  int get_model_n2lines() const { return model_n2lines; }
  int get_convolve() const { return convolve; }
  int get_ni() const { return ni; }
  int get_nj() const { return nj; }
  int get_nr() const { return nr; }
  int get_nv() const { return nv; } 
  double get_x_min() const { return x_min; }
  double get_x_max() const { return x_max; }
  double get_y_min() const { return y_min; }
  double get_y_max() const { return y_max; }
  double get_r_min() const { return r_min; }
  double get_r_max() const { return r_max; }
  double get_dx() const { return dx; }
  double get_dy() const { return dy; }
  double get_dr() const { return dr; }
  double get_db() const { return db; }
  double get_psf_fwhm() const { return psf_fwhm; }
  double get_psf_beta() const { return psf_beta; }
  double get_psf_sigma() const { return psf_sigma; }
  double get_lsf_sigma() const { return lsf_sigma; }
  double get_sigma_cutoff() const { return sigma_cutoff; }
  double get_vsys_gamma() const { return vsys_gamma; }
  double get_vsys_max() const { return vsys_max; }
  double get_vmax_min() const { return vmax_min; }
  double get_vmax_max() const { return vmax_max; }
  double get_fluxmu_min() const { return fluxmu_min; }
  double get_fluxmu_max() const { return fluxmu_max; }
  double get_lnfluxsd_min() const { return lnfluxsd_min; }
  double get_lnfluxsd_max() const { return lnfluxsd_max; }
  double get_vdispmu_min() const { return vdispmu_min; }
  double get_vdispmu_max() const { return vdispmu_max; }
  double get_lnvdispsd_min() const { return lnvdispsd_min; }
  double get_lnvdispsd_max() const { return lnvdispsd_max; }
  double get_qlim_min() const { return qlim_min; }
  double get_sigma_min() const { return sigma_min; }
  double get_sigma_max() const { return sigma_max; }
  int get_x_pad() const { return x_pad; }
  int get_y_pad() const { return y_pad; }
  double get_x_pad_dx() const { return x_pad_dx; }
  double get_y_pad_dy() const { return y_pad_dy; }
  double get_sum_flux() const {return sum_flux; }
  
  int get_sample() const { return sample; }
  double get_dxos() const { return dxos; }
  double get_dyos() const { return dyos; }
  double get_nios() const { return nios; }
  double get_njos() const { return njos; }
  double get_x_pados() const { return x_pados; }
  double get_y_pados() const { return y_pados; }
  double get_x_pad_dxos() const { return x_pad_dxos; }
  double get_y_pad_dyos() const { return y_pad_dyos; }


  const std::vector< std::vector<double> >& get_x_rays() const
  { return x_rays; }
  const std::vector< std::vector<double> >& get_y_rays() const
  { return y_rays; }
  const std::vector<double>& get_r_rays() const
  { return r_rays; }
  const std::vector< std::vector< std::vector<double> > >& get_image() const
  { return image; }
  const std::vector< std::vector< std::vector<double> > >& get_var_cube() const
  { return var_cube; }
  const std::vector< std::vector<int> >& get_valid() const
  { return valid; }

  // Singleton
 private:
  static Data instance;
 public:
  static Data& get_instance() { return instance; }
};

#endif

