#ifndef BLOBBY3D_DATA_H_
#define BLOBBY3D_DATA_H_

#include <cmath>
#include <vector>
#include <string>


class Data
{
 private:

  // files
  std::string metadata_file;
  std::string cube_file;
  std::string var_file;

  // model parameters
  int model = 0;
  int nmax = 300;
  bool nfixed = false;
  int convolve = 0;
  double vsys_gamma = 30.0;
  double vsys_max = 150.0;
  double vmax_min = 40.0;
  double vmax_max = 400.0;
  double fluxmu_min = 1e-3;
  double fluxmu_max = 1e3;
  double lnfluxsd_min = 0.03;
  double lnfluxsd_max = 3.0;
  double vdispmu_min = 1.0;
  double vdispmu_max = 200.0;
  double lnvdispsd_min = 0.01;
  double lnvdispsd_max = 1.0;
  double qlim_min = 0.2;
  double sigma_min = 1e-12;
  double sigma_max = 1e12;
  double gama_inc;

  // Hard-coded parameters
  // TODO: Read in parameters from file
  double radiuslim_max = 30.0;
  double wd_min = 0.03;
  double wd_max = 30.0;
  double vslope_min = 0.03;
  double vslope_max = 30.0;
  double vgamma_min = 1.0;
  double vgamma_max = 100.0;
  double vbeta_min = -0.75;
  double vbeta_max = 0.75;
  int vdisp_order = 1;
  double vdisp0_min = log(1.0);
  double vdisp0_max = log(200.0);
  double vdispn_sigma = 0.2;
  double sigma1_min = 1E-12;
  double sigma1_max = 1E0;
  double Md_min = 1E-3;
  double Md_max = 1E3;
  double wxd_min = 0.3;
  double wxd_max = 30.0;
  double gamma_pos = 2.0;

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
  std::vector<double> psf_amp;
  std::vector<double> psf_fwhm;
  double psf_beta;
  std::vector<double> psf_sigma;
  double sigma_cutoff;
  double sigma_pad;

  // LSF
  double lsf_fwhm;
  double lsf_sigma;

  // Padding due to convolution
  double x_pad, y_pad;
  double x_pad_dx, y_pad_dy;
  std::vector<double> psf_sigma_overdx;
  std::vector<double> psf_sigma_overdy;

  // Spatial oversampling
  double dxos, dyos;
  double x_pados, y_pados;
  double nios, njos;
  double x_pad_dxos, y_pad_dyos;

  // Geometric widths
  double pixel_width;
  double image_width;
  double x_imcentre, y_imcentre;

  // Coordinates of pixel centers
  std::vector< std::vector<double> > x_rays;
  std::vector< std::vector<double> > y_rays;
  std::vector<double> r_rays;

  // Data
  std::vector< std::vector< std::vector<double> > > data;
  std::vector< std::vector< std::vector<double> > > var;

  // Valid spaxels
  std::vector< std::vector<int> > valid;

  // Private functions
  std::vector< std::vector< std::vector<double> > > arr_3d();
  std::vector< std::vector< std::vector<double> > >
    read_cube(std::string filepath);
  void compute_ray_grid();

 public:
  Data();
  void load(const char* moptions_file);

  // Getters
  int get_model() const { return model; }
  int get_nmax() const { return nmax; }
  bool get_nfixed() const { return nfixed; }
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
  std::vector<double> get_psf_amp() const { return psf_amp; }
  std::vector<double> get_psf_fwhm() const { return psf_fwhm; }
  double get_psf_beta() const { return psf_beta; }
  std::vector<double> get_psf_sigma() const { return psf_sigma; }
  double get_lsf_sigma() const { return lsf_sigma; }
  double get_sigma_cutoff() const { return sigma_cutoff; }
  double get_sigma_pad() const { return sigma_pad; }
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
  double get_gama_inc() const { return gama_inc; }
  double get_radiuslim_max() const { return radiuslim_max; }
  double get_wd_min() const { return wd_min; }
  double get_wd_max() const { return wd_max; }
  double get_vslope_min() { return vslope_min; };
  double get_vslope_max() { return vslope_max; }
  double get_vgamma_min() { return vgamma_min; }
  double get_vgamma_max() { return vgamma_max; }
  double get_vbeta_min() { return vbeta_min; }
  double get_vbeta_max() { return vbeta_max; }
  int get_vdisp_order() { return vdisp_order; }
  double get_vdisp0_min() { return vdisp0_min; }
  double get_vdisp0_max() { return vdisp0_max; }
  double get_vdispn_sigma() { return vdispn_sigma; }
  double get_sigma1_min() { return sigma1_min; }
  double get_sigma1_max() { return sigma1_max; }
  double get_Md_min() { return Md_min; }
  double get_Md_max() { return Md_max; }
  double get_wxd_min() { return wxd_min; }
  double get_wxd_max() { return wxd_max; }
  double get_gamma_pos() { return gamma_pos; }

  int get_x_pad() const { return x_pad; }
  int get_y_pad() const { return y_pad; }
  double get_x_pad_dx() const { return x_pad_dx; }
  double get_y_pad_dy() const { return y_pad_dy; }

  double get_pixel_width() const { return pixel_width; }
  double get_image_width() const { return image_width; }
  double get_x_imcentre() const { return x_imcentre; }
  double get_y_imcentre() const { return y_imcentre; }

  std::vector<double> get_psf_sigma_overdx() const { return psf_sigma_overdx; }
  std::vector<double> get_psf_sigma_overdy() const { return psf_sigma_overdy; }

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
  const std::vector< std::vector< std::vector<double> > >& get_data() const
  { return data; }
  const std::vector< std::vector< std::vector<double> > >& get_var() const
  { return var; }
  const std::vector< std::vector<int> >& get_valid() const
  { return valid; }

  // Singleton
 private:
  static Data instance;
 public:
  static Data& get_instance() { return instance; }
};

#endif  // BLOBBY3D_DATA_H_

