#include "Data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>

#include "Constants.h"

Data Data::instance;

Data::Data() {}

void Data::load(const char* moptions_file) {
  // Loading data comment
  std::cout<<"\nLoading data:\n";

  // Load model options
  std::fstream fin(moptions_file, std::ios::in);
  if (!fin)
    std::cerr<<"# ERROR: couldn't open file "<<moptions_file<<"."<<std::endl;

  size_t n;
  std::string line;
  std::string name;
  std::string tmp_str;
  double tmp_double;
  std::vector<double> tmp_vector;
  bool line_flag = false;
  bool lsf_fwhm_flag = false;
  bool psf_amp_flag = false;
  bool psf_fwhm_flag = false;
  bool psf_beta_flag = false;
  bool inc_flag = false;
  bool gamma_pos_flag = false;
  bool radiuslim_min_flag = false;
  while (std::getline(fin, line)) {
    std::istringstream lin(line);
    lin >> name;
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);

    if (line.empty()) {
      continue;
    } else if (name[0] == '#') {
      continue;
    } else if (name == "METADATA_FILE") {
        lin >> metadata_file;
    } else if (name == "DATA_FILE") {
      lin >> data_file;
    } else if (name == "VAR_FILE") {
      lin >> var_file;
    } else if (name == "CONVOLVE_METHOD") {
      lin >> convolve;
    } else if (name == "PSFWEIGHT") {
      while (lin >> tmp_double)
        psf_amp.push_back(tmp_double);
      psf_amp_flag = true;
    } else if (name == "PSFFWHM") {
      while (lin >> tmp_double)
        psf_fwhm.push_back(tmp_double);
      psf_fwhm_flag = true;
    } else if (name == "PSFBETA") {
      lin >> psf_beta;
      psf_beta_flag = true;
    } else if (name == "LSFFWHM") {
      lin >> lsf_fwhm;
      lsf_fwhm_flag = true;
    } else if (name == "NMAX") {
      lin >> nmax;
    } else if (name == "NFIXED") {
      lin >> tmp_str;
      std::transform(
        tmp_str.begin(), tmp_str.end(),
        tmp_str.begin(), ::toupper);
      if ((tmp_str == "FALSE") || (tmp_str == "0")) {
        nfixed = false;
      } else if ((tmp_str == "TRUE") || (tmp_str == "1")) {
        nfixed = true;
      } else {
        std::cerr<<"# ERROR: couldn't determine N_FIXED."<<std::endl;
        exit(0);
      }
    } else if (name == "VSYS_MAX") {
      lin >> vsys_max;
    } else if (name == "VSYS_GAMMA") {
      lin >> vsys_gamma;
    } else if (name == "VC_MIN") {
      lin >> vmax_min;
    } else if (name == "VC_MAX") {
      lin >> vmax_max;
    } else if (name == "MEANFLUX_MIN") {
      lin >> fluxmu_min;
    } else if (name == "MEANFLUX_MAX") {
      lin >> fluxmu_max;
    } else if (name == "LOGSIGMAFLUX_MIN") {
      lin >> lnfluxsd_min;
    } else if (name == "LOGSIGMAFLUX_MAX") {
      lin >> lnfluxsd_max;
    } else if (name == "MEANRATIOFLUX_MIN") {
      lin >> ratiofluxmu_min;
    } else if (name == "MEANRATIOFLUX_MAX") {
      lin >> ratiofluxmu_max;
    } else if (name == "LOGSIGMARATIOFLUX_MIN") {
      lin >> lnratiofluxsd_min;
    } else if (name == "LOGSIGMARATIOFLUX_MAX") {
      lin >> lnratiofluxsd_max;
    } else if (name == "SIGMAV0_MIN") {
      lin >> vdisp0_min;
    } else if (name == "SIGMAV0_MAX") {
      lin >> vdisp0_max;
    } else if (name == "QLIM_MIN") {
      lin >> qlim_min;
    } else if (name == "INC") {
      lin >> inc;
      inc_flag = true;
    } else if (name == "SIGMA0_MIN") {
      lin >> sigma_min;
    } else if (name == "SIGMA0_MAX") {
      lin >> sigma_max;
    } else if (name == "RADIUSLIM_MIN") {
      lin >> radiuslim_min; // Not passed to ConditionalPrior yet
      radiuslim_min_flag = true;
    } else if (name == "RADIUSLIM_MAX") {
      lin >> radiuslim_max;
    } else if (name == "WD_MIN") {
      lin >> wd_min;
    } else if (name == "WD_MAX") {
      lin >> wd_max;
    } else if (name == "VSLOPE_MIN") {
      lin >> vslope_min;
    } else if (name == "VSLOPE_MAX") {
      lin >> vslope_max;
    } else if (name == "VGAMMA_MIN") {
      lin >> vgamma_min;
    } else if (name == "VGAMMA_MAX") {
      lin >> vgamma_max;
    } else if (name == "VBETA_MIN") {
      lin >> vbeta_max;
    } else if (name == "VBETA_MAX") {
      lin >> vbeta_max;
    } else if (name == "VDISP_ORDER") {
      lin >> vdisp_order;
    } else if (name == "LOGVDISP0_MIN") {
      lin >> vdisp0_min;
    } else if (name == "LOGVDISP0_MAX") {
      lin >> vdisp0_max;
    } else if (name == "VDISPN_SIGMA") {
      lin >> vdispn_sigma;
    } else if (name == "SIGMA1_MIN") {
      lin >> sigma1_min;
    } else if (name == "SIGMA1_MAX") {
      lin >> sigma1_max;
    } else if (name == "MD_MIN") {
      lin >> Md_min;
    } else if (name == "MD_MAX") {
      lin >> Md_max;
    } else if (name == "WXD_MIN") {
      lin >> wxd_min;
    } else if (name == "WXD_MAX") {
      lin >> wxd_max;
    } else if (name == "CENTRE_GAMMA") {
      lin >> gamma_pos;
      gamma_pos_flag = true;
    } else if (name == "LINE") {
      tmp_vector.clear();
      while (lin >> tmp_double) {
        tmp_vector.push_back(tmp_double);
      }
      n = em_line.size();
      em_line.resize(n+1);
      em_line[n].resize(tmp_vector.size());
      for (size_t i=0; i<tmp_vector.size(); i++)
        em_line[n][i] = tmp_vector[i];
      line_flag = true;
    } else {
      std::cerr
        <<"Couldn't determine input parameter assignment for keyword: "
        <<name<<"."
        <<std::endl;
        exit(0);
    }
  }
  fin.close();

  // Check required parameters are provided
  if (!line_flag) {
    std::cerr
      <<"# ERROR: Required keyword (LINE) not provided."<<std::endl;
    exit(0);
  }

  if (!lsf_fwhm_flag) {
    std::cerr
      <<"# ERROR: Required keyword (LSFFWHM) not provided."<<std::endl;
    exit(0);
  }

  if (convolve == 0) {
    if (!psf_amp_flag) {
      std::cerr
        <<"# ERROR: Required keyword (PSFWEIGHT) not provided."<<std::endl;
      exit(0);
    }

    if (!psf_fwhm_flag) {
      std::cerr
        <<"# ERROR: Required keyword (PSFFWHM) not provided."<<std::endl;
      exit(0);
    }
  } else if (convolve == 1) {
    if (!psf_beta_flag) {
      std::cerr
        <<"# ERROR: Required keyword (PSFBETA) not provided."<<std::endl;
      exit(0);
    }
  }

  if (!inc_flag) {
    std::cerr<<"# ERROR: Required keyword (INC) not provided."<<std::endl;
    exit(0);
  }

  // sigma cutoff parameter for blobs
  sigma_cutoff = 5.0;

  // PSF convolution method message
  for(size_t i=0; i<psf_fwhm.size(); i++)
    psf_sigma.push_back(psf_fwhm[i]/sqrt(8.0*log(2.0)));

  // LSF in wavelength (needs to be redshift corrected)
  lsf_sigma = lsf_fwhm/sqrt(8.0*log(2.0));

  // Model choice (only for testing)
  model = 0;

  // Override blob parameters for disk model
  if (model == 1) {
    nmax = 0;
    nfixed = true;
  }

  // Spatial sampling of cube
  sample = 1;

  // Read in the metadata
  fin.open(metadata_file, std::ios::in);
  if (!fin)
    std::cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<std::endl;
  fin >> ni >> nj;
  fin >> nr;
  fin >> x_min >> x_max >> y_min >> y_max;
  fin >> r_min >> r_max;
  fin.close();
  std::cout<<"Metadata Loaded..."<<std::endl;

  // Make sure maximum > minimum
  if (x_max <= x_min || y_max <= y_min)
    std::cerr<<"# ERROR: strange input in "<<metadata_file<<"."<<std::endl;

  data = read_cube(data_file);
  std::cout<<"Image Loaded...\n";

  var = read_cube(var_file);
  std::cout<<"Variance Loaded...\n";

  /*
    Determine the valid data pixels
    Considered valid if sigma > 0.0 and there is at least 1 non-zero value.
   */
  valid.assign(1, std::vector<int>(2));
  double tmp_im, tmp_sig;
  std::vector<int> tmp_vec(2);
  nv = 0;
  for (size_t i=0; i<data.size(); i++) {
    for(size_t j=0; j<data[i].size(); j++) {
      tmp_im = 0.0;
      tmp_sig = 0.0;
      for (size_t r=0; r<data[i][j].size(); r++) {
        if (data[i][j][r] != 0.0) { tmp_im = 1.0; }
        tmp_sig += var[i][j][r];
	    }

      // Add valid pixels to array
      if ((tmp_im == 1.0) && (tmp_sig > 0.0)) {
        tmp_vec[0] = i;
        tmp_vec[1] = j;

        if (nv != 0)
          valid.push_back(tmp_vec);
        else
          valid[0] = tmp_vec;

        nv += 1;
      }
    }
  }
  std::cout<<"Valid pixels determined...\n\n";

  // Compute pixel widths
  dx = (x_max - x_min)/nj;
  dy = (y_max - y_min)/ni;
  dr = (r_max - r_min)/nr;
  for (size_t i=0; i<psf_sigma.size(); i++) {
    psf_sigma_overdx.push_back(psf_sigma[i]/dx);
    psf_sigma_overdy.push_back(psf_sigma[i]/dy);
  }

  // Calculate geometric widths
  pixel_width = sqrt(dx*dy);
  image_width = sqrt((x_max - x_min)*(y_max - y_min));

  // rc_max for TruncatedExponential distribution
  rc_max = sqrt(pow(abs(x_max - x_min), 2) + pow(abs(y_max - y_min), 2))/cos(inc);

  // Image centres
  x_imcentre = (x_min + x_max)/2.0;
  y_imcentre = (y_min + y_max)/2.0;

  // Array padding to help edge problems
  x_pad = (int)ceil(sigma_pad*psf_sigma[0]/dx);
  y_pad = (int)ceil(sigma_pad*psf_sigma[0]/dy);
  ni += 2*y_pad;
  nj += 2*x_pad;
  x_pad_dx = x_pad*dx;
  y_pad_dy = y_pad*dy;

  x_min -= x_pad_dx;
  x_max += x_pad_dx;
  y_min -= y_pad_dy;
  y_max += y_pad_dy;

  // Compute spatially oversampled parameters
  dxos = abs(dx)/sample;
  dyos = abs(dy)/sample;
  x_pados = (int)ceil(sigma_pad*psf_sigma[0]/dxos);
  y_pados = (int)ceil(sigma_pad*psf_sigma[0]/dyos);
  nios = sample*ni;
  njos = sample*nj;
  x_pad_dxos = x_pados*dxos;
  y_pad_dyos = y_pados*dyos;


  // Construct defaults that are dependent on data
  if (!radiuslim_min) {
    radiuslim_min = pixel_width;
  }
  if (!gamma_pos_flag) {
    gamma_pos = 0.1*image_width;
  }

  // Compute x, y, r arrays
  compute_ray_grid();

  summarise_model();
}

std::vector< std::vector< std::vector<double> > > Data::arr_3d() {
  std::vector< std::vector< std::vector<double> > > arr;
  // Create 3D array with shape (ni, nj, nr)
  arr.resize(ni);
  for (int i=0; i<ni; i++) {
    arr[i].resize(nj);
    for(int j=0; j<nj; j++) {
      arr[i][j].resize(nr);
    }
  }

  return arr;
}

std::vector< std::vector< std::vector<double> > >
  Data::read_cube (std::string filepath) {
  // Read data file
  std::vector< std::vector< std::vector<double> > > cube = arr_3d();
  std::fstream fin(filepath, std::ios::in);

  if (!fin)
    std::cerr<<"# ERROR: couldn't open file "<<filepath<<"."<<std::endl;

  for (size_t i=0; i<cube.size(); i++)
    for (size_t j=0; j<cube[i].size(); j++)
      for (size_t r=0; r<cube[i][j].size(); r++)
        fin >> cube[i][j][r];
  fin.close();

  return cube;
}

void Data::compute_ray_grid() {
  // Make vectors of the correct size
  x.assign(ni, std::vector<double>(nj));
  y.assign(ni, std::vector<double>(nj));

  for (size_t i=0; i<x.size(); i++) {
    for (size_t j=0; j<x[i].size(); j++) {
      x[i][j] = x_min + (j + 0.5)*dx;
      y[i][j] = y_min + (i + 0.5)*dy; // Assuming origin=lower
    }
  }

  r.assign(nr, 0.0);
  for(size_t k=0; k<r.size(); k++)
    r[k] = r_min + (k + 0.5)*dr;
}

void Data::summarise_model() {
  // Print summarised model to terminal
  std::string dashline =
    "-----------------------------------------------------------";

  std::cout<<dashline<<std::endl;
  std::cout<<"Joint Prior Distribution"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Constants"<<std::endl;
  std::cout<<dashline<<std::endl;

  std::cout<<"Emission Line(s):"<<std::endl;
  for (size_t l=0; l<em_line.size(); l++) {
    std::cout<<"  Line: "<<em_line[l][0]<<std::endl;
    for (size_t i=0; i<(em_line[l].size() - 1)/2; i++) {
      std::cout<<"   Constrained Line: "<<em_line[l][2*i+1];
      std::cout<<", Factor: "<<em_line[l][2*i+2];
      std::cout<<std::endl;
    }
  }

  if (nfixed)
    std::cout<<"N: "<<nmax<<std::endl;
  std::cout<<"PSF Profile: ";
  if (convolve == 0) {
    std::cout<<"Sum of concentric Gaussians."<<std::endl;
    std::cout<<"PSF_WEIGHTS: ";
    for (size_t i=0; i<psf_amp.size(); i++)
      std::cout<<psf_amp[i]<<" ";
    std::cout<<std::endl;
  } else if (convolve == 1) {
    std::cout<<"Moffat (WARNING: Not safe for multi-threading)."<<std::endl;
    std::cout<<"PSF_BETA: "<<psf_beta<<std::endl;
  }
  std::cout<<"PSF_FWHM: ";
  for (size_t i=0; i<psf_fwhm.size(); i++)
    std::cout<<psf_fwhm[i]<<" ";
  std::cout<<std::endl;
  std::cout<<"LSF_FWHM (Gauss Instr. Broadening): "<<lsf_fwhm<<std::endl;
  std::cout<<"i: "<<inc<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Global Parameters"<<std::endl;
  std::cout<<dashline<<std::endl;
  std::cout
    <<"x_c ~ Cauchy("<<x_imcentre<<", "<<gamma_pos<<")"
    <<"T("<<x_min + x_pad_dx<<", "<<x_max - x_pad_dx<<")"
    <<std::endl;
  std::cout
    <<"y_c ~ Cauchy("<<x_imcentre<<", "<<gamma_pos<<")"
    <<"T("<<y_min + y_pad_dy<<", "<<y_max - y_pad_dy<<")"
    <<std::endl;

  std::cout<<"Theta ~ Uniform(0, 2pi)"<<std::endl;

  if (!nfixed)
    std::cout<<"N ~ Loguniform(0, "<<nmax<<")"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Blob hyperparameters"<<std::endl;
  std::cout<<dashline<<std::endl;
  std::cout<<"mu_r ~ Loguniform("<<wd_min<<", "<<wd_max<<")"<<std::endl;
  std::cout<<"mu_F ~ Loguniform("<<fluxmu_min<<", "<<fluxmu_max<<")"<<std::endl;
  std::cout<<"sigma_F ~ Loguniform("<<lnfluxsd_min<<", "<<lnfluxsd_max<<")"<<std::endl;
  std::cout<<"W_max ~ Loguniform("<<radiuslim_min<<", "<<radiuslim_max<<")"<<std::endl;
  std::cout<<"q_min ~ Uniform("<<qlim_min<<", "<<"1)"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Blob parameters"<<std::endl;
  std::cout<<dashline<<std::endl;
  std::cout<<"F_j ~ Lognormal(mu_F, sigma_F^2)"<<std::endl;
  std::cout<<"r_j ~ Exponential(mu_r)"<<std::endl;
  std::cout<<"Theta_j ~ Uniform(0, 2pi)"<<std::endl;
  std::cout<<"w_j ~ Loguniform("<<radiuslim_min<<", W_max)"<<std::endl;
  std::cout<<"q_j ~ Triangular(q_min, 1)"<<std::endl;
  std::cout<<"phi_j ~ Uniform(0, pi)"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Velocity profile parameters"<<std::endl;
  std::cout<<dashline<<std::endl;
  std::cout<<"v_sys ~ Cauchy(0, "<<vsys_max<<")"<<std::endl;
  std::cout<<"w_j ~ Loguniform("<<vmax_min<<", "<<vmax_max<<")"<<std::endl;
  std::cout<<"r_t ~ Loguniform("<<vslope_min<<", "<<vslope_max<<")"<<std::endl;
  std::cout<<"gamma_v ~ Loguniform("<<vgamma_min<<", "<<vgamma_max<<")"<<std::endl;
  std::cout<<"beta_v ~ Uniform("<<vbeta_min<<", "<<vbeta_max<<")"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Velocity dispersion profile parameters"<<std::endl;
  std::cout<<dashline<<std::endl;
  std::cout<<"sigma_v0 ~ Loguniform("<<vdisp0_min<<", "<<vdisp0_max<<")"<<std::endl;
  std::cout<<"sigma_vn ~ Normal(0, "<<vdispn_sigma<<")"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<"Systematic noise parameters"<<std::endl;
  std::cout<<dashline<<std::endl;
  std::cout<<"sigma_0 ~ Loguniform("<<sigma_min<<", "<<sigma_max<<")"<<std::endl;

  std::cout<<dashline<<std::endl;
  std::cout<<std::endl;
}
