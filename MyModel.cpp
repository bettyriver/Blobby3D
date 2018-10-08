#include "MyModel.h"

#include <cmath>

#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "LookupExp.h"
#include "LookupErf.h"
#include "Conv.h"
#include "Constants.h"

// TODO: Remove references to sigma1 throughout code.
// Partial fix: not perturbing.

// TODO: Remove references to inc throughout code.
// Partial fix: not perturbing.

// TODO: Resolve issue regarding 7 vs 6 blob parameters.

/*
  Public
*/
MyModel::MyModel()
    :objects(
      7, Data::get_instance().get_nmax(),
      Data::get_instance().get_nfixed(),
      MyConditionalPrior(
        log(Data::get_instance().get_fluxmu_min()),
        log(Data::get_instance().get_fluxmu_max()),
        Data::get_instance().get_lnfluxsd_min(),
        Data::get_instance().get_lnfluxsd_max(),
        Data::get_instance().get_pixel_width(), 30.0,
        0.03, 30.0,
        0.2,
        Data::get_instance().get_hp_step()
        ),
      DNest4::PriorType::log_uniform) {
  /*
    initialise arrays
  */
  size_t ni = Data::get_instance().get_ni();
  size_t nj = Data::get_instance().get_nj();
  size_t nr = Data::get_instance().get_nr();
  size_t x_pad = Data::get_instance().get_x_pad();
  size_t y_pad = Data::get_instance().get_y_pad();

  // model cube
  image.resize(ni);
  for (size_t i=0; i<ni; ++i) {
    image[i].resize(nj);
    for (size_t j=0; j<nj; ++j) {
      image[i][j].resize(nr);
    }
  }

  // Convolved cube
  convolved.resize(ni - 2*y_pad);
  for (size_t i=0; i<convolved.size(); ++i) {
    convolved[i].resize(nj - 2*x_pad);
    for (size_t j=0; j<convolved[i].size(); ++j) {
        convolved[i][j].resize(nr);
    }
  }

  // Shifted arrays
  x_shft.assign(ni, std::vector<double>(nj));
  y_shft.assign(ni, std::vector<double>(nj));

  // Radius and angle arrays
  rad.assign(ni, std::vector<double>(nj));
  cos_angle.assign(ni, std::vector<double>(nj));

  // Flux profile
  flux.assign(ni, std::vector<double>(nj));

  // relative lambda for velocity
  rel_lambda.assign(ni, std::vector<double>(nj));

  // velocity dispersion
  vdisp.assign(ni, std::vector<double>(nj));
}

void MyModel::from_prior(DNest4::RNG& rng) {
  // Get local variables from data
  const int model = Data::get_instance().get_model();
  const double x_min = Data::get_instance().get_x_min();
  const double x_max = Data::get_instance().get_x_max();
  const double y_min = Data::get_instance().get_y_min();
  const double y_max = Data::get_instance().get_y_max();
  const double x_pad_dx = Data::get_instance().get_x_pad_dx();
  const double y_pad_dy = Data::get_instance().get_y_pad_dy();
  const double sum_flux = Data::get_instance().get_sum_flux();
  const double vsys_max = Data::get_instance().get_vsys_max();
  const double gama_inc = Data::get_instance().get_gama_inc();
  const double image_width = Data::get_instance().get_image_width();

  /*
    Limits: global parameters
  */
  // Error
  sigma0_min = Data::get_instance().get_sigma_min();
  sigma0_width = Data::get_instance().get_sigma_max()/sigma0_min;

  sigma1_min = 1E-12;
  sigma1_width = 1E0/sigma1_min;

  // Position
  gamma_pos = 0.1*image_width;

  // Velocity
  gamma_vsys = Data::get_instance().get_vsys_gamma();

  vmax_min = Data::get_instance().get_vmax_min();
  vmax_width = Data::get_instance().get_vmax_max()/vmax_min;

  vslope_min = 0.03;
  vslope_width = 30.0/vslope_min;

  vgamma_min = 1.0;
  vgamma_width = 100.0/vgamma_min;

  vbeta_min = -0.75;
  vbeta_width = 1.5;

  vdisp0_min = log(1.0);
  vdisp0_width = log(200.0/1.0);

  /*
    Limits: disk parameters
  */
  if (model != 0) {
      Md_min = 1E-2*sum_flux;
      Md_width = 10.0*sum_flux/Md_min;
    }

  /*
    Initialise: global parameters
  */

  // Initialise: Error
  sigma0 = exp(log(sigma0_min) + log(sigma0_width)*rng.rand());
  sigma1 = exp(log(sigma1_min) + log(sigma1_width)*rng.rand());

  // Initialise: Position
  x_imagecentre = (x_max + x_min)/2.0;
  DNest4::Cauchy cauchy_xc(x_imagecentre, gamma_pos);
  do {
    xcd = cauchy_xc.generate(rng);
  } while ((xcd < x_min + x_pad_dx) || (xcd > x_max - x_pad_dx));

  y_imagecentre = (y_max + y_min)/2.0;
  DNest4::Cauchy cauchy_yc(y_imagecentre, gamma_pos);
  do {
    ycd = cauchy_yc.generate(rng);
  } while ((ycd < y_min + y_pad_dy) || (ycd > y_max - y_pad_dy));

  pa = 2.0*M_PI*rng.rand();
  // inc = 0.5*M_PI*rng.rand();
  inc = gama_inc;

  wxd_min = 0.3;
  wxd_width = 30.0/wxd_min;

  /*
    Initialise: Velocity
  */
  DNest4::Cauchy cauchy_vsys(0.0, gamma_vsys);
  do {
    vsys = cauchy_vsys.generate(rng);
  } while ((vsys < -vsys_max) || (vsys >  vsys_max));

  vmax = exp(log(vmax_min) + log(vmax_width)*rng.rand());
  vslope = exp(log(vslope_min) + log(vslope_width)*rng.rand());
  vgamma = exp(log(vgamma_min) + log(vgamma_width)*rng.rand());
  vbeta = vbeta_min + vbeta_width*rng.rand();

  /*
    Initialise: blobs
  */
  if (model == 0) {
    do {
      objects.from_prior(rng);
    } while (objects.get_components().size() == 0);
  } else if (model == 2) {
    objects.from_prior(rng);
  }

  /*
    Initialise: disk
  */
  if ((model == 1) || (model == 2)) {
    Md = exp(log(Md_min) + log(Md_width)*rng.rand());
    wxd = exp(log(wxd_min) + log(wxd_width)*rng.rand());
  } else {
    Md = 0.0;
    wxd = 0.0;
  }

  // Dispersion
  vdisp_order = 1;
  DNest4::Gaussian gaussian_vdisp(0.0, 0.2);
  vdisp_param.assign(vdisp_order+1, 0.0);
  vdisp_param[0] = vdisp0_min + vdisp0_width*rng.rand();
  for(int v=1; v<vdisp_order+1; v++)
      vdisp_param[v] = gaussian_vdisp.generate(rng);

  /*
    Calculate image
  */
  array_perturb = true;
  vel_perturb = true;
  vdisp_perturb = true;
  blob_perturb = true;
  noise_perturb = true;

  if (model == 0)
    disc_flux_perturb = false;
  else
    disc_flux_perturb = true;

  calculate_image();
}

double MyModel::perturb(DNest4::RNG& rng) {
  const int model = Data::get_instance().get_model();
  const double x_min = Data::get_instance().get_x_min();
  const double x_max = Data::get_instance().get_x_max();
  const double y_min = Data::get_instance().get_y_min();
  const double y_max = Data::get_instance().get_y_max();
  const double x_pad_dx = Data::get_instance().get_x_pad_dx();
  const double y_pad_dy = Data::get_instance().get_y_pad_dy();
  const double disc_step = Data::get_instance().get_disc_step();
  const double sigma_step = Data::get_instance().get_sigma_step();

  DNest4::Cauchy cauchy_xc(x_imagecentre, gamma_pos);
  DNest4::Cauchy cauchy_yc(y_imagecentre, gamma_pos);
  DNest4::Cauchy cauchy_vsys(0.0, gamma_vsys);
  DNest4::Gaussian gaussian_vdisp(0.0, 0.2);

  double logH = 0.0;
  double rnd = rng.rand();

  array_perturb = false;
  vel_perturb = false;
  vdisp_perturb = false;
  disc_flux_perturb = false;

  if (rnd <= 0.9) {
    // skip disk perturb for model 0
    if (model == 0)
      rnd = 0.8*rng.rand();
    else if (model == 1)
      rnd = 1.0 - 0.3*rng.rand();
    else
      rnd = rng.rand();

    if (rnd <= 0.7) {
      // Perturb blob parameters
      logH += objects.perturb(rng);
      if ((model == 0) & (objects.get_components().size() == 0))
        return logH = -1E300;
    } else if (rnd <= 0.8) {
      // Perturb disc parameters
      int which = rng.rand_int(10);

      switch (which) {
        case 0:
          logH += cauchy_xc.perturb(xcd, rng);
          if ((xcd < x_min + x_pad_dx) || (xcd > x_max - x_pad_dx))
            return logH = -1E300;
          array_perturb = true;
          break;
        case 1:
          logH += cauchy_yc.perturb(ycd, rng);
          if((ycd < y_min + y_pad_dy) || (ycd > y_max - y_pad_dy))
            return logH = -1E300;
          array_perturb = true;
          break;
        case 2:
          logH += cauchy_vsys.perturb(vsys, rng);
          if((vsys < -5.0*gamma_vsys) || (vsys >  5.0*gamma_vsys))
            return logH = -1E300;
          vel_perturb = true;
          break;
        case 3:
          vmax = log(vmax);
          vmax += disc_step*log(vmax_width)*rng.randh();
          vmax = DNest4::mod(vmax - log(vmax_min), log(vmax_width));
          vmax += log(vmax_min);
          vmax = exp(vmax);
          vel_perturb = true;
          break;
        case 4:
          vslope = log(vslope);
          vslope += disc_step*log(vslope_width)*rng.randh();
          vslope = DNest4::mod(vslope - log(vslope_min), log(vslope_width));
          vslope += log(vslope_min);
          vslope = exp(vslope);
          vel_perturb = true;
          break;
        case 5:
          vgamma = log(vgamma);
          vgamma += disc_step*log(vgamma_width)*rng.randh();
          vgamma = DNest4::mod(vgamma - log(vgamma_min), log(vgamma_width));
          vgamma += log(vgamma_min);
          vgamma = exp(vgamma);
          vel_perturb = true;
          break;
        case 6:
          vbeta += disc_step*vbeta_width*rng.randh();
          vbeta = DNest4::mod(vbeta - vbeta_min, vbeta_width);
          vbeta += vbeta_min;
          vel_perturb = true;
          break;
        case 7:
          pa += disc_step*2.0*M_PI*rng.randh();
          pa = DNest4::mod(pa, 2*M_PI);
          array_perturb = true;
          break;
        case 8:
          vdisp_param[0] += disc_step*vdisp0_width*rng.randh();
          vdisp_param[0] = DNest4::mod(vdisp_param[0] - vdisp0_min, vdisp0_width);
          vdisp_param[0] += vdisp0_min;
          vdisp_perturb = true;
          break;
        case 9:
          which = rng.rand_int(vdisp_order);
          logH += gaussian_vdisp.perturb(vdisp_param[which+1], rng);
          vdisp_perturb = true;
          break;
      }
    } else {
      // Perturb disc flux parameters
      // disk_perturb = true;
      disc_flux_perturb = true;

      int which = rng.rand_int(2);
      switch(which){
        case 0:
          Md = log(Md);
          Md += disc_step*log(Md_width)*rng.randh();
          Md = DNest4::mod(Md - log(Md_min), log(Md_width));
          Md += log(Md_min);
          Md = exp(Md);
          break;
        case 1:
          wxd = log(wxd);
          wxd += disc_step*log(wxd_width)*rng.randh();
          wxd = DNest4::mod(wxd - log(wxd_min), log(wxd_width)) + log(wxd_min);
          wxd = exp(wxd);
          break;
      }
    }

    // Pre-rejection trick
    if(rng.rand() >= exp(logH))
      return -1E300;
    else
      logH = 0.0;

    calculate_image();

  } else {
    int which = rng.rand_int(1);
    switch (which) {
      case 0:
        sigma0 = log(sigma0);
        sigma0 += sigma_step*log(sigma0_width)*rng.randh();
        sigma0 = DNest4::mod(sigma0 - log(sigma0_min), log(sigma0_width));
        sigma0 += log(sigma0_min);
        sigma0 = exp(sigma0);
        break;
      case 1:
        // Currently redundant
        sigma1 = log(sigma1);
        sigma1 += sigma_step*log(sigma1_width)*rng.randh();
        sigma1 = DNest4::mod(sigma1 - log(sigma1_min), log(sigma1_width));
        sigma1 += log(sigma1_min);
        sigma1 = exp(sigma1);
        break;
    }
  }

  return logH;
}

void MyModel::calculate_image() {
  /*
    Calculate image as a function of model parameters.
  */
  const int model = Data::get_instance().get_model();

  bool update;  // Determine if adding blobs
  std::vector< std::vector<double> > components;

  // Calculate position arrays
  if (array_perturb)
    calculate_shifted_arrays();

  // Calculate relative lambda array
  if (vel_perturb || array_perturb)
    calculate_rel_lambda();

  // Calculate velocity dispersion array
  if (vdisp_perturb || array_perturb)
    calculate_vdisp();

  //  Calculate flux map
  switch (model) {
    case 0:
      // Blobs only model
      update = objects.get_removed().size() == 0;
      if (array_perturb || !update) {
        clear_flux_map();
        components = objects.get_components();

      } else {
        components = objects.get_added();
      }
      add_blob_flux(components);
      break;
    case 1:
      // Disc only model
      update = false;
      clear_flux_map();
      if (disc_flux_perturb || array_perturb)
        add_disc_flux();
      break;
    case 2:
      // Disc + blobs model
      update = objects.get_removed().size() == 0;
      if (disc_flux_perturb || array_perturb || !update) {
        clear_flux_map();
        components = objects.get_components();
        add_disc_flux();

      } else {
        components = objects.get_added();
      }
      add_blob_flux(components);
      break;
  }

  construct_cube();

  convolved = conv.apply(image);
}

double MyModel::log_likelihood() const {
  const int model = Data::get_instance().get_model();
  const std::vector< std::vector< std::vector<double> > >& data = Data::get_instance().get_image();
  const std::vector< std::vector< std::vector<double> > >& var_cube = Data::get_instance().get_var_cube();
  const std::vector< std::vector<int> >& valid = Data::get_instance().get_valid();

  // If no blobs return prob = 0
  if ((model == 0) && (objects.get_components().size() == 0))
    return -1E300;

  long double logL = 0.0;
  double sigma0sq = sigma0*sigma0;

  double var;
  int i, j;

  for (size_t h=0; h<valid.size(); h++) {
    i = valid[h][0];
    j = valid[h][1];
    for (size_t r=0; r<data[i][j].size(); r++) {
      if ((data[i][j][r] != 0.0) && (var_cube[i][j][r] != 0.0)) {
        var = var_cube[i][j][r] + sigma0sq;
        logL += -0.5*log(2.0*M_PI*var);
        logL += -0.5*pow(data[i][j][r] - convolved[i][j][r], 2)/var;
      }
    }
  }

  return logL;
}

void MyModel::print(std::ostream& out) const {
  const int x_pad = Data::get_instance().get_x_pad();
  const int y_pad = Data::get_instance().get_y_pad();

  out<<std::setprecision(6);

  // Save deconvolved image
  for(size_t i=y_pad; i<image.size()-y_pad; i++)
      for(size_t j=x_pad; j<image[i].size()-x_pad; j++)
        for(size_t r=0; r<image[i][j].size(); r++)
            out << image[i][j][r] << ' ';

  // Save convolved image
  for(size_t i=0; i<convolved.size(); i++)
    for(size_t j=0; j<convolved[i].size(); j++)
      for(size_t r=0; r<convolved[i][j].size(); r++)
        out << convolved[i][j][r] << ' ';

  // Save components
  objects.print(out); out<<' ';

  // Save global variables
  out<<xcd<<' ';
  out<<ycd<<' ';
  out<<Md<<' ';
  out<<wxd<<' ';
  out<<vsys<<' ';
  out<<vmax<<' ';
  out<<vslope<<' ';
  out<<vgamma<<' ';
  out<<vbeta<<' ';
  for (int i=0; i<vdisp_order+1; i++)
    out<<vdisp_param[i]<<' ';
  out<<inc<<' ';
  out<<pa<<' ';
  out<<sigma0<<' ';
  out<<sigma1<<' ';
}

std::string MyModel::description() const {
  return std::string("objects");
}

/*
  Private
*/
void MyModel::construct_cube() {
  /*
    Create cube from maps.
  */
  const double sigma_lsfsq = pow(Data::get_instance().get_lsf_sigma(), 2);
  const std::vector<double>& wave = Data::get_instance().get_r_rays();
  const double dr = Data::get_instance().get_dr();

  double lambda;
  double sigma_lambda;
  double wlsq, invtwo_wlsq;
  double ha_cdf_min, ha_cdf_max;

  for (size_t i=0; i<image.size(); i++) {
    for (size_t j=0; j<image[i].size(); j++)  {
      if (flux[i][j] > 0.0)  {
        // Calculate mean lambda for lines
        lambda = constants::HA*rel_lambda[i][j];

        // Calculate line width
        sigma_lambda = vdisp[i][j]*constants::HA/constants::C;
        wlsq = sigma_lambda*sigma_lambda + sigma_lsfsq;
        invtwo_wlsq = 1.0/sqrt(2.0*wlsq);

        // Calculate image for 1st wavelength bin
        ha_cdf_min = LookupErf::evaluate((wave[0] - 0.5*dr - lambda)*invtwo_wlsq);
        ha_cdf_max = LookupErf::evaluate((wave[0] + 0.5*dr - lambda)*invtwo_wlsq);
        image[i][j][0] = 0.5*flux[i][j]*(ha_cdf_max - ha_cdf_min);

        // Loop through remaining bins
        for (size_t r=1; r<wave.size(); r++) {
          ha_cdf_min = ha_cdf_max;
          ha_cdf_max = LookupErf::evaluate((wave[r] + 0.5*dr - lambda)*invtwo_wlsq);
          image[i][j][r] = 0.5*flux[i][j]*(ha_cdf_max - ha_cdf_min);
        }
      } else {
        for (size_t r=0; r<wave.size(); r++)
          image[i][j][r] = 0.0;
      }

    }
  }
}

void MyModel::calculate_shifted_arrays() {
  /*
    Calculate arrays shifted by disk parameters.
  */
  const std::vector< std::vector<double> >& x = Data::get_instance().get_x_rays();
  const std::vector< std::vector<double> >& y = Data::get_instance().get_y_rays();

  double sin_pa = sin(pa);
  double cos_pa = cos(pa);
  double invcos_inc = 1.0/cos(inc);

  double xx_rot, yy_rot;

  for (size_t i=0; i<image.size(); i++) {
    for (size_t j=0; j<image[i].size(); j++) {
      // Shift
      x_shft[i][j] = x[i][j] - xcd;
      y_shft[i][j] = y[i][j] - ycd;

      // rotate by pa around z (counter-clockwise, East pa = 0)
      xx_rot = x_shft[i][j]*cos_pa + y_shft[i][j]*sin_pa;
      yy_rot = -x_shft[i][j]*sin_pa + y_shft[i][j]*cos_pa;

      // rotate by inclination around yy_rot
      yy_rot *= invcos_inc;

      // calculate radius
      rad[i][j] = sqrt(xx_rot*xx_rot + yy_rot*yy_rot);

      // calculate angle to receding major axis
      if ((xx_rot != 0.0) || (yy_rot != 0.0))
        cos_angle[i][j] = cos(atan2(yy_rot, xx_rot));
      else
        cos_angle[i][j] = 1.0;
    }
  }
}

void MyModel::add_disc_flux() {
  /*
    Add disc flux component to flux map.
  */
  const double dx = Data::get_instance().get_dx();
  const double dy = Data::get_instance().get_dy();

  double invwxd = invwxd = 1.0/wxd;
  double amp = dx*dy*Md*invwxd;

  for (size_t i=0; i<flux.size(); i++)
    for(size_t j=0; j<flux[i].size(); j++)
      flux[i][j] += amp*LookupExp::evaluate(rad[i][j]*invwxd);
}

void MyModel::add_blob_flux(std::vector< std::vector<double> >& components) {
  /*
    Calculate flux map.
  */
  const double dx = Data::get_instance().get_dx();
  const double dy = Data::get_instance().get_dy();
  const double sigma_cutoffsq = pow(Data::get_instance().get_sigma_cutoff(), 2);
  const double pixel_width = Data::get_instance().get_pixel_width();

  double sin_pa = sin(pa);
  double cos_pa = cos(pa);
  double cos_inc = cos(inc);
  double invcos_inc = 1.0/cos_inc;

  // Blob parameters
  double rc, thetac, M, wx, q, phi;
  double xc, yc;
  double wxsq;
  double qsq, invq, sqrtq;
  double invwxsq;
  double sin_phi, cos_phi;
  double amp;

  // Rotated disk coordinates
  double xd_shft, yd_shft;
  double xb_shft, yb_shft;
  double xxd_rot, yyd_rot;
  double xxb_rot, yyb_rot;
  double rsq;

  // oversampled parameters
  int si;
  double dxfs, dyfs;
  double amps;

  // Blob contribution
  for (size_t k=0; k<components.size(); ++k) {
    // Components
    rc = components[k][0];
    thetac = components[k][1];
    M = components[k][2];
    wx = components[k][3];
    q = components[k][4];
    phi = components[k][5];

    // xc, yc in disk plane
    xc = rc*cos(thetac);
    yc = rc*sin(thetac);

    /*
      Component manipulations
    */
    wxsq = wx*wx;
    qsq = q*q;
    sqrtq = sqrt(q);
    invwxsq = 1.0/(wx*wx);
    invq = 1.0/q;
    sin_phi = sin(phi);
    cos_phi = cos(phi);

    // Oversampling rate
    if (q*wx*cos_inc < 0.5*pixel_width)
      si = 2;
    else if (q*wx*cos_inc < 1.0*pixel_width)
      si = 1;
    else
      si = 0;

    dxfs = dx/(2.0*si + 1.0);
    dyfs = dy/(2.0*si + 1.0);

    // Flux normalised sum
    amp = dxfs*dyfs*M/(2.0*M_PI*wxsq*cos_inc);
    for (size_t i=0; i<image.size(); i++) {
      for (size_t j=0; j<image[i].size(); j++) {
        amps = 0.0;
        for (int is=-si; is<=si; is++) {
          for (int js=-si; js<=si; js++) {
            /*
              Get rotated/inc disk coordinates
            */
            // Shift
            xd_shft = x_shft[i][j] + js*dxfs;
            yd_shft = y_shft[i][j] + is*dyfs;

            // rotate by pa around z (counter-clockwise, East pa = 0)
            xxd_rot = xd_shft*cos_pa + yd_shft*sin_pa;
            yyd_rot = -xd_shft*sin_pa + yd_shft*cos_pa;

            // rotate by inclination around yy_rot
            yyd_rot *= invcos_inc;

            /*
              Get distance wrt centre of blob in
              rotated/inc disk coordinates.
            */
            // Shift
            xb_shft = xxd_rot - xc;
            yb_shft = yyd_rot - yc;

            // Rotate
            xxb_rot = xb_shft*cos_phi + yb_shft*sin_phi;
            yyb_rot = -xb_shft*sin_phi + yb_shft*cos_phi;

            // Calculate normalised squared distance to centre of blob
            rsq = q*pow(xxb_rot, 2) + invq*pow(yyb_rot, 2);
            rsq *= invwxsq;

            if (rsq < sigma_cutoffsq)
              amps += amp*LookupExp::evaluate(0.5*rsq);

          }
        }
        flux[i][j] += amps;
      }
    }
  }
}

void MyModel::calculate_rel_lambda() {
  /*
    Calculate rel_lambda map.
  */
 double sin_inc = sin(inc);

  for (size_t i=0; i<rel_lambda.size(); i++) {
    for (size_t j=0; j<rel_lambda[i].size(); j++) {
      // Calc relative lambda
      if (rad[i][j] == 0.0) {
        rel_lambda[i][j] = 0.0;
      } else {
        rel_lambda[i][j] = vmax*pow(1.0 + vslope/rad[i][j], vbeta);
        rel_lambda[i][j] /= pow(1.0 + pow(vslope/rad[i][j], vgamma), 1.0/vgamma);
        rel_lambda[i][j] *= sin_inc*cos_angle[i][j];
      }
      rel_lambda[i][j] += vsys;
      rel_lambda[i][j] /= constants::C;
      rel_lambda[i][j] += 1.0;
    }
  }
}

void MyModel::calculate_vdisp() {
  /*
    Calculate velocity dispersion map.
  */
  for (size_t i=0; i<vdisp.size(); i++) {
    for (size_t j=0; j<vdisp[i].size(); j++) {
      vdisp[i][j] = vdisp_param[0];
      for (int v=0; v<vdisp_order; v++)
        vdisp[i][j] += vdisp_param[v+1]*pow(rad[i][j], v+1);
      vdisp[i][j] = exp(vdisp[i][j]);
    }
  }
}

void MyModel::clear_flux_map() {
  for (size_t i=0; i<image.size(); i++)
    for (size_t j=0; j<image[i].size(); j++)
      flux[i][j] = 0.0;
}