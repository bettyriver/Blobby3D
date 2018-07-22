#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "LookupExp.h"
#include "LookupErf.h"
#include "Conv.h"
#include "Constants.h"
#include <cmath>

using namespace std;
using namespace DNest4;

MyModel::MyModel()
:objects(7, Data::get_instance().get_nmax(), 
	   Data::get_instance().get_nfixed(), 
	 MyConditionalPrior(Data::get_instance().get_x_min(), 
			    Data::get_instance().get_x_max(),
			    Data::get_instance().get_y_min(), 
			    Data::get_instance().get_y_max(),
			    Data::get_instance().get_r_min(), 
			    Data::get_instance().get_r_max(),
			    Data::get_instance().get_dx(),
			    Data::get_instance().get_dy(),
			    Data::get_instance().get_dr(),
			    Data::get_instance().get_x_pad_dx(),
			    Data::get_instance().get_y_pad_dy(),
			    log(Data::get_instance().get_fluxmu_min()),
			    log(Data::get_instance().get_fluxmu_max())),
	 PriorType::log_uniform)
{

  // initialise arrays
  size_t ni = Data::get_instance().get_ni();
  size_t nj = Data::get_instance().get_nj();
  size_t nr = Data::get_instance().get_nr();
  size_t x_pad = Data::get_instance().get_x_pad();
  size_t y_pad = Data::get_instance().get_y_pad();

  // model cube
  image.resize(ni);
  for(size_t i=0; i<ni; ++i)
      {
	image[i].resize(nj);
	for(size_t j=0; j<nj; ++j)
	  {
	    image[i][j].resize(nr);
	  }
      }

  // convolved cube
  convolved.resize(ni - 2*y_pad);
  for(size_t i=0; i<convolved.size(); ++i)
    {
      convolved[i].resize(nj - 2*x_pad);
      for(size_t j=0; j<convolved[i].size(); ++j)
	{
	  convolved[i][j].resize(nr);
	}
   }

  // Flux profile
  flux.assign(ni, std::vector<double>(nj));

  // relative lambda for velocity
  rel_lambda.assign(ni, std::vector<double>(nj));

  // velocity dispersion
  vdisp.assign(ni, std::vector<double>(nj));

  // radius
  rad.assign(ni, std::vector<double>(nj));

}

void MyModel::from_prior(RNG& rng)
{ 
  // Get local variables from data
  const int model = Data::get_instance().get_model();
  const double x_min = Data::get_instance().get_x_min();
  const double x_max = Data::get_instance().get_x_max();
  const double y_min = Data::get_instance().get_y_min();
  const double y_max = Data::get_instance().get_y_max();
  const double r_min = Data::get_instance().get_r_min();
  const double r_max = Data::get_instance().get_r_max();
  const double dx = Data::get_instance().get_dx();
  const double dy = Data::get_instance().get_dy();
  const double dr = Data::get_instance().get_dr();
  const double db = Data::get_instance().get_db();
  const double x_pad = Data::get_instance().get_x_pad();
  const double y_pad = Data::get_instance().get_y_pad();
  const double x_pad_dx = Data::get_instance().get_x_pad_dx();
  const double y_pad_dy = Data::get_instance().get_y_pad_dy();
  const double sum_flux = Data::get_instance().get_sum_flux();
  const double vsys_max = Data::get_instance().get_vsys_max();

  /*
    Limits: global parameters
  */

  // Error
  sigma0_min = Data::get_instance().get_sigma_min();
  sigma0_width = Data::get_instance().get_sigma_max()/sigma0_min;

  sigma1_min = 1E-4;
  sigma1_width = 1E0/sigma1_min;

  // Position
  gamma_pos = 0.1*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy));

  // Velocity
  gamma_vsys = Data::get_instance().get_vsys_gamma();

  vmax_min = Data::get_instance().get_vmax_min();
  vmax_width = Data::get_instance().get_vmax_max()/vmax_min;

  vslope_min = sqrt(dx*dy);
  vslope_width = 2.0*sqrt((x_max - x_min)*(y_max - y_min))/vslope_min;

  vgamma_min = 1.0;
  vgamma_width = 100.0/vgamma_min;

  vbeta_min = -0.75;
  vbeta_width = 1.5;

  vdisp0_min = log(1.0);
  vdisp0_width = log(200.0/1.0);

  /*
    Limits: disk parameters
  */
  if(model != 0)
    {
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
  do
    {
      xcd = cauchy_xc.generate(rng);
    }while((xcd < x_min + x_pad_dx) || 
	   (xcd > x_max - x_pad_dx));
  
  y_imagecentre = (y_max + y_min)/2.0;
  DNest4::Cauchy cauchy_yc(y_imagecentre, gamma_pos);
  do
    {
      ycd = cauchy_yc.generate(rng);
    }while((ycd < y_min + y_pad_dy) ||
	   (ycd > y_max - y_pad_dy));

  pa = 2.0*M_PI*rng.rand();
  inc = 0.5*M_PI*rng.rand();

  wxd_min = sqrt(dx*dy);
  wxd_width = 3.0*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy))/wxd_min;

  // Initialise: Velocity
  DNest4::Cauchy cauchy_vsys(0.0, gamma_vsys);
  do
    {
      vsys = cauchy_vsys.generate(rng);
    }while((vsys < -vsys_max) ||
	   (vsys >  vsys_max)); 

  vmax = exp(log(vmax_min) + log(vmax_width)*rng.rand());
  vslope = exp(log(vslope_min) + log(vslope_width)*rng.rand());
  vgamma = exp(log(vgamma_min) + log(vgamma_width)*rng.rand());
  vbeta = vbeta_min + vbeta_width*rng.rand();

  // Initialise: blobs
  if(model == 0)
    {
      do{
	objects.from_prior(rng);
      }
      while(objects.get_components().size() == 0);
    }
  else if(model == 2)
    {
      objects.from_prior(rng);
    }

  // Initialise: disk
  if((model == 1) || (model == 2))
    {
      Md = exp(log(Md_min) + log(Md_width)*rng.rand());
    }
  else
    {
      Md = 0.0;
    }

  // Dispersion
  vdisp_order = 1;
  DNest4::Gaussian gaussian_vdisp(0.0, 0.2);
  vdisp_param.assign(vdisp_order+1, 0.0);
  vdisp_param[0] = vdisp0_min + vdisp0_width*rng.rand();
  for(int v=1; v<vdisp_order+1; v++)
      vdisp_param[v] = gaussian_vdisp.generate(rng);

  wxd = exp(log(wxd_min) + log(wxd_width)*rng.rand());

  // Calculate image
  rot_perturb = true;
  if(model == 0)
    disk_perturb = false;
  else
    disk_perturb = true;

  calculate_image();

}



double MyModel::perturb(RNG& rng)
{
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
  disk_perturb = false;
  rot_perturb = false;
  
  if(rnd <= 0.9)
    {

      // skip disk perturb for model 0
      if(model == 0)
	rnd = 0.8*rng.rand();
      else if(model == 1)
	rnd = 1.0 - 0.3*rng.rand();
      else
	rnd = rng.rand();

      if(rnd <= 0.7)
	{
	  // Perturb blob parameters
	  logH += objects.perturb(rng);
	  if((model == 0) & (objects.get_components().size() == 0))
	    return logH = -1E300;
	}
      else if(rnd <= 0.8)
	{
	  // Perturb disc parameters
	  rot_perturb = true;
	  int which = rng.rand_int(12);

	  switch(which)
	    {
	    case 0:
	      logH += cauchy_xc.perturb(xcd, rng);
	      if((xcd < x_min + x_pad_dx) || 
		 (xcd > x_max - x_pad_dx))
		return logH = -1E300;
	      break;
	    case 1:
	      logH += cauchy_yc.perturb(ycd, rng);
	      if((ycd < y_min + y_pad_dy) || 
		 (ycd > y_max - y_pad_dy))
		return logH = -1E300;
	      break;
	    case 2:
	      logH += cauchy_vsys.perturb(vsys, rng);
	      if((vsys < -5.0*gamma_vsys) || 
		 (vsys >  5.0*gamma_vsys))
		return logH = -1E300;
	      break;
	    case 3:
	      vmax = log(vmax);
	      vmax += disc_step*log(vmax_width)*rng.randh();
	      vmax = mod(vmax - log(vmax_min), log(vmax_width)) + log(vmax_min);
	      vmax = exp(vmax);
	      break;
	    case 4:
	      vslope = log(vslope);
	      vslope += disc_step*log(vslope_width)*rng.randh();
	      vslope = mod(vslope - log(vslope_min), log(vslope_width)) + log(vslope_min);
	      vslope = exp(vslope);
	      break;
	    case 5:
	      vgamma = log(vgamma);
	      vgamma += disc_step*log(vgamma_width)*rng.randh();
	      vgamma = mod(vgamma - log(vgamma_min), log(vgamma_width)) + log(vgamma_min);
	      vgamma = exp(vgamma);
	      break;
	    case 6:
	      vbeta += disc_step*vbeta_width*rng.randh();
	      vbeta = mod(vbeta - vbeta_min, vbeta_width) + vbeta_min;
	      break;
	    case 7:
	      pa += disc_step*2.0*M_PI*rng.randh();
	      pa = mod(pa, 2*M_PI);
	      break;
	    case 8:
	      inc += disc_step*0.5*M_PI*rng.randh();
	      inc = mod(inc, 0.5*M_PI);
	      break;
	    case 9:
	      wxd = log(wxd);
	      wxd += disc_step*log(wxd_width)*rng.randh();
	      wxd = mod(wxd - log(wxd_min), log(wxd_width)) + log(wxd_min);
	      wxd = exp(wxd);
	      break;
	    case 10:
	      vdisp_param[0] += disc_step*vdisp0_width*rng.randh();
	      vdisp_param[0] = mod(vdisp_param[0] - vdisp0_min, vdisp0_width);
	      vdisp_param[0] += vdisp0_min;
	      break;
	    case 11:
	      which = rng.rand_int(vdisp_order);
	      logH += gaussian_vdisp.perturb(vdisp_param[which+1], rng);
	      break;
	    }
	}
      else
	{
	  // Perturb disc flux parameters
	  disk_perturb = true;
	  
	  Md = log(Md);
	  Md += disc_step*log(Md_width)*rng.randh();
	  Md = mod(Md - log(Md_min), log(Md_width)) + log(Md_min);
	  Md = exp(Md);
	}
      
      // Pre-rejection trick
      if(rng.rand() >= exp(logH))
	return -1E300;
      else
	logH = 0.0;

      calculate_image();
      
    }
  else
    {
      int which = rng.rand_int(2);
      switch(which)
      {
      case 0:
	sigma0 = log(sigma0);
	sigma0 += sigma_step*log(sigma0_width)*rng.randh();
	sigma0 = mod(sigma0 - log(sigma0_min), log(sigma0_width));
	sigma0 += log(sigma0_min);
	sigma0 = exp(sigma0);
	break;
      case 1:
	sigma1 = log(sigma1);
	sigma1 += sigma_step*log(sigma1_width)*rng.randh();
	sigma1 = mod(sigma1 - log(sigma1_min), log(sigma1_width));
	sigma1 += log(sigma1_min);
	sigma1 = exp(sigma1);
	break;
      }
    }

  return logH;

}

void MyModel::calculate_image()
{

  // Get data
  const int model = Data::get_instance().get_model();
  const vector< vector<double> >& x = Data::get_instance().get_x_rays();
  const vector< vector<double> >& y = Data::get_instance().get_y_rays();
  const vector<double>& wave = Data::get_instance().get_r_rays();
  const double dx = Data::get_instance().get_dx();
  const double dy = Data::get_instance().get_dy();
  const double db = Data::get_instance().get_db();
  const double dr = Data::get_instance().get_dr();
  const int ni = Data::get_instance().get_ni();
  const int nj = Data::get_instance().get_nj();
  const int nr  = Data::get_instance().get_nr();
  const int nv = Data::get_instance().get_nv();
  const double x_min = Data::get_instance().get_x_min();
  const double x_max = Data::get_instance().get_x_max();
  const double y_min = Data::get_instance().get_y_min();
  const double y_max = Data::get_instance().get_y_max();
  const double r_min = Data::get_instance().get_r_min();
  const double r_max = Data::get_instance().get_r_max();
  const double sigma_cutoff = Data::get_instance().get_sigma_cutoff();
  const double lsf_sigma = Data::get_instance().get_lsf_sigma();
  const int convolve = Data::get_instance().get_convolve();
  // const double dxos = Data::get_instance().get_dxos();
  // const double dyos = Data::get_instance().get_dyos();
  
  const double x_pad_dx = Data::get_instance().get_x_pad_dx();
  const double y_pad_dy = Data::get_instance().get_y_pad_dy();
  const double image_width = sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy));


  const double sigma_lsfsq = lsf_sigma*lsf_sigma;

  // Determine if adding blobs
  bool update = false;
  vector< vector<double> > components;  
  if(model != 1)
    {
      update = objects.get_removed().size() == 0;
    }

  /*
    CLEAR CUBE
  */
  if(model != 1)
    {
      if(disk_perturb || rot_perturb || !update)
	{ 
	  // get all components
	  components = objects.get_components();
	  
	  // clear flux map
	  for(int i=0; i<ni; i++)
	    for(int j=0; j<nj; j++)
	      flux[i][j] = 0.0;

	  // clear cube
	  for(int i=0; i<ni; i++)
	    for(int j=0; j<nj; j++)
	      for(int r=0; r<nr; r++)
		image[i][j][r] = 0.0;
	}
      else
	{
	  // clear cube
	  for(int i=0; i<ni; i++)
	    for(int j=0; j<nj; j++)
	      for(int r=0; r<nr; r++)
		image[i][j][r] = 0.0;	  
	  
	  // only get added components
	  components = objects.get_added();
	}
    }
  else
    {
      // clear flux map
      for(int i=0; i<ni; i++)
	for(int j=0; j<nj; j++)
	  flux[i][j] = 0.0;      
      
      // clear cube
      for(int i=0; i<ni; i++)
	for(int j=0; j<nj; j++)
	  for(int r=0; r<nr; r++)
	    image[i][j][r] = 0.0;
    }

  // LSF
  double wlsq, invwlsq;
  double invtwo_wlsq;

  /*
    Calculate exponential disk component
    
    Also calculate lambda array.
    (wavelength corresponding to mean velocity)
  */

  // Disk parameters
  double sin_pa, cos_pa;
  double sin_inc, cos_inc;
  double invcos_inc;
  double invvslope;
  double x_shft, y_shft; 
  double xx_rot, yy_rot, yy_inc;
  double rsq;
  double amp;
  double angle;
  double lambda;
  double lsq;

  // Calculate sin/cosine outside loops
  sin_pa = sin(pa); cos_pa = cos(pa);
  sin_inc = sin(inc); cos_inc = cos(inc);
  invcos_inc = 1.0/cos_inc;

  // Calculate rotational array
  if(rot_perturb)
    {
      invvslope = 1.0/vslope;
      for(int i=0; i<ni; i++)
	{
	  for(int j=0; j<nj; j++)
	    {
	      // Shift
	      x_shft = x[i][j] - xcd;
	      y_shft = y[i][j] - ycd;

	      // rotate by pa around z (counter-clockwise, East pa = 0)
	      xx_rot =  x_shft*cos_pa + y_shft*sin_pa;
	      yy_rot = -x_shft*sin_pa + y_shft*cos_pa;
	  
	      // rotate by inclination around yy_rot
	      yy_rot *= invcos_inc;

	      // calculate radius
	      rad[i][j] = sqrt(xx_rot*xx_rot + yy_rot*yy_rot);

	      // calculate angle to receding major axis
	      if((xx_rot == 0.0) & (yy_rot == 0.0))
		{ angle = 0.0; }
	      else
		{ angle = atan2(yy_rot, xx_rot); }
	      
	      // Calc relative lambda
	      if(rad[i][j] == 0.0)
		{
		  rel_lambda[i][j] = 0.0;
		}
	      else
		{
		  rel_lambda[i][j] = vmax*pow(1.0 + vslope/rad[i][j], vbeta);
		  rel_lambda[i][j] /= pow(1.0 + pow(vslope/rad[i][j], vgamma), 1.0/vgamma);
		  rel_lambda[i][j] *= sin_inc*cos(angle);
		}
	      rel_lambda[i][j] += vsys;
	      rel_lambda[i][j] /= constants::C;
	      rel_lambda[i][j] += 1.0;

	      // Calc dispersion (assuming constant atm)
	      vdisp[i][j] = vdisp_param[0];
	      for(int v=0; v<vdisp_order; v++)
		vdisp[i][j] += vdisp_param[v+1]*pow(rad[i][j], v+1);
	      vdisp[i][j] = exp(vdisp[i][j]);
	    }
	}
    }


  /*
    Disk component
  */
  double sigma_lambda;
  if(model != 0)
    {  
      double invwxd;
      if(disk_perturb || rot_perturb || !update)
	{

	  // disk flux width
	  invwxd = 1.0/wxd;

	  // Flux amplitude
	  amp = dx*dy*Md*invwxd;
	  for(int i=0; i<ni; i++)
	    {
	      for(int j=0; j<nj; j++)
		{
	      
		  flux[i][j] += amp*LookupExp::evaluate(rad[i][j]*invwxd);

		  /*
		  // Carry rel lambda over to each emline
		  lambda = constants::HA*rel_lambda[i][j];

		  // Calculate line width
		  sigma_lambda = vdisp[i][j]*constants::HA/constants::C;
		  wlsq = sigma_lambda*sigma_lambda + sigma_lsfsq;
		  invwlsq = 1.0/wlsq;

		  // Calculate disk contribution to spaxels
		  for(int r=0; r<nr; r++)
		    {

		      // Ha contribution (change this to using erf for gauss in future)
		      lsq = pow(wave[r] - lambda, 2)*invwlsq;
		      image[i][j][r] += amp*LookupExp::evaluate(rad[i][j]*invwxd + 0.5*lsq);

		      // Testing just on flux profile
		      // image[i][j][r] += amp*LookupExp::evaluate(rad[i][j]*invwxd)/30.0;

		    }
		  */
		}
	    }
	}
    }


  // Blob parameters
  double rc, thetac, M, wx, q, phi;
  double xc, yc;
  double xc_inc, yc_inc, xc_rot, yc_rot;
  double wy, wxsq, wysq, corr;
  double qsq, invq, sqrtq;
  double invwxsq;
  double sin_phi, cos_phi;
  double sin_theta, cos_theta;
  double det, invdet;
  int r_ind_min, r_ind_max;
  int i_min, i_max;
  int j_min, j_max;
  const double two_pi = pow(2.0*M_PI, -1.5); 

  // Rotated disk coordinates
  double xd_shft, yd_shft;
  double xxd_rot, yyd_rot;

  // Don't calculate when > sigma_cutoff
  const double sigma_cutoffsq = sigma_cutoff*sigma_cutoff;
  double cutoff_width;

  // oversample for flux
  int si;
  double dxfs, dyfs;
  double amps;

  // calculation for cdf
  double ha_cdf_min, ha_cdf_max;

  // Integration testing
  // double sum_blob; 
  // double sum_blobw;

  // Blob contribution
  if(model != 1)
    {
      for(size_t k=0; k<components.size(); ++k)
	{
	  // Components
	  rc = components[k][0]*wxd; 
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
	  if(q*wx*cos_inc < 0.5*dx)
	    { si = 4; }
	  else if(q*wx*cos_inc < 1.0*dx)
	    { si = 2; }
	  else if(q*wx*cos_inc < 1.5*dx)
	    { si = 1; }
	  else
	    { si = 0; }

	  dxfs = dx/(2.0*si + 1.0);
	  dyfs = dy/(2.0*si + 1.0);

	  // Flux normalised sum
	  amp = dxfs*dyfs*M/(2.0*M_PI*wxsq*cos_inc);

	  // sum_blob = 0.0;
	  for(int i=0; i<ni; i++)
	    {
	    for(int j=0; j<nj; j++)
	      {
		amps = 0.0;
		for(int is=-si; is<=si; is++)
		  {
		    for(int js=-si; js<=si; js++)
		      {

			/*
			  Get rotated/inc disk coordinates
			*/
			// Shift
			xd_shft = x[i][j] + js*dxfs - xcd;
			yd_shft = y[i][j] + is*dyfs - ycd;

			// rotate by pa around z (counter-clockwise, East pa = 0)
			xxd_rot =  xd_shft*cos_pa + yd_shft*sin_pa;
			yyd_rot = -xd_shft*sin_pa + yd_shft*cos_pa;
	  
			// rotate by inclination around yy_rot
			yyd_rot *= invcos_inc;


			/*
			  Get distance wrt centre of blob 
			  in rotated/inc disk coord
			 */
			// Shift
			x_shft = xxd_rot - xc;
			y_shft = yyd_rot - yc;
			
			// Rotate
			xx_rot =  x_shft*cos_phi + y_shft*sin_phi;
			yy_rot = -x_shft*sin_phi + y_shft*cos_phi;
			      
			// Calculate normalised squared distance to centre of blob
			rsq = q*xx_rot*xx_rot + invq*yy_rot*yy_rot;
			rsq *= invwxsq;
			
			if(rsq < sigma_cutoffsq)
			  {
			    amps += amp*LookupExp::evaluate(0.5*rsq);
			  }
		      }
		  }

		flux[i][j] += amps;

		/*
		if(amps > 0.0)
		  {
		    // sum_blobw = 0.0;
		    
		    // Calculate mean lambda for lines
		    lambda = constants::HA*rel_lambda[i][j];

		    // Calculate line width
		    sigma_lambda = vdisp[i][j]*constants::HA/constants::C;
		    wlsq = sigma_lambda*sigma_lambda + sigma_lsfsq;
		    invtwo_wlsq = 1.0/sqrt(2.0*wlsq);

		    for (int r=0; r<nr; r++)
		      {
			// Ha contribution
			if(r == 0)
			  {
			    ha_cdf_min = LookupErf::evaluate((wave[r] - 0.5*dr - lambda)*invtwo_wlsq);
			  }
			else
			  {
			    ha_cdf_min = ha_cdf_max;
			  }
			
			ha_cdf_max = LookupErf::evaluate((wave[r] + 0.5*dr - lambda)*invtwo_wlsq);
		
			image[i][j][r] += 0.5*amps*(ha_cdf_max - ha_cdf_min);
			// sum_blob += image[i][j][r];
		      }
		  }
		*/
	      }
  
	    }

	}
    }

  
  /*
    Create cube
  */
  for(int i=0; i<ni; i++)
    {
      for(int j=0; j<nj; j++)
	{
	  if(flux[i][j] > 0.0)
	    {
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
	      for(int r=1; r<nr; r++)
		{
		  ha_cdf_min = ha_cdf_max;
		  ha_cdf_max = LookupErf::evaluate((wave[r] + 0.5*dr - lambda)*invtwo_wlsq);
		  image[i][j][r] = 0.5*flux[i][j]*(ha_cdf_max - ha_cdf_min);
		}
	    }
	}
    }



  /* 
     Convolve Cube
  */
  if(convolve == 0)
    convolved = conv.brute_gaussian_blur(image);
  else if(convolve == 1)
    convolved = conv.fftw_moffat_blur(image);
 
}


double MyModel::log_likelihood() const
{
  const vector< vector< vector<double> > >& data = Data::get_instance().get_image();
  const vector< vector< vector<double> > >& var_cube = Data::get_instance().get_var_cube();
  
  const int model = Data::get_instance().get_model();
  const vector< vector<int> >& valid = Data::get_instance().get_valid();
  const int nr = Data::get_instance().get_nr();
  const int nv = Data::get_instance().get_nv();
  const int x_pad = Data::get_instance().get_x_pad();
  const int y_pad = Data::get_instance().get_y_pad();
  
  // If no blobs return prob = 0
  if((model == 0) && (objects.get_components().size() == 0)) 
    return -1E300;

  long double logL;
  double var;
  double sigma0sq;

  logL = 0.0;
  sigma0sq = sigma0*sigma0;

  int i, j;  
  for(int h=0; h<nv; h++)
    {
      i = valid[h][0];
      j = valid[h][1];
      
      for(int r=0; r<nr; r++)
	{
	  if((data[i][j][r] != 0.0) && (var_cube[i][j][r] != 0.0))
	    {
	      var = var_cube[i][j][r] + sigma0sq + convolved[i][j][r]*sigma1;
	      // var = var_cube[i][j][r]; // Testing convolution on likelihood first
	      
	      logL += -0.5*log(2.0*M_PI*var)
		-0.5*pow(data[i][j][r] - convolved[i][j][r], 2)/var;
	    }
	}
    }

  return logL;
}

void MyModel::print(std::ostream& out) const
{
  const int x_pad = Data::get_instance().get_x_pad();
  const int y_pad = Data::get_instance().get_y_pad();

  out<<setprecision(6);

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
  for(int i=0; i<vdisp_order+1; i++)
    out<<vdisp_param[i]<<' ';
  out<<inc<<' ';
  out<<pa<<' ';
  out<<sigma0<<' ';
  out<<sigma1<<' ';

}

string MyModel::description() const
{
  return string("objects");
}
