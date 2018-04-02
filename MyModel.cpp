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
:objects(7 + Data::get_instance().get_model_n2lines(), 
	   Data::get_instance().get_nmax(), 
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
			    Data::get_instance().get_fluxmu_min()*Data::get_instance().get_sum_flux(),
			    Data::get_instance().get_fluxmu_max()*Data::get_instance().get_sum_flux(),
			    Data::get_instance().get_fluxmu_min()*Data::get_instance().get_sum_flux(), 
			    Data::get_instance().get_fluxmu_max()*Data::get_instance().get_sum_flux()),
	 PriorType::log_uniform)
{

  // initialise arrays
  size_t ni = Data::get_instance().get_ni();
  size_t nj = Data::get_instance().get_nj();
  size_t nr = Data::get_instance().get_nr();
  size_t x_pad = Data::get_instance().get_x_pad();
  size_t y_pad = Data::get_instance().get_y_pad();
  
  /*
  size_t nios = Data::get_instance().get_nios();
  size_t njos = Data::get_instance().get_njos();
  size_t x_pados = Data::get_instance().get_x_pados();
  size_t y_pados = Data::get_instance().get_y_pados();
  */

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
  
  /*
  // oversampled model cube
  imageos.resize(nios);
  for(size_t i=0; i<nios; ++i)
      {
	imageos[i].resize(njos);
	for(size_t j=0; j < nj; ++j)
	  {
	    imageos[i][j].resize(nr);
	  }
      }
  */

  // convolved cube
  /*
  convolved.resize(ni - 2*y_pad);
  for(size_t i=0; i<convolved.size(); ++i)
    {
      convolved[i].resize(nj - 2*x_pad);
      for(size_t j=0; j<convolved[i].size(); ++j)
	{
	  convolved[i][j].resize(nr);
	}
   }
  */

  // relative lambda for velocity
  rel_lambda.assign(ni, std::vector<double>(nj));

  // radius
  rad.assign(ni, std::vector<double>(nj));

}

void MyModel::from_prior(RNG& rng)
{ 
  // Get local variables from data
  const int model = Data::get_instance().get_model();
  const double prior_inc = Data::get_instance().get_prior_inc();
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

  sigma1_min = 1E-6;
  sigma1_width = 1E6/sigma1_min;


  // Position
  gamma_pos = 0.1*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy));

  // Velocity
  gamma_vsys = Data::get_instance().get_vsys_gamma();

  vmax_min = Data::get_instance().get_vmax_min();
  vmax_width = Data::get_instance().get_vmax_max()/vmax_min;

  vslope_min = sqrt(dx*dy);
  vslope_width = 5.0*sqrt((x_max - x_min)*(y_max - y_min))/vslope_min;

  vgamma_min = 1.0;
  vgamma_width = 50.0/vgamma_min;

  vbeta_min = -0.75;
  vbeta_width = 1.5;

  /*
    Limits: disk parameters
  */
  if(model != 0)
    {
      Md_min = 1E-4*sum_flux;
      Md_width = 10.0*sum_flux/Md_min;

      Mdn2_min = 1E-4*sum_flux;
      Mdn2_width = 10.0*sum_flux/Mdn2_min;

      wxd_min = sqrt(dx*dy);
      wxd_width = 5.0*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy))/wxd_min;

      wxd_n2_min = sqrt(dx*dy);
      wxd_n2_width = 5.0*sqrt((x_max - x_min - 2.0*x_pad_dx)*(y_max - y_min - 2.0*y_pad_dy))/wxd_n2_min;

      sigmad_lambda_min = 0.01*dr;
      sigmad_lambda_width = 0.25*(r_max - r_min)/sigmad_lambda_min;
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
  if(model != 1)
    {
      objects.from_prior(rng);
    }

  // Initialise: disk
  if((model == 1) || (model == 2))
    {
      Md = exp(log(Md_min) + log(Md_width)*rng.rand());
      Mdn2 = exp(log(Mdn2_min) + log(Mdn2_width)*rng.rand());
      wxd = exp(log(wxd_min) + log(wxd_width)*rng.rand());
      wxd_n2 = exp(log(wxd_n2_min) + log(wxd_n2_width)*rng.rand());
      sigmad_lambda = exp(log(sigmad_lambda_min) + log(sigmad_lambda_width)*rng.rand());
    }

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
  const double prior_inc = Data::get_instance().get_prior_inc();
  
  const double x_min = Data::get_instance().get_x_min();
  const double x_max = Data::get_instance().get_x_max();
  const double y_min = Data::get_instance().get_y_min();
  const double y_max = Data::get_instance().get_y_max();
  const double x_pad_dx = Data::get_instance().get_x_pad_dx();
  const double y_pad_dy = Data::get_instance().get_y_pad_dy();

  DNest4::Cauchy cauchy_xc(x_imagecentre, gamma_pos);
  DNest4::Cauchy cauchy_yc(y_imagecentre, gamma_pos);
  DNest4::Cauchy cauchy_vsys(0.0, gamma_vsys);

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
	rnd = 1.0 - 0.4*rng.rand();
      else
	rnd = rng.rand();

      if(rnd <= 0.6)
	{
	  // Perturb blob parameters
	  logH += objects.perturb(rng);
	  if((model == 0) & (objects.get_components().size() == 0))
	    return logH = -1E300;
	}
      else if(rnd <= 0.8)
	{
	  // Perturb disk parameters
	  rot_perturb = true;
	  int which = rng.rand_int(9);

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
	      vmax += log(vmax_width)*rng.randh();
	      vmax = mod(vmax - log(vmax_min), log(vmax_width)) + log(vmax_min);
	      vmax = exp(vmax);
	      break;
	    case 4:
	      vslope = log(vslope);
	      vslope += log(vslope_width)*rng.randh();
	      vslope = mod(vslope - log(vslope_min), log(vslope_width)) + log(vslope_min);
	      vslope = exp(vslope);
	      break;
	    case 5:
	      vgamma = log(vgamma);
	      vgamma += log(vgamma_width)*rng.randh();
	      vgamma = mod(vgamma - log(vgamma_min), log(vgamma_width)) + log(vgamma_min);
	      vgamma = exp(vgamma);
	      break;
	    case 6:
	      vbeta += vbeta_width*rng.randh();
	      vbeta = mod(vbeta - vbeta_min, vbeta_width) + vbeta_min;
	      break;
	    case 7:
	      pa += 2.0*M_PI*rng.randh();
	      pa = mod(pa, 2*M_PI);
	      break;
	    case 8:
	      inc += 0.5*M_PI*rng.randh();
	      inc = mod(inc, 0.5*M_PI);
	      break;
	    }
	}
      else
	{
	  // Perturb disk flux parameters
	  disk_perturb = true;
	  int which = rng.rand_int(5);
	  
	  switch(which)
	    {
	    case 0:
	      Md = log(Md);
	      Md += log(Md_width)*rng.randh();
	      Md = mod(Md - log(Md_min), log(Md_width)) + log(Md_min);
	      Md = exp(Md);
	      break;
	    case 1:
	      Mdn2 = log(Mdn2);
	      Mdn2 += log(Mdn2_width)*rng.randh();
	      Mdn2 = mod(Mdn2 - log(Mdn2_min), log(Mdn2_width)) + log(Mdn2_min);
	      Mdn2 = exp(Mdn2);
	      break;
	    case 2:
	      wxd = log(wxd);
	      wxd += log(wxd_width)*rng.randh();
	      wxd = mod(wxd - log(wxd_min), log(wxd_width)) + log(wxd_min);
	      wxd = exp(wxd);
	      break;
	    case 3:
	      wxd_n2 = log(wxd_n2);
	      wxd_n2 += log(wxd_n2_width)*rng.randh();
	      wxd_n2 = mod(wxd_n2 - log(wxd_n2_min), log(wxd_n2_width)) + log(wxd_n2_min);
	      wxd_n2 = exp(wxd_n2);
	      break;
	    case 4:
	      sigmad_lambda = log(sigmad_lambda);
	      sigmad_lambda += log(sigmad_lambda_width)*rng.randh();
	      sigmad_lambda = mod(sigmad_lambda - log(sigmad_lambda_min), log(sigmad_lambda_width));
	      sigmad_lambda += log(sigmad_lambda_min);
	      sigmad_lambda = exp(sigmad_lambda);
	      break;
	    }
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
	sigma0 += log(sigma0_width)*rng.randh();
	sigma0 = mod(sigma0 - log(sigma0_min), log(sigma0_width));
	sigma0 += log(sigma0_min);
	sigma0 = exp(sigma0);
	break;
      case 1:
	sigma1 = log(sigma1);
	sigma1 += log(sigma1_width)*rng.randh();
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
  const int model_n2lines = Data::get_instance().get_model_n2lines();
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
  // const double dxos = Data::get_instance().get_dxos();
  // const double dyos = Data::get_instance().get_dyos();
  

  // Determine if adding blobs
  bool update = false;
  vector< vector<double> > components;  
  if(model != 1)
    {
      update = objects.get_removed().size() == 0;
    }

  int i, j, r;
  if(model != 1)
    {
      if(disk_perturb || rot_perturb || !update)
	{ 
	  // get all components
	  components = objects.get_components();
	  
	  // clear image
	  for(int i=0; i<ni; i++)
	    for(int j=0; j<nj; j++)
	      for(int r=0; r<nr; r++)
		image[i][j][r] = 0.0;
	}
      else
	{
	  // only get added components
	  components = objects.get_added();
	}
    }
  else
    {
      // clear image
      for(int i=0; i<ni; i++)
	for(int j=0; j<nj; j++)
	  for(int r=0; r<nr; r++)
	    image[i][j][r] = 0.0;
    }

  // LSF
  double wlsq, invwlsq;
  double invtwo_wlsq;
  const double sigma_lsfsq = lsf_sigma*lsf_sigma;
  

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
  /*
  vector<double> lsq; // lookup array for lsq
  lsq.assign(nr, 0.);
  */

  // N2 line parameters
  double n2amp;
  double lambda_n2upp, lambda_n2low;
  int n2low_ind_max, n2upp_ind_min;
  int n2upp_ind_max = nr; // Just to avoid initialisation warning
  int n2low_ind_min = 0; // Avoid initialisation warning
  double n2upp_lsq, n2low_lsq;

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

	      // Determine velocity
	      // rel_lambda[i][j] = 2.0*vmax*atan(rad[i][j]*invvslope)*sin_inc*cos(angle)/M_PI;
	      // rel_lambda[i][j] += vsys;

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
	      
	    }
	}
    }


  /*
    Disk component
  */
  if(model != 0)
    {  

      double invwxd, invwxd_n2;

      if(disk_perturb || rot_perturb || !update)
	{
	  // line width
	  wlsq = sigmad_lambda*sigmad_lambda + sigma_lsfsq;
	  invwlsq = 1.0/wlsq;

	  // disk flux width
	  invwxd = 1.0/wxd;
	  invwxd_n2 = 1.0/wxd_n2;

	  // Flux amplitude
	  amp = db*Md*invwxd; 
	  n2amp = db*Mdn2*invwxd_n2;

	  for(i=0; i<ni; i++)
	    {
	      for(j=0; j<nj; j++)
		{
	      
		  // Carry rel lambda over to each emline
		  lambda = constants::HA*rel_lambda[i][j];
		  lambda_n2upp = constants::N2UPP*rel_lambda[i][j];
		  lambda_n2low = constants::N2LOW*rel_lambda[i][j];
		  
		  for(int r=0; r<nr; r++)
		    {

		      // Ha contribution
		      lsq = pow(wave[r] - lambda, 2)*invwlsq;
		      image[i][j][r] += amp*LookupExp::evaluate(rad[i][j]*invwxd + 0.5*lsq);
		  
		      // N2low contribution (invwlsq isn't exactly right)
		      n2low_lsq = pow(wave[r] - lambda_n2low, 2)*invwlsq;
		      image[i][j][r] += 0.333333333*n2amp*LookupExp::evaluate(rad[i][j]*invwxd_n2 + 0.5*n2low_lsq);
	  
		      // N2upp contribution (invwlsq isn't exactly right)
		      n2upp_lsq = pow(wave[r] - lambda_n2upp, 2)*invwlsq;
		      image[i][j][r] += n2amp*LookupExp::evaluate(rad[i][j]*invwxd_n2 + 0.5*n2upp_lsq);

		    }
		}
	    }
	}
    }


  // Blob parameters
  double rc, thetac, M, wx, q, phi, Mn2;
  double xc, yc;
  double xc_inc, yc_inc, xc_rot, yc_rot;
  double vdisp, sigma_lambda;
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
  double amps, n2amps;

  // calculation for cdf
  double ha_cdf_min, ha_cdf_max;

  // Integration testing
  double sum_blob; 
  // double sum_blobw;

  // Blob contribution
  if(model != 1)
    {
      for(size_t k=0; k<components.size(); ++k)
	{
	  // Components
	  rc = components[k][0]; 
	  thetac = components[k][1];
	  M = components[k][2];
	  wx = components[k][3]; 
	  q = components[k][4];
	  phi = components[k][5];
	  vdisp = components[k][6];
	  if(model_n2lines == 1)
	    Mn2 = components[k][7];
	  
	  // Testing
	  rc = 0.0;
	  wx = 2.0;
	  q = 1.0;

	  // component manipulations
	  sigma_lambda = vdisp*constants::HA/constants::C;
	  
	  // xc, yc in disk plane
	  xc = rc*cos(thetac);
	  yc = rc*sin(thetac);

	  // Variance of instrumentally broadened line
	  wlsq = sigma_lambda*sigma_lambda + sigma_lsfsq;
	  invtwo_wlsq = 1.0/sqrt(2.0*wlsq);

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
	  if(model_n2lines == 1)
	    n2amp = dxfs*dyfs*Mn2/(2.0*M_PI*wxsq*cos_inc);

	  sum_blob = 0.0;
	  for(int i=0; i<ni; i++)
	    {
	    for(int j=0; j<nj; j++)
	      {
		amps = 0.0;
		n2amps = 0.0;
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

			    if(model_n2lines == 1)
			      n2amps += n2amp*LookupExp::evaluate(0.5*rsq);
			  }
		      }
		  }

		if(amps > 0.0)
		  {
		    // sum_blobw = 0.0;
		    
		    // Calculate mean lambda for lines
		    lambda = constants::HA*rel_lambda[i][j];
		    lambda_n2upp = constants::N2UPP*rel_lambda[i][j];
		    lambda_n2low = constants::N2LOW*rel_lambda[i][j];

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
			sum_blob += image[i][j][r];

			if(model_n2lines == 1)
			  {
			    // N2low contribution (invwlsq isn't exactly right)
			    n2low_lsq = pow(wave[r] - lambda_n2low, 2)*invwlsq;
			    image[i][j][r] += 0.333333333*n2amps*LookupExp::evaluate(0.5*n2low_lsq);
	  
			    // N2upp contribution (invwlsq isn't exactly right)
			    n2upp_lsq = pow(wave[r] - lambda_n2upp, 2)*invwlsq;
			    image[i][j][r] += n2amps*LookupExp::evaluate(0.5*n2upp_lsq);
			  }
		      }
		  }
	      } 
  
	    }

	}
    }

  /* 
     Convolve Cube
  */
  const int convolve = Data::get_instance().get_convolve();
  if(convolve == 0)
    convolved = conv.brute_gaussian_blur(image);
  else if(convolve == 1)
    convolved = conv.fftw_moffat_blur(image);

  /*
    Collapse convolved cube
  */
  /*
  for(int r=0; r<nr; r++)
    {
    for(int i=0; i<ni; i++)
      {
      for(int j=0; j<nj; j++)
	{
	convolved[i][j][r] = 0.0;
	for(int si=0; si<sampling; si++)
	  {
	  for(int sj=0; sj<sampling; sj++)
	    convolved[i][j][r] += convolved[sampling*i+si][sampling*j+sj][r];
	  }
	}
      }
    }
  */
	  
 
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

  long double logL = 0.;
  double var;
  double sigma0sq = sigma0*sigma0;

  int i, j;  
  for(int h=0; h<nv; h++)
    {
      i = valid[h][0];
      j = valid[h][1];
      
      for(int r=0; r<nr; r++)
	{
	  if((data[i][j][r] != 0.) && (var_cube[i][j][r] != 0.))
	    {
	      var = var_cube[i][j][r] + sigma0sq + convolved[i][j][r]*sigma1;
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
  out<<Mdn2<<' ';
  out<<wxd<<' ';
  out<<wxd_n2<<' ';
  out<<vsys<<' ';
  out<<vmax<<' ';
  out<<vslope<<' ';
  out<<vgamma<<' ';
  out<<vbeta<<' ';
  out<<sigmad_lambda<<' ';
  out<<inc<<' ';
  out<<pa<<' ';
  out<<sigma0<<' ';
  out<<sigma1<<' ';

}

string MyModel::description() const
{
  return string("objects");
}
