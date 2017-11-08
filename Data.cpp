#include "Data.h"
#include "Constants.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

using namespace std;

Data Data::instance;

Data::Data()
{

}


void Data::load(const char* moptions_file)
{

  /*
    Model parameters:
    model:
     0 -> ONLY blobs
     1 -> ONLY disk
     2 -> disk+blobs

    nmax:
     Maximum number of blobs
    nfixed:
     Whether to fix number of blobs to nmax
    z: galaxy redshift
    
  */
  // std::string metadata_file, cube_file, var_file;
  fstream fin(moptions_file, ios::in);
  if(!fin)
    cerr<<"# ERROR: couldn't open file "<<moptions_file<<"."<<endl;

  while(fin.peek() == '#')
    fin.ignore(1000000, '\n');

  fin>>metadata_file; fin.ignore(1000000, '\n');
  fin>>cube_file; fin.ignore(1000000, '\n');
  fin>>var_file; fin.ignore(1000000, '\n');
  fin>>convolve; fin.ignore(1000000, '\n'); // Gaussian=0, Moffat(FFTW)=1
  fin>>psf_fwhm; fin.ignore(1000000, '\n');
  fin>>psf_beta; fin.ignore(1000000, '\n'); // Beta parameter if using moffat
  fin>>lsf_fwhm; fin.ignore(1000000, '\n');
  fin>>nmax; fin.ignore(1000000, '\n');
  fin>>nfixed; fin.ignore(1000000, '\n');
  fin>>vsys_gamma; fin>>vsys_max; fin.ignore(1000000, '\n');
  fin>>vmax_min; fin>>vmax_max; fin.ignore(1000000, '\n');
  fin>>fluxmu_min; fin>>fluxmu_max; fin.ignore(1000000, '\n');
  fin>>lnfluxsd_min; fin>>lnfluxsd_max; fin.ignore(1000000, '\n');
  fin>>vdispmu_min; fin>>vdispmu_max; fin.ignore(1000000, '\n');
  fin>>lnvdispsd_min; fin>>lnvdispsd_max; fin.ignore(1000000, '\n');
  fin>>qlim_min; fin.ignore(1000000, '\n');
  fin>>sigma_min; fin>>sigma_max; fin.ignore(10000000, '\n');
  fin.close();

  // Print out parameters
  std::cout << "Input Metadata file: "<< metadata_file << std::endl;
  std::cout << "Input Cube file: "<< cube_file << std::endl;
  std::cout << "Input Variance Cube file: "<< var_file << std::endl;


  // PSF convolution method message
  psf_sigma = psf_fwhm/sqrt(8.0*log(2.0));
  sigma_cutoff = 5.0;
  if(convolve == 0)
    {
      std::cout<<"Model will assume Gaussian convolution kernel.\n";
      std::cout<<"PSF FWHM (ASEC): "<<psf_fwhm<<std::endl;
    }
  else if(convolve == 1)
    {
      std::cout<<"Model will assume Moffat convolution kernel.\n";
      std::cout<<"PSF FWHM (ASEC): "<<psf_fwhm<<std::endl;
      std::cout<<"PSF BETA: "<<psf_beta<<std::endl<<std::endl;
    }
  
  // LSF convolution message
  // LSF in wavelength (needs to be redshift corrected)
  lsf_sigma = constants::HA*lsf_fwhm/(constants::C*sqrt(8.0*log(2.0)));
  std::cout<<"Model assumes a gaussian instrumental broadening.\n";
  std::cout<<"LSF FWHM (ANG): "<<constants::HA*lsf_fwhm/constants::C<<std::endl;

  // Print out remaining parameters
  std::cout<<"Nmax blobs: "<<nmax<<std::endl;
  std::cout<<"N fixed to Nmax: "<<nfixed<<std::endl;
  std::cout<<"VSYSgamma, VSYSmax: "<<vsys_gamma<<", "<<vsys_max<<std::endl;
  std::cout<<"VMAXmin, VMAXmax: "<<vmax_min<<", "<<vmax_max<<std::endl;
  std::cout<<"FLUXMUmin, FLUXMUmax: "<<fluxmu_min<<", "<<fluxmu_max<<std::endl;
  std::cout<<"LNFLUXSDmin, LNFLUXSDmax: "<<lnfluxsd_min<<", "<<lnfluxsd_max<<std::endl;
  std::cout<<"VDISPMUmin, VDISPMUmax: "<<vdispmu_min<<", "<<vdispmu_max<<std::endl;
  std::cout<<"LNVDISPSDmin, LNVDISPSDmax: "<<lnvdispsd_min<<", "<<lnvdispsd_max<<std::endl;
  std::cout<<"QLIMmin: "<<qlim_min<<std::endl;
  std::cout<<"SIGMAmin, SIGMAmax: "<<sigma_min<<", "<<sigma_max<<std::endl;

  // Model choice (only for testing)
  model = 0;
  
  // Override blob parameters for disk model
  if(model == 1)
    nmax = 0;
  
  if(model == 1)
    nfixed = true;

  // model n2 lines?
  model_n2lines = 0;


  /*
   * First, read in the metadata
   */
  fin.open(metadata_file, ios::in);
  if(!fin)
    cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
  fin >> ni >> nj;
  fin >> nr;
  fin >> x_min >> x_max >> y_min >> y_max;
  fin >> r_min >> r_max;
  fin.close();

  // Make sure maximum > minimum
  if(x_max <= x_min || y_max <= y_min)
    cerr<<"# ERROR: strange input in "<<metadata_file<<"."<<endl;

  // Compute pixel widths
  dx = (x_max - x_min)/nj;
  dy = (y_max - y_min)/ni;
  dr = (r_max - r_min)/nr;
  db = dx*dy*dr;

  // Check that pixels are square
  if(abs(log(dx/dy)) >= 1E-3)
    cerr<<"# ERROR: pixels aren't square."<<endl;

  // Array needs to be padded due to 
  // convolution and allows centre of blobs 
  // to be given a buffer
  x_pad = (int)ceil(sigma_cutoff*psf_sigma/dx - 0.5*dx/psf_sigma);
  y_pad = (int)ceil(sigma_cutoff*psf_sigma/dy - 0.5*dy/psf_sigma);
  ni += 2*y_pad;
  nj += 2*x_pad;
  x_pad_dx = x_pad*dx;
  y_pad_dy = y_pad*dy;

  x_min -= x_pad_dx;
  x_max += x_pad_dx;
  y_min -= y_pad_dy;
  y_max += y_pad_dy;

  // Comput x, y, r arrays
  compute_ray_grid();

  // Loading data comment
  std::cout<<"\nLoading data:\n";

  /*
   * Now, load the image
  */
  fin.open(cube_file, ios::in);
  if(!fin)
    cerr<<"# ERROR: couldn't open file "<<cube_file<<"."<<endl;
  image = arr_3d();
  for(size_t i=y_pad; i<image.size()-y_pad; i++)
    for(size_t j=x_pad; j<image[i].size()-x_pad; j++)
      for(size_t r=0; r<image[i][j].size(); r++)
	fin >> image[i][j][r];
  fin.close();
  std::cout<<"Image Loaded...\n";

  /*
   * Load the sigma map
   */
  fin.open(var_file, ios::in);
  if(!fin)
    cerr<<"# ERROR: couldn't open file "<<var_file<<"."<<endl;
  // Variance cube from data input
  var_cube = arr_3d();
  for(size_t i=y_pad; i<var_cube.size()-y_pad; i++)
    for(size_t j=x_pad; j<var_cube[i].size()-x_pad; j++)
      for(size_t r=0; r<var_cube[i][j].size(); r++)
	fin >> var_cube[i][j][r];
  fin.close();
  std::cout<<"Variance Loaded...\n";

  /*
   * Determine the valid data pixels
   */
  valid.assign(1, vector<int>(2));
  double tmp_im, tmp_sig;
  vector<int> tmp_vec(2);
  nv = 0;
  sum_flux = 0.0;
  for(int i=y_pad; i<ni-y_pad; i++)
    {
      for(int j=x_pad; j<nj-x_pad; j++)
	{
	  tmp_im = 0.;
	  tmp_sig = 0.;
	  for(int r=0; r<nr; r++)
	    {
	      tmp_im += image[i][j][r];
	      tmp_sig += var_cube[i][j][r];
	    }
	  
	  // Add valid pixels to array
	  if(tmp_im != 0. && tmp_sig != 0.)
	    {
	      tmp_vec[0] = i; tmp_vec[1] = j;
	      sum_flux += tmp_im;
	      if(nv == 0)
		{
		  valid[0] = tmp_vec;
		}
	      else
		{
		  valid.push_back(tmp_vec);
		}
	      nv += 1;
	    }
	}
    }
  std::cout<<"Valid pixels determined...\n\n";
}

// Create desired size array
std::vector< std::vector< std::vector<double> > > Data::arr_3d()
{
  std::vector< std::vector< std::vector<double> > > arr;

  size_t ni = Data::get_instance().get_ni();
  size_t nj = Data::get_instance().get_nj();
  size_t nr = Data::get_instance().get_nr();

  arr.resize(ni);
  for(size_t i=0; i<ni; i++)
    {
      arr[i].resize(nj);
      for(size_t j=0; j<nj; j++)
	{
	  arr[i][j].resize(nr);
	}
    }
  return arr;
}

void Data::compute_ray_grid()
{
  // Make vectors of the correct size
  x_rays.assign(ni, vector<double>(nj));
  y_rays.assign(ni, vector<double>(nj));
  r_rays.assign(nr, 0.0);

  for(size_t i=0; i<x_rays.size(); i++)
    {
      for(size_t j=0; j<x_rays[i].size(); j++)
	{
	  x_rays[i][j] = x_min + (j + 0.5)*dx;
	  y_rays[i][j] = y_max - (i + 0.5)*dy;
	}
    }
	
  for(size_t r=0; r<r_rays.size(); r++)
    r_rays[r] = r_min + (r + 0.5)*dr;
}

