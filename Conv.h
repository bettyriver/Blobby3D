#ifndef BLOBBY3D_CONV_H_
#define BLOBBY3D_CONV_H_

#include <vector>
#include <fftw3.h>

#include "Data.h"

/*
  Class for convolving by the PSF
*/
class Conv {
 private:
  // convolve method
  int convolve;
  
  // PSF
  std::vector<double> psf_amp;
  std::vector<double> psf_fwhm;
  double psf_beta;
  std::vector<double> psf_sigma;
  std::vector<double> psf_sigma_overdx, psf_sigma_overdy;
  int ni, nj, nr;
  double dx, dy;
  double sigma_cutoff;
  int x_pad, y_pad;

  // vectors for separable kernel
  std::vector< std::vector<double> > kernel_x;
  std::vector< std::vector<double> > kernel_y;
  std::vector< std::vector<double> > convolved_tmp_2d;

  // vector for moffat kernel
  std::vector< std::vector<double> > kernel;

  // fftw pointers/plans
  fftw_complex *out, *kernelout, *conv;
  // fftw_complex *out_all;
  // fftw_complex *kernelout_all;
  // fftw_complex *conv_all;
  fftw_plan p, q, k;
  // fftw_plan p_all, q_all, k_all;
  double* in;
  double* in2;
  // double* in_all;
  int Ni, Nj;
  int nik, njk; // size of kernel
  int midik, midjk;

  // Convolved matrix
  std::vector< std::vector< std::vector<double> > > convolved;

  static Conv instance;

 public:
  // Constructor
  Conv();

  // brute force gaussian blur
  std::vector< std::vector< std::vector<double> > > 
    brute_gaussian_blur(
      std::vector< std::vector< std::vector<double> > >& 
	preconvolved);
  
  // fftw moffat blur
  std::vector< std::vector< std::vector<double> > > 
    fftw_moffat_blur(
      std::vector< std::vector< std::vector<double> > >& 
	preconvolved);
};

#endif  // CONV_H_
