#include "Conv.h"

#include <iostream>
#include <cmath>

#include "DNest4/code/DNest4.h"
#include "Data.h"

using namespace std;
using namespace DNest4;

Conv Conv::instance;

Conv::Conv()
    :convolve(Data::get_instance().get_convolve())
    ,psf_amp(Data::get_instance().get_psf_amp())
    ,psf_fwhm(Data::get_instance().get_psf_fwhm())
    ,psf_beta(Data::get_instance().get_psf_beta())
    ,psf_sigma(Data::get_instance().get_psf_sigma())
    ,psf_sigma_overdx(Data::get_instance().get_psf_sigma_overdx())
    ,psf_sigma_overdy(Data::get_instance().get_psf_sigma_overdy())
    ,ni(Data::get_instance().get_ni())
    ,nj(Data::get_instance().get_nj())
    ,nr(Data::get_instance().get_nr())
    ,dx(Data::get_instance().get_dx())
    ,dy(Data::get_instance().get_dy())
    ,x_pad(Data::get_instance().get_x_pad())
    ,y_pad(Data::get_instance().get_y_pad()) {

  // Construct empty convolved matrix
  convolved.resize(ni - 2*y_pad);
  for (size_t i=0; i<convolved.size(); i++) {
      convolved[i].resize(nj - 2*x_pad);
      for (size_t j=0; j<convolved[i].size(); j++) {
        convolved[i][j].resize(nr);
      }
  }

  /*
    Setup convolve method
  */
  double norm;
  if (convolve == 0) {
    // Setup size of kernels due to number of gaussians
    kernel_x.resize(psf_sigma.size());
    kernel_y.resize(psf_sigma.size());

    // Setup each gaussian kernel
    double erf_min, erf_max;
    int szk_x, szk_y;
    for (size_t k=0; k<psf_sigma.size(); k++) {
      szk_x = (int)ceil(5.0*psf_sigma_overdx[k]);
      kernel_x[k].assign(2*szk_x+1, 0.0);
      for (int j=-szk_x; j<=szk_x; j++) {
        erf_min = std::erf((j - 0.5)*dx/(psf_sigma[k]*sqrt(2.0)));
        erf_max = std::erf((j + 0.5)*dx/(psf_sigma[k]*sqrt(2.0)));
        kernel_x[k][j+szk_x] = 0.5*(erf_max - erf_min);
      }

      szk_y = (int)ceil(5.0*psf_sigma_overdy[k]);
      kernel_y[k].assign(2*szk_y+1, 0.0);
      for (int i=-szk_y; i<=szk_y; i++) {
        erf_min = std::erf((i - 0.5)*dy/(psf_sigma[k]*sqrt(2.0)));
        erf_max = std::erf((i + 0.5)*dy/(psf_sigma[k]*sqrt(2.0)));
        kernel_y[k][i+szk_y] = 0.5*(erf_max - erf_min);
      }
    }

    // setup temporary convolved kernel for single separable convolution
    convolved_tmp_2d.resize(ni);
    for(size_t i=0; i<convolved_tmp_2d.size(); i++)
      convolved_tmp_2d[i].resize(nj - 2*y_pad);

  } else if (convolve == 1) {
    /*
      Setup moffat convolve by fftw
    */
    double invalphasq = 1.0/pow(psf_fwhm[0], 2);

    // Sampling
    int sample = 9;

    // Max size of kernel
    // Forced to be odd
    int max_nik = 2*(ni - 2*y_pad) + 1;
    int max_njk = 2*(nj - 2*x_pad) + 1;

    int max_midik = (max_nik - 1)/2;
    int max_midjk = (max_njk - 1)/2;

    // construct temporary moffat kernel
    std::vector< std::vector<double> > kernel_tmp;
    kernel_tmp.assign(max_nik, std::vector<double>(max_njk));

    double rsq;
    int os = (sample - 1)/2;
    double dxfs, dyfs;
    dxfs = dx/sample;
    dyfs = dy/sample;
    norm = dxfs*dyfs*invalphasq*(psf_beta - 1.0)/M_PI;

    for (int i=0; i<max_nik; i++) {
      for (int j=0; j<max_njk; j++) {
        for (int is=-os; is<=os; is++) {
          for(int js=-os; js<=os; js++) {
            rsq = pow((i - max_midik)*dy + is*dyfs, 2);
            rsq += pow((j - max_midjk)*dx + js*dxfs, 2);
            kernel_tmp[i][j] += norm*pow(1.0 + rsq*invalphasq, -psf_beta);
          }
        }
      }
    }

    /*
      Determine required size of kernel
      Sum up square kernel sum and as soon as > 0.997
      don't make kernel any larger.
      0.997 is equivalent to 3-sigma, so seems reasonable.
    */
    int szk = min((max_nik - 1)/2, (max_njk - 1)/2);
    double tl = kernel_tmp[max_midik][max_midjk];
    for (int s=1; s<=szk; s++) {
      for (int j=-s; j<s; j++) {
        tl += 4.0*kernel_tmp[max_midik - s][max_midjk + j];
	    }

      if (tl > 0.997) {
        szk = s;
        break;
      }
    }

    // Overwrite nik, njk using a smaller kernel
    nik = 2*szk + 1;
    njk = 2*szk + 1;

    midik = szk;
    midjk = szk;

    // Size of padded kernel/image
    Ni = ni + nik - 1;
    Nj = nj + njk - 1;

    // fftw arrays
    in = new double[Ni*Nj];
    in2 = new double[Ni*Nj];
    double* kernelin = new double[Ni*Nj];

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*(Nj/2+1));
    kernelout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*(Nj/2+1));
    conv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*(Nj/2+1));

    p = fftw_plan_dft_r2c_2d(Ni, Nj, in, out, FFTW_ESTIMATE);
    k = fftw_plan_dft_r2c_2d(Ni, Nj, kernelin, kernelout, FFTW_ESTIMATE);
    q = fftw_plan_dft_c2r_2d(Ni, Nj, conv, in2, FFTW_ESTIMATE);

    // Clear in arrays
    for (int i=0; i<Ni*Nj; i++) {
      in[i] = 0.0;
      in2[i] = 0.0;
      kernelin[i] = 0.0;
    }

    // Clear out arrays
    for (int i=0; i<Ni*(Nj/2+1); i++) {
      out[i][0] = 0.0; out[i][1] = 0.0;
      kernelout[i][0] = 0.0; kernelout[i][1] = 0.0;
      conv[i][0] = 0.0; conv[i][1] = 0.0;
    }

    // Zero pad kernel
    for (int i=0; i<Ni*(Nj/2+1); i++)
      kernelin[i] = 0.0;

    // add temporary kernel up to the required size
    for (int i=0; i<nik; i++) {
      for (int j=0; j<njk; j++) {
        kernelin[j + Nj*i] = kernel_tmp[(max_nik - nik)/2 + i][(max_njk - njk)/2 + j];
      }
    }

    // transform moffat kernel
    fftw_execute(k);
    fftw_destroy_plan(k);
  }
}

std::vector< std::vector< std::vector<double> > > Conv::brute_gaussian_blur(
    std::vector< std::vector< std::vector<double> > >& preconvolved) {
  /*
    Calculate convolved matrix using a decomposition of concentric Gaussians.

    The below procedure uses a separable convolution, first convolving across
    the columns, then across rows. This approach is only valid for 2D Gaussian
    PSFs.

    It performs the multiple gaussian convolution by convolving by each
    Gaussian kernel in turn. This uses the distributive property of
    convolution.

    This procedure is currently preferred for Moffat convolution compared to
    the fftw_moffat_blur function whenever the code is executed in parallel.
    Reasoning is due to the FFTW3 implementation not being thread-safe at this
    time.
  */
  const std::vector< std::vector<int> >& valid = Data::get_instance().get_valid();

  // Clear convolved matrix
  int i, j;
  for (size_t h=0; h<valid.size(); h++) {
    i = valid[h][0];
    j = valid[h][1];
    for (int r=0; r<nr; r++)
      convolved[i][j][r] = 0.0;
  }

  double norm;
  int szk_x, szk_y;
  double tl_pre, tl_con;
  for (int r=0; r<nr; r++) {
    /*
	    Convolve slice for each Gaussian kernel
    */
    for (size_t k=0; k<psf_sigma.size(); k++) {
      szk_x = (int)ceil(5.0*psf_sigma[k]/dx);
      szk_y = (int)ceil(5.0*psf_sigma[k]/dy);

      // blur across columns
      for (i=0; i<convolved_tmp_2d.size(); i++) {
        for (j=0; j<convolved_tmp_2d[i].size(); j++) {
	        convolved_tmp_2d[i][j] = 0.0;
          norm = 0.0;
          for (int p=-szk_x; p<=szk_x; p++) {
            if ((x_pad + j + p >= 0)
                && (x_pad + j + p < convolved_tmp_2d[i].size())) {
              convolved_tmp_2d[i][j] += preconvolved[i][x_pad+j+p][r]*kernel_x[k][szk_x+p];
              norm += kernel_x[k][szk_x+p];
	          }
	        }
	      }
      }

      // blur across rows for valid pixels
      for (size_t h=0; h<valid.size(); h++) {
        i = valid[h][0];
        j = valid[h][1];
        norm = 0.0;
        for (int p=-szk_y; p<=szk_y; p++) {
          if ((y_pad + i + p >= 0)
              && (y_pad + i + p < convolved_tmp_2d.size())) {
            convolved[i][j][r] += psf_amp[k]*convolved_tmp_2d[y_pad+i+p][j]*kernel_y[k][szk_y+p];
            norm += kernel_y[k][szk_y+p];
          }
        }
      }
    }
  }

  return convolved;
}

std::vector< std::vector< std::vector<double> > > Conv::fftw_moffat_blur(
    std::vector< std::vector< std::vector<double> > >& preconvolved) {
  /*
    Calculate convolved matrix using a Moffat profile.

    This method uses a non-thread safe implementation of FFTW3. As such, the
    brute_gaussian_blur method is usually preferred.
  */
  for (int r=0; r<nr; r++) {
    // put wavelength slice vector into fftw double
    for (int i=0; i<ni; i++)
      for (int j=0; j<nj; j++)
        in[j + Nj*i] = preconvolved[i][j][r];

    // transform slice
    fftw_execute(p);

    // convolve
    for (int i=0; i<Ni*(Nj/2+1); i++) {
      conv[i][0] = out[i][0]*kernelout[i][0] - out[i][1]*kernelout[i][1];
      conv[i][1] = out[i][0]*kernelout[i][1] + out[i][1]*kernelout[i][0];
    }

    // backwards transform to slice
    fftw_execute(q);

    // put renormalised convolved slice in convolved matrix
    const double invnorm = 1.0/(Ni*Nj);
    for (int i=0; i<ni-2*y_pad; i++) {
      for (int j=0; j<nj-2*x_pad; j++) {
        convolved[i][j][r] = in2[midjk + x_pad + j + Nj*(i + midik + y_pad)];
	      convolved[i][j][r] *= invnorm;
      }
    }
  }

  return convolved;
}
