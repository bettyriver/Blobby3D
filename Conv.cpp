#include <iostream>
#include <cmath>
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "Conv.h"

using namespace std;
using namespace DNest4;


Conv Conv::instance;

Conv::Conv()
  :convolve(Data::get_instance().get_convolve())
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
  ,y_pad(Data::get_instance().get_y_pad())
{

  // Construct empty convolved matrix
  convolved.resize(ni - 2*y_pad);
  for(size_t i=0; i<convolved.size(); i++)
    {
      convolved[i].resize(nj - 2*x_pad);
      for(size_t j=0; j<convolved[i].size(); j++)
	{
	  convolved[i][j].resize(nr);
	}
    }

  /*
    Setup convolve method
  */
  double norm;
  if(convolve == 0)
    {
      /*
	Setup gaussian convolution
      */
      double dist;
      double erf_min, erf_max;

      int szk_x = (int)ceil(5.0*psf_sigma_overdx);
      kernel_x.assign(2*szk_x+1, 0.0);
      for(int j=-szk_x; j<=szk_x; j++)
	{
	  erf_min = std::erf((j - 0.5)*dx/(psf_sigma*sqrt(2.0)));
	  erf_max = std::erf((j + 0.5)*dx/(psf_sigma*sqrt(2.0)));
	  kernel_x[j+szk_x] = 0.5*(erf_max - erf_min);
	}
      
      int szk_y = (int)ceil(5.0*psf_sigma_overdy);
      kernel_y.assign(2*szk_y+1, 0.0);
      for(int i=-szk_y; i<=szk_y; i++)
	{
	  erf_min = std::erf((i - 0.5)*dy/(psf_sigma*sqrt(2.0)));
	  erf_max = std::erf((i + 0.5)*dy/(psf_sigma*sqrt(2.0)));
	  kernel_y[i+szk_y] = 0.5*(erf_max - erf_min);
	}
      
      // setup temporary convolved kernel for single separable convolution
      convolved_tmp_2d.resize(ni);
      for(size_t i=0; i<convolved_tmp_2d.size(); i++)
	convolved_tmp_2d[i].resize(nj - 2*y_pad);

    }
  else if(convolve == 1)
    {
      /*
	Setup moffat convolve by fftw
      */
      double invalphasq = 1.0/pow(psf_fwhm, 2);

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
      
      for(int i=0; i<max_nik; i++)
	{
	  for(int j=0; j<max_njk; j++)
	    {
	      for(int is=-os; is<=os; is++)
		{
		  for(int js=-os; js<=os; js++)
		    {
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
      for(int s=1; s<=szk; s++)
	{
	  for(int j=-s; j<s; j++)
	    {
	      tl += 4.0*kernel_tmp[max_midik - s][max_midjk + j];
	    }
	  
	  if(tl > 0.997)
	    {
	      szk = s;
	      break;
	    }
	}
      std::cout<<szk<<std::endl;


      // overwrite nik, njk using a smaller kernel
      nik = 2*szk + 1;
      njk = 2*szk + 1;
      
      midik = szk;
      midjk = szk;

      // Size of padded kernel/image
      Ni = ni + nik - 1;
      Nj = nj + njk - 1;

      // variables for FFTW advanced interface
      int rank = 2;
      int Nin[] = {Ni, Nj};
      int Nout[] = {Ni, Nj/2 + 1};
      int howmany = 1;
      int idist = 0;
      int odist = 0;
      int istride = 1;
      int ostride = 1;
      int *inembed = Nin;
      int *onembed = Nout;

      // fftw arrays
      in2 = new double[Ni*Nj];
      double* kernelin = new double[Ni*Nj];
      in = new double[Ni*Nj];

      // in_all = new double[Ni*Nj*nr];
      // double* kernelin_all = new double[Ni*Nj*nr];

      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*(Nj/2+1));
      kernelout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*(Nj/2+1));
      conv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*(Nj/2+1));

      // out_all = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*Nj*(nr/2+1));
      // kernelout_all = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Ni*Nj*(nr/2+1));

      // create plans
      /*
      p = fftw_plan_many_dft_r2c(rank, Nin, howmany,
				 in, inembed, istride, idist,
				 out, onembed, ostride, odist, 
				 FFTW_PATIENT);
      */
      p = fftw_plan_dft_r2c_2d(Ni, Nj, in, out, FFTW_PATIENT);
      
      /*
      p_all = fftw_plan_many_dft_r2c(rank, Nin, nr,
				     in_all, inembed, nr, 1,
				     out_all, onembed, nr, 1,
				     FFTW_PATIENT);
      */

      k = fftw_plan_dft_r2c_2d(Ni, Nj, kernelin, kernelout, FFTW_ESTIMATE);
      
      /*
      k_all = fftw_plan_many_dft_r2c(rank, Nin, nr, 
				     kernelin_all, inembed, nr, 1,
				     kernelout_all, onembed, nr, 1,
				     FFTW_ESTIMATE);
      */
      
      /*
      q = fftw_plan_many_dft_c2r(rank, Nout, howmany,
				 conv, onembed, ostride, odist,
				 in2, inembed, istride, idist,
				 FFTW_PATIENT);
      */
      
      q = fftw_plan_dft_c2r_2d(Ni, Nj, conv, in2, FFTW_PATIENT);
      
      /*
	Construct padded kernel
      */

      // zero pad
      for(int i=0; i<Ni*(Nj/2+1); i++)
	kernelin[i] = 0.0;
      
      // add temporary kernel up to the required size
      for(int i=0; i<nik; i++)
	{
	  for(int j=0; j<njk; j++)
	    {
	    kernelin[j + Nj*i] = kernel_tmp[(max_nik - nik)/2 + i][(max_njk - njk)/2 + j];
	    }
	}
      
      /*
      // kernel all
      std::cout<<"start test\n";
      for(int i=0; i<Ni*Nj*nr; i++)
	kernelin_all[i] = 0.0;      
      std::cout<<"kernel all zeroed\n";

      for(int i=0; i<nik; i++)
	for(int j=0; j<njk; j++)
	  for(int r=0; r<nr; r++)
	    kernelin_all[r + nr*j + i*nr*Nj] = kernelin[j + Nj*i];
      std::cout<<"kernelin_all constructed\n";
      */
      
      // transform moffat kernel
      fftw_execute(k);
      fftw_destroy_plan(k);

      /*
      fftw_execute(k_all);
      fftw_destroy_plan(k_all);
      */
      
      // std::cout<<"kernel transformed\n";

    }

}


std::vector< std::vector< std::vector<double> > > 
Conv::brute_gaussian_blur(std::vector< std::vector< std::vector<double> > >& 
			  preconvolved)
{
  const std::vector< std::vector<int> >& valid = Data::get_instance().get_valid();
  const double szk_x = (int)ceil(5.0*psf_sigma/dx);
  const double szk_y = (int)ceil(5.0*psf_sigma/dy);

  /*
    Construct convolved_tmp 2d matrix.
    Used for blurring across columns.
  */

  /*
    Calculate convolved matrix
  
    The below procedure uses a separable convolution,
    first convolving across the columns, then across rows.

    Only valid for 2d gaussian PSFs.
  */
  int i, j;
  double norm;
  for(int r=0; r<nr; r++)
    {

      // blur across columns
      for(i=0; i<convolved_tmp_2d.size(); i++)
	{
	  for(j=0; j<convolved_tmp_2d[i].size(); j++)
	    {
	      convolved_tmp_2d[i][j] = 0.0;
	      norm = 0.0;
	      for(int p=-szk_x; p<=szk_x; p++)
		{
		  if((x_pad + j + p >= 0) & (x_pad + j + p < preconvolved[i].size()))
		    {
		     convolved_tmp_2d[i][j] += preconvolved[i][x_pad+j+p][r]*kernel_x[szk_x+p];
		     norm += kernel_x[szk_x+p];
		    }
		}
	       // renormalise
	       convolved_tmp_2d[i][j] /= norm;
	    }
	}

      // blur across rows for valid pixels
      for(size_t h=0; h<valid.size(); h++)
	{
	  i = valid[h][0];
	  j = valid[h][1];

	  convolved[i][j][r] = 0.0;
	  norm = 0.0;
	  for(int p=-szk_y; p<=szk_y; p++)
	    {
	      if((y_pad + i + p >= 0) && (y_pad + i + p < convolved_tmp_2d.size()))
		{
		  convolved[i][j][r] += convolved_tmp_2d[y_pad+i+p][j]*kernel_y[szk_y+p];
		  norm += kernel_y[szk_y+p];
		}
	    }
	  // renormalise
	  convolved[i][j][r] /= norm;
	}



      /*
      for(size_t i=0; i<nios; i++)
	{
	for(size_t j=0; j<njos; j++)
	  {
	    convolved[i][j][r] = 0.0;
	    for(int p=0; p<2*x_pados+1; p++)
	      convolved[i][j][r] += convolved_tmp_2d[i+p][j]*kernel_y[p];
	  }
	}
      */
	  



    }

  return convolved;
}


std::vector< std::vector< std::vector<double> > > 
Conv::fftw_moffat_blur(std::vector< std::vector< std::vector<double> > >& 
			  preconvolved)
{

  /*
    Convolving cube at same time
  */
  
  /*
  // put model into fftw double
  std::cout<<"in all"<<std::endl;
  for(int i=0; i<ni; i++)
    {
    for(int j=0; j<nj; j++)
      {
      for(int r=0; r<nr; r++)
	{
	  std::cout<<r<<" "<<j<<" "<<i<<std::endl;
	  std::cout<<"Sz: "<<nr<<" "<<nj<<" "<<ni<<std::endl;
	  std::cout<<"I: "<<r + j*nr + i*nr*Nj<<" "<<Ni*Nj*nr<<std::endl;
	  std::cout<<"PRES: "<<preconvolved.size()<<" "
		   <<preconvolved[0].size()<<" "
		   <<preconvolved[0][0].size()<<std::endl;
	  in_all[r + j*nr + i*nr*Nj] = preconvolved[i][j][r];
	}
      }
    }

  // transform model
  fftw_execute(p_all);

  // perform convolution
  std::cout<<"Conv_all\n";
  for(int i=0; i<Ni*Nj*(nr/2+1); i++)
	{
	  conv_all[i][0] = out_all[i][0]*kernelout_all[i][0] - out_all[i][1]*kernelout_all[i][1];
	  conv_all[i][1] = out_all[i][0]*kernelout_all[i][1] + out_all[i][1]*kernelout_all[i][0];
	}


  // inverse transform

  // put back into a vector

  */

  for(int r=0; r<nr; r++)
    {
      
      // put wavelength slice vector into fftw double
      // double tls = 0.;
      for(int i=0; i<ni; i++)
	{
	  for(int j=0; j<nj; j++)
	    {
	      in[j + Nj*i] = preconvolved[i][j][r];
	      // tls += in[j + Nj*i];
	    }
	}
      // std::cout<<"TLS: "<<tls<<std::endl;

      // transform slice
      fftw_execute(p);   


      // convolve
      for(int i=0; i<Ni*(Nj/2+1); i++)
	{
	  conv[i][0] = out[i][0]*kernelout[i][0] - out[i][1]*kernelout[i][1];
	  conv[i][1] = out[i][0]*kernelout[i][1] + out[i][1]*kernelout[i][0];
	  /*
	  if(conv[i][0] != conv[i][0])
	    {
	      std::cout<<"REAL: "<<conv[i][0]<<" "<<i<<std::endl;
	      std::cout<<"INPUTS: "
		       <<in[i]<<" "
		       <<out[i][0]<<" "<<out[i][1]
		       <<kernelout[i][0]<<" "<<kernelout[i][1]<<std::endl;;
	    }
	  if(conv[i][1] != conv[i][1])
	    std::cout<<"IMAG: "<<conv[i][1]<<" "<<i<<std::endl;
	  */
	}

      // backwards transform to slice
      fftw_execute(q);

      // put renormalised convolved slice in convolved matrix
      const double invnorm = 1.0/(Ni*Nj);
      for(int i=0; i<ni-2*y_pad; i++)
	{
	for(int j=0; j<nj-2*x_pad; j++)
	  {
	    convolved[i][j][r] = in2[midjk + x_pad + j + Nj*(i + midik + y_pad)];
	    convolved[i][j][r] *= invnorm;
	  }
	}
    }

  return convolved;

}

