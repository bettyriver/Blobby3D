# Blobby3d

Blobby3D is designed to perform bayesian inference for gas kinematics on the H-alpha/N2 complex for IFU galaxy observations. It is robust even for galaxies that have clumpy H2 regions.

## Getting Started

These instructions will get you a copy of the project up and running on your machine for testing purposes. If you have any trouble, feel free to contact Mathew Varidel at mathew.varidel@sydney.edu.au.

### Prerequisites

Blobby3D is dependent on the following:

C++11
[DNest4](https://github.com/eggplantbren/DNest4)
[FFTW3](http://www.fftw.org)
[Python 2.7](https://www.python.org) with [Cython](http://cython.org) package (for postprocessing)

[FFTW3](http://www.fftw.org) and [Boost](http://www.boost.org) will typically be available on scientific machines/clusters, so you will often either have them or can get them relatively easily by contacting your system administrator.

[DNest4](https://github.com/eggplantbreen/DNest4) is a bespoke nested sampling algorithm written by Brendon Brewer, and will need to be cloned from github in the usual manner. It is only available for Unix like machines.

### Installing

Once you have downloaded/installed [FFTW3](www.fftw.org) and [Boost](www.boost.org), I would suggest setting up DNest4 and Blobby3D within the same directory. So, go to the root directory where you would like to install the packages, and create a directory:

mkdir Blobby3D
cd Blobby3D

#### Installing DNest4

Then install DNest4 using the available instructions:

git clone https://github.com/eggplantbren/DNest4
cd DNest4/code
make
cd ../python
python setup.py install

#### Installing Blobby3D

Then download/install Blobby3D:

cd ../..
git clone https://github.com/SpaceOdyssey/Blobby3D
cd Blobby3D
make

### Running Blobby3D on Examples

To get a feel for Blobby3D you will want to start with an example. You can find a few examples in the Blobby3D/Examples folder. There are 3 example toy models provided. The toy models are constructed in accordance with the paper. The parent folder describes the PSF FWHM (psfa), PSF Moffat beta parameter (psfb) and the signal-noise (sn). The subfolder describes the velocity dispersion (vdisp) and the v(R_opt) (vropt) value.

Navigate to one of the subfolders. The code can be run as an executable from within the subdirectory as follows:

../../Blobby3D -t 1 -f MODEL_OPTIONS

You may want to setup a terminal alias such that you can run the executable from anywhere on your machine. The -t 1 parameter sets the number of threads -- I recommend using multiple threads. The -f MODEL_OPTIONS tells the model to use the input hyperparameters/parameters required to defined the priors. Note the the DNest4 parameters are defined in the OPTIONS file.

There will be 3 output files created as the executable runs:

sample_info.txt
sample.txt
sampler_state.txt

[Run posterior scripts...]

### Running Blobby3D on Your Own Data

To run Blobby3D on your own data you need to setup several files.
