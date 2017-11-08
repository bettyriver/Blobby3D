# Blobby3d

Blobby3D is designed to perform bayesian inference for gas kinematics on the H-alpha/N2 complex for IFU galaxy observations. It is robust even for galaxies that have clumpy H2 regions.

## Getting Started

These instructions will get you a copy of the project up and running on your machine for testing purposes. If you have any trouble, feel free to contact Mathew Varidel at mathew.varidel@sydney.edu.au. 

### Prerequisites

Blobby3D is dependent on the following:

[DNest4](https://github.com/eggplantbren/DNest4)  
[FFTW3](https://www.fftw.org)  
[Boost C++ Library](www.boost.org)  
Python with Cython package  

[FFTW3](www.fftw.org) and [Boost](www.boost.org) will typically be available on scientific machines/clusters, so you will often either have them or can get them relatively easily by contacting your system administrator. 

[DNest4](https://github.com/eggplantbreen/DNest4) is a bespoke nested sampling algorithm written by Brendon Brewer, and will need to be cloned from github in the usual manner. It is only available for Unix like machines.

### Installing

Once you have downloaded/installed [FFTW3](www.fftw.org) and [Boost](www.boost.org), I would suggest setting up DNest4 and Blobby3D within the same directory. So, go to the root directory where you would like to install the packages, and create a directory

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

Then install Blobby3D:

cd ../.. # return to root directory
git clone https://github.com/SpaceOdyssey/Blobby3D
cd Blobby3D
make