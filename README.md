# Blobby3D

Blobby3D is designed to perform Bayesian inference for gas kinematics on emission lines observations of IFS galaxy observations. The code was designed with two main goals in mind:

 - To robustly infer gas kinematics for regularly rotating galaxies even if the gas profiles have significant substructure.
 - Infer gas kinematic properties free from the effects of beam smearing. Where beam smearing is the effect of the observation seeing blurring the gas profiles spatially. This has significant effects on the observed gas kinematic properties - particularly the observed velocity dispersion.

A paper on the technique implemented in this code is forthcoming.

## Getting Started

These instructions will get you a copy of the project up and running on your machine for testing purposes. If you have any trouble, feel free to contact Mathew Varidel at mathew.varidel@sydney.edu.au.

### Prerequisites

Blobby3D is dependent on the following:

C++11  
[DNest4](https://github.com/eggplantbren/DNest4)  
[FFTW3](http://www.fftw.org)  
[Python 3](https://www.python.org) (for postprocessing)  

[FFTW3](http://www.fftw.org) will typically be available on scientific machines/clusters, so you will often either have them or can get them relatively easily by contacting your system administrator.

[DNest4](https://github.com/eggplantbren/DNest4) is a bespoke nested sampling algorithm written by Brendon Brewer, and will need to be cloned from github in the usual manner. It is only available for Unix like machines.

### Installing

To install Blobby3D you need to clone this repository. Then go to the cloned repository directory and make using:

make DNEST4_PATH=/path/to/DNest4

Where /path/to/DNest4 is the root directory for DNest4 on your machine.

### Running Blobby3D on Examples

To get a feel for Blobby3D you will want to start with an example. You can find a few examples in the Blobby3D/Examples folder. There are 3 example toy models provided. The toy models are constructed in accordance with the paper. The parent folder describes the PSF FWHM (psfa), PSF Moffat beta parameter (psfb) and the signal-noise (sn). The subfolder describes the velocity dispersion (vdisp) and the v(R_opt) (vropt) value.

Navigate to one of the folders in the Examples where the MODEL_OPTIONS file resides. The code can be run as an executable from within the subdirectory as follows:

../../Blobby3D -t 1 -f MODEL_OPTIONS

You may want to setup a terminal alias such that you can run the executable from anywhere on your machine. The -t 1 parameter sets the number of threads -- I recommend using multiple threads. The -f MODEL_OPTIONS tells the model to use the input hyperparameters/parameters that define the priors. Note the the DNest4 parameters are defined in the OPTIONS file. Option files are provided in the examples.

There will be 3 output files created as Blobby3D runs:

sample_info.txt  
sample.txt  
sampler_state.txt  

At any time during or after the run, you can perform postprocessing of the current Blobby3D output. Example postprocessing scripts are available in the Examples folders labeled post.py.

<!--
### Running Blobby3D on Your Own Data

To run Blobby3D on your own data you need to setup several files.
[TO COME...]
-->

