# Blobby3D

Blobby3D is designed to perform Bayesian inference for gas kinematics on emission lines observations of IFS galaxy observations. The code was designed with two main goals in mind:

 - To robustly infer gas kinematics for regularly rotating galaxies even if the gas profiles have significant substructure.
 - Infer gas kinematic properties free from the effects of beam smearing. Where beam smearing is the effect of the observational seeing spatially blurring the gas profiles. This has significant effects on the observed gas kinematic properties - particularly the observed velocity dispersion.

You can find the motivation and first implementation of Blobby3D in this [paper](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4024V/abstract). Please cite this paper if you use the code in your research.

## Getting Started

These instructions will get you a copy of the project up and running on your machine for testing purposes. More specific instructions will be provided at a later date. If you have any trouble, feel free to contact Mathew Varidel at mathew.varidel@sydney.edu.au.

### Prerequisites

Blobby3D is dependent on the following:

C++11  
[DNest4](https://github.com/eggplantbren/DNest4)  
[FFTW3](http://www.fftw.org)  
[Python 3](https://www.python.org) (for postprocessing)  

[FFTW3](http://www.fftw.org) will typically be available on scientific machines/clusters, so you will often either have them or can get them relatively easily by contacting your system administrator.

[DNest4](https://github.com/eggplantbren/DNest4) is a bespoke nested sampling algorithm written by Brendon Brewer, and will need to be cloned from github in the usual manner. It is only available for Unix like machines.

### Compilling

To install Blobby3D you need to clone this repository. You then need to compile the C++ code and install the python module.

To compile the C++ code type 'make' in the root directory. Note that this assumes that you have setup a DNEST4_PATH environment variable as suggested in the [DNEST4 paper](https://arxiv.org/abs/1606.03757). If you have not setup a DNEST4_PATH environment variable, then simply make using the following:

make DNEST4_PATH=/path/to/DNest4

Where /path/to/DNest4 is the root directory for DNest4 on your machine.

To install the python module you can just do the following:

pip install pyblobby3d

The python module is then available using:

import pyblobby3d

### Running Blobby3D on Examples

To get a feel for Blobby3D you will want to start with an example. There is an example using GAMA 485885 in the examples folder. The code can be run as an executable from within the subdirectory as follows:

../../Blobby3D -t 1 -f MODEL_OPTIONS

You may want to setup a terminal alias such that you can run the executable from anywhere on your machine. The -t 1 parameter sets the number of threads. Blobby3D is computationally expensive, so you should run it using multiple threads. The -f MODEL_OPTIONS tells the model to use the input hyperparameters that define the priors. Note the DNest4 parameters are defined in the OPTIONS file (see [DNest4](https://github.com/eggplantbren/DNest4) documentation for that). Option files are provided in the examples.

There will be 3 output files created as Blobby3D runs:

sample_info.txt  
sample.txt  
sampler_state.txt  

At any time during or after the run, you can perform postprocessing of the current Blobby3D output. An example postprocessing script is available in the examples folders labeled post.py.
