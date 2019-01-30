# GAMA 485885

This is an example for data set using data from the [SAMI Galaxy Survey](https://sami-survey.org/) of [GAMA 485885](http://www.gama-survey.org/dr3/tools/sov.php).

## Key Files

MODEL_OPTIONS:  
Blobby3D model options.

OPTIONS:  
DNEST4 options file.

data.txt:  
Continuum subtracted data cube around the Halpha and surrounding [NII] lines at 6583.1 A and 6548.1 A. Data is in row-major order.

var.txt:  
Variance cube in row-major order.

metadata.txt:  
Metadata for the data specifying the number of spaxels in each direction and the (min, max) pairs. In order: Ni, Nj, Nk, x_min, x_max, y_min, y_max, wave_min, wave_max. Note that (min, max) values are specified by the edge of the lower/upper-most bins (ie. NOT the centre of the bin).

post:  
An example postprocessing file.