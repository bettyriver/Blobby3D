#include <iostream>
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "MyModel.h"


using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{ 
  // Get command line options
  CommandLineOptions options(argc, argv);
  
  // Get model specific options file
  const char* moptions_file;
  if(options.get_config_file() == "")
    {
    cerr<<"# ERROR: No model options file provided. "
	<<"Specify model options file as a command line option "
	<<"using -f <MODEL_OPTIONS_FILE>."<<std::endl;
    exit(0);
    }
  else
    moptions_file = options.get_config_file().c_str();

  // Load data
  Data::get_instance().load(moptions_file);

  // Setup and run sampler
  Sampler<MyModel> sampler = setup<MyModel>(options);
  sampler.run();

  return 0;
}

