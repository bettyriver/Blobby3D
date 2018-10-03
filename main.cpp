#include <iostream>

#include "DNest4/code/DNest4.h"

#include "Data.h"
#include "MyModel.h"

int main(int argc, char** argv) {
  // Get command line options
  DNest4::CommandLineOptions options(argc, argv);

  // Get model specific options file
  const char* moptions_file;
  if (options.get_config_file() == "") {
    std::cerr
      <<"# ERROR: No model options file provided. "
	    <<"Specify model options file as a command line option "
	    <<"using -f <MODEL_OPTIONS_FILE>."
      <<std::endl;
    exit(0);
  } else {
    moptions_file = options.get_config_file().c_str();
  }

  // Load data
  Data::get_instance().load(moptions_file);

  // Setup and run sampler
  DNest4::Sampler<MyModel> sampler = DNest4::setup<MyModel>(options);
  sampler.run();

  return 0;
}

