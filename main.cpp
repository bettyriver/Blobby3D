#include <iostream>
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include "MyModel.h"


using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{

  Data::get_instance().load("MODEL_OPTIONS");

  Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
  sampler.run();

  return 0;
}

