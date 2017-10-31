#ifndef DNest4_RJObject_GalaxyField_Lookup_CDF
#define DNest4_RJObject_GalaxyField_Lookup_CDF

#include <vector>

/*
* Lookup tables for speeding things up
* Singleton pattern
*/

class Lookup_CDF
{
	private:
		int num;
		double xMin, xMax, dx, one_over_dx;
		std::vector<double> _cdf; // exp(-x) for x >= 0

		Lookup_CDF();
		Lookup_CDF(const Lookup_CDF& other);

		static Lookup_CDF instance;

	public:
		static double evaluate(double x);

};

#endif

