#include "Lookup_CDF.h"
#include <cmath>
#include <cassert>

Lookup_CDF Lookup_CDF::instance;

Lookup_CDF::Lookup_CDF()
:num(5001)
,xMin(-5.0)
,xMax(5.0)
,dx((xMax - xMin)/(num-1))
,one_over_dx(1.0/dx)
,_cdf(num)
{
  for(int i=0; i<num; i++)
    _cdf[i] = erf(xMin + i*dx);
}

double Lookup_CDF::evaluate(double x)
{
  int i = static_cast<int>(x/instance.dx);
  double frac = (x - i*instance.dx)*instance.one_over_dx;

  if(i < 0 || i >= (static_cast<int>(instance._cdf.size()) - 1))
    return 0.0;
  
  return frac*instance._cdf[i+1] + (1.0 - frac)*instance._cdf[i];
}
