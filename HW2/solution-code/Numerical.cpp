#include <stdio.h>
#include <math.h>
#include "Numerical.h"

using namespace std;

double DE_product_transform(Integrand& f, Integrand& g, double t){
  double c_constant =M_PI /2.0;
  double sinht = sinh(t);
  double x = sinh(c_constant *sinht);
  double dx_dt = c_constant * cosh(t) * cosh(c_constant*sinht);
  double fxgx = f.eval(x) * g.eval(x) ;
  return fxgx*dx_dt;
}


double DEInfinity(Integrand& f1, Integrand& f2, const int nmupoints_oneside, const double stepsize){
  double Sum = DE_product_transform(f1, f2, 0.0 );
  for(int k = 1; k <= nmupoints_oneside; k++){
    Sum += DE_product_transform(f1, f2, k*stepsize);
    Sum += DE_product_transform(f1, f2,  -k*stepsize );
  }
  return Sum * stepsize;
}


double Gussian::eval(double x)
{
  double xd = x - x0;  
  double y = -xd*xd * alpha;
  double l_part = pow(xd , l);

  return l_part * exp(y);
}


