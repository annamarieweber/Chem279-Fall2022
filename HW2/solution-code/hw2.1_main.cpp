#include <stdlib.h>
#include <stdio.h>
#include "Numerical.h"


int main(int argc, char* argv[])
{
  //Create 1d Gaussain functions, the input is x_coor, alpha, l
  Gussian s1(0.0, 1.0 , 0);
  Gussian p1(0.0, 1.0 , 1);
  Gussian s2(1.0, 1.0 , 0);
  Gussian p2(1.0, 1.0 , 1);
  //I use DE rule to calculate the integral, the input is func, func, numpoints, stepsize
  const double s1s1 = DEInfinity(s1, s1, 3000, 0.001);
  const double s1p1 = DEInfinity(s1, p1, 3000, 0.001);
  const double s1s2 = DEInfinity(s1, s2, 3000, 0.001);
  const double s1p2 = DEInfinity(s1, p2, 3000, 0.001);
  printf("Numerical Integral %1.17e for s1s1\n", s1s1);
  printf("Numerical Integral %1.17e for s1p1\n", s1p1);
  printf("Numerical Integral %1.17e for s1s2\n", s1s2);
  printf("Numerical Integral %1.17e for s1p2\n", s1p2);
  return EXIT_SUCCESS;
}
