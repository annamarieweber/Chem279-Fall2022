/*
  UC Berkeley - MSSE Program
  Chem 279
  Fall 2022
  This file, filereader.h, contains an API to all combinatorics functions
  being used in Chem 279
*/
#ifndef COMBO
#define COMBO
#include <string>
#include <armadillo>

using arma::vec;

int factorial(int n);

// calculates the a factorial where the output is the product of all integers from 1 to n where x%a == n%a
int factorial(int n, int a);

int ncr(int n, int r);

double calc_binomial(int m, int n);

vec vec_factorial(int n, int a);

vec vec_factorial(int n);

vec calc_binomial(vec m, vec n);
#endif
