/*
  UC Berkeley - MSSE Program
  Chem 279
  Fall 2022
  This file, combinatorics.cxx, contains utilities for calculating combinations and factorials

*/

#include <iostream>

int factorial(int n)
{
    int product = n;
    for (int i = n - 1; i > 0; i--)
    {
        product *= i;
    }
    if (product <= 1)
    {
        return 1;
    }
    return product;
}

// calculates the a factorial where the output is the product of all integers from 1 to n where x%a == n%a
int factorial(int n, int a)
{
    int product = n;
    for (int i = n; i > 0; i -= a)
    {
        product *= i;
    }
    if (product < 1)
    {
        return 1;
    }
    return product;
}

int ncr(int n, int r)
{
    if (n == r)
        return 1;
    if (r == 0 && n != 0)
        return 1;
    else
    {
        return factorial(n) / (factorial(r) * factorial(n - r));
    }
}

double calc_binomial(int m, int n)
{
    return 1.0 * factorial(m) / (1.0 * factorial(n) * factorial(m - n));
}
