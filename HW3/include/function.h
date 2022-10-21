#include <armadillo>
#include <iostream>
#include <cmath>
#include <stdio.h>
using arma::mat;
using arma::vec;
using std::string;

class Function
{
public:
    template <typename T>
    T operator()(T x)
    {
        return 1 * x;
    }

    template <typename T>
    T first_deriv(T x)
    {
        return 1 * x;
    }

    template <typename T>
    T second_deriv(T x)
    {
        return 1 * x;
    }
};
