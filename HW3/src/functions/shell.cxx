#include <armadillo>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include "combinatorics.h"
#include "shell.h"
#include "constants.h"

using arma::mat;
using arma::vec;
using std::string;
using namespace constants;

mat Shell::angularMomentum(int l)
{
    mat amMatrix(ncr(l + 3 - 1, 3 - 1), 3);

    int idx = 0;
    for (int i = 0; i <= l; i++)
    {
        for (int j = 0; j <= l - i; j++)
        {
            amMatrix(idx, 0) = i;
            amMatrix(idx, 1) = j;
            amMatrix(idx, 2) = l - i - j;
            idx++;
        }
    }

    return amMatrix;
}

Shell::Shell(){};

Shell::Shell(int atomicNum, vec r_a, int l)
{
    _r_a = r_a;
    _l_a = angularMomentum(l);

    mat coeffs(ATOM_COEFS_MAP[atomicNum - 1]);
    _d_k = coeffs.col(1 + l);
    _alpha = coeffs.col(0);
};

void Shell::printShellMatrix()
{
    std::cout << "r_a: " << std::endl;
    std::cout << _r_a << std::endl;
    std::cout << "l_a: " << std::endl;
    std::cout << _l_a << std::endl;
    std::cout << "d_k: " << std::endl;
    std::cout << _d_k << std::endl;
    std::cout << "alpha: " << std::endl;
    std::cout << _alpha << std::endl;
};

mat Shell::operator()(vec r)
{
    mat term1(_r_a.n_elem, _l_a.n_rows);
    for (int i = 0; i < _l_a.n_rows; i++)
    {
        term1(i) = prod(pow(r - _r_a, _l_a(i)));
    }
    vec term2 = exp((-1.0 * _alpha) % pow((r - _r_a), 2));
    return term1 % term2;
}

vec Shell::r_a()
{
    return _r_a;
}

mat Shell::l_a()
{
    return _l_a;
}

vec Shell::d_k()
{
    return _d_k;
}

vec Shell::alpha()
{
    return _alpha;
}
