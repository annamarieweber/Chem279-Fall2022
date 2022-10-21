#ifndef SHELL
#define SHELL
#include <armadillo>
#include <iostream>
#include <cmath>
#include <stdio.h>
using arma::mat;
using arma::vec;
using std::string;

class Shell
{
private:
    vec _r_a; // center (x,y,z, ...)
    mat _l_a; // matrix containing all possible (l,m,n....) combinations
    vec _d_k;
    vec _alpha;
    mat angularMomentum(int l);

public:
    Shell();

    Shell(int atomicNum, vec r_a, int l);

    mat operator()(vec r);

    vec r_a();

    mat l_a();

    vec alpha();

    vec d_k();

    void printShellMatrix();
};
#endif