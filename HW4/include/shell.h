#ifndef SHELL
#define SHELL
#include <armadillo>
#include <iostream>
#include <cmath>
#include <stdio.h>
using arma::mat;
using arma::rowvec;
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

    Shell(int atomicNum, vec r_a, vec d_k, vec alpha, rowvec l);

    mat operator()(vec r);

    vec r_a();

    vec l_a(int i);

    int num_quantum_arrangements();

    vec alpha();

    mat alphaMat();

    vec d_k();

    void printShellMatrix();
};
#endif