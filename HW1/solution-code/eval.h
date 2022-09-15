#ifndef EVAL_H
#define EVAL_H

#include "readfile.h"
#include <armadillo>

double E_LJ(const std::vector<Atom> &Atoms);

double E_LJ(const arma::mat &coords);

void F_LJ_forward(arma::mat &force, const std::vector<Atom> &Atoms, double stepsize);

void F_LJ_central(arma::mat &force, const std::vector<Atom> &Atoms, double stepsize);

void F_LJ_analytic(arma::mat &force, const std::vector<Atom> &Atoms);

void Steepest_descend( std::vector<Atom> &opt_Atoms, const std::vector<Atom> &Atoms, double fdstepsize, double thresh);

void Steepest_descend_line_search( std::vector<Atom> &opt_Atoms, const std::vector<Atom> &Atoms,
double fdstepsize, double thresh);
#endif // EVAL_H