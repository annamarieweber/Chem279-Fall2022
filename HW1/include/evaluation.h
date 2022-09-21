#ifndef EVALUATE
#define EVALUATE
#include "cluster.h"
#include <sciplot/Vec.hpp>

using namespace sciplot;
using arma::mat;
using arma::vec;

double calcTruncationError(Cluster c, double h);

void evaluateFDApproximation(Cluster c);

#endif
