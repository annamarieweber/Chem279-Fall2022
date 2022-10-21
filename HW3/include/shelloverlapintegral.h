#ifndef SHELL_OVERLAP_INTEGRAL
#define SHELL_OVERLAP_INTEGRAL
#include "shell.h"

class ShellOverlapIntegral
{
private:
    Shell _s_a;
    Shell _s_b;

    double alphas_product(int k);

    double alphas_sum(int k);

    double dim_dist_sqr(int dim);

    double exponential_prefactor(int dim, int k);

    double root_term(int k);

    double overlap_summation(double x_p, int l_pair_a, int l_pair_b, int dim, int k);

    double product_center(int dim, int k);

public:
    ShellOverlapIntegral();

    ShellOverlapIntegral(Shell s_a, Shell s_b);

    mat operator()(int k);
};
#endif