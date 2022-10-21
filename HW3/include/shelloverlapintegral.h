#ifndef SHELL_OVERLAP_INTEGRAL
#define SHELL_OVERLAP_INTEGRAL
#include "shell.h"

class ShellOverlapIntegral
{
private:
    Shell _s_a;
    Shell _s_b;

    vec alphas_product();

    vec alphas_sum();

    double dim_dist_sqr(int dim);

    double exponential_prefactor(int dim);

    vec root_term();

    double overlap_summation(double x_p, int l_pair_a, int l_pair_b, int dim);

    double product_center(int dim);

public:
    ShellOverlapIntegral();

    ShellOverlapIntegral(Shell s_a, Shell s_b);

    mat operator()();
};
#endif