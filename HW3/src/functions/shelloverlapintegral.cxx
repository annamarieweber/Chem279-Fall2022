#include "shelloverlapintegral.h"
#include "combinatorics.h"
#include <armadillo>

using arma::mat;
using arma::vec;

double ShellOverlapIntegral::alphas_product(int k)
{
    return _s_a.alpha()(k) * _s_b.alpha()(k);
}

double ShellOverlapIntegral::alphas_sum(int k)
{
    return _s_a.alpha()(k) + _s_b.alpha()(k);
}

double ShellOverlapIntegral::dim_dist_sqr(int dim)
{
    return pow(_s_a.r_a()(dim) - _s_b.r_a()(dim), 2);
}

double ShellOverlapIntegral::exponential_prefactor(int dim, int k)
{
    return exp(-1.0 * ((alphas_product(k) * dim_dist_sqr(dim)) / alphas_sum(k)));
}

double ShellOverlapIntegral::root_term(int k)
{
    return sqrt(M_PI / alphas_sum(k));
}

double ShellOverlapIntegral::overlap_summation(double x_p, int l_pair_a, int l_pair_b, int dim, int k)
{
    double summation = 0.0;
    for (int i = 0; i <= l_pair_a; i++)
    {
        for (int j = 0; j <= l_pair_b; j++)
        {
            if ((i + j) % 2 == 0)
            {
                double binomial_term = calc_binomial(l_pair_a, i) * calc_binomial(l_pair_b, j);
                double factorial_term = factorial(i + j - 1, 2);
                double a_term = pow(x_p - _s_a.r_a()(dim), l_pair_a - i);
                double b_term = pow(x_p - _s_b.r_a()(dim), l_pair_b - j);
                double denominator = pow(2.0 * alphas_sum(k), (i + j) / 2.0);
                double step = (binomial_term * ((factorial_term * a_term * b_term) / denominator));
                summation += step;
            }
        }
    }
    return summation;
}

double ShellOverlapIntegral::product_center(int dim, int k)
{
    return (_s_a.r_a()(dim) * _s_a.alpha()(k) + _s_b.r_a()(dim) * _s_b.alpha()(k)) / (alphas_sum(k));
}

ShellOverlapIntegral::ShellOverlapIntegral() {}

ShellOverlapIntegral::ShellOverlapIntegral(Shell s_a, Shell s_b)
{
    _s_a = s_a;
    _s_b = s_b;
};

mat ShellOverlapIntegral::operator()(int k)
{
    mat result(_s_a.l_a().n_rows, _s_b.l_a().n_rows, arma::fill::ones);
    for (int i = 0; i < _s_a.r_a().n_elem; i++)
    {
        double x_p = product_center(i, k);
        for (int k = 0; k < _s_a.l_a().n_rows; k++)
        {
            for (int l = 0; l < _s_b.l_a().n_rows; l++)
            {
                double z = exponential_prefactor(i, k) * root_term(k) * overlap_summation(x_p, _s_a.l_a()(k, i), _s_b.l_a()(l, i), i, k);
                std::cout << "res: " << z << std::endl;
                result(k, l) *= z;
            }
        }
    }
    return result;
}
