#include "shelloverlapintegral.h"
#include "combinatorics.h"
#include <armadillo>

using arma::cube;
using arma::mat;
using arma::vec;

vec ShellOverlapIntegral::alphas_product()
{
    return _s_a.alpha() % _s_b.alpha();
}

vec ShellOverlapIntegral::alphas_sum()
{
    return _s_a.alpha() + _s_b.alpha();
}

vec ShellOverlapIntegral::dim_dist_sqr()
{
    return pow(_s_a.r_a() - _s_b.r_a(), 2);
}

vec ShellOverlapIntegral::exponential_prefactor()
{
    return exp(-1.0 * ((alphas_product() % dim_dist_sqr()) / alphas_sum()));
}

vec ShellOverlapIntegral::root_term()
{
    return sqrt(M_PI / alphas_sum());
}

vec ShellOverlapIntegral::overlap_summation(vec x_p, int l_pair_a, int l_pair_b)
{
    vec summation(x_p.n_elem);
    for (int i = 0; i <= l_pair_a; i++)
    {
        for (int j = 0; j <= l_pair_b; j++)
        {
            if ((i + j) % 2 == 0)
            {
                double binomial_term = calc_binomial(l_pair_a, i) * calc_binomial(l_pair_b, j);
                double factorial_term = factorial(i + j - 1, 2);
                vec a_term = pow(x_p - _s_a.r_a(), l_pair_a - i);
                vec b_term = pow(x_p - _s_b.r_a(), l_pair_b - j);
                vec denominator = pow(2.0 * alphas_sum(), (i + j) / 2.0);
                vec step = (binomial_term * ((factorial_term * a_term % b_term) / denominator));
                summation += step;
            }
        }
    }
    return summation;
}

vec ShellOverlapIntegral::product_center()
{
    return ((_s_a.alpha() % _s_a.r_a()) + (_s_b.alpha() % _s_b.r_a())) / (alphas_sum());
}

ShellOverlapIntegral::ShellOverlapIntegral() {}

ShellOverlapIntegral::ShellOverlapIntegral(Shell s_a, Shell s_b)
{
    _s_a = s_a;
    _s_b = s_b;
};

mat ShellOverlapIntegral::operator()()
{
    vec x_p = product_center();
    mat x_cent = _s_a.alphaMat() * _s_a.r_a() + _s_b.alphaMat() * _s_b.r_a();

    mat result(_s_a.l_a(0).n_elem, _s_b.l_a(0).n_elem, arma::fill::ones);
    for (int i = 0; i < product_center().n_elem; i++)
    {

        int elem = 0;
        for (int k = 0; k < _s_a.num_quantum_arrangements(); k++)
        {
            vec ktotal(x_p.n_elem, arma::fill::zeros);

            for (int l = 0; l < _s_b.num_quantum_arrangements(); l++)
            {
                ktotal += (exponential_prefactor() % root_term() % overlap_summation(x_p, _s_a.l_a(k)(i), _s_b.l_a(l)(i)));
            }
            result.col(k) %= ktotal;
        }
    }

    return result;
}