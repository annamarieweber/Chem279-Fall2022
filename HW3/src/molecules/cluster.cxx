#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "cluster.h"
#include "baseexception.h"
#include "constants.h"
#include "combinatorics.h"
#include "shell.h"
#include "shelloverlapintegral.h"

using arma::cube;
using arma::field;
using arma::mat;
using arma::rowvec;
using arma::span;
using arma::sum;
using arma::vec;
using constants::ATOM_BASIS_FN_MAP;
using std::string;
using namespace constants;
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

class Gaussian
{
private:
  double _x_a;
  double _l_a;
  double _alpha;

public:
  /**
   * copy constructor for a gaussian
   **/
  void operator=(const Gaussian &g)
  {
    _x_a = g._x_a;
    _l_a = g._l_a;
    _alpha = g._alpha;
  }

  Gaussian(){};
  /**
   * constructor for a gaussian
   **/
  Gaussian(double x_a, double l_a, double alpha)
  {
    _x_a = x_a;
    _l_a = l_a;
    _alpha = alpha;
  }

  double alpha()
  {
    return _alpha;
  }

  double l_a()
  {
    return _l_a;
  }

  double x_a()
  {
    return _x_a;
  }

  /**
   * evaluation of the Gaussian Function = (x − XA)
   **/
  template <typename T>
  T operator()(T x)
  {
    T res = pow(x - _x_a, _l_a) * exp(-1.0 * _alpha * pow(x - _x_a, 2));
    return res;
  };

  /**
   * first derivative of the gaussian
   * e^(-alpha (x - _x_a)^2) (x - _x_a)^(_l_a - 1) (_l_a - 2 alpha (x - _x_a)^2)
   **/
  template <typename T>
  T first_deriv(T x)
  {
    return exp(-1.0 * _alpha * pow(x - _x_a, 2)) * pow(x - _x_a, _l_a - 1) * (_l_a - 2 * _alpha * pow(x - _x_a, 2));
  }

  /**
   * second derivative of the gaussian
   * e^(-_alpha (-_x_a + x)^2) (-_x_a + x)^(_l_a - 2) (4 _alpha^2 (_x_a - x)^4 - 2 _alpha (2 _l_a + 1) (-_x_a + x)^2 + _l_a (-1 + _l_a))
   **/
  template <typename T>
  T second_deriv(T x)
  {
    return exp(-1.0 * _alpha * pow(-1.0 * (_x_a + x), 2) * pow(-1.0 * _x_a + x, _l_a - 2)) * (4 * pow(_alpha, 2) * pow(_x_a - x, 4) - 2 * _alpha * (2 * _l_a + 1) * (-1 * pow(_x_a + x, 2)) + _l_a * (-1 + _l_a));
  }
};

template <typename FA, typename FB>
class FunctionProduct
{
private:
  FA _fn_a;
  FB _fn_b;

public:
  /**
   * copy constructor for function product
   **/
  void operator=(const FunctionProduct &f)
  {
    _fn_a = f._fn_a;
    _fn_b = f._fn_b;
  };

  FunctionProduct(){};

  /**
   * constructor for product of functions
   **/
  FunctionProduct(FA fn_a, FB fn_b)
  {
    _fn_a = fn_a;
    _fn_b = fn_b;
  }

  /**
   * computes product of 2 functions fn_a and fn_b at x
   **/
  template <typename T>
  T operator()(T x)
  {

    T res = _fn_a(x) * _fn_b(x);
    return res;
  }

  /**
   * first derivitive using the product rule
   **/
  template <typename T>
  T first_deriv(T x)
  {
    return _fn_a.first_deriv(x) * _fn_b(x) + _fn_b.first_deriv(x) * _fn_a(x);
  }

  /**
   * second derivative using the product rule: d^2/dx^2(f(x) g(x)) = g(x) f''(x) + 2 f'(x) g'(x) + f(x) g''(x)
   **/
  template <typename T>
  T second_deriv(T x)
  {
    return _fn_b(x) * _fn_a.second_deriv(x) + 2 * _fn_a.first_deriv(x) * _fn_b.first_deriv(x) + _fn_a(x) * _fn_b.second_deriv(x);
  }
};

class OverlapIntegral
{
private:
  Gaussian _s_a;
  Gaussian _s_b;

public:
  OverlapIntegral()
  {
  }

  /**
   * constructor for integral
   **/
  OverlapIntegral(Gaussian s_a, Gaussian s_b)
  {
    _s_a = s_a;
    _s_b = s_b;
  }

  /*
   * calculates the integeral using the extended trapazoidal rule
   **/

  double operator()()
  {
    double alpha_prod = _s_a.alpha() * _s_b.alpha();
    double alpha_sum = _s_a.alpha() * _s_b.alpha();
    double dim_dist_sqr = pow(_s_a.x_a() - _s_b.x_a(), 2);
    double exponential_prefactor = exp(-1.0 * ((alpha_prod * dim_dist_sqr) / alpha_sum));
    double x_p = (_s_a.x_a() * _s_a.alpha() + _s_b.x_a() * _s_b.alpha()) / (alpha_sum);
    double l_pair_a = _s_a.l_a();
    double l_pair_b = _s_b.l_a();
    double summation = 0.0;
    for (int i = 0; i <= l_pair_a; i++)
    {
      for (int j = 0; j <= l_pair_b; j++)
      {
        if ((i + j) % 2 == 0)
        {
          double binomial_term = calc_binomial(l_pair_a, i) * calc_binomial(l_pair_b, j);
          double factorial_term = factorial(i + j - 1, 2);
          double a_term = pow(x_p - _s_a.x_a(), l_pair_a - i);
          double b_term = pow(x_p - _s_b.x_a(), l_pair_b - j);
          double denominator = pow(2.0 * alpha_sum, (i + j) / 2.0);
          double step = (binomial_term * ((factorial_term * a_term * b_term) / denominator));
          summation += step;
        }
      }
    }

    double result = exponential_prefactor * (sqrt(M_PI / alpha_sum)) * summation;
    return result;
  }
};

Cluster::Cluster(){};

Cluster::Cluster(int numAtoms)
{
  atomMatrix = mat(numAtoms, 4);

  // todo get rid of these member vectors and be able to look up params for given atom on the fly (could store in constants)
  epsilons = vec(numAtoms);
  sigmas = vec(numAtoms);
}

bool Cluster::allAtomsMatch(int num)
{
  return all(atomMatrix.col(0) == num);
}

int Cluster::countBasisFunctions()
{
  int sum = 0;
  vec atoms = atomMatrix.col(0);
  for (int i = 0; i < atoms.size(); i++)
  {
    sum += ATOM_BASIS_FN_MAP[atoms[i] - 1];
  }
  return sum;
}

int Cluster::countElectronPairs()
{
  int n = countBasisFunctions();
  if (n / 2 * 2 != n)
  {
    throw BaseException("InvalidConfig: number of electron pairs mus be an integer");
  }
  else
  {
    return n / 2;
  }
}

// void Cluster::printBasisFuncs()
// {

//   for (int i = 0; i < atomMatrix.n_rows; i++)
//   {
//     vec l_vals(ATOM_ANGULAR_MOMENTUM_MAP[atomMatrix(i, 0) - 1]);
//     for (int j = 0; j < l_vals.n_elem; j++)
//     {
//       Shell s(atomMatrix(i, 0), atomMatrix.row(i).cols(1, 3).t(), l_vals(j));
//       s.printShellMatrix();
//     }
//   }
// }

void Cluster::evalHw3()
{
  int K = 3;
  int b = countBasisFunctions();
  mat r(b, K);
  vec l(b);
  vector<int> atom_num(b);
  vector<mat> n_k;
  mat basisfns(b, 5 * K + 1);

  int fn = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    vec l_vals(ATOM_ANGULAR_MOMENTUM_MAP[atomMatrix(i, 0) - 1]);
    for (int j = 0; j < l_vals.n_elem; j++)
    {
      // for (int k = 0; k < SECONDARY_QUANTUM_NUMBER_COMBINATIONS[l_vals(j)]; k++)
      // {
      //       }

      // for (int k = 0; k < K; k++)
      // {
      // l(fn) = l_vals(j);
      // r.row(fn) = atomMatrix.row(i).cols(1, 3);
      // atom_num[fn] = atomMatrix(i, 0);

      Shell s(atomMatrix(i, 0), r.row(fn).t(), l_vals(j));
      ShellOverlapIntegral s_aa(s, s);
      mat norm_consts = 1.0 / sqrt(s_aa());
      for (int p = 0; p < s.num_quantum_arrangements(); p++)
      {
        basisfns.row(fn).cols(0, K) = atomMatrix.row(i).cols(0, K);
        basisfns.row(fn).cols(K + 1, 2 * K) = s.l_a(p).t();
        basisfns.row(fn).cols(2 * K + 1, 3 * K) = s.alpha().t();
        basisfns.row(fn).cols(3 * K + 1, 4 * K) = s.d_k().t();
        basisfns.row(fn).cols(4 * K + 1, 5 * K) = norm_consts.col(p).t();
        fn++;
      }

      // }
    }
  }

  std::cout << basisfns << std::endl;

  mat result(b, b);
  for (int m = 0; m < b; m++)
  {
    for (int v = 0; v < b; v++)
    {

      Shell s_m(basisfns(m, 0), basisfns.row(m).cols(1, K).t(), basisfns.row(m).cols(3 * K + 1, 4 * K).t(), basisfns.row(m).cols(2 * K + 1, 3 * K).t(), basisfns.row(m).cols(K + 1, 2 * K));
      Shell s_v(basisfns(v, 0), basisfns.row(v).cols(1, K).t(), basisfns.row(v).cols(3 * K + 1, 4 * K).t(), basisfns.row(v).cols(2 * K + 1, 3 * K).t(), basisfns.row(v).cols(K + 1, 2 * K));
      ShellOverlapIntegral s(s_m, s_v);
      // std::cout << s_a.d_k() << std::endl;
      // std::cout << s_b.d_k() << std::endl;
      // std::cout << n_k[e] << std::endl;

      // std::cout << n_k[f] << std::endl;
      // std::cout << n_k(e) << std::endl;
      // std::cout << n_k(f) << std::endl;
      // std::cout << o() << std::endl;
      // std::cout << "------------------------------------------------" << std::endl;
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      rowvec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K);
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      rowvec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K);
      // std::cout << (dk_m * dk_v) % (nk_m * nk_v) << std::endl;
      // std::cout << (dk_m * dk_v) % (nk_m * nk_v) % s() << std::endl;
      // std::cout << (dk_m * dk_v) % (nk_m * nk_v) * s() << std::endl;
      // std::cout << (dk_m * dk_v) * (nk_m * nk_v) * s() << std::endl;
      // std::cout << (dk_m % dk_v.t()) * (nk_m % nk_v.t()).t() * s() << std::endl;
      // std::cout << (dk_m % dk_v.t()) % (nk_m % nk_v.t()) * s().t() << std::endl;
      // std::cout << "d_m" << std::endl;
      // std::cout << dk_m << std::endl;
      // std::cout << "d_v" << std::endl;
      // std::cout << dk_v << std::endl;
      // std::cout << (dk_m % nk_m) * (dk_v % nk_v) << std::endl;
      // std::cout << (dk_m % nk_m % s()) % (dk_v % nk_v % s().t()).t() << std::endl;
      // std::cout << (dk_m % nk_m % s()) * (dk_v % nk_v % s().t()) << std::endl;
      // std::cout << sum(sum((dk_m % nk_m % s()) * (dk_v % nk_v % s().t()))) << std::endl;
      // std::cout << sum(sum((dk_m % nk_m) * (dk_v % nk_v))) * s() << std::endl;
      // std::cout << sum(sum((dk_m % nk_m) * (dk_v % nk_v))) * (3 * s()) << std::endl;
      // result(e, f) = sum(sum((s_a.d_k() % s_b.d_k()) % (n_k[e] % n_k[f]) % o()));

      // std::cout << "---n_k-d_k-crossproduct-------------------------" << std::endl
      //           << (s_a.d_k() * s_b.d_k().t()) % (n_k[e] * n_k[f].t()) << std::endl
      //           << "-------------------------------------------------" << std::endl;
      // std::cout << "---n_k-d_k-crossproduct--------------------------" << std::endl
      //           << (s_a.d_k() * s_b.d_k().t()) * (n_k[e] * n_k[f].t()) << std::endl
      //           << "-------------------------------------------------" << std::endl;
      // std::cout << "---d_k-matrix-----------------------------------" << std::endl
      //           << (s_a.d_k() * s_b.d_k().t()) << std::endl
      //           << "-------------------------------------------------" << std::endl;
      // std::cout << "----n_k-matrix-----------------------------------" << std::endl
      //           << (n_k[e] * n_k[f].t()) << std::endl
      //           << "-------------------------------------------------" << std::endl;
      // std::cout << (n_k[e] % n_k[f] * prod(o(), 0)) << std::endl;
      // std::cout << "------------------------------------------------" << std::endl;
      // std::cout << (prod(n_k[e] % n_k[f] % o(), 0)) << std::endl;
      // std::cout << "------------------------------------------------" << std::endl;

      for (int k = 0; k < K; k++)
      {
        for (int l = 0; l < K; l++)
        {
          Gaussian s_k(basisfns(m, 1 + k), basisfns(m, K + 1 + k), basisfns(m, 2 * K + 1 + k));
          Gaussian s_l(basisfns(v, 1 + l), basisfns(v, K + 1 + l), basisfns(v, 2 * K + 1 + l));
          OverlapIntegral s_kl(s_k, s_l);
          // std::cout << "ra" << std::endl
          //           << basisfns.row(m).cols(1 + k, 1 + k).t() << std::endl;
          // std::cout << "dk" << std::endl
          //           << basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t() << std::endl;
          // std::cout << "alpha" << std::endl
          //           << basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t() << std::endl;
          // std::cout << "l" << std::endl
          //           << basisfns.row(m).cols(K + 1 + k, K + 1 + k) << std::endl;
          std::cout
              << "SKL" << std::endl
              << s_kl() << std::endl;

          std::cout << dk_m(k) << std::endl;
          std::cout << dk_v(l) << std::endl;
          std::cout << nk_m(k) << std::endl;
          std::cout << nk_v(l) << std::endl;
          std::cout << "kl total: " << std::endl
                    << s_kl() << std::endl;
          std::cout << "other: " << std::endl
                    << dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l)
                    << std::endl;

          result(m, v) += s_kl();
        }
      }
      std::cout << dk_m << std::endl;
      std::cout << dk_v << std::endl;
      std::cout << nk_m << std::endl;
      std::cout << nk_v << std::endl;
    }
  }

  std::cout << result << std::endl;

  // std::cout << "printing self shell overlaps: " << std::endl;
  // std::cout << n_k << std::endl;
  // std::cout << n_k(1) << std::endl;
  // for (int z = 0; z < shellOverlaps.size(); z += 3)
  // {
  //   double test = 1.0 / sqrt(norm(shellOverlaps[z]));
  //   mat overlapVal_a = 1.0 / sqrt(shellOverlaps[z]);
  //   mat overlapVal_b = 1.0 / sqrt(shellOverlaps[z + 1]);
  //   mat overlapVal_c = 1.0 / sqrt(shellOverlaps[z + 2]);
  //   std::cout << test << std::endl;
  //   std::cout << overlapVal_a << std::endl;
  //   std::cout << overlapVal_b << std::endl;
  //   std::cout << overlapVal_c << std::endl;
  // }
}

double Cluster::calcDistance(mat a1, mat a2)
{

  double distance = 0;
  rowvec coordDists = square(a2.row(0) - a1.row(0));
  return sqrt(sum(coordDists));
}

double Cluster::sigma_ij(int i, int j)
{
  return sqrt(sigmas(i) * sigmas(j));
}

double Cluster::epsilon_ij(int i, int j)
{
  return sqrt(epsilons(i) * epsilons(j));
}

double Cluster::calcLJ(double r_ij, int i, int j)
{
  if (r_ij == 0)
  {
    throw BaseException("InvalidAtomPosition: Distance between atoms must be > 0. Check coordinates to ensure coordinates are valid.");
  }
  double sigma_over_r_ij = sigma_ij(i, j) / r_ij;
  double term_r6 = pow(sigma_over_r_ij, 6.0);
  double term_r12 = term_r6 * term_r6;
  return epsilon_ij(i, j) * (term_r12 - (2 * term_r6));
}

double Cluster::calcLJEnergy()
{
  double total_energy = 0.0;

  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    for (int j = i + 1; j < atomMatrix.n_rows; j++)
    {
      double dist_ij = calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3));
      total_energy += calcLJ(dist_ij, i, j);
    }
  };

  return total_energy;
}

double Cluster::calcLJPrime(double r_ij, int i, int j)
{
  // \todo: store literal values in definitions for readability
  double term_r6 = pow(sigma_ij(i, j), 6) / pow(r_ij, 7);
  double term_r12 = -1 * pow(sigma_ij(i, j), 12) / pow(r_ij, 13);
  return epsilon_ij(i, j) * 12 * (term_r12 + term_r6);
}

mat Cluster::calcAnalyticalForce()
{
  mat analytical_force(atomMatrix.n_rows, 3);
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    for (int j = 0; j < atomMatrix.n_rows; j++)
    {
      if (i != j)
      {
        double dist_ij = calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3));
        analytical_force(i, 0) += -1 * calcLJPrime(dist_ij, i, j) * ((atomMatrix(i, 1) - atomMatrix(j, 1)) / dist_ij);
        analytical_force(i, 1) += -1 * calcLJPrime(dist_ij, i, j) * ((atomMatrix(i, 2) - atomMatrix(j, 2)) / dist_ij);
        analytical_force(i, 2) += -1 * calcLJPrime(dist_ij, i, j) * ((atomMatrix(i, 3) - atomMatrix(j, 3)) / dist_ij);
      }
    }
  }
  return analytical_force;
}

mat Cluster::calcForwardFdForce(double h)
{
  mat energy_h(atomMatrix.n_rows, 3);
  mat energy(atomMatrix.n_rows, 3);
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {

    for (int j = 0; j < atomMatrix.n_rows; j++)
    {
      if (i != j)
      {
        mat shift_x = atomMatrix.row(j).cols(1, 3);
        shift_x(0, 0) += h;
        mat shift_y = atomMatrix.row(j).cols(1, 3);
        shift_y(0, 1) += h;
        mat shift_z = atomMatrix.row(j).cols(1, 3);
        shift_z(0, 2) += h;

        energy_h(j, 0) += (calcLJ(calcDistance(shift_x, atomMatrix.row(i).cols(1, 3)), i, j));
        energy(j, 0) += (calcLJ(calcDistance(atomMatrix.row(j).cols(1, 3), atomMatrix.row(i).cols(1, 3)), i, j));
        energy_h(j, 1) += (calcLJ(calcDistance(shift_y, atomMatrix.row(i).cols(1, 3)), i, j));
        energy(j, 1) += (calcLJ(calcDistance(atomMatrix.row(j).cols(1, 3), atomMatrix.row(i).cols(1, 3)), i, j));
        energy_h(j, 2) += (calcLJ(calcDistance(shift_z, atomMatrix.row(i).cols(1, 3)), i, j));
        energy(j, 2) += (calcLJ(calcDistance(atomMatrix.row(j).cols(1, 3), atomMatrix.row(i).cols(1, 3)), i, j));
      }
    }
  }

  return (energy_h - energy) / (-1 * h);
}

mat Cluster::calcCentralFdForce(double h)
{
  mat energy_h(atomMatrix.n_rows, 3);
  mat energy(atomMatrix.n_rows, 3);
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    for (int j = 0; j < atomMatrix.n_rows; j++)
    {
      if (i != j)
      {

        mat shift_x_1 = atomMatrix.row(i).cols(1, 3);
        mat shift_x_2 = atomMatrix.row(i).cols(1, 3);
        shift_x_1(0, 0) += h;
        shift_x_2(0, 0) -= h;
        mat shift_y_1 = atomMatrix.row(i).cols(1, 3);
        shift_y_1(0, 1) += h;
        mat shift_y_2 = atomMatrix.row(i).cols(1, 3);
        shift_y_2(0, 1) -= h;
        mat shift_z_1 = atomMatrix.row(i).cols(1, 3);
        shift_z_1(0, 2) += h;
        mat shift_z_2 = atomMatrix.row(i).cols(1, 3);
        shift_z_2(0, 2) -= h;

        energy_h(i, 0) += (calcLJ(calcDistance(shift_x_1, atomMatrix.row(j).cols(1, 3)), i, j));
        energy_h(i, 0) -= (calcLJ(calcDistance(shift_x_2, atomMatrix.row(j).cols(1, 3)), i, j));
        energy_h(i, 1) += (calcLJ(calcDistance(shift_y_1, atomMatrix.row(j).cols(1, 3)), i, j));
        energy_h(i, 1) -= (calcLJ(calcDistance(shift_y_2, atomMatrix.row(j).cols(1, 3)), i, j));
        energy_h(i, 2) += (calcLJ(calcDistance(shift_z_1, atomMatrix.row(j).cols(1, 3)), i, j));
        energy_h(i, 2) -= (calcLJ(calcDistance(shift_z_2, atomMatrix.row(j).cols(1, 3)), i, j));
      }
    }
  }

  return energy_h / (-2 * h);
}

void Cluster::addAtom(int index, int atomNum, double x, double y, double z, double e, double s)
{
  atomMatrix(index, 0) = atomNum;
  atomMatrix(index, 1) = x;
  atomMatrix(index, 2) = y;
  atomMatrix(index, 3) = z;
  epsilons(index) = e;
  sigmas(index) = s;
}

std::ostream &operator<<(std::ostream &os, const Cluster &c)
{
  os << "Atom Matrix: \n"
     << c.atomMatrix << "\n";

  return os;
}
