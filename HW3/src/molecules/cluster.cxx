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

using std::string;
using namespace constants;
using namespace arma;

Cluster::Cluster(){};

Cluster::Cluster(int numAtoms)
{
  atomMatrix = mat(numAtoms, 4);
  K = 3;
  epsilons = vec(numAtoms);
  sigmas = vec(numAtoms);
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

mat Cluster::basisFunctions()
{
  int b = countBasisFunctions();
  vec l(b);

  mat basisfns(b, 5 * K + 2);

  int fn = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    vec l_vals(ATOM_ANGULAR_MOMENTUM_MAP[atomMatrix(i, 0) - 1]);
    for (int j = 0; j < l_vals.n_elem; j++)
    {

      Shell s(atomMatrix(i, 0), atomMatrix.row(i).cols(1, atomMatrix.n_cols - 1).t(), l_vals(j));
      ShellOverlapIntegral s_aa(s, s);
      mat norm_consts = 1.0 / sqrt(s_aa());
      for (int p = 0; p < s.num_quantum_arrangements(); p++)
      {
        basisfns.row(fn).cols(0, K) = atomMatrix.row(i).cols(0, K);
        basisfns.row(fn).cols(K + 1, 2 * K) = s.l_a(p).t();
        basisfns.row(fn).cols(2 * K + 1, 3 * K) = s.alpha().t();
        basisfns.row(fn).cols(3 * K + 1, 4 * K) = s.d_k().t();
        basisfns.row(fn).cols(4 * K + 1, 5 * K) = (l_vals(j) + 1) * norm_consts.col(p).t();
        basisfns.row(fn).col(5 * K + 1) = j;
        fn++;
      }
    }
  }
  return basisfns;
}

mat Cluster::overlapMatrix()
{
  mat basisfns = basisFunctions();
  int b = countBasisFunctions();

  mat result(b, b);
  for (int m = 0; m < b; m++)
  {
    for (int v = 0; v < b; v++)
    {
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      rowvec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K);
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      rowvec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K);

      for (int k = 0; k < K; k++)
      {
        for (int l = 0; l < K; l++)
        {
          Shell s_kx(basisfns(m, 0), basisfns.row(m).cols(1, 1).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 1, K + 1));
          Shell s_lx(basisfns(v, 0), basisfns.row(v).cols(1, 1).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 1, K + 1));
          ShellOverlapIntegral s_klx(s_kx, s_lx);

          Shell s_ky(basisfns(m, 0), basisfns.row(m).cols(2, 2).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 2, K + 2));
          Shell s_ly(basisfns(v, 0), basisfns.row(v).cols(2, 2).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 2, K + 2));
          ShellOverlapIntegral s_kly(s_ky, s_ly);

          Shell s_kz(basisfns(m, 0), basisfns.row(m).cols(3, 3).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 3, K + 3));
          Shell s_lz(basisfns(v, 0), basisfns.row(v).cols(3, 3).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 3, K + 3));
          ShellOverlapIntegral s_klz(s_kz, s_lz);

          result(m, v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx() * s_kly() * s_klz()));
        }
      }
    }
  }

  return result;
}

mat Cluster::extendedHuckelHamiltonian()
{
  mat basisfns = basisFunctions();
  mat ov_mat = overlapMatrix();
  int b = countBasisFunctions();

  mat result(b, b);
  for (int m = 0; m < b; m++)
  {
    for (int v = 0; v < b; v++)
    {
      int a_m = basisfns(m, 0);
      int l_m = basisfns(m, 5 * K + 1);
      int a_v = basisfns(v, 0);
      int l_v = basisfns(v, 5 * K + 1);
      string h_m_key = ATOM_SYMBOLS[a_m - 1] + VALENCE_SHELL_CONFIGS[a_m - 1][l_m];
      string h_v_key = ATOM_SYMBOLS[a_v - 1] + VALENCE_SHELL_CONFIGS[a_v - 1][l_v];
      if (m == v)
      {
        result(m, v) = EXTENDED_HUCKEL_MAP.at(h_m_key);
      }
      else
      {
        result(m, v) = 0.5 * K_VAL * (EXTENDED_HUCKEL_MAP.at(h_m_key) + EXTENDED_HUCKEL_MAP.at(h_v_key)) * ov_mat(m, v);
      }
    }
  }
  return result;
}

mat Cluster::molecularOrbitalCoefficients()
{
  vec eigval;
  mat eigvec;
  mat overlap_matrix = overlapMatrix();
  mat hamiltonian = extendedHuckelHamiltonian();

  eig_sym(eigval, eigvec, overlap_matrix);

  // make the orthogonalization transformation
  mat x_diag = arma::diagmat(arma::pow(eigval, -0.5));
  mat x = eigvec * x_diag * arma::trans(eigvec);
  // form the hamiltonian in the orthogonalized basis
  mat h_p = arma::trans(x) * hamiltonian * x;
  // diagonalize
  vec e;
  mat C_p;
  eig_sym(e, C_p, h_p);
  // Form the MO coefficients:
  mat C = x * C_p;
  mat mo_overlap = arma::trans(C) * overlap_matrix * C;

  return mo_overlap;
}

vec Cluster::eigenvalues()
{
  vec eigval;
  mat eigvec;
  mat overlap_matrix = overlapMatrix();
  mat hamiltonian = extendedHuckelHamiltonian();

  eig_sym(eigval, eigvec, overlap_matrix);

  // make the orthogonalization transformation
  mat x_diag = arma::diagmat(arma::pow(eigval, -0.5));
  mat x = eigvec * x_diag * arma::trans(eigvec);
  // form the hamiltonian in the orthogonalized basis
  mat h_p = arma::trans(x) * hamiltonian * x;
  // diagonalize
  vec e;
  mat C_p;
  eig_sym(e, C_p, h_p);

  return e;
}

double Cluster::molecularEnergy()
{
  return 2 * sum(eigenvalues().subvec(0, countElectronPairs() - 1));
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
