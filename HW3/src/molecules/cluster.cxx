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

using arma::field;
using arma::mat;
using arma::rowvec;
using arma::span;
using arma::vec;
using constants::ATOM_BASIS_FN_MAP;
using std::string;
using namespace constants;

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
  field<mat> n_k(b);
  mat r(b, K);
  vec l(b);
  vector<int> atom_num(b);

  int fn = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    vec l_vals(ATOM_ANGULAR_MOMENTUM_MAP[atomMatrix(i, 0) - 1]);
    for (int j = 0; j < l_vals.n_elem; j++)
    {
      // for (int k = 0; k < SECONDARY_QUANTUM_NUMBER_COMBINATIONS[l_vals(j)]; k++)
      // {
      //       }

      for (int k = 0; k < K; k++)
      {
        l(fn) = l_vals(j);
        r.row(fn) = atomMatrix.row(i).cols(1, 3);
        atom_num[fn] = atomMatrix(i, 0);
        Shell s(atomMatrix(i, 0), r.row(fn).t(), l_vals(j));
        ShellOverlapIntegral s_aa(s, s);

        n_k(fn) = 1.0 / sqrt(s_aa(k));

        fn++;
      }
    }
  }

  mat result(b, b);
  for (int e = 0; e < b; e++)
  {
    for (int f = 0; f < b; f++)
    {
      double tot = 0.0;

      Shell s_a(atom_num[e], r.row(e).t(), l(e));
      std::cout << r.row(e).t() << std::endl;
      Shell s_b(atom_num[f], r.row(f).t(), l(f));
      ShellOverlapIntegral o(s_a, s_b);
      std::cout << s_a.d_k() << std::endl;
      std::cout << s_b.d_k() << std::endl;
      std::cout << n_k(e) << std::endl;
      std::cout << n_k(f) << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << o(0) << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      result(e, f) += sum(sum((s_a.d_k() % s_b.d_k()) * (n_k(e) % n_k(f)) * o(0)));

      std::cout << (s_a.d_k() % s_b.d_k()) * (n_k(e) % n_k(f)) * o(0) << std::endl;
    }
  }
  std::cout << result << std::endl;

  std::cout << "printing self shell overlaps: " << std::endl;
  std::cout << n_k << std::endl;
  std::cout << n_k(1) << std::endl;
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
