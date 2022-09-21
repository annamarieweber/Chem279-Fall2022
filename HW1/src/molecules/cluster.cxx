#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "cluster.h"
#include "baseexception.h"

using arma::mat;
using arma::rowvec;
using arma::vec;
using std::string;

// struct lj {
//   double operator()( mat x ) const {
//   //lj equation applied to matrix here
//   }
// };

// struct lj_prime {
//   double operator()(mat x) const{
//     //partial derivative of lj equation
//   }
// }

// steepest_descent(lj fn,lj_prime gradient,double l,mat x,int count, double threshold){
//     double x1 = x - l*gradient(x);
//     while(abs(x1[0] - x[0]) < threshold){
//       if(fn(x1) < fn(x)){
// 	l = 1.2*l;
// 	x=x1;
// 	count+=1;
//       }else{
// 	l=0.5*l;
// 	cout+=1;
//       }
//     }
//     return x;
// }

Cluster::Cluster(){};

Cluster::Cluster(int numAtoms)
{
  atomMatrix = mat(numAtoms, 4);

  // \todo get rid of these member vectors and be able to look up params for given atom on the fly (could store in constants)
  epsilons = vec(numAtoms);
  sigmas = vec(numAtoms);
}

bool Cluster::allAtomsMatch(int num)
{
  return all(atomMatrix.col(0) == num);
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
  // \todo: store literal values in definitions for readabiliy
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
