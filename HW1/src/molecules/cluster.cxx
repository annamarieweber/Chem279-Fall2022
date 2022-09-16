#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "cluster.h"

using std::string;
using arma::mat;
using arma::vec;
using arma::rowvec;

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

Cluster::Cluster(int numAtoms){
  atomMatrix = mat(numAtoms,4);

  // \todo get rid of these member vectors and be able to look up params for given atom on the fly (could store in constants)
  epsilons = vec(numAtoms);
  sigmas = vec(numAtoms);
}

bool Cluster::allAtomsMatch(int num){
  return all(atomMatrix.col(0) == num);
}


double Cluster::calcDistance(mat a1, mat a2)
{

    double distance = 0;
    rowvec coordDists = square(a2.row(0)- a1.row(0));
    return sqrt(sum(coordDists));
}

double Cluster::sigma_ij(int i, int j){
  return sqrt(sigmas(i)*sigmas(j));
}

double Cluster::epsilon_ij(int i, int j){
  return sqrt(epsilons(i)*epsilons(j));
}

double Cluster::calculateLJ(double r_ij, int i, int j)
{
  double sigma_over_r_ij =  sigma_ij(i,j)/ r_ij;
  double term_r6 = pow(sigma_over_r_ij, 6.0);
  double term_r12 = term_r6 * term_r6;
  return epsilon_ij(i,j)*(term_r12-(2*term_r6));
}

double Cluster::calcLJEnergy(){

  double total_energy = 0.0;

  for(int i = 0; i < atomMatrix.n_rows; i++){
    for(int j = i + 1; j < atomMatrix.n_rows; j++){
      double dist_ij = calcDistance(atomMatrix.row(i).cols(1,3), atomMatrix.row(j).cols(1,3));
      total_energy += calculateLJ(dist_ij,i,j);
    }
  };

  return total_energy;
}

double Cluster::calcLJPrime(double r_ij, int i, int j){
  // \todo: store literal values in definitions for readabiliy
  double term_r6 =  pow(sigma_ij(i,j),6)/pow(r_ij,7);
  double term_r12 = -1 * pow(sigma_ij(i,j),12)/pow(r_ij,13);
  return epsilon_ij(i,j)*12*(term_r12 + term_r6);
}

mat Cluster::calcAnalyticalForce(){
  mat analytical_force(atomMatrix.n_rows,3);
  for(int i = 0; i < atomMatrix.n_rows; i++){
    for(int j = 0; j < atomMatrix.n_rows; j++){
      if(i != j){
	double dist_ij = calcDistance(atomMatrix.row(i).cols(1,3), atomMatrix.row(j).cols(1,3));
	analytical_force(i,0) += calcLJPrime(dist_ij, i, j)*((atomMatrix(i,1)-atomMatrix(j,1))/dist_ij);
	analytical_force(i,1) += calcLJPrime(dist_ij, i, j)*((atomMatrix(i,2)-atomMatrix(j,2))/dist_ij);
	analytical_force(i,2) += calcLJPrime(dist_ij, i, j)*((atomMatrix(i,3)-atomMatrix(j,3))/dist_ij);
      }
    }
  }
  return analytical_force;
}

void Cluster::addAtom(int index, int atomNum, double x, double y, double z, double e, double s){
  atomMatrix(index, 0) = atomNum;
  atomMatrix(index, 1) = x;
  atomMatrix(index, 2) = y;
  atomMatrix(index, 3) = z;
  epsilons(index) = e;
  sigmas(index) = s;
}

std::ostream& operator<<(std::ostream& os, const Cluster& c){
  os << "Atom Matrix: \n"<< c.atomMatrix << "\n";
  
  return os;
}

