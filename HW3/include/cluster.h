#ifndef CLUSTER
#define CLUSTER
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>

using arma::mat;
using arma::vec;
using std::string;

class Cluster
{
private:
  mat atomMatrix;
  vec epsilons;
  vec sigmas;

public:
  Cluster();

  /**
   * @brief Cluster constructor
   * @detail Creates a cluster of atoms represented as a matrix with numAtoms rows and 4 columns
   * where the first column is the atomic number and the remaining there columns are the x, y, and z coordinates
   * of the given atom respectively.
   * @param numAtoms int: the number of atoms to be included in the cluster
   **/
  Cluster(int numAtoms);

  /**
   * @brief Check if all atoms in cluster match provided atomic number
   * @param num int the atomic number to match
   * @return bool
   **/
  bool allAtomsMatch(int num);

  /**
   * @brief Calculate the number of basis functions N for the described cluster
   * @return int
   */
  int countBasisFunctions();

  /**
   * @brief Calculate the number of electron pairs n for the described cluster
   * @return int
   */
  int countElectronPairs();

  void printBasisFuncs();

  void evalHw3();

  /**
   * @brief calculates the distance between two atoms represented by 1d matrices with 3 values for x,y,and z respectively
   * @param a1 (mat): matrix representing the first atoms
   * @param a2 (mat): matrix representing the second value
   * @return double: the distance between the two atoms described as having positions a1 and a2
   **/
  double calcDistance(mat a1, mat a2);

  /**
   * @brief calculates the geometric mean of the sigma values for the ith and jth atoms in the cluster
   * @param i (int): i index
   * @param j (int): j index
   * @return double: geometric mean of sigmas(i) and sigmas(j)
   **/
  double sigma_ij(int i, int j);

  /**
   * @brief calculates the geometric mean of the epsilon value for the ith and jth atoms in the cluster
   * @param i (int): i index
   * @param j (int): j index
   * @return double: geometric mean of epsilons(i) and epsilons(j)
   **/
  double epsilon_ij(int i, int j);

  /**
   * @brief calculates the pairwise LJ potential energy between atom i and j
   * @param r_ij (double): the distance between a1 and a2
   * @param i (int): the i index
   * @param j (int): the j index
   * @returns double: the LJ potential energy between atom i and j
   * @throws BaseException (InvalidAtomPosition) - thrown when provided r_ij == 0;
   **/
  double calcLJ(double r_ij, int i, int j);

  /**
   * @brief calculates the Lenard jones potential energy of the cluster
   * @return double the total potential energy
   **/
  double calcLJEnergy();

  /**
   * @brief calculates the analytical force between two atoms
   * @param r_ij (double): the distance between a1 and a2
   * @param i (int): the i index
   * @param j (int): the j index
   * @returns double: the analytical force between a1 and a2
   **/
  double calcLJPrime(double r_ij, int i, int j);

  /**
   * @brief calculates the analytical force within a cluster
   * @return mat: a matrix containing the total analytical force in the x y and z directions respectively for each atom in the cluster
   **/
  mat calcAnalyticalForce();

  /**
   * @brief calculates the forward finite differences approximation of the force within a cluster
   * @param h (double): the stepsize to use
   * @return mat: a matrix containing the approximate force in the x y and z directions respectively for each atom in the cluster
   */
  mat calcForwardFdForce(double h);

  /**
   * @brief calculates the forward finite differences approximation of the force within a cluster
   * @param h (double): the stepsize to use
   * @return mat: a matrix containing the approximate force in the x y and z directions respectively for each atom in the cluster
   */
  mat calcCentralFdForce(double h);

  /**
   * @brief inserts an atom into the cluster
   * @param index (int): the index where the atom should be inserted
   * @param atomNum (int): the atomic number representing the atom type
   * @param x (double): the x coordinate of the atom
   * @param y (double): the y coordinate of the atom
   * @param z (double): the z coordinate of the atom
   * @param e (double): the epsilon corresponding to the atom at the specified index
   * @param s (double): the sigma corresponding to the atom at the specified index
   **/
  void addAtom(int index, int atomNum, double x, double y, double z, double e, double s);

  /**
   * @brief overloading the << operator
   **/
  friend std::ostream &operator<<(std::ostream &os, const Cluster &c);
};
#endif
