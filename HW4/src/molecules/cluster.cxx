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
#include <cmath>

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
  int numElectrons = countBasisFunctions();
  _p = numElectrons / 2 + numElectrons % 2;
  _q = numElectrons / 2;
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
    throw BaseException("InvalidConfig: number of electron pairs must be an integer");
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

mat Cluster::sBasisFunctions()
{

  mat basisfns(atomMatrix.n_rows, 5 * K + 2);
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    rowvec l_vals("0 0 0");
    int atomNum = atomMatrix(i, 0);
    mat coeffs(ATOM_COEFS_MAP[atomNum - 1]);
    Shell s(atomMatrix(i, 0), atomMatrix.row(i).cols(1, atomMatrix.n_cols - 1).t(), coeffs.col(1), coeffs.col(0), l_vals);
    ShellOverlapIntegral s_aa(s, s);
    mat norm_consts = 1.0 / sqrt(s_aa());
    basisfns.row(i).cols(0, K) = atomMatrix.row(i).cols(0, K);
    basisfns.row(i).cols(K + 1, 2 * K) = s.l_a(0).t();
    basisfns.row(i).cols(2 * K + 1, 3 * K) = s.alpha().t();
    basisfns.row(i).cols(3 * K + 1, 4 * K) = s.d_k().t();
    basisfns.row(i).cols(4 * K + 1, 5 * K) = (0 + 1) * norm_consts.col(0).t();
    basisfns.row(i).col(5 * K + 1) = 0;
  }
  return basisfns;
}

mat Cluster::gammaMatrix()
{
  mat basisfns = sBasisFunctions();
  int atomCount = atomMatrix.n_rows;
  mat gammas(atomCount, atomCount);
  mat result(atomCount, atomCount);

  for (int m = 0; m < atomCount; m++)
  {
    for (int v = 0; v < atomCount; v++)
    {
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      vec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K).t();
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      vec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K).t();

      vec dk_m_prime = dk_m % nk_m;
      vec dk_v_prime = dk_v % nk_v;

      // std::cout << "dk m" << std::endl;
      // std::cout << dk_m << std::endl;

      // std::cout << "dk v" << std::endl;
      // std::cout << dk_v << std::endl;

      // std::cout << "nk m" << std::endl;
      // std::cout << nk_m << std::endl;

      // std::cout << "nk v" << std::endl;
      // std::cout << nk_v << std::endl;

      // std::cout << "d'K m" << std::endl;
      // std::cout << dk_m_prime << std::endl;
      // std::cout << "d'k v" << std::endl;
      // std::cout << dk_v_prime << std::endl;
      // std::cout << "d'K m elem mult d'k v" << std::endl;
      // std::cout << dk_m_prime % dk_v_prime << std::endl;

      mat dk_m_prime_combos = dk_m_prime * dk_m_prime.t();
      mat dk_v_prime_combos = dk_v_prime * dk_v_prime.t();

      vec dk_m_prime_combos_vec = arma::vectorise(dk_m_prime_combos.t());
      vec dk_v_prime_combos_vec = arma::vectorise(dk_v_prime_combos.t());

      float r_12 = calcDistance(atomMatrix.row(m).cols(1, 3), atomMatrix.row(v).cols(1, 3));

      mat gamma_summation_matrix = dk_m_prime_combos_vec * dk_v_prime_combos_vec.t();

      vec alphas = basisfns.row(m).cols(2 * K + 1, 3 * K).t();
      mat alpha_copy(alphas.n_elem, alphas.n_elem, arma::fill::ones);
      alpha_copy.each_col() %= alphas;

      mat alpha_prime = alpha_copy + alpha_copy.t();

      vec betas = basisfns.row(v).cols(2 * K + 1, 3 * K).t();
      mat beta_copy(betas.n_elem, betas.n_elem, arma::fill::ones);
      beta_copy.each_col() %= betas;

      mat beta_prime = beta_copy + beta_copy.t();

      vec alpha_prime_vec = arma::vectorise(alpha_prime, 1).t();
      vec beta_prime_vec = arma::vectorise(beta_prime, 1).t();

      vec sigma_a = arma::pow(alpha_prime_vec, -1.0);
      mat sigma_a_copy(sigma_a.n_elem, sigma_a.n_elem, arma::fill::ones);
      sigma_a_copy.each_col() %= sigma_a;
      vec u_a = arma::pow(sigma_a * M_PI, 3.0 / 2.0);
      vec sigma_b = arma::pow(beta_prime_vec, -1.0);
      mat sigma_b_copy(sigma_a.n_elem, sigma_a.n_elem, arma::fill::ones);
      sigma_b_copy.each_col() %= sigma_b;
      vec u_b = arma::pow(sigma_b * M_PI, 3.0 / 2.0);
      mat v_sq = arma::pow(sigma_a_copy + sigma_b_copy.t(), -1.0);
      mat u = u_a * u_b.t();
      float r_ab_dist = calcDistance(atomMatrix.row(m).cols(1, 3), atomMatrix.row(v).cols(1, 3));
      mat t = v_sq * pow(r_ab_dist, 2);
      t.for_each([](mat::elem_type &val)
                 { val = erf(sqrt(val)); });

      if (m == v)
      {
        // std::cout << "u" << std::endl;
        // std::cout << u << std::endl;
        // std::cout << "circ sig b " << std::endl;
        // std::cout << circ_toeplitz(sigma_b) << std::endl;
        // std::cout << "vsq " << std::endl;
        // std::cout << 1 / v_sq << std::endl;
        // std::cout << "sqrt(2v^2)" << std::endl;
        // std::cout << sqrt(2.0 * v_sq) << std::endl;
        // std::cout << "sqrt(2/pi)" << std::endl;
        // std::cout << sqrt(2.0 / M_PI) << std::endl;
        // std::cout << "[0]^(0)" << std::endl;
        // std::cout << u % sqrt(2.0 * v_sq) * sqrt(2.0 / M_PI) << std::endl;
        // std::cout << "d_k * d_k' * d_l * d_l'" << std::endl;
        // std::cout << gamma_summation_matrix << std::endl;
        // std::cout << "d_k * d_k' * d_l * d_l * [0]^(0)'" << std::endl;
        // std::cout << gamma_summation_matrix % (u % sqrt(2.0 * v_sq) * sqrt(2.0 / M_PI)) << std::endl;
        // std::cout << "sum(d_k * d_k' * d_l * d_l * [0]^(0))" << std::endl;
        // std::cout << sum(sum(gamma_summation_matrix % (u % sqrt(2.0 * v_sq) * sqrt(2.0 / M_PI)))) << std::endl;
        // std::cout << "sum(d_k * d_k' * d_l * d_l) * sum([0]^(0))" << std::endl;
        // std::cout << sum(sum(gamma_summation_matrix)) * sum(sum((u % sqrt(2.0 * v_sq) * sqrt(2.0 / M_PI)))) << std::endl;
        gammas(m, v) = CONVERSION_FACTOR * sum(sum(gamma_summation_matrix % (u % sqrt(2.0 * v_sq) * sqrt(2.0 / M_PI))));
      }
      else
      {
        // std::cout << "alpha prime" << std::endl;
        // std::cout << alpha_prime << std::endl;
        // std::cout << "beta prime" << std::endl;
        // std::cout << beta_prime.t() << std::endl;
        // std::cout << "alpha prime vec" << std::endl;
        // std::cout << alpha_prime_vec << std::endl;
        // std::cout << "beta prime vec" << std::endl;
        // std::cout << beta_prime_vec << std::endl;
        // std::cout << "vsq " << std::endl;
        // std::cout << v_sq << std::endl;
        // std::cout << "1/vsq" << std::endl;
        // std::cout << sigma_a_copy + circ_toeplitz(sigma_b) << std::endl;
        // std::cout << "sig a cp" << std::endl;
        // std::cout << sigma_a_copy << std::endl;
        // std::cout << "circ sig b " << std::endl;
        // std::cout << circ_toeplitz(sigma_b) << std::endl;
        // std::cout << "u_a" << std::endl;
        // std::cout << u_a << std::endl;
        // std::cout << "u_b" << std::endl;
        // std::cout << u_b << std::endl;
        // std::cout << "u" << std::endl;
        // std::cout << u << std::endl;
        // std::cout << "sqrt(1/(r_a-r_b))" << std::endl;
        // std::cout << sqrt(1.0 / pow(r_ab_dist, 2)) << std::endl;
        // std::cout << "sqrt(t)" << std::endl;
        // std::cout << t << std::endl;
        // std::cout << "[0]^(0)" << std::endl;
        // std::cout << u % t * sqrt(1.0 / pow(r_ab_dist, 2)) << std::endl;
        // std::cout << "d_k * d_k' * d_l * d_l'" << std::endl;
        // std::cout << gamma_summation_matrix << std::endl;
        // std::cout << "d_k * d_k' * d_l * d_l * [0]^(0)'" << std::endl;
        // std::cout << gamma_summation_matrix % (u % t * sqrt(1.0 / pow(r_ab_dist, 2))) << std::endl;
        // std::cout << "sum(d_k * d_k' * d_l * d_l * [0]^(0))" << std::endl;
        // std::cout << sum(sum(gamma_summation_matrix % (u % t * sqrt(1.0 / pow(r_ab_dist, 2))))) << std::endl;
        // std::cout << "sum(d_k * d_k' * d_l * d_l) * sum([0]^(0))" << std::endl;
        // std::cout << sum(sum(gamma_summation_matrix)) * sum(sum((u % t * sqrt(1.0 / pow(r_ab_dist, 2))))) << std::endl;

        gammas(m, v) += CONVERSION_FACTOR * sum(sum(gamma_summation_matrix % (u % t * sqrt(1.0 / pow(r_ab_dist, 2)))));
      }
    }
  }

  return gammas;
}

vec Cluster::z_vals()
{
  mat valenceElectrons = atomMatrix.col(0);
  valenceElectrons.for_each([](vec::elem_type &val)
                            { val = VALENCE_ATOMIC_NUM[val - 1]; });
  return valenceElectrons;
}

mat Cluster::hCore()
{
  mat gamma = gammaMatrix();
  mat z = z_vals();
  int electron_count = sum(sum(z));

  // std::cout << "gamma: " << std::endl;
  // std::cout << gamma << std::endl;
  // std::cout << "z: " << std::endl;
  // std::cout << z << std::endl;
  vec gamma_z_tot = (gamma - diagmat(diagvec(gamma))) * z;
  // std::cout << "gamma * z : " << std::endl;
  std::cout << (gamma - diagmat(diagvec(gamma))) * z << std::endl;
  // std::cout << "gamma z tot : " << std::endl;
  // std::cout << gamma_z_tot << std::endl;

  vec diag_gamma_z = (z.col(0) - (vec(z.n_rows, fill::ones) / 2.0)) % diagvec(gamma);
  // std::cout << "diag gamma z: " << std::endl;
  // std::cout << diag_gamma_z << std::endl;

  vec expanded_diag_gamma_z(electron_count);
  vec expanded_gamma_z_tot(electron_count);
  vec ion_energy_electron_affinity(electron_count);
  vec atoms = atomMatrix.col(0);
  mat bonding_params(electron_count, electron_count);
  mat h_off_diag(electron_count, electron_count, fill::ones);
  h_off_diag -= diagmat(vec(electron_count, fill::ones));

  int k = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    ion_energy_electron_affinity.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) = vec(ATOM_TO_IONIZATION_ENERGY_ELECTRON_AFFINITY_PARAMS_MAPPING[atoms(i) - 1]);
    expanded_diag_gamma_z.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(diag_gamma_z(i));
    expanded_gamma_z_tot.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(gamma_z_tot(i));
    bonding_params.rows(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], electron_count).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
    bonding_params.cols(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(electron_count, VALENCE_ATOMIC_NUM[atoms(i) - 1]).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
    k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }

  h_off_diag %= (-1.0 / 2.0 * bonding_params % overlapMatrix());

  vec h_diag = -1 * ion_energy_electron_affinity - expanded_diag_gamma_z - expanded_gamma_z_tot;
  return diagmat(h_diag) + h_off_diag;
}

mat Cluster::cndo2FockMatrix()
{
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
