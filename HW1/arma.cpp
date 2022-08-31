#include <iostream>
#include <armadillo>

int main() {
	arma::mat A = arma::randu(6,5);
	std::cout << "A:\n" << A << std::endl;
	arma::vec b = arma::randu(5);
	std::cout << "b:\n" << b << std::endl;
	std::cout << "b:\n" << b.rows(1,3) << std::endl;
	arma::vec v(3);
	v[0] = 3; v[1] = 6; v[2] = 9;
	for (int i = 0; i < 3; i++) {
		std::cout << "v[" << i << "] = " << v[i] << std::endl;
	}
	arma::mat zero; 
	zero.zeros(6,4);
	zero.print("zero");
	arma::mat init = {{1/20., 0},{49./180., 0.5*(1-std::sqrt(3./7.))}};
	init.print("init matrix");
	return 0;
}

