#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <string>
#include "cluster.h"
#include "filereader.h"
#include "baseexception.h"
using std::string;
using std::vector;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "Please input the file name!" << std::endl;
	}

	try
	{
		string testpa = "5.9699e-02   1.7844e-01  -1.5560e-01  -1.4706e-16  -5.3459e-16  -1.8099e-03   8.7444e-03  -1.8575e-16   1.1526e-16   2.3566e-03;1.7844e-01   7.6666e-01  -2.1533e-01  -1.0842e-16  -6.6102e-16   2.0008e-01  -2.4628e-01   3.3755e-16   5.4159e-16  -1.8099e-03;-1.5560e-01  -2.1533e-01   6.7364e-01  -2.0508e-16   9.3512e-16   2.4628e-01  -2.9756e-01   4.7761e-16  -1.5131e-15  -8.7444e-03;-1.4706e-16  -1.0842e-16  -2.0508e-16   5.0000e-01   3.1225e-16  -1.4312e-16   7.3565e-16   5.0000e-01   6.7307e-16   1.2653e-16;-5.3459e-16  -6.6102e-16   9.3512e-16   3.1225e-16   5.0000e-01  -8.4017e-17  -1.0038e-15  -6.0368e-16   5.0000e-01  -2.2517e-16;-1.8099e-03   2.0008e-01   2.4628e-01  -1.4312e-16  -8.4017e-17   7.6666e-01   2.1533e-01  -1.6960e-16   2.5516e-16   1.7844e-01;8.7444e-03  -2.4628e-01  -2.9756e-01   7.3565e-16  -1.0038e-15   2.1533e-01   6.7364e-01  -3.3195e-16   7.4095e-16   1.5560e-01;-1.8575e-16   3.3755e-16   4.7761e-16   5.0000e-01  -6.0368e-16  -1.6960e-16  -3.3195e-16   5.0000e-01  -2.4286e-16  -6.2465e-17;1.1526e-16   5.4159e-16  -1.5131e-15   6.7307e-16   5.0000e-01   2.5516e-16   7.4095e-16  -2.4286e-16   5.0000e-01   1.4999e-16;2.3566e-03  -1.8099e-03  -8.7444e-03   1.2653e-16  -2.2517e-16   1.7844e-01   1.5560e-01  -6.2465e-17   1.4999e-16   5.9699e-02";
		string testpb = "5.9699e-02   1.7844e-01  -1.5560e-01  -1.4706e-16  -5.3459e-16  -1.8099e-03   8.7444e-03  -1.8575e-16   1.1526e-16   2.3566e-03;1.7844e-01   7.6666e-01  -2.1533e-01  -1.0842e-16  -6.6102e-16   2.0008e-01  -2.4628e-01   3.3755e-16   5.4159e-16  -1.8099e-03;-1.5560e-01  -2.1533e-01   6.7364e-01  -2.0508e-16   9.3512e-16   2.4628e-01  -2.9756e-01   4.7761e-16  -1.5131e-15  -8.7444e-03;-1.4706e-16  -1.0842e-16  -2.0508e-16   5.0000e-01   3.1225e-16  -1.4312e-16   7.3565e-16   5.0000e-01   6.7307e-16   1.2653e-16;-5.3459e-16  -6.6102e-16   9.3512e-16   3.1225e-16   5.0000e-01  -8.4017e-17  -1.0038e-15  -6.0368e-16   5.0000e-01  -2.2517e-16;-1.8099e-03   2.0008e-01   2.4628e-01  -1.4312e-16  -8.4017e-17   7.6666e-01   2.1533e-01  -1.6960e-16   2.5516e-16   1.7844e-01;8.7444e-03  -2.4628e-01  -2.9756e-01   7.3565e-16  -1.0038e-15   2.1533e-01   6.7364e-01  -3.3195e-16   7.4095e-16   1.5560e-01;-1.8575e-16   3.3755e-16   4.7761e-16   5.0000e-01  -6.0368e-16  -1.6960e-16  -3.3195e-16   5.0000e-01  -2.4286e-16  -6.2465e-17;1.1526e-16   5.4159e-16  -1.5131e-15   6.7307e-16   5.0000e-01   2.5516e-16   7.4095e-16  -2.4286e-16   5.0000e-01   1.4999e-16;2.3566e-03  -1.8099e-03  -8.7444e-03   1.2653e-16  -2.2517e-16   1.7844e-01   1.5560e-01  -6.2465e-17   1.4999e-16   5.9699e-02";
		string filename = argv[1];
		Cluster cluster = readfile(filename);
		std::cout << cluster << std::endl;

		std::cout << "numBasisFunctions: " << cluster.countBasisFunctions() << std::endl;
		std::cout << "numElectronPairs: " << cluster.countElectronPairs() << std::endl;
		std::cout << "BasisFunctions:" << std::endl
				  << cluster.basisFunctions() << std::endl;
		std::cout << "OV_mat: " << std::endl
				  << cluster.overlapMatrix() << std::endl;
		std::cout << "H_mat: " << std::endl
				  << cluster.extendedHuckelHamiltonian() << std::endl;
		std::cout << "MO Coeffs:" << std::endl
				  << cluster.molecularOrbitalCoefficients() << std::endl;
		std::cout << "The molecule in file " << filename << " has energy: " << cluster.molecularEnergy() << std::endl;
		std::cout << "sBasisFunctions:" << std::endl
				  << cluster.sBasisFunctions() << std::endl;
		std::cout << "Gamma: " << std::endl
				  << cluster.gammaMatrix() << std::endl;
		std::cout << "FockMatrix: " << std::endl
				  << cluster.cndo2FockMatrix(mat(testpa), mat(testpb)) << std::endl;

		cluster.calcSCFEnergy(10e-6);
		double h = 1e-4;
	}
	catch (std::invalid_argument &e)
	{
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (BaseException &e)
	{
		e.displayError();
		return EXIT_FAILURE;
	}
	catch (const std::exception &ex)
	{
		std::cout << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cout << "An unknown Error occurred" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
