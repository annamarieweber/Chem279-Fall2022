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
		string filename = argv[1];
		Cluster cluster = readfile(filename);
		std::cout << cluster << std::endl;

		std::cout << "numBasisFunctions: " << cluster.countBasisFunctions() << std::endl;
		std::cout << "numElectrons: " << cluster.countElectronPairs() << std::endl;

		cluster.evalHw3();
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
