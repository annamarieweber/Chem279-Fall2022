#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <string>
#include "atominfos.h"
#include "cluster.h"
#include "filereader.h"
#include "clustermap.h"
#include "baseexception.h"
#include "evaluation.h"
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
		double h = 1e-4;

		evaluateFDApproximation(cluster);
		double energyLJ = cluster.calcLJEnergy();
		mat analyticForce = cluster.calcAnalyticalForce();
		mat forwardFdForce = cluster.calcForwardFdForce(h);
		mat centralFdForce = cluster.calcCentralFdForce(h);
		std::cout << "Forward FD Force for Cluster: \n"
				  << forwardFdForce << std::endl;
		std::cout << "Forward FD Force Norm: " << norm(forwardFdForce) << std::endl;
		std::cout << "Central FD Force for Cluster: \n"
				  << centralFdForce << std::endl;
		std::cout << "Central FD Force Norm: " << norm(centralFdForce) << std::endl;
		std::cout << "Analytic Force for Cluster: \n"
				  << analyticForce << std::endl;
		std::cout << "Analytic Force Norm: " << norm(analyticForce) << std::endl;
		std::cout << "LJ Potential Energy of Cluster: " << energyLJ << std::endl;
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
		std::cout << "An unkknown Error occured" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
