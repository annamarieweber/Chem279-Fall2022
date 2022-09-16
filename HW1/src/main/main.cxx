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

using std::string;
using std::vector;


int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cout << "Please input the file name!" << std::endl;
	}

	try {  
	  string filename = argv[1];
	  Cluster cluster = readfile(filename);
	  std::cout << cluster << std::endl;
	  double energyLJ = cluster.calcLJEnergy();
	  double analyticForce = cluster.calcAnalyticalForce();
	  std::cout << "LJ Potential Energy of Cluster: " << energyLJ << std::endl;
	  std::cout << "Analytic Force for Cluster: " << analyticForce << std::endl;
	}
	catch (std::invalid_argument &e) {
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	//catch (...) {
	//	std::cout << "Something wrong happened!" << std::endl;
	//	return EXIT_FAILURE;
	//}
	
	return EXIT_SUCCESS;
}



