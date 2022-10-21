/*
  UC Berkeley - MSSE Program
  Chem 279
  Fall 2022
  This file, util.cxx, contains utilities for reading molecule input

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "cluster.h"
#include "filereader.h"
#include "baseexception.h"
using std::string;
using std::vector;

/**
 * @brief reads files with molecule data
 *
 * @detail reads molecule files where the first line is the number of molecules and
 *         the remaining lines follow the format below:
 *         atomic_number x_coord y_coord z_coord
 *
 *         atomic_numbers must be > 0 and < MAX_ATOMIC_NUM (Set by MAX_ATOMIC_NUM environment variable)
 *         atomic_numbers can be restricted with ALLOWED_ATOMIC_NUM environment variable which is a list of allowed atomic_numbers
 *
 *         data is stored in a ClusterMap that stores the collective informaiton about all atoms of the specified type within the cluster
 *         described by the file being read.
 *
 *         the number of atoms described by the file must match the specified number of atoms in the first line
 *
 * @param atoms ClusterMap the information about all atom types in the cluster
 * @param string filename the name of the file containing molecule data for the cluster
 *
 **/
Cluster readfile(string &filename)
{
	std::ifstream infile(filename);
	if (infile.is_open())
	{
		int atomNum;
		double x_coord;
		double y_coord;
		double z_coord;
		// read first line
		int expectedAtoms;
		int miscnum;
		infile >> expectedAtoms;
		infile >> miscnum;
		int atomsCounted = 0;
		Cluster cluster(expectedAtoms);

		while (infile >> atomNum >> x_coord >> y_coord >> z_coord)
		{
			if (atomsCounted > expectedAtoms)
			{
				throw BaseException("AtomCountMismatch: Atom Count provided in input File must match number of atoms listed");
			}
			cluster.addAtom(atomsCounted, atomNum, x_coord, y_coord, z_coord, 0, 0);
			atomsCounted++;
		}

		infile.close();

		if (atomsCounted > expectedAtoms)
		{
			throw BaseException("AtomCountMismatch: Atom Count provided in input File must match number of atoms listed.");
		}
		return cluster;
	}
	else
	{
		throw std::invalid_argument("Can't open file to read.");
	}
}
