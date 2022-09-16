#include "readfile.h"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream> //to use istringstream

using namespace std;

void Readatomsfromfile(vector<Atom> &Atoms, string &fname)
{
	ifstream in;
	in.open(fname, ios::in);
	// cout << fname;
	string line;
	getline(in, line);
	int num_atoms = stod(line);
	// cout<<num_atoms;
	while (getline(in, line))
	{
		istringstream iss(line);
		Atom readatom;
		if (!(iss >> readatom.AtomicNumber >> readatom.r[0] >> readatom.r[1] >> readatom.r[2]))
		{
			throw invalid_argument("There is some problem with Atom format.");
		}
		if (readatom.AtomicNumber != 79)
		{
			throw invalid_argument("There are atoms other than Au atom.");
		}
		// cout << readatom << std::endl;
		Atoms.push_back(readatom);
	}
	if(Atoms.size() != num_atoms){
		throw invalid_argument("Number of atoms are not consistent ");
	}
	in.close();
}
