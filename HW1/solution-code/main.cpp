#include <stdexcept>
#include "readfile.h"
#include "eval.h"

using namespace std;

int main(int argc, char* argv[])
{
	if (argc !=2)
	{
		printf("usage hw1 filename, for example hw1 example.txt");
		return EXIT_FAILURE;
	}
	string fname(argv[1]);
	vector<Atom> Atoms;
	try
	{
		Readatomsfromfile(Atoms, fname);
	}
	catch (invalid_argument &e)
	{
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}

	// cout << E_LJ(Atoms)<< endl;
	// arma::mat Force(3, 3);
	// F_LJ_fd(Force, Atoms, 0.0001);
	// Force.print("Force");
	// F_LJ_fd(Force, Atoms, 0.00001);
	// Force.print("Force");

	vector<Atom> opt_Atoms = Atoms;
	// normal steepest descent
	// Steepest_descend(opt_Atoms, Atoms, 1e-4, 1e-2);
	// steepest descent with golden section line search
	Steepest_descend_line_search(opt_Atoms, Atoms, 1e-4, 1e-2);
	for (auto x : opt_Atoms)
		cout << x << std::endl;

	return EXIT_SUCCESS;
}
