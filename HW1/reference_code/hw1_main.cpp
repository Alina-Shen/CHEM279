#include <iostream>
#include <vector>
#include <stdexcept>
#include <armadillo>
#include <string>
#include <math.h>

#include "Atom.h"
#include "eval.h"

using namespace std;

int main(int argc, char* argv[]) {
   
    if (argc !=2)
    {
    	printf("Usage: hw1 <filename>, for example hw1 example.txt\n");
    	return EXIT_FAILURE;
    }
	
    string fname = argv[1];
    vector<Atom> atoms;
   
    //read in the atomic number and the xyz coordinates from file 
    try {
        read_atoms_from_file(atoms, fname);
    } catch (invalid_argument &e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    
    for (int i = 0; i < atoms.size(); i++) {
        cout << atoms[i] << endl;
    }

    //write the xyz values to an output file
    string ofname = "output.txt";
    write_atoms_to_file(atoms, ofname);
    double E = E_LJ(atoms);
    cout << "E_LJ = " << E << endl;
	
    //Calculate analytical force
    arma::mat F;
    F_LJ_analytical(F, atoms);
    F.print("F_LJ_analytical");

    //Calculate force with finite difference method with different step size
    for (int i = 1;i < 5;i++){
	double fdstepsize = pow(10,-i);
	cout << "Stepsize for finite difference:" << fdstepsize << endl;
    	F_LJ_forward_difference(F, atoms, fdstepsize);
    	F.print("F_LJ_forward_difference");
    	F_LJ_central_difference(F, atoms, fdstepsize);
    	F.print("F_LJ_central_difference");
    }
   
    vector<Atom> opt_atoms = atoms;
		
    //standard steepest descent
    cout << "start steepest descent" << endl;
    Steepest_descent(opt_atoms, atoms, 1e-4, 0.3, 1e-2);
    cout << "Optimized structure:" << endl;
    for (int i = 0; i < opt_atoms.size(); i++) {
        cout << opt_atoms[i] << endl;
    }
    
    //steepest descent with golden section line search
    cout << "start steepest descent with golden section line search" << endl;
    Steepest_descent_line_search(opt_atoms, atoms, 1e-4, 0.3, 1e-2, 1e-2);
    cout << "Optimized structure:" << endl;
    for (int i = 0; i < opt_atoms.size(); i++) {
        cout << opt_atoms[i] << endl;
    }

    return 0;
}
