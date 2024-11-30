#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "CNDO.h"

using namespace std;

int main(int argc, char *argv[])
{
  arma::mat H_STO3G, C_STO3G;
  H_STO3G.load("H_STO3G.txt");
  C_STO3G.load("C_STO3G.txt");
  // C_STO3G.print("C_STO3G");

  if (argc != 2)
  {
    printf("usage hw5 filename, for example ./hw5 example.txt\n");
    return EXIT_FAILURE;
  }
  string fname(argv[1]);
  try
  {
    Molecule_basis mol(fname, H_STO3G, C_STO3G);
    // mol.PrintAtoms();

    // // Check fdiff gamma and OV
    // int dim = mol.getnum_atoms();
    // arma::mat gamma1st_mat(3, dim *dim);
    // mol.eval_gamma1stmat(gamma1st_mat);
    // gamma1st_mat.print("gamma1st_mat");
    // arma::mat gamma_mat(dim, dim);
    // mol.eval_gammamat(gamma_mat);
    // gamma_mat.print("gamma_mat");
    // arma::vec R0 = mol.mAtoms[0].mAOs[0].get_R0();
    // R0(0) += 1e-5;
    // mol.mAtoms[0].set_R0(R0);
    // mol.eval_gammamat(gamma_mat);
    // arma::mat gammafd_mat = gamma_mat;
    // R0(0) -= 2e-5;
    // mol.mAtoms[0].set_R0(R0);
    // mol.eval_gammamat(gamma_mat);
    // gammafd_mat -=gamma_mat;
    // gammafd_mat /= 2.e-5;
    // gammafd_mat.print("gammafd_mat");
    // // mol.PrintAtoms();


    CNDO ourSCF(mol, 50, 1e-5);
    int ok = ourSCF.init();
    if(ok != 0) return EXIT_FAILURE;
    ok = ourSCF.run();
    if(ok != 0) return EXIT_FAILURE;
    double Energy = ourSCF.getEnergy();
    arma::mat gradient = ourSCF.getGradient();
    gradient.print("gradient");
   
   
    //Debug with numerical gradient
    printf("\nIf I use finite difference, change the x coordinate of second atom:\n");

    printf("The molecule in file %s has energy %f\n", argv[1], Energy);
    arma::vec R0 = mol.mAtoms[1].mAOs[0].get_R0();
    R0(0) += 1e-6;
    mol.mAtoms[1].set_R0(R0);
    CNDO ourSCF2(mol, 50, 1e-6);
    ourSCF2.init();
    ourSCF2.run();
    double Energy2 = ourSCF2.getEnergy();
    printf("The molecule in file %s has energy %f\n", argv[1], Energy2);
    R0(0) -= 2e-6;
    mol.mAtoms[1].set_R0(R0);
    CNDO ourSCF3(mol, 50, 1e-6);
    ourSCF3.init();
    ourSCF3.run();
    double Energy3 = ourSCF3.getEnergy();
    printf("The molecule in file %s has energy %f\n", argv[1], Energy3);
    printf("The molecule has gradient %f\n", (Energy2 - Energy3)/2.e-6);

  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

