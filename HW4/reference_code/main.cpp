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
    printf("usage hw4 filename, for example hw4 example.txt\n");
    return EXIT_FAILURE;
  }
  string fname(argv[1]);
  try
  {
    Molecule_basis mol(fname, H_STO3G, C_STO3G);
    // AO H1s = mol.mAtoms[0].mAOs[0];
    // AO H2s = mol.mAtoms[1].mAOs[0];
    // double temp2 = Eval_2eI_sAO(H1s, H2s);
    // printf("Overlap of H1s and H2s %1.10f\n", temp2);
    // int dimatoms = mol.getnumAOs();
    // arma::mat gamma_mat(dimatoms, dimatoms);
    // mol.eval_gammamat(gamma_mat);
    // gamma_mat.print("gamma_mat");
    // int dim = mol.getnumAOs();
    // arma::mat OV_mat(dim, dim);
    // mol.eval_OVmat(OV_mat);
    // OV_mat.print("OV_mat");
    // mol.PrintAtoms();
    CNDO ourSCF(mol, 50, 1e-5);
    int ok = ourSCF.init();
    if(ok != 0) return EXIT_FAILURE;
    ok = ourSCF.run();
    if(ok != 0) return EXIT_FAILURE;
    double Energy = ourSCF.getEnergy();
    printf("The molecule in file %s has energy %f\n", argv[1], Energy);

  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

