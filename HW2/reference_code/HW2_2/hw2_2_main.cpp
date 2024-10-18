#include <stdlib.h>
#include <stdexcept>
#include <armadillo>
#include <stdio.h>
#include "Shell.h"

using namespace std;

int main(int argc, char* argv[])
{
  if (argc !=2)
  {
  printf("usage hw2_2 filename, for example hw2_2 example.txt\n");
  printf(" example.txt should have two lines and each line should have x y z alpha l \n");
  printf(" Printed matrix is dim_a * dim_b, each dim follow the alphabetical order\n");
  return EXIT_FAILURE;
  }
  string fname(argv[1]);
  Shell sh1,sh2;
  try
  {
    ReadShellparameter(sh1, sh2, fname);
  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  /*
  const double s1s1 = Overlap_onedim(0.0, 0.0, 1.0, 1.0, 0, 0);
  printf("Analytical Integral %1.17e for s1s1\n", s1s1);
  const double s1p1 = Overlap_onedim(0.0, 0.0, 1.0, 1.0, 0, 1);
  printf("Analytical Integral %1.17e for s1p1\n", s1p1);
  const double s1s2 = Overlap_onedim(1.0, 0.0, 1.0, 1.0, 0, 0);
  printf("Analytical Integral %1.17e for s1s2\n", s1s2);
  const double s1p2 = Overlap_onedim(0.0, 1.0, 1.0, 1.0, 0, 1);
  printf("Analytical Integral %1.17e for s1p2\n", s1p2);
  */
  printf("Shell 1 has %d functions.\n", sh1.dim_func());
  sh1.printinfo();
  printf("Shell 2 has %d functions.\n", sh2.dim_func());
  sh2.printinfo();

  int dim1 = sh1.dim_func(), dim2 = sh2.dim_func();
  arma::mat Overlap_matrix(dim1, dim2,arma::fill::zeros);
  Eval_Ov(Overlap_matrix, sh1, sh2);
  Overlap_matrix.print("Overlap integral between Shell 1 and Shell 2");
  

  return EXIT_SUCCESS;
}

