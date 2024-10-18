#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>
#include <iostream>
#include "Numerical.h"

int main(int argc, char* argv[])
{
  if (argc !=2)
  {
    printf("usage hw2_1 filename, for example hw2_1 example.txt\n");
    printf(" example.txt should have two lines and each line should have x alpha l \n");
    return EXIT_FAILURE;
  }
  //Create 1d Gaussain functions, the input is x_coor, alpha, l
  Gaussian g1, g2;
  std::string fname(argv[1]);
  try
  {
    ReadGaussianparameter(g1, g2, fname);
  }
  catch (std::invalid_argument &e)
  {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  //I use DE rule to calculate the integral, the input is func, func, numpoints, stepsize
  const double g1g2_de = DEInfinity(g1, g2, 3000, 0.001);

  //Trapezoidal rule for comparison
  //The numpoints and stepsize should be chosen based on the exponents of two gaussian functions
  const double g1g2_trapzd = Trapzd(g1, g2, 3000, 0.01);

  printf("1d numerical overlap integral (DE rule) between Gaussian functions is %1.17e\n", g1g2_de);
  printf("1d numerical overlap integral (Trapezoidal rule) between Gaussian functions is %1.17e\n", g1g2_trapzd);
  return EXIT_SUCCESS;
}
