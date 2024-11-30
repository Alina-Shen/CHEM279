#include <stdexcept>
#include <math.h>
#include "util.h"
#include "AO.h"
#include <cassert>
#include <iostream>

using namespace std;

double Combination(int n, int k)
{
  if (n < k || k < 0)
  {
    throw invalid_argument("Comination: the number of elements should be bigger than selection numbers AND two numbers should be positive\n");
    return EXIT_FAILURE;
  }
  double result = 1e308;
  if (pow(result, 1.0 / k) < n)
  {
    throw invalid_argument("The Combination number may be bigger than the maxium of double precision\n");
    return EXIT_FAILURE;
  }
  double n_d = (double)n;
  result = 1.0;
  int k_small = min(k, n - k);
  for (double j = (double)k_small; j > 0; j--)
  {
    result /= j;
    result *= n_d;
    n_d--;
  }
  return result;
}

double DoubleFactorial(int n)
{
  if (n < -1)
    throw invalid_argument("DoubleFactorial: Input should be 1 at least\n");
  if (n == 0 || n == -1)
    return 1;
  double result = 1;
  while (n > 1)
  {
    result *= n;
    n -= 2;
  }
  return result;
}

double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb)
{
  // double x = Combination(3, 1);
  // double y = DoubleFactorial(4);
  double prefactor = exp(-alphaa * alphab * (xa - xb) * (xa - xb) / (alphaa + alphab)) * sqrt(M_PI / (alphaa + alphab));
  double xP = (alphaa * xa + alphab * xb) / (alphaa + alphab);

  double result = 0.0;
  for (int i_index = 0; i_index <= la; i_index++)
    for (int j_index = 0; j_index <= lb; j_index++)
    {
      if ((i_index + j_index) % 2 == 1)
        continue;
      double C_part = Combination(la, i_index) * Combination(lb, j_index);
      double DF_part = DoubleFactorial(i_index + j_index - 1);
      double numerator = pow(xP - xa, la - i_index) * pow(xP - xb, lb - j_index);
      // Caution: convert i_index + j_index to float!
      double dominator = pow(2 * (alphaa + alphab), double(i_index + j_index) / 2.0);
      double temp = C_part * DF_part * numerator / dominator;
      result += temp;
      // printf("%f %f %f %f  %f\n", C_part, DF_part, numerator, dominator, temp);
    }

  result *= prefactor;
  return result;
}

double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb)
{
  double Overlap = Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, lmna(0), lmnb(0)) *
                   Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, lmna(1), lmnb(1)) *
                   Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, lmna(2), lmnb(2));
  return Overlap;
}


double I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB){
  double U =  pow( M_PI*M_PI* sigmaA * sigmaB, 1.5);
  double V2 =  1.0/ (sigmaA +sigmaB);
  // cout << sigmaA << "  " << sigmaB << endl;
  // Ra.print("Ra");
  // Rb.print("Rb");
  double Rd= arma::norm(Ra -Rb, 2);
  if(Rd == 0.0)
    return U * 2.0 * sqrt(V2/M_PI);
  double srT = sqrt(V2) * Rd;;
  // cout << U << "  " << V2 << endl;
  // cout << Rd << "  " << srT << endl;
  double result = U /Rd * erf(srT);
  return result; 
}


double Overlap1st_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb)
{
  // double x = Combination(3, 1);
  // double y = DoubleFactorial(4);
  double First = 0, Second;

  if(la != 0)
    First = -la * Overlap_onedim(xa, xb, alphaa, alphab, la -1, lb);

  Second = 2* alphaa * Overlap_onedim(xa, xb, alphaa, alphab, la + 1, lb);

  return First + Second;
}

void Overlap1st_3d(arma::vec &S_Ra, arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb){
  assert(S_Ra.n_elem == 3);
  S_Ra(0) = Overlap1st_onedim(Ra(0), Rb(0), alphaa, alphab, lmna(0), lmnb(0)) *
            Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, lmna(1), lmnb(1)) *
            Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, lmna(2), lmnb(2));
  S_Ra(1) = Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, lmna(0), lmnb(0)) *
            Overlap1st_onedim(Ra(1), Rb(1), alphaa, alphab, lmna(1), lmnb(1)) *
            Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, lmna(2), lmnb(2));
  S_Ra(2) = Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, lmna(0), lmnb(0)) *
            Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, lmna(1), lmnb(1)) *
            Overlap1st_onedim(Ra(2), Rb(2), alphaa, alphab, lmna(2), lmnb(2));
  // S_Ra.print();
}

void I2e_pG_1st(arma::vec &gammaAB_Ra, arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB){
  double Rd= arma::norm(Ra -Rb, 2);
  if(Rd == 0.0)
    gammaAB_Ra.zeros();
  double U =  pow( M_PI*M_PI* sigmaA * sigmaB, 1.5);
  double V2 =  1.0/ (sigmaA +sigmaB);
  double T = V2 * Rd *Rd;;
  double result = U /Rd /Rd;
  result *=( 2 *sqrt(V2/M_PI) *exp(-T) - erf(sqrt(T)) /Rd) ;
  gammaAB_Ra = Ra -Rb;
  gammaAB_Ra *= result;
}

