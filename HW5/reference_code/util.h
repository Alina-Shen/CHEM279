#ifndef UTIL_H
#define UTIL_H

#include <armadillo>
double Combination(int n, int k);
double DoubleFactorial(int n);

double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb);
double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb);

// 2 electron integral of two primitive Gaussians
double I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaa, double sigmab);


double Overlap1st_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb);
void Overlap1st_3d(arma::vec &S_Ra, arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb);
// 1st derivative of 2 electron integral of two primitive Gaussians
void I2e_pG_1st(arma::vec &gammaAB_Ra, arma::vec &Ra, arma::vec &Rb, double sigmaa, double sigmab);

#endif // UTIL_H
